"""Numerical diagnostics, not an unfolding implementation. Inputs: bin integrals."""
import numpy as np

REPORT_BINS = ((40., 80.), (80., 120.), (120., 160.))


def finite(x, name='array'):
    x = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError(name + ' contains nonfinite values')
    return x


def edges(x):
    x = finite(x, 'edges')
    if x.ndim != 1 or len(x) < 2 or np.any(np.diff(x) <= 0):
        raise ValueError('Edges must be strictly increasing')
    return x


def physical_index(available, value):
    """Float_t serialization equivalence only; no nearest-edge approximation."""
    e = edges(available).astype(np.float32)
    if np.any(np.diff(e) <= 0):
        raise ValueError('Serialization-degenerate edges')
    found = np.flatnonzero(e == np.float32(value))
    if len(found) != 1:
        raise ValueError('Unsupported physical edge: ' + str(value))
    return int(found[0])


def interval_index(available, interval):
    first, last = [physical_index(available, v) for v in interval]
    if last != first + 1:
        raise ValueError('Requested interval must be a native bin; reproject explicitly')
    return first


def flattened_layout(flat_edges, jet_edges):
    flat, jets = edges(flat_edges), edges(jet_edges)
    blocks = len(jets)-1
    if (len(flat)-1) % blocks:
        raise ValueError('Flattened dimension disagrees with jet binning')
    n = (len(flat)-1)//blocks
    dr = flat[:n+1].copy()
    if n < 1 or dr[0] < 0:
        raise ValueError('Invalid native DeltaR block')
    for i in range(blocks):
        block = flat[i*n:(i+1)*n+1]
        if not np.allclose(block-block[0], dr-dr[0], rtol=2e-6, atol=2e-8):
            raise ValueError('DeltaR blocks are not repeated consistently')
    return n, dr


def covariance(x, n=None, check_psd=False):
    x = finite(x, 'covariance')
    if x.ndim != 2 or x.shape[0] != x.shape[1] or (n is not None and x.shape != (n,n)):
        raise ValueError('Wrong covariance dimensions/space')
    scale = max(float(np.max(np.abs(x), initial=0)), np.finfo(float).tiny)
    if not np.allclose(x,x.T,rtol=1e-7,atol=scale*1e-10):
        raise ValueError('Nonsymmetric covariance')
    x = (x+x.T)/2
    if np.min(np.diag(x), initial=0) < -scale*1e-10:
        raise ValueError('Negative covariance diagonal')
    if check_psd and np.min(np.linalg.eigvalsh(x), initial=0) < -scale*1e-8:
        raise ValueError('Materially indefinite covariance')
    return x


def response_operator(c, m, t, tolerance=5e-6):
    c, m, t = finite(c), finite(m), finite(t)
    if c.ndim != 2 or m.ndim != 1 or t.ndim != 1 or c.shape != (len(m),len(t)):
        raise ValueError('Response must use reco rows, truth columns')
    if any(np.any(v < 0) for v in (c,m,t)):
        raise ValueError('Signed response needs separate validated treatment')
    miss, fake = t-c.sum(axis=0), m-c.sum(axis=1)
    if np.any(miss < -tolerance*np.maximum(t,c.sum(axis=0))) or np.any(fake < -tolerance*np.maximum(m,c.sum(axis=1))):
        raise ValueError('Inconsistent response/marginal weights or bookkeeping')
    k = np.divide(c,t[None,:],out=np.zeros_like(c),where=t[None,:]>0)
    return k, miss, fake  # Preserve small negative roundoff residuals, never clip.


def diagnostics(k):
    k = finite(k)
    singular = np.linalg.svd(k, compute_uv=False)
    tol = max(k.shape)*np.finfo(float).eps*singular[0] if len(singular) else 0.
    rank = int(np.count_nonzero(singular > tol))
    nullity = k.shape[1]-rank
    condition = float(singular[0]/singular[-1]) if len(singular) and not nullity and singular[-1]>tol else 'infinite'
    return dict(singular_values=singular.tolist(), rank=rank, truth_nullity=nullity,
                condition_number=condition, rank_tolerance=tol)


def forward_fold(k, u, fake, exposure):
    """Explicit unmatched-jet exposure, not pair-background subtraction."""
    k, u, fake = finite(k), finite(u), finite(fake)
    if k.shape != (len(fake),len(u)) or not np.isfinite(exposure) or exposure < 0:
        raise ValueError('Invalid refolding dimensions/exposure')
    return k@u + exposure*fake


def normalization_bins(dr, region):
    dr = edges(dr)
    lo, hi = map(float, region)
    if not np.isfinite(lo+hi) or lo <= 0 or hi <= lo or lo < dr[0] or hi >= dr[-1]:
        raise ValueError('Unsupported positive normalization interval')
    # Same complete-bin endpoint convention as the AN macro's FindBin calls.
    a, b = np.searchsorted(dr,[lo,hi],side='right')-1
    mask = np.zeros(len(dr)-1,dtype=bool); mask[a:b+1] = True
    return mask, [float(dr[a]),float(dr[b+1])]


def normalize(y, v, widths, mask):
    y, widths = finite(y), finite(widths)
    v = covariance(v,len(y))
    mask = np.asarray(mask,dtype=bool)
    if y.ndim != 1 or widths.shape != y.shape or mask.shape != y.shape or np.any(widths<=0):
        raise ValueError('Invalid normalization arrays')
    total = y[mask].sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError('Nonpositive normalization')
    jac = np.diag(1/(total*widths))-np.outer(y/widths,mask.astype(float))/total**2
    return y/(total*widths), covariance(jac@v@jac.T), jac


def ratio(u, vu, t, vt, cross=None):
    u, t = finite(u), finite(t)
    vu, vt = covariance(vu,len(u)), covariance(vt,len(t))
    if u.shape != t.shape or np.any(t<=0):
        raise ValueError('Invalid ratio or nonpositive denominator')
    a, b = np.diag(1/t), np.diag(-u/t**2)
    v = a@vu@a.T + b@vt@b.T
    if cross is not None:
        cross = finite(cross)
        if cross.shape != vu.shape:
            raise ValueError('Wrong joint covariance dimensions')
        v += a@cross@b.T + b@cross.T@a.T
    return u/t, covariance(v)


def summarize_scan(data, jets, dr, region):
    """AN-style diagnostic only: truth Sumw2 proxy, no joint cross covariance."""
    n = len(dr)-1; widths = np.diff(dr)
    norm, effective = normalization_bins(dr,region)
    chi = np.zeros((2,20)); rows = []; masks = {}
    for bounds in REPORT_BINS:
        j = interval_index(jets,bounds); sl = slice(j*n,(j+1)*n)
        valid = norm.copy()
        for s in (0,1): valid &= data[f'truth_{s}'][sl] > 0
        ix = np.flatnonzero(valid)
        if not len(ix): raise ValueError('No common truth support in '+str(bounds))
        values = {}
        for s in (0,1):
            t, vt, _ = normalize(data[f'truth_{s}'][sl],np.diag(data[f'truth_var_{s}'][sl]),widths,norm)
            for k in range(1,21):
                u, vu, _ = normalize(data[f'u_{s}_{k}'][sl],data[f'cov_{s}_{k}'][sl,sl],widths,norm)
                values[s,k] = ratio(u[ix],vu[np.ix_(ix,ix)],t[ix],vt[np.ix_(ix,ix)])
        fixed = np.ones(len(ix),bool)
        for r,v in values.values(): fixed &= np.diag(v)>0
        if not fixed.any(): raise ValueError('No fixed positive-variance bin mask')
        masks[str(j)] = ix[fixed].tolist()
        for (s,k),(r,v) in values.items():
            statistic = float(np.sum((r[fixed]-1)**2/np.diag(v)[fixed]))
            chi[s,k-1] += statistic
            display, error = np.full(n,np.nan), np.full(n,np.nan)
            display[ix[fixed]] = r[fixed]; error[ix[fixed]] = np.sqrt(np.diag(v)[fixed])
            data[f'ratio_{s}_{k}_{j}'], data[f'error_{s}_{k}_{j}'] = display,error
            rows.append(dict(split=s+1,iteration=k,jet_low=bounds[0],jet_high=bounds[1],
                             chi2=statistic,n_bins=int(fixed.sum()),statistic='DIAGONAL_UNITY_DIAGNOSTIC'))
    data['chi2'], data['chi2_average'] = chi,chi.mean(axis=0)
    best = int(np.argmin(data['chi2_average']))+1
    candidate = dict(iterations=best,status='DIAGNOSTIC_MINIMUM_NOT_VALIDATED',frozen=False,
        boundary_minimum=best in (1,20),fixed_bin_masks=masks,effective_region=effective,
        criterion='Minimum two-fold mean of the sum of three per-jet-bin diagonal unity chi2 values',
        limitations=['Truth covariance uses diagonal Sumw2 proxy','No unfolded/truth cross covariance',
                     'No response-MC statistical uncertainty','No untouched final validation sample'])
    return rows,candidate
