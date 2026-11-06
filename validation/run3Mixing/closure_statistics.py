"""Joint event statistics for the OO GenGen closure pilot.

Histogram values come from the existing analyzer / histogram manager. This
module only combines their categories and reduces the retained event vectors.
Errors condition on the scanned donor pool and fixed partner assignment. They
are not the final covariance including resampling the reused donor pool.
"""
import math
from collections import defaultdict

SAME, JM1, M11, JM2, M12, M22 = 0, 3, 4, 5, 6, 7
JP1, P11, JP2, P12, P22 = 8, 9, 10, 11, 12
PP, PH, HP, HH, ALL = 0, 1, 2, 3, 4
JET_BINS = ((40., 80.), (80., 120.), (120., 160.))
SELECTIONS = (("selected_integrated", -1., 100.), ("C0_30", 0., 30.),
              ("C30_50", 30., 50.), ("C50_90", 50., 90.))
UNCERTAINTY = ("Poisson(1) resampling of observed signal-background fingerprint clusters "
               "(physical signal events for empty/unavailable fingerprints), conditional on "
               "the scanned fixed donor pool and deterministic partner assignment; "
               "joint methods/truth covariance retained; donor-pool resampling, "
               "unobserved shared-background lineage and final chain covariance NOT_RUN")


def cluster_vectors(events, cluster_keys):
    """Retain observed repeated signal backgrounds in the same resampling unit."""
    result = defaultdict(lambda: defaultdict(float))
    for identity, vector in events.items():
        cluster = cluster_keys.get(identity, "physical:" + identity)
        for name, value in vector.items():
            result[cluster][name] += value
    return result



def coefficients():
    """Coefficients act on manager projections: projected PH already merges HP."""
    def terms(*items):
        result = defaultdict(float)
        for coefficient, pairing, origin in items:
            result[pairing, origin] += coefficient
        return dict(result)
    total = terms((1., SAME, ALL))
    pp, ph, hh = (terms((1., SAME, origin)) for origin in (PP, PH, HH))
    def background(pairings):
        return terms(*[(c, p, ALL) for c, p in zip((.5, .5, .5, .5, -1.), pairings)])
    mixed = background((JM1, JM2, M11, M22, M12))
    perp = background((JP1, JP2, P11, P22, P12))
    def combine(*weighted):
        result = defaultdict(float)
        for scale, mapping in weighted:
            for key, value in mapping.items():
                result[key] += scale * value
        return dict(result)
    result = dict(T=total, PP=pp, PH=ph, HH=hh,
                  Btrue=combine((1., ph), (1., hh)), Bmix=mixed, Bperp=perp,
                  Smix=combine((1., total), (-1., mixed)),
                  Sperp=combine((1., total), (-1., perp)),
                  truth_sum_residual=combine((1., total), (-1., pp), (-1., ph), (-1., hh)))
    for prefix, first, second, cross, same1, same2 in (
            ("mix", JM1, JM2, M12, M11, M22),
            ("perp", JP1, JP2, P12, P11, P22)):
        result[prefix + "_PMavg"] = terms((.5, first, PH), (.5, second, PH))
        result[prefix + "_HMavg_minus_cross"] = terms((.5, first, HH), (.5, second, HH), (-1., cross, ALL))
        result[prefix + "_same_avg"] = terms((.5, same1, ALL), (.5, same2, ALL))
    return result


COEFFICIENTS = coefficients()
METRICS = tuple(COEFFICIENTS)


def project_origin(raw, pairing, origin):
    """Apply the manager's origin merge exactly once to raw ordered categories."""
    if origin == ALL:
        return sum(raw.get((pairing, value), 0.) for value in (PP, PH, HP, HH))
    if origin in (PH, HP):
        return raw.get((pairing, PH), 0.) + raw.get((pairing, HP), 0.)
    return raw.get((pairing, origin), 0.)


def derive(raw):
    return {name: sum(scale * project_origin(raw, pairing, origin)
                      for (pairing, origin), scale in mapping.items())
            for name, mapping in COEFFICIENTS.items()}


def joint_moments(vectors):
    """Same physical-event outer products as run3Workflow/joint_statistics.py."""
    total = {name: 0. for name in METRICS}
    cov = {(first, second): 0. for i, first in enumerate(METRICS) for second in METRICS[i:]}
    for vector in vectors:
        for name in METRICS:
            value = float(vector.get(name, 0.))
            if not math.isfinite(value):
                raise ValueError("Nonfinite event contribution")
            total[name] += value
        for first, second in cov:
            cov[first, second] += vector.get(first, 0.) * vector.get(second, 0.)
    return total, cov


def covariance(cov, first, second):
    return cov.get((first, second), cov.get((second, first), 0.))


def comparison(total, cov, estimate, truth):
    """Joint residual and delta-method ratio; undefined ratios stay null."""
    numerator, denominator = total[estimate], total[truth]
    vnum, vden = covariance(cov, estimate, estimate), covariance(cov, truth, truth)
    cross = covariance(cov, estimate, truth)
    residual_variance = max(0., vnum + vden - 2. * cross)
    valid = math.isfinite(denominator) and denominator > 0.
    ratio = numerator / denominator if valid else None
    ratio_variance = ((vnum + ratio * ratio * vden - 2. * ratio * cross) /
                      denominator ** 2) if valid else None
    return dict(estimate=estimate, truth=truth, residual=numerator - denominator,
                residual_error=math.sqrt(residual_variance), valid_ratio=valid,
                ratio=ratio, ratio_error=math.sqrt(max(0., ratio_variance)) if valid else None)


def weighted_pairs(first, second=None, power=1):
    """Small combinatoric oracle used only by validation tests, never analysis."""
    if second is None:
        return sum((a * b) ** power for i, a in enumerate(first) for b in first[i + 1:])
    return sum((a * b) ** power for a in first for b in second)
