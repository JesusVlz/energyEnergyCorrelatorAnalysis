#!/usr/bin/env python3
"""Render descriptive MC response diagnostics from already extracted arrays.

The public entry point is ``render(data, output_dir)``. Matrices always have
reconstructed rows and truth columns. This module does not access ROOT files,
fit models, propagate uncertainties, or decide whether physics closure passes.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.patches import Rectangle


FOLD_COLORS = {'train': '#0072B2', 'validation': '#D55E00'}
FOLD_LABELS = {'train': 'Training fold', 'validation': 'Validation fold'}
BASE_SCOPE = 'OO MC diagnostic  |  R = 0.4  |  particle pT > 2 GeV  |  inclusive activity'
MATRIX_SCOPE = ('Reconstructed rows, truth columns; regular stored bins only. '
                'No unfolding or efficiency inference.')
EEC_SCOPE = ('Generator-particle EEC with reconstructed/truth jets. '
             'Raw weighted sums; no correction or unfolding applied here.')
SHAPE_SCOPE = ('Each spectrum is normalized by its own sum and angular-bin width. '
               'Ratios are descriptive; no precision errors or closure claim.')
MOMENT_SCOPE = ('Uncentered within-jet second moments; this is not a centered '
                'event covariance or an uncertainty matrix.')


def _divide(numerator, denominator):
    numerator, denominator = np.broadcast_arrays(numerator, denominator)
    result = np.full(numerator.shape, np.nan, dtype=float)
    np.divide(numerator, denominator, out=result,
              where=np.isfinite(numerator) & np.isfinite(denominator) & (denominator != 0))
    return result


def _positive_norm(arrays):
    values = np.concatenate([np.asarray(a, dtype=float).ravel() for a in arrays])
    finite = values[np.isfinite(values)]
    nonzero = np.abs(finite[finite != 0])
    if not len(nonzero):
        return Normalize(0, 1)
    if np.any(finite < 0):
        bound = float(np.max(np.abs(finite)))
        return SymLogNorm(linthresh=max(float(np.min(nonzero)), bound * 1e-6),
                          vmin=-bound, vmax=bound, base=10)
    lower, upper = float(np.min(nonzero)), float(np.max(nonzero))
    if lower == upper:
        lower = upper / 10.
    return LogNorm(vmin=lower, vmax=upper)


def _fraction_norm(arrays):
    maxima = [float(np.nanmax(a)) for a in arrays if np.any(np.isfinite(a))]
    return Normalize(0., max(maxima or [1.]) or 1.)


def _cmap(name):
    result = plt.get_cmap(name).copy()
    result.set_bad('#eeeeee')
    return result


def _interval_labels(edges):
    return ['{:g}–{:g}'.format(low, high) for low, high in zip(edges[:-1], edges[1:])]


def _figure(title, rows=1, cols=1, figsize=None, footer=MATRIX_SCOPE, **kwargs):
    if figsize is None:
        figsize = (8.8 * cols, 6.8 * rows)
    fig, axes = plt.subplots(rows, cols, figsize=figsize, squeeze=False, **kwargs)
    fig.suptitle(title, fontsize=16, fontweight='semibold', y=.975)
    fig.text(.5, .935, BASE_SCOPE, ha='center', va='center', fontsize=10, color='#404040')
    fig.text(.5, .024, footer, ha='center', va='bottom', fontsize=9, color='#404040')
    return fig, axes


def _matrix(ax, array, edges, title, colorbar_label, norm=None, full=False,
            cmap='viridis', highlight=True):
    array = np.asarray(array, dtype=float)
    values = np.ma.masked_invalid(array)
    if isinstance(norm, LogNorm):
        values = np.ma.masked_less_equal(values, 0)
    artist = ax.imshow(values, origin='lower', interpolation='nearest', aspect='equal',
                       extent=(0, array.shape[1], 0, array.shape[0]),
                       cmap=_cmap(cmap), norm=norm)
    blocks = len(edges) - 1
    step = array.shape[0] // blocks if full else 1
    centers = (np.arange(blocks) + .5) * step
    labels = _interval_labels(edges)
    ax.set_xticks(centers)
    ax.set_yticks(centers)
    ax.set_xticklabels(labels, rotation=55, ha='right', fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    if full:
        for boundary in range(step, array.shape[0], step):
            ax.axvline(boundary, color='#999999', alpha=.5, lw=.55)
            ax.axhline(boundary, color='#999999', alpha=.5, lw=.55)
        ax.set_xlabel('Truth jet pT interval [GeV]\nΔR bins increase within each block', fontsize=10)
        ax.set_ylabel('Reco jet pT interval [GeV]\nΔR bins increase within each block', fontsize=10)
    else:
        ax.set_xlabel('Truth jet pT interval [GeV]')
        ax.set_ylabel('Reco jet pT interval [GeV]')
    if highlight:
        ax.add_patch(Rectangle((3 * step, 3 * step), 3 * step, 3 * step,
                               fill=False, edgecolor='#202020', linestyle='--', lw=1.3))
    ax.set_title(title, fontsize=12, pad=11)
    bar = ax.figure.colorbar(artist, ax=ax, fraction=.045, pad=.025)
    bar.set_label(colorbar_label, fontsize=9)
    bar.ax.tick_params(labelsize=8)
    return artist


def _shape(raw, jet_bin, dr_edges):
    block = np.asarray(raw, dtype=float).reshape(-1, len(dr_edges) - 1)[jet_bin]
    total = np.sum(block)
    if not np.isfinite(total) or total <= 0:
        return np.full(block.shape, np.nan)
    return block / total / np.diff(dr_edges)


def _step(ax, edges, values, **kwargs):
    return ax.step(edges, np.r_[values, values[-1]], where='post', **kwargs)


def _angular_axes(ax, edges, density=False):
    if np.all(edges > 0):
        ax.set_xscale('log')
    ax.set_xlim(edges[0], edges[-1])
    if density:
        ax.set_yscale('log')
        ax.set_ylabel('Unit-area EEC density / ΔR')
    else:
        ax.set_xlabel('ΔR')
        ax.axhline(1., color='#555555', linestyle=':', lw=1)
    ax.grid(True, which='major', color='#dddddd', linewidth=.55, alpha=.75)
    ax.tick_params(labelsize=9)


def _ratio_limits(axes, arrays):
    finite = np.concatenate([np.asarray(a).ravel() for a in arrays])
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        lower, upper = .5, 1.5
    else:
        lower = min(.9, float(np.min(finite)))
        upper = max(1.1, float(np.max(finite)))
        padding = .08 * (upper - lower)
        lower, upper = max(0., lower - padding), upper + padding
    for ax in axes:
        ax.set_ylim(lower, upper)


def _shape_grid(title, edges, reported, jet_edges):
    fig, axes = _figure(title, rows=2, cols=len(reported), figsize=(15.2, 8.7),
                        footer=SHAPE_SCOPE, sharex='col',
                        gridspec_kw={'height_ratios': [3, 1.4]})
    for column, jet_bin in enumerate(reported):
        _angular_axes(axes[0, column], edges, density=True)
        _angular_axes(axes[1, column], edges)
        axes[0, column].set_title('{:g} ≤ jet pT < {:g} GeV'.format(
            jet_edges[jet_bin], jet_edges[jet_bin + 1]), fontsize=12, pad=10)
        if column:
            axes[0, column].set_ylabel('')
    fig.subplots_adjust(left=.075, right=.985, bottom=.12, top=.845, hspace=.10, wspace=.27)
    return fig, axes


def render(data, output_dir):
    """Write PDF/PNG figures and a combined PDF; return the figure catalog.

    ``data`` follows the extraction/metrics contract: integer power keys,
    reco-row/truth-column matrices, and NaN for undefined conditional metrics.
    Inputs are read-only and never modified.
    """
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    jet_edges = np.asarray(data['jet_edges'], dtype=float)
    dr_edges = np.asarray(data['dr_edges'], dtype=float)
    reported = list(data['reported_bins'])
    folds, powers = list(data['folds']), list(data['powers'])
    records = []
    style = {'font.family': 'DejaVu Sans', 'font.size': 10,
             'axes.spines.top': False, 'axes.spines.right': False,
             'savefig.facecolor': 'white', 'figure.facecolor': 'white',
             'pdf.fonttype': 42, 'ps.fonttype': 42}
    with plt.rc_context(style), PdfPages(str(output / 'response_diagnostics.pdf')) as combined:
        def save(fig, stem, title, family):
            pdf, png = output / (stem + '.pdf'), output / (stem + '.png')
            fig.savefig(str(pdf))
            fig.savefig(str(png), dpi=155)
            combined.savefig(fig)
            plt.close(fig)
            records.append({'stem': stem, 'title': title, 'family': family,
                            'pdf': str(pdf), 'png': str(png)})

        # Raw and matched-conditional jet migration, with a common fold scale.
        jet_norm = _positive_norm([data['jet'][fold]['response'] for fold in folds])
        jet_fraction = _fraction_norm([data['metrics'][fold]['jet']['conditional'] for fold in folds])
        for fold in folds:
            title = 'Jet response — ' + FOLD_LABELS[fold]
            fig, axes = _figure(title, cols=2, figsize=(15.5, 8.5))
            _matrix(axes[0, 0], data['jet'][fold]['response'], jet_edges,
                    'Raw matched jet migration', 'Raw weighted bin sum', norm=jet_norm)
            _matrix(axes[0, 1], data['metrics'][fold]['jet']['conditional'], jet_edges,
                    'Conditional reconstructed bin at fixed truth bin',
                    'Fraction of matched truth-bin weight', norm=jet_fraction)
            fig.subplots_adjust(left=.07, right=.95, top=.85, bottom=.20, wspace=.36)
            save(fig, 'jet_response_' + fold, title, 'jet_response')

            title = 'Jet spectra — ' + FOLD_LABELS[fold]
            footer = ('Separate reconstructed and truth populations; raw weighted bin sums, '
                      'no bin-width scaling or efficiency interpretation.')
            fig, axes = _figure(title, figsize=(11.4, 6.9), footer=footer)
            ax = axes[0, 0]
            x = np.arange(len(jet_edges) - 1)
            for field, color, marker, label in [('measured', '#0072B2', 'o', 'Reconstructed jets'),
                                                ('truth', '#D55E00', 's', 'Truth jets')]:
                values = np.asarray(data['jet'][fold][field])
                ax.plot(x, values, color=color, marker=marker, ms=5, lw=1.6, label=label)
            ax.set_yscale('log')
            ax.set_xticks(x)
            ax.set_xticklabels(_interval_labels(jet_edges), rotation=40, ha='right')
            ax.set_xlabel('Jet pT interval [GeV]')
            ax.set_ylabel('Raw weighted bin sum')
            ax.axvspan(2.5, 5.5, color='#999999', alpha=.12)
            ax.grid(True, axis='y', alpha=.25)
            ax.legend(frameon=False)
            fig.subplots_adjust(left=.10, right=.97, top=.85, bottom=.23)
            save(fig, 'jet_spectra_' + fold, title, 'jet_spectra')

        # Full angular-and-jet migration, and its angularly summed jet blocks.
        for power in powers:
            norm = _positive_norm([data['eec'][fold][power]['response'] for fold in folds])
            fraction = _fraction_norm([data['metrics'][fold][power]['conditional'] for fold in folds])
            for fold in folds:
                title = 'EEC response, N = {} — {}'.format(power, FOLD_LABELS[fold])
                fig, axes = _figure(title, cols=2, figsize=(17., 9.), footer=EEC_SCOPE)
                _matrix(axes[0, 0], data['eec'][fold][power]['response'], jet_edges,
                        'Raw weighted migration', 'Raw weighted bin sum', norm=norm, full=True)
                _matrix(axes[0, 1], data['metrics'][fold][power]['conditional'], jet_edges,
                        'Conditional reconstructed (jet pT, ΔR) at fixed truth bin',
                        'Fraction of matched truth-bin weight', norm=fraction, full=True)
                fig.subplots_adjust(left=.075, right=.96, top=.85, bottom=.21, wspace=.33)
                save(fig, 'eec_response_n{}_{}'.format(power, fold), title, 'eec_response')
            title = 'EEC jet migration summed over ΔR, N = {}'.format(power)
            block_norm = _positive_norm([data['metrics'][fold][power]['jet_migration'] for fold in folds])
            fig, axes = _figure(title, cols=2, figsize=(15.5, 8.5), footer=EEC_SCOPE)
            for ax, fold in zip(axes[0], folds):
                _matrix(ax, data['metrics'][fold][power]['jet_migration'], jet_edges,
                        FOLD_LABELS[fold], 'Raw EEC-weighted sum over angular bins', norm=block_norm)
            fig.subplots_adjust(left=.07, right=.95, top=.85, bottom=.20, wspace=.36)
            save(fig, 'eec_jet_migration_n{}'.format(power), title, 'eec_jet_migration')

        quantities = [('jet', 'Jet', 'jet')] + [(p, 'EEC N = {}'.format(p), 'eec_n{}'.format(p)) for p in powers]
        for key, label, stem in quantities:
            title = label + ' matched-bin purity and stability'
            footer = ('Matched-entry fractions, not reconstruction efficiencies. '
                      'Shading marks the reported 40–160 GeV region; guard bins are retained.')
            fig, axes = _figure(title, cols=2, figsize=(14.5, 6.9), footer=footer)
            x = np.arange(len(jet_edges) - 1)
            for ax, field, panel_title, ylabel in zip(axes[0], ('purity', 'stability'),
                    ('Purity: fixed reconstructed jet bin', 'Stability: fixed truth jet bin'),
                    ('Diagonal / matched reco-bin sum', 'Diagonal / matched truth-bin sum')):
                ax.axvspan(2.5, 5.5, color='#999999', alpha=.14)
                for fold in folds:
                    ax.plot(x, data['metrics'][fold][key][field], marker='o', ms=5, lw=1.6,
                            color=FOLD_COLORS[fold], label=FOLD_LABELS[fold])
                ax.set_xticks(x)
                ax.set_xticklabels(_interval_labels(jet_edges), rotation=50, ha='right', fontsize=9)
                ax.set_ylim(0, 1.04)
                ax.set_xlabel('Jet pT interval [GeV]')
                ax.set_ylabel(ylabel)
                ax.set_title(panel_title, fontsize=12)
                ax.grid(True, axis='y', alpha=.25)
                ax.legend(frameon=False, fontsize=9)
            fig.subplots_adjust(left=.075, right=.985, top=.85, bottom=.25, wspace=.29)
            save(fig, stem + '_purity_stability', title, 'purity_stability')

            title = label + ' migration between reported and guard jet bins'
            footer = ('Fractions are conditional on matched entries in stored regular bins. '
                      'Blank cells are undefined; these are not inefficiency estimates.')
            fig, axes = _figure(title, rows=len(folds), cols=2, figsize=(13., 9.3), footer=footer)
            for row, fold in enumerate(folds):
                for column, field in enumerate(('feed_in', 'feed_out')):
                    ax = axes[row, column]
                    values = np.asarray(data['metrics'][fold][key][field])
                    artist = ax.imshow(np.ma.masked_invalid(values), vmin=0, vmax=1,
                                       cmap=_cmap('cividis'), origin='upper', aspect='auto')
                    ax.set_xticks(range(3))
                    ax.set_xticklabels(['< 40 GeV', '40–160 GeV', '≥ 160 GeV'])
                    ax.set_yticks(range(len(reported)))
                    ax.set_yticklabels([_interval_labels(jet_edges)[b] for b in reported])
                    ax.set_xlabel(('Truth' if field == 'feed_in' else 'Reco') + ' jet pT group')
                    ax.set_ylabel(('Reco' if field == 'feed_in' else 'Truth') + ' jet pT interval [GeV]')
                    ax.set_title(FOLD_LABELS[fold] + (' — sources of reco bins' if field == 'feed_in'
                                                      else ' — destinations of truth bins'), fontsize=11)
                    for i in range(values.shape[0]):
                        for j in range(values.shape[1]):
                            value = values[i, j]
                            if np.isfinite(value):
                                ax.text(j, i, '{:.1%}'.format(value), ha='center', va='center',
                                        color='white' if value < .5 else '#111111', fontsize=12)
                            else:
                                ax.text(j, i, '—', ha='center', va='center', color='#777777', fontsize=12)
                    bar = fig.colorbar(artist, ax=ax, fraction=.045, pad=.025)
                    bar.set_label('Matched-weight fraction', fontsize=9)
            fig.subplots_adjust(left=.10, right=.95, top=.85, bottom=.12, hspace=.50, wspace=.33)
            save(fig, stem + '_guard_migration', title, 'guard_migration')

            difference = 100. * (np.asarray(data['metrics']['validation'][key]['conditional'])
                                 - np.asarray(data['metrics']['train'][key]['conditional']))
            finite = difference[np.isfinite(difference)]
            bound = float(np.max(np.abs(finite))) if len(finite) else 1.
            bound = bound or 1.
            title = label + ' conditional matrix difference: validation − training'
            fig, axes = _figure(title, figsize=(10.5, 9.2),
                                footer='Descriptive fold comparison in percentage points; undefined matched columns remain blank.')
            _matrix(axes[0, 0], difference, jet_edges, 'Validation − training',
                    'Difference [percentage points]', norm=Normalize(-bound, bound),
                    full=key != 'jet', cmap='RdBu_r')
            fig.subplots_adjust(left=.16, right=.90, top=.85, bottom=.22)
            save(fig, stem + '_fold_conditional_difference', title, 'fold_comparison')

        # Independently normalized EEC shapes: reco/truth and fold comparisons.
        for power in powers:
            title = 'Reconstructed-jet and truth-jet EEC shapes, N = {}'.format(power)
            fig, axes = _shape_grid(title, dr_edges, reported, jet_edges)
            ratios = []
            for column, jet_bin in enumerate(reported):
                for fold in folds:
                    measured = _shape(data['eec'][fold][power]['measured'], jet_bin, dr_edges)
                    truth = _shape(data['eec'][fold][power]['truth'], jet_bin, dr_edges)
                    for values, linestyle, kind in ((measured, '-', 'reco jet'), (truth, '--', 'truth jet')):
                        _step(axes[0, column], dr_edges, values, lw=1.6, color=FOLD_COLORS[fold],
                              linestyle=linestyle, label=FOLD_LABELS[fold] + ', ' + kind)
                    ratio = _divide(measured, truth)
                    ratios.append(ratio)
                    _step(axes[1, column], dr_edges, ratio, lw=1.5, color=FOLD_COLORS[fold])
                axes[1, column].set_ylabel('Reco / truth' if column == 0 else '')
            _ratio_limits(axes[1], ratios)
            axes[0, 0].legend(frameon=False, fontsize=8, loc='best')
            save(fig, 'eec_reco_truth_shapes_n{}'.format(power), title, 'eec_shapes')

            for field, population in (('measured', 'Reconstructed-jet'), ('truth', 'Truth-jet')):
                title = '{} EEC shape comparison between folds, N = {}'.format(population, power)
                fig, axes = _shape_grid(title, dr_edges, reported, jet_edges)
                ratios = []
                for column, jet_bin in enumerate(reported):
                    shapes = {}
                    for fold in folds:
                        shapes[fold] = _shape(data['eec'][fold][power][field], jet_bin, dr_edges)
                        _step(axes[0, column], dr_edges, shapes[fold], lw=1.7,
                              color=FOLD_COLORS[fold], label=FOLD_LABELS[fold])
                    ratio = _divide(shapes['validation'], shapes['train'])
                    ratios.append(ratio)
                    _step(axes[1, column], dr_edges, ratio, lw=1.6, color='#333333')
                    axes[1, column].set_ylabel('Validation / training' if column == 0 else '')
                _ratio_limits(axes[1], ratios)
                axes[0, 0].legend(frameon=False, fontsize=9)
                save(fig, 'eec_{}_fold_shapes_n{}'.format(field, power), title, 'fold_comparison')

            for fold in folds:
                title = 'Within-jet EEC moment correlation, N = {} — {}'.format(power, FOLD_LABELS[fold])
                fig, axes = _figure(title, figsize=(10.5, 9.2), footer=MOMENT_SCOPE)
                _matrix(axes[0, 0], data['metrics'][fold][power]['correlation'], jet_edges,
                        'Diagonal-normalized within-jet second moments',
                        r'$M_{ij}/\sqrt{M_{ii} M_{jj}}$', norm=Normalize(-1, 1), full=True,
                        cmap='RdBu_r', highlight=False)
                # A moment matrix uses reconstructed EEC coordinates on both axes.
                axes[0, 0].set_xlabel('Reco jet pT interval [GeV]\nΔR bins increase within each block')
                axes[0, 0].set_ylabel('Reco jet pT interval [GeV]\nΔR bins increase within each block')
                fig.subplots_adjust(left=.16, right=.90, top=.85, bottom=.22)
                save(fig, 'eec_moment_correlation_n{}_{}'.format(power, fold), title, 'within_jet_moments')
    return records
