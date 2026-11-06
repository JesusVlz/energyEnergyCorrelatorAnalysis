"""ROOT-free OO MC background-study figures on the input histogram scale.

The caller supplies histogram arrays and optional diagonal variances.  Every
ratio masks a nonpositive denominator.  Variances are propagated as independent
histogram uncertainties; no covariance or event-level significance is inferred.
"""
from pathlib import Path
import subprocess
import tempfile
import re
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


_DATASET_LABEL = "Pythia+Hijing OO 5.36 TeV"


_COLORS = {
    "T": "#555555", "PP": "#171717", "PH": "#ad7413", "HH": "#8c49a3",
    "Btrue": "#171717", "Bmix": "#2463ab", "Bperp": "#c43c35",
    "Smix": "#2463ab", "Sperp": "#c43c35",
}
_LABELS = {
    "T": "Same-jet total", "PP": "Pythia–Pythia", "PH": "Pythia–Hijing",
    "HH": "Hijing–Hijing", "Btrue": "True background",
    "Bmix": "Mixed-cone estimate", "Bperp": "Perpendicular estimate",
    "Smix": "Mixed-cone extraction", "Sperp": "Perpendicular extraction",
}
_FOOTER = (
    "Stored EEC scale; no separate curve normalization. Error bars: histogram statistical errors; "
    "sum/ratio errors use an uncorrelated approximation."
)
_STYLE = {
    "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": "white", "figure.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42,
}


def _array(channel, name):
    """Return a detached value array or None for an unavailable component."""
    value = channel.get("values", {}).get(name)
    if value is None:
        return None
    return np.asarray(value, dtype=float)


def _variance(channel, name):
    value = _array(channel, name)
    variance = channel.get("variances", {}).get(name)
    if value is None:
        return None
    if variance is None:
        return np.zeros_like(value)
    variance = np.asarray(variance, dtype=float)
    return np.where(np.isfinite(variance), np.maximum(variance, 0.0), np.nan)


def _ratio_arrays(numerator, denominator, numerator_variance=None,
                  denominator_variance=None):
    if numerator is None or denominator is None:
        return None, None
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0.0)
    result = np.full_like(numerator, np.nan)
    variance = np.full_like(numerator, np.nan)
    result[valid] = numerator[valid] / denominator[valid]
    numerator_variance = (np.zeros_like(numerator) if numerator_variance is None
                          else np.asarray(numerator_variance, dtype=float))
    denominator_variance = (np.zeros_like(denominator) if denominator_variance is None
                            else np.asarray(denominator_variance, dtype=float))
    variance[valid] = (numerator_variance[valid] / denominator[valid] ** 2
                       + numerator[valid] ** 2 * denominator_variance[valid]
                       / denominator[valid] ** 4)
    return result, variance


def _ratio(channel, numerator, denominator):
    return _ratio_arrays(_array(channel, numerator), _array(channel, denominator),
                         _variance(channel, numerator), _variance(channel, denominator))


def _centers(channel):
    edges = np.asarray(channel["edges"], dtype=float)
    return (edges[:-1] + edges[1:]) * 0.5


def _domain(channel):
    edges = np.asarray(channel["edges"], dtype=float)
    centers = _centers(channel)
    radius = float(channel.get("radius", 0.4))
    positive = centers[np.isfinite(centers) & (centers > 0) & (centers <= 2 * radius)]
    if not len(positive):
        return 1.e-4, 2 * radius
    lower_edges = edges[np.isfinite(edges) & (edges > 0)]
    low = min(positive.min() * 0.65, lower_edges.min()) if len(lower_edges) else positive.min() * 0.65
    high = 2 * radius
    display_range = channel.get("display_range")
    if display_range is not None:
        requested_low, requested_high = map(float, display_range)
        stored_low = lower_edges.min() if len(lower_edges) else positive.min()
        low = max(low, stored_low, requested_low)
        high = min(high, float(edges[-1]), requested_high)
        if not np.isfinite(low) or not np.isfinite(high) or high <= low:
            raise ValueError("display_range must overlap the positive stored Delta-r range")
    return low, high


def _axes(ax, channel, ylabel=None, xlabel=True):
    low, high = _domain(channel)
    ax.set_xscale("log")
    ax.set_xlim(low, high)
    ax.axvline(float(channel.get("radius", 0.4)), color="0.55", lw=0.8,
               ls=":", zorder=0)
    ax.grid(False, which="both", axis="both")
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(r"$\Delta r$")
    else:
        ax.tick_params(labelbottom=False)
    ax.tick_params(direction="in", top=False, right=False)


def _draw(ax, channel, values, variances=None, label=None, color=None,
          linestyle="-", marker=".", alpha=1.0, coefficient=1.0):
    if values is None:
        return False
    x = _centers(channel)
    y = np.asarray(values, dtype=float) * coefficient
    low, high = _domain(channel)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (x >= low) & (x <= high)
    if not np.any(valid):
        return False
    errors = None
    if variances is not None:
        errors = np.sqrt(np.maximum(0.0, np.asarray(variances, dtype=float))) * abs(coefficient)
        valid &= np.isfinite(errors)
        errors = errors[valid]
    ax.errorbar(x[valid], y[valid], yerr=errors, label=label, color=color,
                linestyle=linestyle, marker=marker, markersize=2.3, lw=1.05,
                elinewidth=0.5, capsize=0, alpha=alpha)
    return True


def _curve(ax, channel, name, **kwargs):
    kwargs.setdefault("label", _LABELS.get(name, name))
    kwargs.setdefault("color", _COLORS.get(name))
    return _draw(ax, channel, _array(channel, name), _variance(channel, name), **kwargs)


def _ratio_curve(ax, channel, numerator, denominator, **kwargs):
    values, variances = _ratio(channel, numerator, denominator)
    kwargs.setdefault("label", _LABELS.get(numerator, numerator))
    kwargs.setdefault("color", _COLORS.get(numerator))
    return _draw(ax, channel, values, variances, **kwargs)


def _yscale(ax, channel, names, force_linear=False):
    samples = []
    low, high = _domain(channel)
    centers = _centers(channel)
    for name in names:
        values = _array(channel, name)
        if values is None:
            continue
        valid = np.isfinite(values) & (centers >= low) & (centers <= high)
        samples.extend(values[valid])
    values = np.asarray(samples, dtype=float)
    nonzero = np.abs(values[values != 0])
    if force_linear or not len(nonzero):
        ax.set_yscale("linear")
    elif np.any(values < 0):
        # Keep signed subtraction results; a modest linear core avoids hiding zeros.
        threshold = max(np.quantile(nonzero, 0.12), nonzero.max() * 1.e-5)
        ax.set_yscale("symlog", linthresh=threshold, linscale=0.5)
        ax.axhline(0, color="0.6", lw=0.65)
    else:
        ax.set_yscale("log")
        # ROOT density histograms can contain zero bins; zero is retained in data,
        # while the log drawing naturally omits it.
    ax.margins(y=0.13)


def _legend(ax, ncol=1, fontsize=8):
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="lower left" if getattr(ax, "_oo_required_term_note", False) else "best", frameon=True, framealpha=0.88,
                  edgecolor="0.9", ncol=ncol, fontsize=fontsize,
                  title=_DATASET_LABEL, title_fontsize=fontsize)


def _unity(ax, value=1):
    ax.axhline(value, color="0.35", ls="--", lw=0.8, zorder=0)


def _split(fig, slot, channel, title, lower_ylabel):
    grid = slot.subgridspec(2, 1, height_ratios=[2.3, 1], hspace=0.045)
    upper = fig.add_subplot(grid[0])
    lower = fig.add_subplot(grid[1], sharex=upper)
    upper.set_title(title, loc="left", pad=5, fontweight="medium")
    _axes(upper, channel, "EEC (stored scale)", xlabel=False)
    _axes(lower, channel, lower_ylabel)
    return upper, lower


def _jet_selection(channel):
    """Keep complete caller labels; add a quantity to bare numeric ranges."""
    jet = str(channel.get("jet_label") or "%g–%g GeV" % tuple(channel["jet"]))
    plain = re.sub(r"[\s${}_\\]", "", jet).lower()
    if "jet" in plain or "pt" in plain:
        return jet
    return "Jet $p_T$: " + jet


def _centrality_selection(channel):
    label = str(channel.get("label", ""))
    return label if "centrality" in label.lower() else "Centrality: " + label


def _header(channel, page):
    inclusive = channel.get("scope") == "integrated"
    dataset = "Inclusive OO MC GenGen" if inclusive else "OO MC GenGen"
    selections = []
    if not inclusive:
        selections.append(_centrality_selection(channel))
    selections.extend([
        _jet_selection(channel),
        f"$p_T^{{ch}} > {channel['track_threshold']:g}$ GeV",
        f"$n={channel.get('exponent', 1):g}$, $R={channel.get('radius', 0.4):g}$",
    ])
    header = f"{dataset} — {page}\n" + "   |   ".join(selections)
    if channel.get("display_range") is not None:
        low, high = _domain(channel)
        header += f"\nCore view: ${low:g} < \\Delta r < {high:g}$"
    return header


def _required_term_messages(channel, prefix):
    """Identify absent or all-zero stored estimator inputs in the shown range."""
    required = [("J1", "JC₁"), ("J2", "JC₂"), ("same1", "C₁C₁"),
                ("same2", "C₂C₂"), ("cross", "C₁C₂")]
    centers = _centers(channel)
    low, high = _domain(channel)
    shown = np.isfinite(centers) & (centers > 0) & (centers >= low) & (centers <= high)
    zero, absent = [], []
    for key, label in required:
        values = _array(channel, prefix + "_" + key)
        if values is None or not np.any(shown & np.isfinite(values)):
            absent.append(label)
        elif np.all(values[shown & np.isfinite(values)] == 0):
            zero.append(label)
    messages = []
    if len(zero) == len(required):
        messages.append("All required cone terms are zero in the shown range.")
    elif zero:
        messages.append("Zero required terms: " + ", ".join(zero) + ".")
    if absent:
        messages.append("Unavailable required terms: " + ", ".join(absent) + ".")
    return messages


def _annotate_required_terms(ax, channel, prefixes=("mix", "perp")):
    lines = []
    for prefix in prefixes:
        method = "Mixed" if prefix == "mix" else "Perpendicular"
        for message in _required_term_messages(channel, prefix):
            lines.extend(textwrap.wrap(method + ": " + message, width=67))
    if lines:
        ax._oo_required_term_note = True
        ax.text(0.975, 0.975, "\n".join(lines), transform=ax.transAxes,
                fontsize=7.2, ha="right", va="top", color="#7f330b",
                bbox={"facecolor": "white", "edgecolor": "#ddbd9b",
                      "alpha": 0.94, "boxstyle": "round,pad=0.3"}, zorder=20)


def _page(channel, page):
    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.055, h_pad=0.06, wspace=0.06, hspace=0.10)
    fig.suptitle(_header(channel, page), fontsize=13)
    fig.supxlabel(_FOOTER + "  Dotted line: jet R.", fontsize=7.4)
    return fig, fig.add_gridspec(2, 2)


def _primary_page(channel):
    fig, grid = _page(channel, "truth and subtraction closure")
    top, bottom = _split(fig, grid[0, 0], channel, "Same-jet subevent decomposition", "Fraction / total")
    for name in ("T", "PP", "PH", "HH"):
        _curve(top, channel, name)
    for name in ("PP", "PH", "HH"):
        _ratio_curve(bottom, channel, name, "T")
    _yscale(top, channel, ("T", "PP", "PH", "HH"))
    bottom.set_ylim(-0.04, 1.08)
    _legend(top, ncol=2)
    residual = channel.get("decomposition_max_relative_residual")
    if residual is not None and np.isfinite(residual):
        top.text(0.03, 0.03, f"Max. decomposition residual: {residual:.2g}",
                 transform=top.transAxes, fontsize=7.5, color="0.35")

    top, bottom = _split(fig, grid[0, 1], channel, "Background closure", "Estimate / truth")
    for name in ("Btrue", "Bmix", "Bperp"):
        _curve(top, channel, name)
    for name in ("Bmix", "Bperp"):
        _ratio_curve(bottom, channel, name, "Btrue")
    _yscale(top, channel, ("Btrue", "Bmix", "Bperp"))
    _annotate_required_terms(top, channel)
    _unity(bottom)
    _legend(top)

    top, bottom = _split(fig, grid[1, 0], channel, "Signal extraction closure", "Extraction / truth")
    for name in ("PP", "T", "Smix", "Sperp"):
        _curve(top, channel, name, alpha=0.6 if name == "T" else 1.0)
    for name in ("Smix", "Sperp"):
        _ratio_curve(bottom, channel, name, "PP")
    _yscale(top, channel, ("PP", "T", "Smix", "Sperp"))
    _unity(bottom)
    _legend(top, ncol=2)

    top, bottom = _split(fig, grid[1, 1], channel, "Cone symmetry and method comparison", "Mixed / perpendicular")
    top.set_ylabel("Second cone / first cone")
    for prefix, color, label in (("mix", "#2463ab", "Mixed"), ("perp", "#c43c35", "Perpendicular")):
        _ratio_curve(top, channel, prefix + "_J2", prefix + "_J1",
                     color=color, label=label + ": jet–cone")
        _ratio_curve(top, channel, prefix + "_same2", prefix + "_same1",
                     color=color, label=label + ": same-cone", linestyle="--", marker=None)
    _unity(top)
    _ratio_curve(bottom, channel, "Bmix", "Bperp", color="#6a4592", label="Background")
    _ratio_curve(bottom, channel, "Smix", "Sperp", color="#15765f", label="Extracted signal", linestyle="--")
    _unity(bottom)
    _legend(top, ncol=2, fontsize=7.5)
    _legend(bottom, ncol=2, fontsize=7.2)
    return fig


def _component_page(channel):
    fig, grid = _page(channel, "background estimator components")
    for column, prefix, title, estimate in ((0, "mix", "Mixed cones", "Bmix"),
                                            (1, "perp", "Perpendicular cones", "Bperp")):
        top, bottom = _split(fig, grid[0, column], channel, title + ": signed estimator terms", "Term / true background")
        # Plot the signed terms actually summed, rather than silently rescaling shapes.
        specifications = [
            (prefix + "_jet_avg", 1.0, "½ (JC₁ + JC₂)", "#2463ab"),
            (prefix + "_cross", -1.0, "− C₁C₂", "#ad7413"),
            (prefix + "_same_avg", 1.0, "½ (C₁C₁ + C₂C₂)", "#8c49a3"),
        ]
        all_values = []
        for name, coefficient, label, color in specifications:
            _curve(top, channel, name, coefficient=coefficient, label=label, color=color)
            values, variance = _ratio(channel, name, "Btrue")
            _draw(bottom, channel, values, variance, coefficient=coefficient, label=label, color=color)
            source = _array(channel, name)
            if source is not None:
                all_values.append(source * coefficient)
        _curve(top, channel, estimate, label="Full estimate", color="#171717", linestyle="--", marker=None)
        # Negative cross-cone contribution is integral to the estimator.
        if all_values:
            nonzero = np.abs(np.concatenate(all_values))
            nonzero = nonzero[np.isfinite(nonzero) & (nonzero > 0)]
            if len(nonzero):
                top.set_yscale("symlog", linthresh=max(np.quantile(nonzero, 0.12), nonzero.max() * 1.e-5), linscale=0.5)
        top.axhline(0, color="0.6", lw=0.65)
        _unity(bottom, 0)
        _annotate_required_terms(top, channel, (prefix,))
        _legend(top, ncol=2, fontsize=7.6)

    top, bottom = _split(fig, grid[1, 0], channel, "Pythia–Hijing component closure", "Estimate / Pythia–Hijing")
    _curve(top, channel, "PH")
    for prefix, color, label in (("mix", "#2463ab", "Mixed"), ("perp", "#c43c35", "Perpendicular")):
        key = prefix + "_PMavg"
        _curve(top, channel, key, color=color, label=label + ": ½ (Pythia–C₁ + Pythia–C₂)")
        _ratio_curve(bottom, channel, key, "PH", color=color, label=label)
    _yscale(top, channel, ("PH", "mix_PMavg", "perp_PMavg"))
    _unity(bottom)
    _legend(top, fontsize=7.8)

    top, bottom = _split(fig, grid[1, 1], channel, "Hijing–Hijing component and cancellation", "Estimate / Hijing–Hijing")
    _curve(top, channel, "HH")
    for prefix, color, label in (("mix", "#2463ab", "Mixed"), ("perp", "#c43c35", "Perpendicular")):
        same = prefix + "_same_avg"
        residual = prefix + "_HMavg_minus_cross"
        _curve(top, channel, same, color=color, label=label + ": mean same-cone")
        _curve(top, channel, residual, color=color, label=label + ": Hijing–cone − cross-cone", linestyle="--", marker=None, alpha=0.7)
        _ratio_curve(bottom, channel, same, "HH", color=color, label=label + ": same-cone")
        # The full HH estimator contains its same-cone term plus this cancellation residual.
        if _array(channel, same) is not None and _array(channel, residual) is not None:
            numerator = _array(channel, same) + _array(channel, residual)
            variance = _variance(channel, same) + _variance(channel, residual)
            values, ratio_variance = _ratio_arrays(numerator, _array(channel, "HH"), variance, _variance(channel, "HH"))
            _draw(bottom, channel, values, ratio_variance, color=color,
                  label=label + ": full HH estimator", linestyle="--", marker=None)
    _yscale(top, channel, ("HH", "mix_same_avg", "perp_same_avg", "mix_HMavg_minus_cross", "perp_HMavg_minus_cross"))
    _unity(bottom)
    _legend(top, ncol=1, fontsize=7.2)
    return fig


def _save_pages(figures, output_base, atlas=None, component_png=True):
    """Render PNGs from the finalized PDF so backend layout cannot clip titles."""
    base = Path(output_base)
    base.parent.mkdir(parents=True, exist_ok=True)
    artifacts = {"pdf": str(base) + ".pdf", "png": str(base) + ".png", "pages": len(figures)}
    try:
        with PdfPages(artifacts["pdf"]) as pdf:
            for fig in figures:
                pdf.savefig(fig)
                if atlas is not None:
                    atlas.savefig(fig)
        targets = [(1, Path(artifacts["png"]))]
        if len(figures) > 1 and component_png:
            artifacts["components_png"] = str(base) + "_components.png"
            targets.append((2, Path(artifacts["components_png"])))
        with tempfile.TemporaryDirectory(prefix=".preview-", dir=base.parent) as scratch:
            for page, target in targets:
                prefix = Path(scratch) / str(page)
                subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page),
                                "-r", "135", "-png", "-singlefile", artifacts["pdf"],
                                str(prefix)], check=True, capture_output=True)
                prefix.with_suffix(".png").replace(target)
    finally:
        for fig in figures:
            plt.close(fig)
    return artifacts


def draw_channel(channel, output_base, atlas=None):
    """Save the two-page truth/closure and component dashboard.

    ``atlas`` may be a matplotlib.backends.backend_pdf.PdfPages instance. Return
    a dictionary containing PDF/PNG paths and page count. Empty channels get a
    single explicit status page if the caller chooses to include them.
    """
    with plt.rc_context(_STYLE):
        if channel.get("status") == "empty":
            fig, grid = _page(channel, "empty selection")
            ax = fig.add_subplot(grid[:, :])
            ax.axis("off")
            ax.text(0.5, 0.5, "No same-jet entries in this stored selection", ha="center", va="center", fontsize=18)
            return _save_pages([fig], output_base, atlas)
        return _save_pages([_primary_page(channel), _component_page(channel)], output_base, atlas)


def _group_label(channel, group):
    if group == "centrality":
        return str(channel.get("label", "%g–%g%%" % tuple(channel["centrality"])))
    if group == "jet":
        return str(channel.get("jet_label", "%g–%g GeV" % tuple(channel["jet"])))
    if group == "track":
        return "$p_T^{ch}>%g$ GeV" % channel["track_threshold"]
    raise ValueError("group must be 'jet', 'track', or 'centrality'")


def draw_overlays(channels, output_base, group, atlas=None):
    """Compare raw closure ratios as one selection dimension is varied.

    Four panels show background and extracted-signal closure for the mixed and
    perpendicular estimators. No normalization is applied. The caller groups
    compatible channels (same exponent and unchanged selection dimensions).
    Empty input returns a zero-page artifact dictionary without writing files.
    """
    if group not in ("jet", "track", "centrality"):
        raise ValueError("group must be 'jet', 'track', or 'centrality'")
    selected = [channel for channel in channels if channel.get("status", "populated") != "empty"]
    if not selected:
        return {"pages": 0}
    with plt.rc_context(_STYLE):
        first = selected[0]
        fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
        fixed = []
        inclusive = first.get("scope") == "integrated"
        if group != "centrality" and not inclusive:
            fixed.append(_centrality_selection(first))
        if group != "jet":
            fixed.append(_jet_selection(first))
        if group != "track":
            fixed.append("$p_T^{ch}>%g$ GeV" % first["track_threshold"])
        fixed.append("$n=%g$, $R=%g$" % (first.get("exponent", 1), first.get("radius", 0.4)))
        dataset = "Inclusive OO MC GenGen" if inclusive else "OO MC GenGen"
        comparison = {"jet": "jet momentum", "track": "track threshold", "centrality": "centrality"}[group]
        fig.suptitle(dataset + ": closure versus " + comparison + "\n" + "   |   ".join(fixed), fontsize=13)
        fig.supxlabel(_FOOTER + "  Dotted line: jet R.", fontsize=7.4)
        specifications = [
            (axes[0, 0], "Bmix", "Btrue", "Mixed-cone background", "Estimate / true background"),
            (axes[0, 1], "Bperp", "Btrue", "Perpendicular background", "Estimate / true background"),
            (axes[1, 0], "Smix", "PP", "Mixed-cone extracted signal", "Extraction / true signal"),
            (axes[1, 1], "Sperp", "PP", "Perpendicular extracted signal", "Extraction / true signal"),
        ]
        palette = plt.get_cmap("tab10")
        for ax, numerator, denominator, title, ylabel in specifications:
            _axes(ax, first, ylabel)
            ax.set_title(title, loc="left")
            for index, channel in enumerate(selected):
                _ratio_curve(ax, channel, numerator, denominator,
                             label=_group_label(channel, group), color=palette(index % 10),
                             linestyle=("-", "--", "-.", ":")[(index // 10) % 4], marker=".")
            _unity(ax)
            _legend(ax, ncol=2 if len(selected) > 5 else 1, fontsize=7.5)
        return _save_pages([fig], output_base, atlas, component_png=False)
