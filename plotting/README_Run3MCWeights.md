# Run-3 MC event weights: derivation and Jussi-style QA

## Files and installation

The public `plotting/deriveMonteCarloWeights.C` entry point now calls the
Run-3 implementation in `deriveMonteCarloWeightsRun3.C`. The previous Run-2
implementation remains callable as `deriveMonteCarloWeightsRun2.C`.
The macro uses ROOT and the existing `JDrawer.h`; it does not require
`libEECDrawing.so`, a separate drawing package, or Python.
The synthetic ROOT test is `testDeriveMonteCarloWeightsRun3.C`.
Running the macro never modifies the analyzer, existing event weights, or input files.

## Run from plotting/

```bash
root -l -b -q 'deriveMonteCarloWeightsRun3.C("OO_DATA.root","OO_MC.root","PP_DATA.root","PP_MC.root")'
```

Arguments, in order:

```cpp
deriveMonteCarloWeightsRun3(
  const char* ooData = "",
  const char* ooMC = "",
  const char* ppData = "",
  const char* ppMC = "",
  const char* outputDirectory = "figures/mcWeightsRun3",
  const char* mode = "both",             // stored, derive, both
  double maxCentrality = 90.0,           // percent, NOT hiBin
  int vertexPlotRebin = 1,               // display only
  int centralityPlotRebin = 1);           // display only
```

Supply analyzer output or processed ROOT files containing event-information
histograms, NOT `HiForestMiniAOD_*.root` forests. Data and MC must describe the
intended compatible event selections. The macro prints the input card's trigger,
vertex, multiplicity, and producer-hash information into its audit file. It does
not establish trigger equivalence, luminosity normalization, or cross-section
normalization from these histograms.

The following filenames are recorded in the existing repository weight studies.
This command assumes they are in the repository root; their presence on lxplus
has not been checked here. Replace them with your new merged files as appropriate.

```bash
root -l -b -q 'deriveMonteCarloWeightsRun3.C("../merged_eecOO2025_processed_121026.root","../eecAnalysis_OOMC_eventWeights_merged_processed.root","../merged_eecpp2024_run3_v7_processed.root","../merged_eecppMC2024_GenGen_processed.root","figures/mcWeightsRun3","both")'
```

OO only, audit the weights actually stored in the production file:

```bash
root -l -b -q 'deriveMonteCarloWeightsRun3.C("OO_DATA.root","OO_MC.root","","","figures/mcWeightsRun3_stored","stored")'
```

pp only, derive a candidate vertex weight:

```bash
root -l -b -q 'deriveMonteCarloWeightsRun3.C("","","PP_DATA.root","PP_MC.root","figures/mcWeightsRun3_pp","derive")'
```

Use `maxCentrality=80` to show hiBin 0-160 as in the example. The default is 90%
(hiBin 0-180), matching the currently implemented OO weight domain. A fit derived
only to 80% must NOT be used in the 80-90% range. The generated parameter comments
state the required hiBin guard. Vertex plots remain inclusive in the selected
event activity; cropping a 1D centrality histogram cannot impose that cut on vz.

## What is already active in the inspected framework

Audited repository: `JesusVlz/energyEnergyCorrelatorAnalysis`.
Pinned revision: `6d4800b8c57747e1db761fc4fecf72a38a1a00c3`.

`src/EECAnalyzer.cxx` contains:

- OO MC vertex `pol6` in the `fIsOOData` constructor branch, starting with
  `0.838856911780552, -0.00776989024519781, 0.0058636391429985266, ...`.
- OO MC centrality: two `pol6` functions with a common value at 30%, defined for
  0-30% and 30-90%. `GetCentralityWeight` evaluates them at `hiBin/2.0` and returns
  unity outside `0 <= hiBin < 180` for OO MC.
- ppRef MC vertex `pol6`, with coefficients
  `0.901743, -0.0100148, 0.00349298, -5.14891e-05, 6.59866e-06, 7.97589e-08, -5.10072e-08`.
  ppRef has no centrality weight.
- `GetVzWeight` returns unity for real data and evaluates the function for MC.
- In the regular centrality-based Run-3 path, the event loop forms
  `fTotalEventWeight = fVzWeight * fCentralityWeight * fPtHatWeight`.
  For OO MC and ppRef MC, `fPtHatWeight` is explicitly set to **1.0**. Generator/
  cross-section weighting is therefore a separate outstanding question, not
  something this shape-reweighting macro derives or enables.

Existing studies are recorded in `validation/ooEventWeights/summary.json` and
`validation/ppEventWeights/summary.json`. The historical input files described
there had identical raw and weighted histogram contents. Active source code
alone does not establish that an old ROOT output contains those weights.
The macro checks the supplied files rather than assuming they were regenerated.

## Exact input histogram semantics

| Purpose | Vertex z key | Centrality key |
|---|---|---|
| Data, including any trigger-combination event weights | `vertexZweighted` | `centralityWeighted` |
| MC baseline, without vertex/centrality factors | `vertexZ` | `centrality` |
| MC as actually weighted in that production | `vertexZweighted` | `centralityWeighted` |

In the inspected analyzer, baseline MC histograms use `fPtHatWeight`; the weighted
ones use `fTotalEventWeight`. Thus the stored OO vertex histogram also includes
the centrality factor, and the stored OO centrality histogram also includes the
vertex factor. A `stored` plot is a genuine readout of that production's full
weight product; it is NOT a new application of the weights in the plotting macro.

The macro warns if raw and stored-weighted contents are identical. An identical
shape by itself cannot identify exactly which factors were applied; inspect the
producer hash and regenerate with the intended analyzer when necessary.

## Derivation and important distinctions

The target is the unit-area ratio `data / baseline MC` in the stated domain.
Candidate weights are full REPLACEMENTS for the corresponding old event-shape
functions, not residual multipliers on top of them.

Vertex fits use degree six and the native histogram bin centers. OO centrality
fits use two degree-six polynomials fitted jointly with a shared value at 30%.
The coordinate is the discrete `c = hiBin/2.0`, at the lower edge of the native
0.5%-wide centrality bins, not the bin center (which would shift it by 0.25%).
The plotting x coordinate alone is converted to hiBin.

Weighted linear least squares is solved with ROOT's `TDecompSVD` in scaled
coordinates, then converted to the ordinary `pol6` coefficients used by the
analyzer. This changes the numerical basis, not the polynomial degree or physical
coordinate. The new centrality fit imposes continuity during the joint fit;
it does not merely adjust the peripheral intercept afterward.

Fits use the original fine bins. The optional rebin factors affect drawings
only and must divide the number of selected bins. Missing MC support under
populated data, wrong systems, incompatible native binning, rank-deficient fits,
and non-positive polynomial excursions fail explicitly rather than being hidden
by clipping or substituted weights.

`*_refit` green curves are deterministic, single-marginal, same-sample previews:
MC bin contents AND their errors are multiplied by the fitted weight and the
result is renormalized. They are not independent closure and do not validate
simultaneously applying both OO factors when vz and centrality are correlated.
Use a new analyzer pass (or suitable joint event-level information) to check the
combined product. Do not divide an old weighted 1D projection by a marginal
function to attempt to recover the joint event distribution.

Fit chi2/covariance and displayed ratio errors use diagonal input errors.
Normalization-induced correlations and fitted-weight uncertainty are not included
in the preview error bars. These are QA/fit diagnostics, not final uncertainties.

## Figures and outputs

The style uses the repository's own `JDrawer::SetDefaultAppearanceSplitCanvas`,
logarithmic top panels, linear MC/Data panels, a dashed unity line, red data
squares, blue baseline-MC circles, and green after-weighting triangles. Headers
are CMS Preliminary and 5.36 TeV. `P+H` is retained as the OO generator shorthand
from the example, without assuming an expansion; change `kOOLabel` if desired.

The default `both` mode creates:

```text
figures/mcWeightsRun3/
  OO_eventWeights_stored.pdf/.png      # hiBin left, vz right
  OO_eventWeights_refit.pdf/.png       # same layout, candidate 1D checks
  OO_hiBin_stored.pdf/.png
  OO_hiBin_refit.pdf/.png
  OO_vz_stored.pdf/.png
  OO_vz_refit.pdf/.png
  pp_vz_stored.pdf/.png
  pp_vz_refit.pdf/.png
  run3_mc_weights.root
  weight_audit.txt
  candidate_parameters.txt
```

The fixed displayed ratio range is 0.5-1.5, matching the example. The audit reports
any off-scale points and full ratio graphs are saved in ROOT. Adjust `kRatioMin`
and `kRatioMax` to inspect them visually.

ROOT subdirectories are `OO/vz`, `OO/centrality`, and `pp/vz`. They contain
normalized input histograms, full ratios, canvases, candidate predictions, scaled
fit parameters/covariance and candidate TF1 objects. Relevant functions are:

```text
OO/vz/weight_vz
OO/centrality/weight_centrality_0_30
OO/centrality/weight_centrality_30_max
pp/vz/weight_vz
```

The root-level `TParameter<int>("complete")` must be 1. It remains 0 after a
failure; do not use a partially written candidate file. The text file provides
`SetParameters(...)` lines for the appropriate existing constructor branches.
Nothing is deployed automatically. Reusing an output directory replaces the
macro's previous output files; use a new directory to retain comparisons.

## Tests and execution status

The portable polynomial conversion/continuity routines were compiled with
`clang++ -std=c++17 -Wall -Wextra -Werror` and tested locally. The test passed
200 randomized parameter sets, both 0-80% and 0-90% centrality domains, vertex
round trips, and the invalid-domain guard. Maximum conversion round-trip
absolute error was about `1.33e-15` in those tests.

ROOT/PyROOT and the actual data/MC ROOT inputs were not available in the authoring
environment. For the repository integration, the ROOT smoke test is run by CI; consult its
exact run result. Actual data/MC fits have not been performed by this source update. A ROOT smoke test is included for lxplus:

```bash
# From plotting/, after loading your usual CMSSW/ROOT environment:
root -l -b -q 'testDeriveMonteCarloWeightsRun3.C()'

# Optional ACLiC compilation of the actual entry point (prints usage):
root -l -b -q 'deriveMonteCarloWeightsRun3.C+()'
```

The smoke test creates synthetic inputs in a temporary directory and checks
polynomial recovery, continuity, ROOT outputs, figure creation, and rejection of
a data-populated bin with no MC support. It is a software test, not physics
closure. The standalone download also carried portable-math tests; the repository uses
the ROOT smoke test above in its integration checks.

## Source locations inspected

At the pinned repository revision:

- `plotting/deriveMonteCarloWeights.C`: existing Run-2-only derivation, polynomial
  fits and weighted-data/unweighted-MC key choices.
- `plotting/JDrawer.h`: group style, split-pad interface and drawing defaults.
- `src/EECAnalyzer.cxx`: constructor weight parameters; `GetVzWeight`,
  `GetCentralityWeight`; Run-3 event-weight product; event-information fills.
- `src/EECHistograms.cxx`: 80 vertex bins on [-20,20] and 202 centrality bins on
  [-1,100], i.e. 0.5 cm and 0.5% native bin widths.
- `validation/ooEventWeights/summary.json` and
  `validation/ppEventWeights/summary.json`: historical input provenance and limits
  of single-marginal checks.

ROOT APIs checked against the official `TDecompSVD` and `TObject` references.
