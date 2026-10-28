# Run-3 residual pair diagnostic tools

These standalone tools implement a candidate derivation and independent weighted
holdout checks. **They do not create a production payload or validate physics.**
No existing analysis correction, coordinate convention, calibration or estimator
is changed. Native R4 and R8 are handled in parallel when available; absent native
R8 trees are recorded as blockers and never replaced by wider AK4 cones.

`export_forest.py` reads a bounded real MC forest using PyROOT and writes JSONL
physical events with native reconstructed WTA jets, their original reconstructed
tracks and generator particles in the **same reconstructed jet cone**. Each row
retains campaign namespace plus run/lumi/event, track indices, original detector
eta/phi, generator weight and provisional stored hiBin. It mirrors the current
Run-3 nominal track cuts (strict measurement pT>2 GeV, |eta|<2.4). Original
low-pT truth completeness must be checked separately, for example with the
preserved input audit. The input minimum of 0.5 GeV and the measurement minimum
of 2 GeV are different requirements.

The exporter explicitly reproduces the inherited EEC single-track lookup using
track eta relative to the WTA jet axis. It retains the physical eta as separate
information; this diagnostic does not resolve the coordinate convention or
single-track calibration applicability. OO lookups retain hiBin and the existing
6.1 GeV high-pT cap. Invalid lookups fail instead of silently becoming unity.
The exported event extra weight is one; the stored generator weight is used once.
No old 0–90% MC matching function is extrapolated. Stored forest `jtpt`, rather
than a new applicable JEC evaluation, is used for this diagnostic; full event
filter equivalence to the analyzer and payload applicability remain explicit
physics blockers. Unsupported OO hiBin events are counted and excluded, so the
output cannot be labeled a 0–100% measurement.

`pair_counts.py` uses differential higher/lower pT, reconstructed jet pT,
provisional hiBin activity and pair separation bins. Training has no EEC energy
factor: generator/event weights remain once and reconstructed counts include the
two nominal single-track factors once. For each differential cell, the candidate
is `C = alpha * Ngen / Nreco_single`, with `alpha` obtained from the large-angle
control within that same cell. The inherited 0.1–0.4 interval is a starting
control, **not a validated control or a final reporting range**. Use `--control
0.15 0.3` to study an alternative aligned to stored bin edges. Both angular
supports extend through 1.6; no hard-coded 0.4 unity rule is introduced. Sparse,
nonpositive or missing cells have null corrections and explicit masked yields;
there is no smoothing, edge clamping or fallback.

Physical events are split by the same namespace/run/lumi/event FNV1a64 key as
`src/OOEventIdentity.h`; folds 0/1 are training/holdout. File, entry, radius and
power never enter this split. Duplicate physical events are rejected. On the
holdout, `(pT1*pT2)^n` is evaluated independently for n=1 and n=2. Single-track,
pair, event and generator factors are each used once; the map is never squared.
Raw, single-track, supported pair-corrected and truth yields are saved separately
from diagnostic full-support unit-area shapes. Any missing map contribution
makes the corrected bin null instead of presenting a partially corrected result.

Common holdout event contribution vectors retain all radii and powers for later
joint replicas. They are sufficient statistics for holdout counting covariance;
finite-MC uncertainty in the derived map additionally requires rederivation in
shared training-event replicas. The supplied sums of squared pair weights do not
claim independent-Poisson pair errors. No covariance or uncertainty is invented.

Run from the repository root; choose an explicit namespace that is identical
across radius streams from the same campaign and distinct for independent MC. The example namespace below identifies the audited
single-file local test only; no campaign provenance has been established for
this local copy, and it must not be assumed to define a full production sample:

```bash
python3 -m unittest discover -s validation/run3Pairs -p 'test_*.py' -v
python3 validation/run3Pairs/export_forest.py \
  /eos/user/j/jvelazqu/Foresting/OOMC/unembedded_L1HLT_objects_MC_1.root \
  /tmp/OO_pairs.jsonl --system OO \
  --sample-namespace OO_local_unembedded_L1HLT_MC_1 --max-events 1000 \
  --payload trackCorrectionTables/OO2025/Eff_OO_2025_Hijing_MB_Centrality_fromHihfpf_NoPU_3D_Nominal_Official_18Nov2025.root
python3 validation/run3Pairs/pair_counts.py /tmp/OO_pairs.jsonl /tmp/OO_pair_candidate.json --min-pairs 10
```

The exporter saves SHA256 hashes of input, single-track payload and exported
JSONL in an adjacent `.metadata.json`; derivation saves its input/tool hashes.
Event ranges and counts are explicit. The bound is 1–10,000 events and is not a
production launcher. Actual read-only input audit and bounded pilot evidence for
this task are under `/tmp/jvelazqu/oo-eec-audit/`; these temporary artifacts are not
repository data or validation certificates.

Production integration remains gated on applicable era/reconstruction/cuts and
jet calibration, OO expert guidance, sufficient independent statistics across
all twelve channels in each system, normalization/binning/mixture sensitivity,
R8 range selection, detector topology/background interaction and joint finite-MC
covariance. Candidate output is always marked `BLOCKED_NOT_A_PRODUCTION_PAYLOAD`.
