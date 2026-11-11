# Run-3 R4 track/pair preparation

Status: framework preparation. No new MC input is assumed, no production is
submitted, and no new detector correction or physics closure is claimed here.
Use AN2023_081_v13 Secs. 6.1-6.2 and 7.5.1; the correction procedure is unchanged.
The software check results are in the repair branch's GitHub Actions run.

## Run-3 selection and payloads

Both OO2025 and ppRef2024 require highPurity. For all three selections the
momentum-error cut is sigma(pT)/pT <= 0.1 only for pT > 10 GeV. At or below
10 GeV require a finite, nonnegative uncertainty, without that upper cut.
The inclusive dz/dzError and dxy/dxyError significance limits are 3, 5 and 2
for Nominal, Loose and Tight respectively. Undefined significances fail.
PbPb2018/pp2017 selection and correction code remains on its legacy branch.

Set both keys consistently, for example:

```
Run3TrackSelection Tight
MaxTrackPtRelativeError 0.1
VertexMaxDistance 2
```

The analyzer selects the corresponding existing Nominal/Loose/Tight table for
the collision system. It rejects conflicting old 0.05/0.15 Run-2 combinations.
The formula is (1-fake)*(1-secondary)/(efficiency*(1+multiple)); no component is
applied twice. OO uses eta,pT,hiBin; pp uses eta,pT. Physical detector eta and
the source event's activity are used before translating cone coordinates.
The existing OO 6.1 GeV lookup cap is retained, explicitly recorded, and is not
newly justified by this repair. Binary table contents are not modified.

The table audit checks all six files, dimensions, common axes, checksums and
invalid/empty bins. Empty bins are reported, never filled with invented values.
An unavailable correction in a selected EEC pair fails rather than contributing
zero silently. Unsupported full-event auxiliary counts keep their existing
explicit masking treatment.

Physical reconstructed charge is now read from AOD or MiniAOD trkCharge and
checked across chained files. Missing charge is an input error, never +1 by
default for Run-3. This repairs the exported charge-matching diagnostic.

## Ready configuration cards

Twelve cards are under `validation/run3Pairs/cards/`: OO and pp, each with
Nominal/Loose/Tight and native RecoReco/RecoGen configurations. They use common
reconstructed R4 WTA jets, |etaJet|<1.6, 40-80/80-120/120-160 GeV, both energy
powers, 0.5 GeV input retention and a strict 2 GeV EEC-pair threshold.
OO has 0-30/30-50/50-90/90-100 activity input bins; the existing pair code
provides the integrated selected-MB and three reported split selections.
HF shifting, background estimation and residual pair application are OFF during
derivation. These are pair-preparation cards, not changes to unfolding guards.

Regenerate or inspect the cards without processing MC:

```bash
python3 validation/run3Pairs/prepare_inputs.py cards --output /tmp/run3-pair-cards
python3 validation/run3Pairs/prepare_inputs.py audit-tables --output /tmp/tracking_tables.json
```

## When the MC inputs are supplied

Provide separate OO and pp text lists, one local or EOS-mounted ROOT pathname
per line. Relative paths resolve beside the list. Stage root:// inputs first;
this preparation command does not fetch remote files or invent paths. Provide
one stable physical campaign namespace per system, identical across its files
and nominal/loose/tight runs. Include the sample production information needed
to identify upstream generator-particle retention and any reused embedding
sources. The preparation preserves the analyzer's existing event-weight policy;
it does not establish a new cross-section normalization.

This creates a bounded pilot plan but does NOT execute it:

```bash
python3 validation/run3Pairs/prepare_inputs.py plan \
  --oo-list OO_MC.txt --oo-namespace OO_CAMPAIGN \
  --pp-list pp_MC.txt --pp-namespace PP_CAMPAIGN \
  --selections Nominal --max-events 10000 \
  --output /path/to/new/campaign
```

The command checks input trees, physical IDs and charge, freezes file identities,
and writes cards, task manifests and `run_commands.sh`. A single-system pilot
is supported. Use a fresh directory for each plan. Nominal should run first;
`--selections Nominal Loose Tight` prepares all existing AN selection variations.
`--all-events` replaces the pilot entry limit with the exact input-file counts.
Do not call a bounded pilot the full campaign.

When execution is authorized:

```bash
bash /path/to/new/campaign/run_commands.sh
```

The script runs existing native exports, compares native/exported EEC sums on
the first file, retains inclusive track QA and representative native responses,
then uses the existing mergeable batch machinery. Each shard saves unnormalized
pair counts and a holdout cache. Global merging precedes the fixed
0.1<DeltaR<0.4 normalization; never average per-file correction maps. Exact raw
identity ledgers prevent duplicate physical events and training/holdout leakage.
The reduced pipeline avoids retaining the full all-MC event-vector result in
one Python process. Native JSONL/ROOT outputs are retained; allocate scratch
space accordingly. No automatic Condor submission or deletion is performed.

## Correction and closure products

The existing AN machinery writes jet-dependent and full-event jet-independent
maps, counts/control/support masks, and truth/raw/single/pair closure plots for
all 30 R4 channels. First show 40-80 GeV, then the same plots for the full grid.
Detector factors enter once for n=1 and n=2. A supported jet-independent cell is
the existing sparse fallback; unsupported cells remain unavailable.
The small-angle reporting boundary and any nonclosure remain MC conclusions.

New exports carry `run3-tracking-2026-11-v1`; new pair payloads carry
`run3-an751-v2`, system/selection tags and the single-track payload identity.
Old pre-repair exports/results are preserved but must be regenerated before
use with the new selection. The normal analyzer rejects missing or conflicting
system/selection tags. Historical plots can still be drawn without being
accepted as new correction payloads.

The same-source background self-pair/cross-cone treatment already present in
the analyzer is retained. Pair correction does not arise from an artificially
small translated cross-event separation. Coordinate with Background mixing
through its existing corrected closure, not a new estimator.

## Activation is a later, explicit step

Default cards keep `DisableTrackPairEfficiencyCorrection 1`. After deriving
supported maps and completing the AN closure, use the actual new ROOT payload
with `Run3TrackPairCorrectionFile` and switch to zero on a compatible card.
The initial maps cover 40 <= reco jet pT < 160 GeV and tracks above 2 GeV;
do not enable them on a wider filled-jet or lower-track domain without matching
support. Retain appropriate unfolding guard ranges in their own task/cards.

The remaining MC-dependent work is the map itself, standard single-track and
pair closure, representative angular/momentum-weight assessment, and the AN
tracking/pair uncertainty evaluation. Empty-cell support or zero stored map
errors are NOT a measured zero uncertainty. Run-2 uncertainty percentages are
not promoted to Run-3 inputs by this preparation.
