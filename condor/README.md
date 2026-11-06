# O–O and pp MC on CERN HTCondor

Prepare and submit from AFS. The analysis source and input forests remain on EOS;
workers use a frozen payload in local scratch and copy validated ROOT outputs to
EOS with XRootD.

The common AFS site is `/afs/cern.ch/user/j/jvelazqu/condor`, with production
campaigns under `production/<sample>/<run-tag>/` and prepared merge campaigns
under `merge/<sample>/<run-tag>/`. Each campaign keeps its own submission files,
manifests, and logs; production campaigns also freeze their cards and payload.

From lxplus, prepare a production with the shared entry point:

```bash
/afs/cern.ch/user/j/jvelazqu/condor/production/prepare.sh --sample OOMC --tag RUN_TAG
```

The wrapper runs the repository's `condor/prepareOOMC.py`, which can also be
called directly. Omit `--tag` to generate a timestamped tag.

Defaults:

- Submission files and logs: `/afs/cern.ch/user/j/jvelazqu/condor/production/OOMC/<run-tag>/`.
- Inputs: `/eos/user/j/jvelazqu/Foresting/OOMC/unembedded_L1HLT_objects_MC_*.root`.
- Results: `/eos/user/j/jvelazqu/HIN/eec-outputs/OOMC/<run-tag>/<mode>/`.
- Runtime: the CVMFS `CMSSW_15_0_11` release, `el9_amd64_gcc12`.

The preparation command discovers the available input files, builds the current
source once in a temporary directory, and records the inputs, configuration
cards, source hashes, and payload checksum in the AFS workspace. It includes
uncommitted source changes. It does not build in the shared checkout or submit
jobs. Existing prepared workspaces are protected against replacement; use a new
`--tag` when preparing another production. An explicit `--submit-dir` remains
available for a custom AFS workspace.

The default queue contains one job per input file and configuration:

| Mode | McCorrelationType | MatchJets | FilledHistograms |
|---|---:|---:|---:|
| RecoReco | 0 | 0 | 31 |
| RecoGen | 1 | 0 | 31 |
| GenReco | 2 | 0 | 31 |
| GenGen | 3 | 0 | 31 |
| Response | 3 | 1 | 479 |

All cards use `DataType 10`, `TriggerSelection 0`, AK4 WTA jets, nominal jet
corrections, both energy-weight exponents from the base MC card, and disabled
track-pair efficiency corrections. Offline event, jet, and track cuts remain
active. The four regular modes retain the MC baseline `BackgroundMethods 0`.
Use `--background-methods 3` when preparing a separate production for reflected
and perpendicular cone closure; the response card always uses zero.
Add `--match-eec-binning` to use the EEC jet and track pT edges for the
corresponding sparse axes in every mode, including Response. This also aligns
`TrackPtBinEdges` with `TrackPtBinEdgesEEC` in the saved cards. Existing EEC edges
are retained; the option is recorded as `MatchEECBinning 1` and checked in worker
outputs before copying them to EOS.

Response mask `479` includes the regular distributions, jet closure, jet/EEC
response, and track–particle matching. These are response inputs for the existing
unfolding workflow; they do not perform the final unfolding.

Review `cards/`, `inputs.txt`, `jobs.tsv`, and `manifest.json`, then validate
the submission description and submit from AFS:

```bash
cd /afs/cern.ch/user/j/jvelazqu/condor
./production/submit.sh OOMC RUN_TAG --dry-run
./production/submit.sh OOMC RUN_TAG
condor_q
```

The ROOT results are copied to EOS by each worker. Condor transfers only the
small payload, executable wrapper, and job logs through AFS. Workers check that
events pass selection, the saved card matches the requested mode, and the output
histograms are readable before copying the result. Sparse response histograms
can have zero entries for a file without qualifying jets; their counts are
reported in the job log.

The submission requests one CPU, 4 GB RAM for regular modes, 8 GB RAM for the
response, 6 GB scratch disk, and the `workday` flavour. Failed jobs retry twice.
The worker uses CERN's supplied batch credentials rather than copying a login
shell's credential environment.

Follow jobs with `condor_q`; inspect `logs/` for the worker output, errors, and
one event log per cluster. ROOT files are analysis histogram outputs; merge and
project each mode separately using the repository's normal postprocessing
workflow. Keep the response and four correlation modes separate.

Submit an already prepared merge campaign through the common merge entry point:

```bash
cd /afs/cern.ch/user/j/jvelazqu/condor
./merge/submit.sh OOMC RUN_TAG --dry-run
./merge/submit.sh OOMC RUN_TAG
```

Both submission wrappers accept `SAMPLE TAG [--dry-run] [--file jobs.sub]`.
They select the campaign under their respective site; preparation and review
must be completed before submission. Existing campaigns retain their original
manifests and logs; `relocation.json` records their previous and current paths.

Preparing new jobs after changing code or cards requires a new payload. Editing
the checkout after submission does not change queued jobs.

References: [CERN job submission](https://batchdocs.web.cern.ch/local/submit.html)
and [HTCondor file transfer](https://htcondor.readthedocs.io/en/24.0/users-manual/file-transfer.html).

## pp MC

Use the same five configurations for the 5.36 TeV pp MC forests with:

```bash
/afs/cern.ch/user/j/jvelazqu/condor/production/prepare.sh --sample ppMC --tag RUN_TAG
cd /afs/cern.ch/user/j/jvelazqu/condor
./production/submit.sh ppMC RUN_TAG --dry-run
./production/submit.sh ppMC RUN_TAG
```

This uses `cardEECpp5v36TeVMC.input`, `DataType 12`, and inputs matching
`/eos/user/j/jvelazqu/Foresting/ppMC/HiForestMiniAOD_MC_*.root`. Results go to
`/eos/user/j/jvelazqu/HIN/eec-outputs/ppMC/<run-tag>/<mode>/`. The pp card
retains its own centrality binning and tracking corrections. The same payload,
retry, and output-copy workflow applies.


## Explicit OO GenGen mixed/perpendicular campaigns

Use a fresh campaign for each requested centrality selection. `--base-card`
selects a dedicated card (relative paths resolve against `--repo`); its original
bytes, path and checksum are frozen alongside the generated mode card. The
normal five-mode defaults remain available for existing workflows.

The OO GenGen mixed/perpendicular configuration is explicitly selected with
`--modes GenGen --background-methods 6`. It requires the fixed matching contract,
strict particle-pair pT above 2 GeV, seed 20261030, and
`OOClosurePartitionPolicy 1`. Every job reserves the identical first 100,000
audited pure-HIJING donor events and rejects signal events sharing a reserved
physical ID or observed nonempty background fingerprint. All source signal
entries are read; the schema audit verifies that the card cap covers the
largest complete input. This does not establish unseen generator lineage or
freeze the background prescription.

Required preparation inputs are `--mixing-pool-list`,
`--mixing-pool-manifest`, `--signal-input-manifest` and `--pool-events 100000`.
The source list remains unchanged. The frozen worker list contains only the
four donor files covering the audited 100,000-event prefix; their source
counts, used entry ranges, UUIDs and sizes accompany both provenance audits.
The ROOT files are downloaded with XRootD into worker scratch. Workers rewrite
only the local mixing list and closure-output prefix, preserving the frozen
physics settings.

These jobs request 8 GB RAM, 16 GB scratch and the `workday` flavour. Their
isolated payload compiles every analysis object with `g++ -O2`, retaining ROOT
compiler flags and avoiding fast-math. `build.log`, the campaign manifest and
`campaign_provenance.json` record the exact compiler/make command and source
hashes. Other preparations retain their existing compiler configuration.

Alongside each `<stem>.root`, the worker publishes `<stem>.partners.csv`,
`<stem>.pool.csv`, `<stem>.jets.csv`, `<stem>.contributions.csv` and
`<stem>.run.json`. Keep these identity and contribution ledgers with the ROOT
outputs: the methods share signal events and reuse reserved donors across
jobs, so those repetitions are correlated exposure. An execution/validation
pass is separate from physics closure and final covariance validation.

Enabled jet-veto ROOT maps and their pinned metadata sidecars are included in
all new payloads that reference them. Missing maps/sidecars or mismatched
checksums fail before a campaign directory is created.

## Merge complete OO GenGen mixed/perpendicular productions

Prepare each centrality campaign independently after all production jobs finish:

```bash
python3 condor/prepareOOMerge.py \
  --production /afs/cern.ch/user/j/jvelazqu/condor/production/OOMC/GenGen_R04_Run3_NH_20261029_split
python3 condor/prepareOOMerge.py \
  --production /afs/cern.ch/user/j/jvelazqu/condor/production/OOMC/GenGen_R04_Run3_NH_20261029_integrated
```

This workflow supports the explicit OO GenGen reserved-donor closure contract
(BackgroundMethods 6). Preparation requires all manifest inputs, successful
per-input run records, frozen campaign/source/pool provenance, and every recorded
artifact at its expected size. It rejects duplicate, extra or missing ROOT
outputs. `--audit-only` performs these checks without creating a campaign.

The new campaigns live under the shared AFS `merge/OOMC/<production-tag>/`
site. Review their `manifest.json`, then use the existing shared submission
wrapper with `--dry-run` followed by submission. Each campaign queues one
8 GB, one-CPU merge worker; scratch is calculated from the frozen artifact
sizes. Production campaigns and their outputs are retained.

Workers stage every ROOT file and all four CSV ledgers with XRootD, verify
the SHA256 of each artifact and its original run record, and fail if any
input is inconsistent. They merge every expected ROOT input without skipping
bad files, preserve the original physics card, and verify compatible axes,
configuration and histogram totals. The complete original per-job ledgers
and run records are archived without concatenating repeated donor pools.

Results are published into `<production-output-dir>/merged/`:

- `merged_OOMC_GenGen.root`: merged histograms.
- `merged_OOMC_GenGen.ledgers.tar.gz`: original CSVs/run records and frozen
  merge manifest.
- `merged_OOMC_GenGen.run.json`: final success record with artifact sizes,
  SHA256 values and merge validation.

Consume the ROOT result only with the final successful merge record and verify
its checksums. Existing merge destinations/campaigns are protected; partial
publication requires review before retrying. Keep split/integrated results
separate because they repeat physical signals. A successful merge does not
establish physics closure or donor/cross-job covariance.


## Merge an available OO/pp CRAB data snapshot

Data can be merged while CRAB processing continues. Freeze the currently stable
outputs first; later files require a new snapshot. To select OO parts 01–05,
including both chunks, and all submitted ppRef tasks:

```bash
python3 condor/inventoryDataOutputs.py --max-oo-part 5 \
  --output /tmp/eec-data-inventory.json
python3 condor/prepareDataMerge.py \
  --inventory /tmp/eec-data-inventory.json \
  --submit-dir /afs/cern.ch/user/j/jvelazqu/condor/merge/Run3Data/NEW_TAG \
  --output-dir /eos/cms/store/group/phys_heavyions/jvelazqu/eec/merged/NEW_TAG \
  --tag NEW_TAG
cd /afs/cern.ch/user/j/jvelazqu/condor/merge/Run3Data/NEW_TAG
./submit.sh --dry-run
./submit.sh
```

Use unused inventory and campaign paths. Omit `--max-oo-part` for all authoritative
OO parts. OO part02/chunk02 uses its v3 replacement; the killed unsplit v1 and
refused v2 submission are excluded. Files must be nonzero, at least two minutes
old, and unchanged across two metadata checks. The inventory records the selected
scope, full campaign input-file count, missing jobs, and withheld files. A part
cutoff represents a share of planned inputs, not necessarily the same share of
completed events.

The campaign-local submission script handles the data dependency graph; the
shared `merge/submit.sh` remains the MC entry point. Batches feed one result per
CRAB part/chunk under `OO/parts/` or `ppRef/parts/`, then separate `merged_OO.root`
and `merged_ppRef.root` totals. Condor runs at most 20 nodes at once. Each ROOT
product has a JSON report published last, with its checksum, original-source
membership, histogram validation and any corrupt-input exclusions. Read that
report before using a result. `PASS` certifies the frozen inputs, not completion
of the full CRAB campaign.

Workers preserve the analysis configuration, validate every histogram bin and
stored error, and retain exact event counts using a double-precision `nEvents`
histogram. Parent reports and checksums prevent overlapping reductions. Reduction
scratch sizes are calculated from reports returned to AFS, without requiring an
EOS mount on the scheduler. Existing outputs are never overwritten; inspect any
failed publication before submitting a rescue DAG.


## Prepare pair and response production from explicit inputs

Use a frozen full-MC inventory and a verified native runtime snapshot:

```bash
python3 condor/preparePairResponse.py \
  --tag NEW_PAIR_RESPONSE_TAG \
  --state-dir /path/to/new-preparation-state \
  --build-state /path/to/build_state.json \
  --inventory /path/to/sample_manifest.json \
  --condor-root /path/to/condor
```

The build-state JSON must contain `status: "BUILT_BATCH_RUNTIME"`, `snapshot`,
`binary_sha256`, `cmssw_base` and `scram_arch`. The snapshot must contain the
compiled `eecAnalysis` and its runtime payloads; the binary checksum is checked
before preparation. Generate the full-file inventory with
`validation/run3Pairs/all_mc_campaign.py prepare`. Historical date-specific
defaults refer to local evidence and are not included in a fresh checkout.
`--reference-oomc` and `--reference-ppmc` optionally record existing campaign
provenance; they are not prerequisites. Preparation writes queues and receipts
without submitting them. Output storage still follows the CERN layout described
above.

The alternative `validation/run3Campaign/prepare_campaign.py` DAG accepts the
same explicit build-state, state-directory and inventory options. Its worker is
pinned to CMSSW_15_0_11 with `el9_amd64_gcc12`; incompatible build states are
rejected before preparation. The plain pair/response queue uses the runtime
specified in its validated build state.

Data and pair/response workers check frozen source size and modification time
before and after staging. XRootD exposes modification times at second precision;
these checks are metadata validation, not a content checksum of the original
inventory. Recorded output checksums attest the staged and produced bytes.
Supplying an independently staged pair input requires its frozen SHA256.
