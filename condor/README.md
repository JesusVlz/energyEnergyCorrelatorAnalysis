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
