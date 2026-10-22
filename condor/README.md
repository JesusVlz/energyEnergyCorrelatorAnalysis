# O–O MC on CERN HTCondor

Prepare and submit from AFS. The analysis source and input forests remain on EOS;
workers use a frozen payload in local scratch and copy validated ROOT outputs to
EOS with XRootD.

From lxplus, prepare the submission with:

```bash
python3 /eos/user/j/jvelazqu/HIN/CMSSW_15_0_11/src/energyEnergyCorrelatorAnalysis/condor/prepareOOMC.py
```

Defaults:

- Submission files and logs: `/afs/cern.ch/user/j/jvelazqu/condor/OOMC`.
- Inputs: `/eos/user/j/jvelazqu/Foresting/OOMC/unembedded_L1HLT_objects_MC_*.root`.
- Results: `/eos/user/j/jvelazqu/HIN/eec-outputs/OOMC/<run-tag>/<mode>/`.
- Runtime: the CVMFS `CMSSW_15_0_11` release, `el9_amd64_gcc12`.

The preparation command discovers the available input files, builds the current
source once in a temporary directory, and records the inputs, configuration
cards, source hashes, and payload checksum in the AFS workspace. It includes
uncommitted source changes. It does not build in the shared checkout or submit
jobs. Existing prepared workspaces are protected against replacement; pass a new
`--submit-dir` and `--tag` when preparing another production.

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

Response mask `479` includes the regular distributions, jet closure, jet/EEC
response, and track–particle matching. These are response inputs for the existing
unfolding workflow; they do not perform the final unfolding.

Review `cards/`, `inputs.txt`, `jobs.tsv`, and `manifest.json`, then validate
the submission description and submit from AFS:

```bash
cd /afs/cern.ch/user/j/jvelazqu/condor/OOMC
condor_submit -dry-run dryrun.ads jobs.sub
./submit.sh
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

Preparing new jobs after changing code or cards requires a new payload. Editing
the checkout after submission does not change queued jobs.

References: [CERN job submission](https://batchdocs.web.cern.ch/local/submit.html)
and [HTCondor file transfer](https://htcondor.readthedocs.io/en/24.0/users-manual/file-transfer.html).

## pp MC

Use the same five configurations for the 5.36 TeV pp MC forests with:

```bash
python3 condor/prepareOOMC.py --sample ppMC
cd /afs/cern.ch/user/j/jvelazqu/condor/ppMC
./submit.sh
```

This uses `cardEECpp5v36TeVMC.input`, `DataType 12`, and inputs matching
`/eos/user/j/jvelazqu/Foresting/ppMC/HiForestMiniAOD_MC_*.root`. Results go to
`/eos/user/j/jvelazqu/HIN/eec-outputs/ppMC/<run-tag>/<mode>/`. The pp card
retains its own centrality binning and tracking corrections. The same payload,
retry, and output-copy workflow applies.
