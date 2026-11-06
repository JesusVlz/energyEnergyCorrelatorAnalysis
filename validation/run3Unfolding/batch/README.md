# Full-MC response matrices

The worker runs the unchanged native response on synchronized train and validation forests. Every physical event is assigned by the existing `pair_counts.fold` function: namespaced FNV1a64 over run/lumi/event, modulo two; zero trains and one validates. The split is independent of filename, input entry, radius and power. Latest ROOT key cycles are used, matching the native readers. All event TTrees must match the event-header entry count; only HiForestInfo metadata trees are copied unchanged. The compressed ledger retains each original input entry and its fold entry.

Run from the packaged analysis checkout under the established CMSSW runtime:

```sh
python3 validation/run3Unfolding/batch/response_worker.py INPUT.root OUTPUT_DIRECTORY --system OO --source-input ORIGINAL_INPUT.root
```

`--analysis-repo` can select a separate frozen runtime checkout. Default entry count is the complete input forest. `--max-events` is only for bounded verification. The two output files are `train_raw.root` and `validation_raw.root`; both contain n=1 and n=2 and every configured activity bin. The receipt, identity ledger and two native logs accompany them. Native physics and calibration sources are unchanged. The cards have truth support5–5020GeV, reco support40–5020GeV, report40–80/80–120/120–160GeV and retain the5,10,20,40 low guards. Track acceptance is pT>2GeV, |eta|<2.4, R4. Detector pair corrections are disabled in the pure generator-particle response.

OO inclusive means all otherwise-selected MB events, with no activity event cut. Stored activity bins[-1,0,30,50,90,100] retain the requested splits and90–100 bookkeeping. The inclusive projection is the exact union[-1,100], after verifying no accepted response content lies in activity flow. Flow causes a merge failure requiring explicit accounting; no event is silently rejected and no unity content is invented. pp uses[-1,100]. `SourceJCard` preserves the full native card; `InputJCard` preserves the actual projection input, and `JCard` describes the projected product.

After every frozen input task has a complete receipt, create a JSON list of local worker output directories and run:

```sh
python3 validation/run3Unfolding/batch/merge_response.py INPUTS.json OUTPUT_DIRECTORY --system OO --manifest FROZEN_ALL_MC_MANIFEST.json
```

The merge checks exact full-file coverage, input sizes and namespaces, coherent native/helper/card hashes, file checksums, full JCard equality and globally unique physical identities using a disk-backed SQLite index. It merges raw weighted histogram counts separately for training and validation before projecting. It writes raw files and n=1/n=2 standard projected files for both splits and inclusive OO; pp gets inclusive products. No normalized matrices are averaged. No train+validation combination is substituted for independent validation. `--no-project` retains only the validated raw matrices. `--allow-bounded` explicitly bypasses full-file coverage only for engineering checks and is not a production option.

Projection uses the unchanged standard macro through `project_response.py`, mask108803. The legacy fixed-edge closure matrix is excluded. The standard native one-dimensional response, EEC response and within-jet covariance are retained. ROOT updates occur in local scratch, with one checked copy to final storage. These are engineering response products: the inherited generator factor1 is still unvalidated for cross-section normalization; within-jet covariance is not certified full event covariance. This workflow performs no unfolding, iteration choice or new calibration payload derivation.

Verification evidence: `../results/2026-11-03-batch-response-checks/`. Seven fixture regressions cover synchronized exact partition, old partial cycles, malformed trees, duplicate identities, missing manifest coverage and changed outputs. Real1000-entry OO/pp split outputs reproduce all seven unsplit native response/covariance entry counts; cell contents agree within1.52e-6 from THnSparseF accumulation order. Real inclusive OO N2 projection preserves the response integral exactly and retains all three card views.
