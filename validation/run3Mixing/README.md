# Run-3 two-partner mixing preflight

`preflight.py` measures pool feasibility using native jets and stored measured
`hiHF_pf` and vertex coordinates. It never fills or subtracts background EECs.
The HF/vz tolerances are explicit **study settings**, not an approved matching
prescription. It applies no HF shift, centrality replacement, matching weights,
pair correction, component normalization or clipping.

The bounded scan saves every selected event, including zero-jet events, to an
adjacent `.events.jsonl` ledger. It retains real hiBin (including unsupported
values), signed generator weights, namespace/run/lumi/event identity, source
entry, native radius and jet bin. Source/accepted 0–90, 90–100 and invalid hiBin
counts are recorded independently of the matching coordinate. Invalid hiBin in
pp is bookkeeping, not an OO percentile acceptance diagnosis.

Required MC filters follow the current Run-3 reader: OO primary vertex,
`OOphfCoincFilterPF2Th4`, pileup; pp primary vertex and `ppNoScrapingFilter`;
abs(vz)<=15. Missing required branches fail. No data trigger is inferred. Source
skim losses and MB representativeness cannot be established from selected MC.
Native pp R8 absence is explicit and never replaced by AK4 jets. Stored jtpt is
provisional while applicable calibration remains unresolved.

Each signal physical event gets two distinct other physical events in the same
sample namespace and FNV training/holdout fold. Vertex and HF mismatches must
satisfy both requested tolerances; sparse pools fail without relaxing them.
The stored HF range and distinct-value count are audited. A constant or empty
coordinate is explicitly BLOCKED_DEGENERATE_COORDINATE; feasible partner counts
then describe vertex-only feasibility and cannot establish activity matching.
Candidates are ranked by normalized L1 mismatch, then stable event identity.
Partners are shared between radii/powers of the same signal event. Reuse counts,
failed signals, jet-bin coverage and signed weights remain visible. These records
support later covariance work; they are not a covariance calculation. In
particular, partner reuse is not independent exposure.

Example (choose a new output path; existing outputs are protected):

```bash
python3 -m unittest discover -s validation/run3Mixing -p 'test_*.py' -v
python3 validation/run3Mixing/preflight.py \
  /eos/user/j/jvelazqu/Foresting/OOMC/unembedded_L1HLT_objects_MC_1.root \
  /tmp/OO_pool_study.json --system OO \
  --sample-namespace OO_local_unembedded_L1HLT_MC_1 \
  --events 10000 --vz-tolerance 1 --hf-tolerance 20
```

Use the same verified physical-sample namespace as the pair exporter. The
example namespace identifies the local audited file, not certified campaign
provenance. Study tolerance sensitivity by importing `match` and rerunning it on
the saved ledger. Do not optimize tolerances on an OO/pp physics ratio.

Promotion requires independent partner-source applicability, matching and jet-HF
bias studies, adequate channel coverage, measured data/MC compatibility, component
exposures/correlations, R8 overlap treatment and extracted-signal MC closure.
The 2 November method comparison and the inherited final correction order remain
unchanged. No production campaign or nominal estimator is authorized by this
preflight.
