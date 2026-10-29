# Archived Run-2 local cards

These seven cards moved from the repository root without changing any byte.
Their historical physics settings are retained, including settings that differ
from their operational CRAB counterparts. This archive does not certify those
settings or change the analysis code.

| Original root path | Current path from repository root | `DataType` |
|---|---|---:|
| `cardEECPbPb.input` | `cards/archive/run2/cardEECPbPb.input` | 1 |
| `cardEECPbPbMC.input` | `cards/archive/run2/cardEECPbPbMC.input` | 3 |
| `cardEECpp.input` | `cards/archive/run2/cardEECpp.input` | 0 |
| `cardEECppMC.input` | `cards/archive/run2/cardEECppMC.input` | 2 |
| `cardEECpPb.input` | `cards/archive/run2/cardEECpPb.input` | 4 |
| `cardEECpPb5TeV.input` | `cards/archive/run2/cardEECpPb5TeV.input` | 6 |
| `cardEECpPbMC.input` | `cards/archive/run2/cardEECpPbMC.input` | 8 |

For a local legacy run, work from the repository root and supply the archived
card path. For example, after supplying a compatible legacy forest:

```bash
./eecAnalysis /path/to/legacy_PbPb_forest.root cards/archive/run2/cardEECPbPb.input legacy_PbPb.root 1 true
```

The forest and output paths above are placeholders; this command was not run.
Use the corresponding archived card for legacy pp, pPb, or MC. Existing
`crab/*.input` cards and CRAB configurations retain their original paths.
They remain independent operational copies, not links to this archive. Current
Run-3 Jet ID/veto settings are described in the [main README](../../../README.md).

See [the complete card inventory](../../README.md) and
[the baseline hashes and parameters](../../inventory.json).
