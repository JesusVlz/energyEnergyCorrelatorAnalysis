# Analysis card layout

The root contains the existing OO2025 and ppRef2024 5.36 TeV local cards.
The seven legacy Run-2 local cards are retained in
[archive/run2/](archive/run2/README.md). Operational copies in `crab/` and
bounded diagnostic cards in `validation/` keep their original paths.

The seven Run-2 card moves preserve every byte, including selection, binning,
trigger and correction settings. Current Run-3 settings are documented in the
[main README](../README.md). The existing
R=0.8 workstream remains separate; its presence here does not establish
calibration availability or override the existing analysis guards.

## Root Run-3 cards

| Family | Data card | MC card |
|---|---|---|
| OO, existing split R=0.4 | [cardEECOO.input](../cardEECOO.input) | [cardEECOOMC.input](../cardEECOOMC.input) |
| OO, inclusive R=0.4 | [cardEECOOInclusiveR04.input](../cardEECOOInclusiveR04.input) | [cardEECOOMCInclusiveR04.input](../cardEECOOMCInclusiveR04.input) |
| ppRef, base R=0.4 | [cardEECpp5v36TeV.input](../cardEECpp5v36TeV.input) | [cardEECpp5v36TeVMC.input](../cardEECpp5v36TeVMC.input) |
| ppRef, dedicated R=0.4 | [cardEECpp5v36TeVR04.input](../cardEECpp5v36TeVR04.input) | [cardEECpp5v36TeVMCR04.input](../cardEECpp5v36TeVMCR04.input) |
| OO, separate R=0.8 | [cardEECOOInclusiveR08.input](../cardEECOOInclusiveR08.input) | [cardEECOOMCInclusiveR08.input](../cardEECOOMCInclusiveR08.input) |
| ppRef, separate R=0.8 | [cardEECpp5v36TeVR08.input](../cardEECpp5v36TeVR08.input) | [cardEECpp5v36TeVMCR08.input](../cardEECpp5v36TeVMCR08.input) |

## Complete inventory snapshot

[inventory.json](inventory.json) records the exact HEAD and branch at inventory,
all 37 old and new paths, SHA-256 and Git blob hashes, and every parsed card
parameter. It describes the working tree on 28 October 2026, including
pre-existing local edits. It is a baseline snapshot, not a generated list that
automatically tracks later edits.

| Location after reorganization | Cards | Role |
|---|---:|---|
| Repository root | 12 | OO2025/ppRef2024 local cards |
| `cards/archive/run2/` | 7 | Legacy PbPb/pp/pPb local cards |
| `crab/` | 10 | Independent operational CRAB cards, including two response variations |
| `validation/run3Workflow/cards/` | 8 | Existing diagnostic fixtures |
| `skimmer/` | 0 | No `.input` files found |

The inventory command, run from the repository root, is:

```bash
rg --files --hidden --no-ignore -g '*.input' -g '!.git/**'
```

All CRAB card references retain their original directory semantics. The main
executable reads the card path supplied as its second argument; it has no
hard-coded legacy default. The existing tarball script packages the analysis
code and corrections, while CRAB configurations stage their own card through
`config.JobType.inputFiles`. No launcher or packaging path needed changing
for these seven moves.
