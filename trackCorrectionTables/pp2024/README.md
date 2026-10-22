# Run-3 ppRef tracking input

The analyzer uses
`Eff_ppref_2024_Pythia_QCD_pThat15_NopU_2D_Nominal_Official.root`
for data types 11 and 12. This nominal file is explicitly exempted from the
repository's `*.root` ignore rule so it can be committed with the code, as the
legacy pp/PbPb tracking inputs are. Loose and Tight files remain unused.

The exact existing nominal file is pinned here: 39,482 bytes, SHA-256
`f92495a0c0bfc0738e3fdcd06e385eaf643b86b077b1430bede17239c7138ee8`.
It contains four `TH2F` histograms, `hEff_2D`, `hFak_2D`, `hSec_2D`, and
`hMul_2D`. All four have identical variable-width axes: 14 eta bins on X
covering [-3, 3), and 50 track-pT bins on Y covering [0, 500) GeV.

`TrkEff2024ppRef` loads this nominal file for both ppRef data and MC and uses
the bin containing `(eta, pT)` without interpolation. Its correction is
`(1 - fake) * (1 - secondary) / (efficiency * (1 + multipleReco))`.
The centrality argument is unused. Values outside the axis ranges, or a
denominator at or below 0.001, return zero. The deployment update changes
neither this formula nor the file content; Loose and Tight tables remain
ignored and unselected.

The file was already present locally and in `crab/eecOO.tar.gz` at the time of
the technical audit. No authoritative download location or derivation record
was recorded in this repository. The checksum establishes which file is
deployed; it does not establish independent provenance or MC closure.

From the repository root, verify the nominal input checksum with:

```sh
sha256sum trackCorrectionTables/pp2024/Eff_ppref_2024_Pythia_QCD_pThat15_NopU_2D_Nominal_Official.root
```

Include this exact file in any newly built CRAB payload. Changing the source
checkout does not refresh an existing `crab/eecOO.tar.gz` archive.
