# ppRef2024 R4 preliminary pair application

These six Data/RecoReco cards enable the original 16 November 2026 Nominal,
Loose and Tight residual pair candidates, with JetID ON and L1FastJet OFF.
They are application tests, not a new detector-closure result.

## Install the actual ROOT files first

A card or manifest is NOT an installed ROOT table. From the repository root,
in the analysis ROOT/PyROOT environment, install the three original binaries:

```sh
python3 cards/install_ppref_pair_tables.py \
  --archive /path/to/ChatGPT_tracking_closure_20261116_151819_UTC.tar.gz --write
python3 cards/install_ppref_pair_tables.py --check
```

The small `ppref2024_pair_payloads_20261116.tar.gz` supplied with the installation
handoff is also accepted: it contains exactly the same three regular archive
members. SHA256 values in `trackCorrectionTables/pp2024/pair/20261116/manifest.json`
pin every original byte and the corresponding single-track payload. The installer
also checks ROOT tags, axes and validity masks. It never executes archive code,
decodes embedded code, contacts a server, launches MC jobs, commits or pushes.

## Apply on reconstructed tracks

```sh
make -j2
./eecAnalysis /path/to/ppRef_MC.root \
  cards/ppRef2024/R4/ppRef_MC_Nominal_pairOn.input \
  ppRef_MC_pairOn_preliminary.root 0 true
```

Replace MC with Data, and Nominal with Loose or Tight, for the other cards.
Each uses the matching single-track and pair working point. The ordinary EEC
weights contain each detector factor once, for both n=1 and n=2.

These bounded cards use 1000 events per file, reconstructed R4 WTA jets in
15-200 GeV, |etaJet|<1.6, 0.5 GeV track input retention and strict pair pT>2 GeV.
The 15-20 GeV interval uses the stored 10-20 GeV correction bin, not a separately
derived map. Background, systematic and particle/unfolding response filling are
OFF. Response guard settings are not changed. The data trigger is inherited;
check its applicability to the actual dataset before interpreting results.

## Explicit activation of the existing pp cards

After installation, enable the six nominal R4 compatibility cards:

```sh
python3 cards/install_ppref_pair_tables.py --check --write --activate-base-cards
```

This updates the four root-level ppRef cards and the two ppRef CRAB cards.
It preserves their event/jet selections, background settings, reporting bins and
response guards, while enforcing JetID ON, L1FastJet OFF and MinPairTrackPt 2.
It rejects response/systematic or non-RecoReco configurations. OO is unchanged.
Compatibility cards remain pair-OFF until this explicit installation/activation;
the separate generated `cards/R4` response catalog retains its original policy.

Publish the installed binary files and activated compatibility cards normally:

```sh
git add trackCorrectionTables/pp2024/pair/20261116 \
  cardEECpp5v36TeV.input cardEECpp5v36TeVMC.input \
  cardEECpp5v36TeVR04.input cardEECpp5v36TeVMCR04.input \
  crab/cardEECpp5v36TeV.input crab/cardEECpp5v36TeVMC.input
git diff --cached --stat
git commit -m "Install original pp2024 pair ROOT tables and enable R4 cards"
git push origin main
```

The three exact binary paths are allowed by the local `.gitignore`. Rebuild worker
archives after installation; a repository-relative path must also exist inside
the worker package. No worker archive is rebuilt or submitted by these helpers.

## Physics and software limits

Original derivation: JetID OFF, L1FastJet OFF. This provenance is preserved.
JetID-on independent-event closure and the tracking/pair uncertainty conclusions
remain pending. Full support in the earlier no-JetID holdout is not a new closure.
No OO map is enabled, no unsupported bin is assigned unity, and no uncertainty
is fabricated. The installed tables are not approved final physics payloads.

The generic native consumer does not yet enforce all SHA/version constraints;
run `install_ppref_pair_tables.py --check` before application and packaging.
Keep broad particle-response diagnostics below 2 GeV separate from pair-corrected
EEC application. The C++ `check_pp2024_installed_tables.cxx` fixture tests lookup
inversion, supported fallback and domain rejection, not detector closure.
