# ppRef2024 residual pair candidates: 16 November 2026

Status: PRELIMINARY_NOT_PHYSICS_QUALIFIED.

Install the three original ROOT binaries with `cards/install_ppref_pair_tables.py`
using the supplied campaign archive, then run `--check`. The installation handoff
also provides a small tar.gz containing only those original ROOT files.
`manifest.json` identifies exact paths, sizes, SHA256, source archive members,
selection names, pair/track versions and matching single-track table checksums.
A manifest alone is not an installed ROOT payload.

Each file contains normalized reco/gen histograms and validity/count/control
information. The native reader applies the inverse. Installation preserves every
original byte. Derivation used JetID OFF and L1FastJet OFF; requested application
cards keep JetID ON, and closure for that selection is still pending.

See `cards/ppRef2024/R4/README.md` for activation, ordinary git publication and
worker packaging. OO and R8 corrections are not enabled.
