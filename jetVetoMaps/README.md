# Run-3 jet-veto map deployment contract

The owner-uploaded OO2025 and ppRef2024 ROOT maps match the pinned source bytes
recorded in their metadata sidecars. The R=0.4 local and CRAB cards enable them
for production tests with `jetvetomap_all`, stored nominal PAT/HiForest `jtpt`
and an owner-declared 15 GeV input-retention bound. Their deployment status is
`ENABLED_FOR_PRODUCTION_TESTS`. The tarball script includes this directory.

See [the research report](../validation/run3Jets/vetoMapSources.md),
[the candidate inventory](../validation/run3Jets/vetoMapCandidates.json) and
[the analyzer configuration](../validation/run3Jets/IMPLEMENTATION.md).
The research artifacts preserve the earlier candidate-stage findings.

Every deployed ROOT file has a sidecar satisfying `metadata.schema.json`.
The sidecar pins the source commit and SHA-256, histogram key and detector
coordinates, and retains validation and owner-authorization evidence.
Full numeric run coverage, precise collection/radius applicability and upstream
retention completeness remain to be established for these campaigns.
`APPROVED_FOR_PRODUCTION` requires those fields to be resolved and
`exact_production_applicability_established` to be true.

The analyzer uses the explicit card file and histogram key; it never searches
these directories or infers an era. Bounded MC tests verify software behavior
with the selected inputs, without establishing campaign coverage or physics
closure. Use the metadata sidecars for current deployment status.
