# OO2025 maps

The owner-uploaded `Summer24Prompt25_RunCDEFG.root` matches the immutable JECDatabase source
recorded in its `.metadata.json` sidecar. Its current deployment status is
`ENABLED_FOR_PRODUCTION_TESTS`. R=0.4 local and CRAB cards select `jetvetomap_all`
and stored nominal PAT/HiForest `jtpt`, with a declared 15 GeV input bound.
The tarball script includes the map and metadata.

Numerical run coverage, exact campaign applicability and full upstream retention
coverage remain to be established. See [the deployment contract](../README.md),
[the historical inventory](../../validation/run3Jets/vetoMapCandidates.json) and
[the analyzer configuration](../../validation/run3Jets/IMPLEMENTATION.md).
