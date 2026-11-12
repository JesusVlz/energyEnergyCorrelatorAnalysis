# Combined Run-3 integration testing

Updated: 2026-11-11. Main promotion is explicitly authorized by the user.
See [OverallCoordination.md](OverallCoordination.md) for the consolidated scope.
This is a software integration, not a qualified correction payload, production
campaign, or final physics result.

## Integrated revisions

| Area | Source | Integration |
| --- | --- | --- |
| Base including merged jet calibration/ID/QA and Nebraska input handling | `main` at `adac8d1256e30ff09e7975fea89b29f01b332b15` | Fast-forward from `fd575b3`; main not changed |
| Tracking/pair preparation | PR #2, `c77d46006e0e5ed9b00f1d6880ba36e353db5c2b` | `30410ff94cd2a4d3ad5b93f25cf7b5fc27349a09` |
| R4 unfolding preparation | PR #3, `ae66819997dda3072127271342cf926c9e1da78e` | `a9ad2e30e14f582874fe426245abc0f14875672f` |
| Legacy-style Run-3 closure plotting | PR #4, `fdc6e02b024d5ea346b4c0f0c1daf6dd07da34a3` | `f88c0f672b98213875cbfea0a345f0296ead7a3b` |
| Slim Run-3 track inputs and exact PF multiplicities | PR #5, `b2d61d8840614fe8a966b1e53501131545ba1278` | `7b91e6ddc7aca9fa77312e3362cb554113f1af40` |

PRs #2-4 were retargeted from main to oo-integration before merging. PR #5
also targets oo-integration. Source branches are preserved. Jet work was already
in main through #1 and is not applied twice. Merge commits preserve ancestry.

## Review findings and retained behavior

- Combined HighForestReader retains physical reconstructed charge, disabled
  PbPb-only fields for slim Run-3 inputs, and exact PF multiplicity handling.
  Missing required charge and unresolved Jet ID multiplicities remain errors.
- R4 JEC, detector-coordinate lookup, additive trackMax/rawPt and Jet ID,
  whole-event veto, and QA remain from the jet update. pp L1FastJet stays OFF
  in current cards; no rho or jet-area substitute is introduced.
- Tracking uses the Nominal/Loose/Tight prescription and matching existing
  tables. New exports/pair payloads carry their updated selection versions.
  Pre-repair pair maps must not be relabeled as compatible with new selections.
- Residual pair corrections remain disabled in baseline cards until supported
  maps and independent closure qualify them. No tables or uncertainty sizes
  were derived by merging these branches.
- The unfolding adapter reuses the native classes and RooUnfoldBayes, keeps the
  two folds and 1-20 scan, and retrieves Eunfold/Ereco. The legacy
  plotting/unfoldEEChistograms.C covariance path is NOT repaired in-place.
- Closure plots preserve the legacy arithmetic by default. Symmetric same-cone
  averaging is an explicit option. Shared-signal approximate errors are NOT
  full covariance, and cannot be used for unrelated normalized-shape ratios.
- R8 remains raw exploratory scope. This merge does not supply an AK8 calibration
  or claim that the separate raw-only execution gap has been resolved.
- No event weighting, background matching or optional HF-shift policy is changed
  merely to make the integration tests pass. Integrated OO remains normal
  centrality binning/projection, not a separate InclusiveMode.

## Read-only branch CI

`.github/workflows/oo-integration-validation.yml` runs only for oo-integration
and main pushes (or an explicit dispatch). Its token has contents:read, checkout does not
persist credentials, and it cannot commit, push, merge or submit analysis jobs.
It builds the analyzer and drawing library and runs tracking/table/card, jet,
PF-count, unfolding numerical, and closure-macro compilation checks.

The installed ROOT wheel is a software-test environment, not a production
CMSSW prescription. The real RooUnfold test may skip if its library is absent;
that skip is recorded explicitly and is NOT a passed runtime or MC test.
Inspect each job/step and logs; a conflict-free merge does not imply a clean build.

## Test in a separate checkout on the analysis host

From an existing repository, preserve local work and use a new worktree:

```sh
git fetch origin
git worktree add --detach ../eec-oo-integration-test origin/oo-integration
cd ../eec-oo-integration-test
git rev-parse HEAD
# Use the established compatible ROOT/CMSSW environment.
make -j2
make -C plotting -j2
python3 -m unittest discover -s validation/run3Pairs -p 'test_*.py' -v
python3 validation/run3Pairs/prepare_inputs.py audit-tables --output /tmp/integration-tracking-tables.json
python3 validation/run3Jets/run_selection_tests.py --output /tmp/integration-jets-NEW
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover -s validation/run3Unfolding -p 'test_study*.py' -v
EEC_REQUIRE_ROOT_TESTS=1 python3 -m unittest discover -s validation/run3Unfolding -p 'test_study_root_runtime.py' -v
```

Use unused worktree/output paths. The EEC_REQUIRE_ROOT_TESTS command is required
on the analysis host and must fail rather than skip when RooUnfold is unavailable.
Consult the tracking, unfolding and plotting READMEs for actual input contracts.

## Physics qualification remaining after software promotion

1. Combined analyzer/drawing builds and software regressions pass. Verify slim
   and legacy forest schemas, real charge, PF counts and ROOT chain transitions.
2. Run bounded current OO and pp MC on the same new revision, preserving the
   existing event identities, both powers and integrated/split R4 definitions.
   Check calibration-stage and selection QA, source-coordinate tracking and
   native/exported pair parity. Do not launch full production from this README.
3. Run the real RooUnfold covariance regression, then at least one compatible
   OO and pp two-fold scan. Check metadata and version compatibility of old
   response/pair products before reuse; old files do not inherit new corrections.
4. Compile and execute the two new Run-3 closure macros on compatible projected
   inputs. Compare central values and explicitly chosen error conventions;
   run checkRun3Closure.C on its documented local fixture.
5. Record actual commands, input/card/payload hashes, outputs and failures.
   Distinguish software PASS, runtime PASS and physics closure. Only then
   qualify the corrections/responses for the measured result.

No real MC production or closure is claimed by this integration record.
