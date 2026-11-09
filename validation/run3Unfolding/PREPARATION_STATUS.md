# R4 Run-3 preparation: completed work and remaining validation

Reviewed base: `792d76c0cd043ab1871b3217ac988f6455f2d9f0`.
Scope: MC-only preparation for the HIN-style response, iteration and closure plots.

## Implemented in this change

1. Manifest-driven Run-3 input and physical activity/power/bin selection.
2. Explicit 40-80, 80-120, 120-160 GeV output, retaining native migration guards.
3. Reuse of the existing train/validation products and worker identity ledgers.
4. Both independent directions with 1-20 RooUnfoldBayes iterations.
5. Required measured per-jet covariance and correct result covariance getter.
6. Truth-space covariance checks and normalized-shape Jacobian propagation.
7. Separate response-defined and optional ordinary RecoGen/GenGen closure.
8. K/Kprime conditioning, misses/fakes, support and explicit-flow checks.
9. Fixed-exposure efficiency-aware refolding helper, without guessing exposure.
10. Split-1/split-2/average unity-chi2 diagnostic and one candidate per response.
11. Response, iteration and three-panel shape/yield closure plotting code.
12. Input/source hashes, channel status, command records, examples and tests.

Native response production, selections, calibrations, background code, R8 and
Run-2 entrypoints are unchanged. In particular, the old macro's covariance-output
bug is NOT fixed in-place; the new Run-3 path bypasses it. Do not call that old
path repaired or use it as the validated Run-3 entrypoint.

## Actually executed in the preparation environment

- 43 numerical, metadata, API-mock, ledger and synthetic PNG renderer tests passed.
- Python syntax compilation passed for all new modules.
- Both command-line help entrypoints executed successfully.
- One real RooUnfold test class was skipped because ROOT/PyROOT is absent.

Commands:

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover \
  -s validation/run3Unfolding -p 'test_study*.py' -v
python -m compileall -q validation/run3Unfolding
python validation/run3Unfolding/run3_study.py --help
python validation/run3Unfolding/draw_run3_study.py --help
```

Synthetic plots are software tests, not physics evidence. No OO/pp MC unfolding,
closure plots, iteration determination or full CMS-checklist completion is claimed.

## Next required execution

On the established ROOT host, run the real regression with
`EEC_REQUIRE_ROOT_TESTS=1`, bind one supported OO channel and pp to their existing
projected folds and production evidence, then run check, scan and drawing.
Inspect 40-80 GeV support and both-fold closure before extending the channel grid.

Still needed for physics validation: retention/weight/calibration and guard tests;
joint closure covariance and response statistics; prior/model/unregularized,
refolding/bottom-line/coverage studies; full corrected/background-subtracted closure.
The diagonal closure chi2 is diagnostic only, not a calibrated p-value or freeze.
