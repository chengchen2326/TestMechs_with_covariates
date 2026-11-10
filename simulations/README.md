# Monte Carlo for the covariate extension

These scripts are **not part of the package**. They call it. Keeping them out
of `R/` is deliberate: `R CMD check` would treat simulation scripts as package
code, and the new estimators need to be testable on their own.

## Running

The scripts `source()` the package files rather than requiring the package to
be installed, because `TestMechs` imports `HonestDiD` (GitHub only), `fixest`
and `osqp`, none of which the simulation touches. `lb_frac_affected` itself
needs only `purrr`, `Rglpk` and base/stats.

Sourcing needs the `TestMechs:::` prefixes stripped, since there is no loaded
namespace. Generate the stripped copies into `src/` (git-ignored):

```sh
mkdir -p src
for f in remove_missing_from_df lb_frac_affected simulate_data_binaryM \
         simulate_data_binaryM_x lb_frac_affected_x; do
  sed 's/TestMechs:::/ /g' R/$f.R > src/src_$f.R
done
```

Then, from `simulations/`:

```sh
TM_SRC=../src TM_REPS=1000 Rscript design_A.R
```

Output goes to `results/design_A.csv` (git-ignored). `TM_REPS=3` for a smoke
test.

Installing the package normally and using `library(TestMechs)` works too, and
is preferable once the dependencies are available; the `src/` route exists only
so the simulation can run in a bare environment.

## Designs

| script | design |
|---|---|
| `design_A.R` | Binary Y, binary M, one binary covariate. |

Planned: B (many cells, then continuous X, where binning is infeasible),
C (multivalued M calibrated to Baranov), D (continuous Y, on Baranov).

## Arms

| arm | what it is |
|---|---|
| `ours` | `lb_frac_affected_x(method = "split")` — stage 1 on one fold, ATE on the other |
| `pooled` | `method = "pooled_split"` — the same inference ignoring covariates, i.e. a confidence statement for the pooled bound |
| `plugin_x` | `method = "plugin"` with covariates — the supremum taken on the data used to report |
| `plugin_0` | `method = "plugin"` without covariates — current practice |

`plugin_*` arms report a point estimate only. A normal confidence statement
would not be valid for them, which is the reason the comparison is here.

## Calibrations in `design_A.R`

- **`sign`** — the sign-changing design. The binding outcome value differs
  across cells, so the pooled bound is exactly 0 while the conditional bound is
  0.06. Any power in the `ours` arm is attributable to conditioning, and the
  `pooled` arm should have none at any sample size.
- **`bursztyn`** — the empirical Bursztyn distribution split by `college_deg`.
  The cell bounds are 0.1333 and 0.0130 but the average is 0.0681 either way,
  so the gap is 0.0000. This is the honest calibration and it is expected to
  show no gain in the average bound. It is included so the simulation reports
  the case that looks like real data, not only the case that flatters the
  method.

## Reference values

Regression targets, for checking any change to the estimators:

| quantity | value |
|---|---|
| `lb_frac_affected` on Bursztyn, `at_group = 0` | 0.1067776457 |
| `lb_frac_affected_x(x = const, method = "plugin")`, same data | must agree to ~1e-16 |
| `population_bound_binaryM_x` on the `sign` calibration | cells 0.06/0.06, `E[eta(X)]` 0.06, pooled 0, `nu` 0.15 |

The first is Kwon & Roth's published figure ("at least 11 percent of
never-takers"). The second is the test that the ATE reformulation is right: with
a constant covariate the new estimator must collapse onto the original.

## Standalone verification

These scripts check a claim in the paper rather than estimating anything. They
need no data and no `src/` setup — `Rscript <file>` is enough.

| script | claim checked |
|---|---|
| `verify_dual.R` | Proposition 1 is the Kantorovich dual restricted to one direction: `(-g_k, g_k)` is dual-feasible for every `A`, its dual objective is `psi_k(A)`, and `psi_k(A)` never exceeds the primal transport value. Includes a falsification check — drop monotonicity and feasibility fails in 87% of draws. |
| `verify_defiers.R` | The defier LP and its dual agree, and the normalised duals are free of the box. |
| `verify_C.R`, `verify_C_exact.R` | The Design C population bound matches the realised law. |
| `nu_check.R` | The share scale via Charnes--Cooper reproduces `lb_frac_affected`. |

`verify_dual.R` runs in about two minutes and prints four tolerances; all four
should be at machine precision, and the last line should report roughly 87%.
