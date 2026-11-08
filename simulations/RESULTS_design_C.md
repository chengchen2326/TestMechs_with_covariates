# Design C — multivalued mediator, calibrated to Baranov

Run 2026-11-08. 300 replications per cell, R 4.4.3, branch
`covariates-simulation`. Code: `R/simulate_data_multiM_x.R`,
`simulations/calib_C.R`, `simulations/design_C.R`,
`simulations/verify_C_exact.R`.

$K = 5$ mediator categories from Baranov's `relationship_husb`, $Q = 5$ outcome
bins from `motherfinancial`. Focal always-taker groups $k = 5$ and $k = 3$.

## 0. Verification

| check | result |
|---|---|
| general-$K$ population bound vs `lb_frac_affected`, population passed as weights | exact: $7.8\times10^{-15}$, $0$, $6.3\times10^{-16}$, $2.1\times10^{-16}$, $5.6\times10^{-17}$ across the five mediator values |
| bound $\le$ true $\theta_{kk}\nu_k$ | holds in every configuration |
| slack under the maximal coupling | exactly $0$ at $k=1$ and $k=5$, positive in between |

The slack pattern is what the identified set implies: at $k=1$ (the bottom) and
at $k=5$ (where $\theta_{55} = 0.337 \approx P(M=5 \mid D=1)$) almost no
compliers enter, so KR's bound is attained; in between, compliers can absorb
the affected mass and the bound is strictly below the truth.

**Coverage, which this design can check against the truth and not only the
bound:** 0.940–0.997 against the bound and $\ge 0.94$ against the truth
everywhere, at a nominal 0.95. Validity holds for multivalued $M$.

## 1. The central claim, now for a multivalued mediator

`signflip` at $k=5$: the two cells' treated outcome laws are perturbed in
opposite directions, so $\eta_5^{\text{pooled}} = 0$ **exactly** while
$\mathbb{E}[\eta_5(X)] = 0.0674$.

| arm | $n=1000$ | $n=4000$ |
|---|---|---|
| ours | **0.67** | **1.00** |
| pooled | 0.06 | 0.06 |

Pooled sits at size at both sample sizes. This is cleaner than the binary
Design A, where ours reached only 0.42 at $n=1000$.

## 2. The first case where conditioning helps on a real-data calibration

`split` — cells from Baranov's `grandmother_baseline`, everything else
empirical — at $k=3$:

| | conditional | pooled |
|---|---|---|
| bound | **0.0289** | **0.0142** |
| power, $n=1000$ | 0.16 | 0.08 |
| power, $n=4000$ | **0.51** | 0.25 |

**The conditional bound is twice the pooled bound, and power doubles.** Same
under both couplings (0.49 against 0.26), so it is not an artifact of how the
always-takers' potential outcomes were coupled.

Everything before today showed no gain: binary $M$ on Bursztyn, binary $M$ on
Baranov, and the identification calculation on both. This is the first
exception, and it is on a real covariate with real marginals.

**Where the gain comes from matters.** At $k=5$ the conditional bound is 0.0741
against a pooled 0.0710 — a 4% gain and no power difference. The gain is
specific to the middle of the mediator's support, where
$\tilde\theta_{kk}^{\min}(x)$ differs across cells and the pooled computation
averages that structure away. So **the gain operates through the compliance
structure of a multivalued mediator, not through the outcome distribution** —
which is the mechanism the binary designs could not exhibit, because with
$K=2$ there is no middle.

That reorients the paper's positive case: it lives in multivalued $M$, which is
also where KR's own contribution lies.

## 3. What conditioning costs when it gains nothing

`homog` — identical cells — at $k=5$:

| arm | $n=1000$ | $n=4000$ |
|---|---|---|
| ours | 0.58 | 0.99 |
| pooled | 0.80 | 0.99 |

A real cost at $n=1000$ (0.58 against 0.80), gone by $n=4000$. The same at
$k=3$: 0.18 against 0.30 at $n=4000$. So conditioning is not free, and the
paper should report the cost next to the gain.

## 4. How lossy the identified set is

Under the independent coupling the truth is far above the bound:

| $k$ | bound | truth | bound / truth |
|---|---|---|---|
| 5 | 0.0776 | 0.2470 | 0.31 |
| 3 | 0.0170 | 0.1640 | **0.10** |

At $k=3$ the bound recovers a tenth of the true affected mass. This is a
property of KR's identified set, not of any estimator, and it is measurable
only because this DGP knows the truth. Worth stating in the paper: the bound is
informative about whether the mechanism story fails, and weakly informative
about by how much — more weakly in the middle of the mediator's support than at
its edges.

## 5. Plug-in bias, proportionally worse for small bounds

At $n=1000$, mean plug-in estimate against its own population bound:

| $k$ | arm | bound | plug-in | bias |
|---|---|---|---|---|
| 5 | `plugin_x` | 0.0776 | 0.0844 | +9% |
| 5 | `plugin_0` | 0.0776 | 0.0810 | +4% |
| 3 | `plugin_x` | 0.0170 | 0.0243 | **+43%** |
| 3 | `plugin_0` | 0.0170 | 0.0193 | +14% |

Consistent with Design A: the winner's curse bites hardest where the bound is
small, which is where a practitioner most needs to know the answer is "not
much".

## Open after Design C

- The `signflip` perturbation is applied at one mediator value, so the $k=3$
  rows of that configuration are effectively a no-complier homogeneous design.
  A configuration whose sign flip is spread across mediator values would test
  whether the gains compound.
- $K = 10$ is untested; Baranov's mediator has only five categories, so it
  needs either a different mediator or a synthetic design.
- The $k=3$ gain should be traced analytically to
  $\tilde\theta_{kk}^{\min}(x)$ rather than left as a numerical observation.
  That is the one piece of theory this result asks for, and it would make the
  paper's positive case an argument rather than a table.
