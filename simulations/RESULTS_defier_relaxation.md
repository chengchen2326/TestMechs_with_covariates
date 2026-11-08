# The defier relaxation, and how to impose it

Run 2026-11-08. Code: `R/theta_kk_min_defiers.R`, `R/eta_cond_joint.R`,
`simulations/verify_defiers.R`, `simulations/defier_gain.R`,
`simulations/compare_imposition.R`.

## Why it was needed

Baranov's five-valued mediator `relationship_husb` — KR's own main
specification — violates monotonicity in sample, with a minimum defier share of
0.0154. Design C handled this by projecting the marginals onto the monotone
region. That projection is not innocuous: it is what produced the $k=3$ gain
reported for Design C. The raw-data question could only be answered by
implementing the relaxation.

## The LP and its dual

With a defier cap, KR's closed form for $\tilde\theta_{kk}^{\min}$ no longer
applies and the value comes from
$$\min \theta_{kk} \ \text{ s.t. } \textstyle\sum_j \theta_{lj} = a_l,\ \sum_l \theta_{lj} = b_j,\ \sum_{l>j}\theta_{lj} \le \bar d,\ \theta \ge 0.$$
Size $K^2$, no outcome anywhere: the separation result survives the relaxation.

The dual is the point:
$$\max_{u,v,w\ge0}\ a'u + b'v - \bar d w \quad \text{s.t. } u_l + v_j - w\mathbf 1\{l>j\} \le \mathbf 1\{(l,j)=(k,k)\}.$$
Weak duality gives $a'u + b'v - \bar d w \le \tilde\theta_{kk}^{\min}$ for every
dual-feasible triple, and $\tilde\theta_{kk}^{\min}$ enters the bound with a
positive sign, so substituting a triple chosen on one fold can only lower the
bound — conservative, never invalid. The dual objective is linear in the
marginals, so the second stage stays a linear functional. **Exactly the same
device as fixing the set $A$.**

**Verification.** 1,451 random tests: at $\bar d = 0$ the LP reproduces KR's
eq. (8) to $2.8\times10^{-16}$ and the duality gap is $\le 2.2\times10^{-16}$.
At $\bar d > 0$ the implied $\nu_k$ matches
`lb_frac_affected(max_defiers_share = ·)` exactly at every $(\bar d, k)$ tried.

## The N2 question, answered

Two ways to impose the cap. **Cell by cell**, each cell getting its own
$\bar d$, raised to its own minimum feasible share; or at the **population
level**, one cap on $\sum_x w(x)\sum_{l>j}\theta^{(x)}_{lj}$, which couples the
cells into one program:
$$\min_{\{\theta^{(x)}\}} \sum_x w(x)\Big(S(x) - b^{(x)}_k + \theta^{(x)}_{kk}\Big)_+ .$$
The positive part takes one auxiliary variable per cell, so this is an LP with
$JK^2 + J$ variables — still free of the outcome support.

**Dominance holds for the population version, by construction.** Given any
jointly feasible $\{\theta^{(x)}\}$, the aggregate
$\bar\theta = \sum_x w(x)\theta^{(x)}$ has the pooled marginals and defier mass
at most $\bar d$, so it is pooled-feasible; convexity of $(\cdot)_+$ together
with $S_{\text{pool}} \le \sum_x w(x) S(x)$ makes the pooled objective at
$\bar\theta$ no larger than the conditional objective at $\{\theta^{(x)}\}$.
Minimising the right side gives $\eta^{\text{pooled}} \le \mathbb{E}[\eta(X)]$.

**It fails for the cell-by-cell version**, because a cell whose own marginals
need more defiers than the pooled marginals do faces a strictly larger
identified set than the pooled problem. Measured, on Baranov at $k=5$:

| covariate | $\bar d$ | pooled | cond (cell) | gain |
|---|---|---|---|---|
| `grandmother_baseline` | 0.020 | 0.069955 | 0.059593 | **−0.0104** |
| `MIL` | 0.020 | 0.072772 | 0.060454 | **−0.0123** |

The population version shows no such failure anywhere in the grid.

**And it is also tighter**, at $\bar d = 0.03$:

| covariate | $k$ | pooled | cond (cell) | cond (pop) |
|---|---|---|---|---|
| `grandmother_baseline` | 3 | **0.0000** | 0.0092 | **0.0146** |
| `grandmother_baseline` | 4 | **0.0000** | 0.0115 | **0.0169** |
| `MIL` | 4 | **0.0000** | 0.0088 | **0.0167** |

47–90% tighter than the cell version. So there is no trade-off on this axis:
**the population-level imposition dominates, and the N2 decision should be
reversed.**

The one cost is feasibility. The joint program's minimum feasible cap is the
$w$-weighted average of the cells' minima — $0.5 \times 0 + 0.5 \times 0.0407 =
0.0204$ for `grandmother_baseline` — so it needs a slightly larger $\bar d$ than
the cell version. That weighted average is itself a useful diagnostic and should
be reported alongside any conditional bound.

**Checks.** On synthetic monotone cells at $\bar d = 0$ the joint program equals
the cell-by-cell computation to $8.3\times10^{-17}$, as it must, since the cap
then forces every cell's defier mass to zero.

## The empirical headline

On raw Baranov, five-valued mediator, no projection, $\bar d = 0.03$: the pooled
bound is **exactly zero** at $k=3$ and $k=4$, while the conditional bound is
positive — 0.0146 and 0.0169 with `grandmother_baseline`, 0.0167 at $k=4$ with
`MIL`. Pooling finds nothing; conditioning finds something.

This is the result the project was looking for, and unlike the Design C version
it does not rest on a monotone projection.

## A trap that has now bitten twice

The pooled comparison must be built by mixing the cells with the common weights,
$\sum_x w(x) a^{(x)}$ and $\sum_x w(x) b^{(x)}$ — not from the raw pooled
tabulation, which weights $p_d(\cdot \mid x)$ by $P(X = x \mid D = d)$ and
differs across arms in finite samples even under randomisation. Using the raw
tabulation made the population-level program appear to violate dominance at
$k=5$. The same trap appeared in the morning's identification calculation. It
needs a prominent warning in the paper, not a footnote.

## Next

The dual-variable route is built and verified but not yet wired into
`lb_frac_affected_x`, so all the numbers above are identification-level. Wiring
it in gives inference on the headline result, and that is the immediate task.

---

# Addendum — inference, and why the headline is not detectable

## What was built and verified

`R/lb_frac_affected_x_defiers.R`: stage 1 solves the dual of the joint program
on one fold, stage 2 evaluates the resulting linear functional on the other.

| check | result |
|---|---|
| strong duality, joint program | 2,359 tests, max $\|$primal $-$ dual$\|$ = $1.4\times10^{-16}$ |
| stage-2 transformed-outcome form vs the dual objective | identical; with an oracle stage 1 the estimate equals the bound exactly in every configuration |
| boxed dual still attains the primal | 1,821 tests, loose in 0 cases at a box of $\pm 2$ |

The normalised dual constraints involve only $\tilde p_x \in [0,1]$ and
$q \ge 0$, never the marginals. So the marginals enter the objective alone, any
projection of the fold-1 marginals is harmless for validity, and a point solved
with estimated cell weights stays feasible for the population program.

## Two things I got wrong, and the checks that caught them

**$\bar d$ must not be estimated.** I first raised $\bar d$ to the fold-1
estimate of the minimum feasible defier share. That share is a maximum of
survivor-function shortfalls — a supremum of noisy estimates — so it is biased
**upward**, the mirror image of the winner's curse in $\sup_A \Delta$. On a 30%
fold it came out at 0.045–0.052 against a full-sample 0.020, and at 0.05 the
bound is zero. The symmetry is worth stating in the paper:

| quantity | bias | effect on the bound | consequence |
|---|---|---|---|
| $\sup_A \Delta$ | up | too large | anticonservative |
| minimum feasible $\bar d$ | up | too small | conservative, powerless |

$\bar d$ is a sensitivity parameter, as in KR's `max_defiers_share`. It is taken
as given and the full-sample minimum is reported separately.

**The dual box is a pure scaling, not a regularisation.** The dual variables
$u, v$ are free in the textbook dual, so a vertex can carry huge components; I
boxed them expecting a bias–variance trade-off. There is none. Across
$\text{box} \in \{0.1, 0.25, 0.5, 1, 2\}$ the bound and its standard error scale
in exact proportion and **the $t$-ratio is invariant**:

| covariate | $k$ | box 0.1 | 0.25 | 0.5 | 1 | 2 |
|---|---|---|---|---|---|---|
| `const` | 5 | 2.52 | 2.52 | 2.52 | 2.52 | 2.52 |
| `grandmother_baseline` | 3 | 0.66 | 0.66 | 0.66 | 0.66 | 0.66 |
| `grandmother_baseline` | 4 | 0.63 | 0.63 | 0.63 | 0.63 | 0.63 |

The box stops binding at 0.5, beyond which the dual attains the primal. So
there is no free power in the choice of dual point within this family, and the
$t$-ratio is a property of the data.

## The answer

Oracle stage 1 (dual point from the full sample), $\bar d = 0.021$:

| covariate | $k$ | bound | SE | LCB | bound/SE |
|---|---|---|---|---|---|
| `grandmother_baseline` | 3 | 0.0236 | 0.0355 | 0 | 0.66 |
| `grandmother_baseline` | 4 | 0.0259 | 0.0409 | 0 | 0.63 |
| `const` (pooled) | 5 | 0.0735 | 0.0291 | **0.0256** | 2.52 |
| `grandmother_baseline` | 5 | 0.0690 | 0.0293 | **0.0208** | 2.35 |

**The identification result at $k=3$ and $k=4$ is real and not detectable at
$n = 568$.** The bound is smaller than one standard error. This is not a defect
of the estimator: the oracle stage 1 produces the bound exactly and carries the
same standard error.

**Where there is power, conditioning does not help.** At $k=5$ the bound is
comfortably significant, but the conditional bound (0.0690) is *below* the
pooled one (0.0735) — the cell-versus-population effect again.

**Sample splitting is not viable here.** With $n = 568$, $K = 5$ and two cells,
a 30% first fold has about 170 observations to estimate ten marginals and
$\hat A_k$, and the split-sample point estimates come out systematically
negative. The oracle numbers above are therefore an upper bound on what any
split-sample procedure achieves on this data.

## Where that leaves the project

The methodological results stand and are verified: the ATE reformulation, the
two-channel decomposition, the defier LP and its dual, and — the main one — that
a population-level defier cap preserves dominance while a cell-level cap does
not. The empirical section cannot carry a positive finding on Baranov at
$n = 568$ with a five-valued mediator.

Three routes, in order of how much they would change the paper:

1. **Report the power calculation as the finding.** "The conditional bound is
   twice the pooled bound at $k=3,4$; detecting it at 80% power needs
   $n \approx (z/0.65)^2 \times 568 \approx 3{,}400$." That is a useful,
   honest contribution and it costs nothing more to produce.
2. **Drop sample splitting for the applications** and report the plug-in bound
   with its bias characterised by the simulations, as KR do. Loses the
   confidence statement, which was the second contribution.
3. **Find a larger experiment** with a multivalued mediator. This is the only
   route to a positive empirical finding, and it is the question deferred to
   Jon.
