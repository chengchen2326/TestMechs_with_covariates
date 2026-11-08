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
