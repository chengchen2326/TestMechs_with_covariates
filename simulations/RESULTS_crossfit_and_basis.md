# Cross-fitting and the fixed-basis question

Run 2026-11-09. Code: `R/lb_frac_affected_x_cf.R`,
`simulations/crossfit_mc.R`, `simulations/c2_basis.R`.

## C2 — Jon's question, answered

*With a fixed basis $b(y)$, is stage 1 a linear program in $\beta(x)$?* Yes, and
in a way that subsumes KR's discretisation.

1. **The relaxation from sets to functions is exact.** $\sup_A \Delta_k(A)$ over
   subsets equals the supremum over functions $a : \mathcal{Y} \to [0,1]$,
   because the objective is linear in $a$ and the extreme points of
   $[0,1]^{\mathcal{Y}}$ are the indicators. Nothing is lost.
2. **With a partition-of-unity basis the LP has a closed form.** Writing
   $a(y) = \beta'b(y)$, the constraint $0 \le \beta'b(y) \le 1$ collapses to
   $\beta \in [0,1]^J$, a box, so the optimum is
   $\beta_j = \mathbf 1\{\text{coefficient}_j > 0\}$.
3. **Indicator bins are such a basis**, and the closed form is then exactly the
   discretisation KR use. So *choice of bins* is a special case of *choice of
   basis*, and smoother bases are available at the same computational cost.

Verified against brute-force enumeration of all $2^Q$ subsets for $Q \le 10$ and
against the LP over $[0,1]^Q$: **zero discrepancy across 800 draws.**

**A caveat KR do not state.** Only *nested* refinements are ordered. The nested
chain of bin widths $1 \subset 2 \subset 4 \subset 12 \subset 24$ gives a
decreasing $0.429, 0.288, 0.202, 0.045, 0$; but widths 2 and 3 do not nest and
give $0.288$ against $0.290$. More bins is not always better — it depends on
where the boundaries fall relative to the sign pattern of $p_1 - p_0$.

**Where this leads.** With covariates, $\beta$ becomes $\beta(x)$. Letting each
cell have its own $\beta$ separates the problem into one box-LP per cell;
modelling $\beta(x) = \Gamma' c(x)$ for a covariate basis $c$ makes stage 1 one
joint LP in $\Gamma$, which shares strength across cells and is a way to
condition without binning the covariates at all. That is worth trying.

## Cross-fitting — closed, negatively, and the reason matters

**The point estimate is valid with nothing needing to converge.** Each fold
selects for the other; each half-estimate is conditionally unbiased for
$\eta_k^{\hat A_f}$; any fixed $A$ gives a valid lower bound; and the average of
two valid lower bounds is a valid lower bound.

**The variance is the problem**, because the two halves are dependent. Two
standard errors were compared: the Cauchy–Schwarz bound
$(\mathrm{sd}_a + \mathrm{sd}_b)/2$, which never relies on the nuisance
settling down, and the usual zero-covariance one.

### The usual standard error does not hold size

Sign calibration, 600 replications:

| $t$ | $n$ | empirical sd of $\hat\eta_{cf}$ | conservative SE | zero-cov SE |
|---|---|---|---|---|
| 0 | 1000 | 0.02937 | 0.03259 | **0.02368** |
| 0 | 4000 | 0.01370 | 0.01629 | **0.01182** |
| 1 | 1000 | 0.03391 | 0.03443 | **0.02453** |
| 1 | 4000 | 0.01405 | 0.01774 | **0.01256** |

The truth sits between the two: the conservative bound is conservative, and the
zero-covariance standard error is **too small by 10–25%**. That shows directly
in size and coverage at $t = 0$:

| estimator | size ($n = 1000, 4000$) | coverage |
|---|---|---|
| `cf_cons` | 0.020, 0.012 | 0.980, 0.988 |
| `cf_pool` | **0.078, 0.067** | **0.922, 0.933** |

So the covariance between the two halves is positive and not negligible,
exactly as the absence of a convergence guarantee for $\hat A_k$ warns.

### With a valid standard error, cross-fitting is a wash

| $t$ | $n$ | single split | `cf_cons` |
|---|---|---|---|
| 0.5 | 1000 | **0.138** | 0.087 |
| 0.5 | 4000 | **0.453** | 0.373 |
| 1.0 | 1000 | 0.410 | 0.392 |
| 1.0 | 4000 | 0.947 | **0.972** |

I expected the conservative standard error to equal the single split's. It does
not, and the reason is simple: the single split uses `split_frac = 0.3`, so
stage 2 gets 70% of the data, while cross-fitting is forced to 50/50 and each
half gets 50%. The cross-fitted **point estimate is genuinely better** — 0.0458
against 0.0369 at $t = 1$, $n = 1000$, on a bound of 0.060, because both halves
select on 50% rather than 30% — but the larger standard error eats the gain.

**The deeper reason is visible in the table above.** The empirical dispersion
(0.0339) is almost exactly the Cauchy–Schwarz bound (0.0344): the inequality is
nearly tight, so the two halves are strongly positively correlated. They are the
same data with the roles exchanged, so there is very little variance to
recover.

### What this settles

Cross-fitting does not close the gap to KR's test. **The power we lose to
sample splitting is structural, not a tuning failure**, and the right procedure
is a single split at `split_frac` around 0.3. An open item is now a result, and
a negative result here is more useful than an unresolved one: it tells a reader
not to try this.

## A note on method

Four mechanism hypotheses were proposed today and all four were refuted by
checks rather than by reading output: that the multivalued gain would come from
the complier channel $G_B$; that zeroing non-monotone cells was the right fix;
that the dual box was a bias–variance trade-off; and that cross-fitting would
recover the lost power. Each was stated as a testable prediction first. The
checks that caught them — comparison against KR's own function, against an
earlier run, and against brute-force enumeration — belong in the replication
package, not only in the working process.

---

# The standard error is conditional on $\hat A$, and that is the point

Two things were flagged as problems today and both were my own measurement
errors. Recording them because the resolution is a remark the paper needs.

## False alarm 1 — the Design C size anomaly

Design C at 400 replications showed size 0.075 and 0.083 in the `signflip`
configuration at $n = 4000$, against a nominal 0.05. I called 3 standard errors
"not noise" without accounting for having scanned 24 cells, where roughly one
excursion past two standard errors is expected by chance. At 2000 replications
the same cells give 0.0485 and 0.0530, and all twelve `signflip` cells lie in
0.0475–0.0580.

## False alarm 2 — the apparent variance bug

An audit comparing the reported standard error to the empirical dispersion of
$\hat\eta$ found the ratio deteriorating with $n$ in the `homog` configuration
at an interior mediator value: 0.965, 0.828, 0.741 at $n = 1000, 4000, 16000$.
A bootstrap agreed, putting the influence-function standard error at 0.57 of
the bootstrap one at $k = 3$.

Both comparisons were mis-specified. The influence function estimates the
standard deviation **conditional on $\hat A$**, while the Monte Carlo and the
bootstrap both re-estimate $\hat A$ on every draw. When $\hat A$ moves, the
estimand $\eta_k^{\hat A}$ moves with it, so the extra dispersion is not
estimation error that needs covering.

Holding fold 1 fixed and redrawing only fold 2 settles it:

| $k$ | $n_2$ | IF SE | conditional empirical sd | ratio | conditional coverage |
|---|---|---|---|---|---|
| 3 | 2800 | 0.013520 | 0.013777 | **0.981** | 1.000 |
| 3 | 11200 | 0.004938 | 0.004901 | **1.008** | 1.000 |
| 5 | 2800 | 0.015432 | 0.015661 | 0.985 | 0.935 |
| 5 | 11200 | 0.007516 | 0.007442 | 1.010 | 0.950 |

## The remark the paper needs

A reader will make the same objection: the reported standard error is smaller
than the sampling dispersion of the estimate. It is supposed to be. The
argument is:

1. Conditional on $\hat A$, the lower confidence bound covers
   $\eta_k^{\hat A}$ with probability at least $1 - \alpha$.
2. $\eta_k^{\hat A} \le \eta_k$ for **every** $A$, by dual feasibility — this
   is the same weak-duality step that licenses the whole procedure.
3. By iterated expectations the bound covers $\eta_k$ unconditionally.

So the variation in $\hat A$ never has to be paid for. That is precisely what
sample splitting buys, and without the remark a reader will read the gap
between the standard error and the dispersion as a defect rather than as the
mechanism working.
