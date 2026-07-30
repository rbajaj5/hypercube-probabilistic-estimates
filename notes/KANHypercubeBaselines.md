# KAN baselines on the Boolean hypercube

## 1. Why the KAN survey is relevant, and where it is not

Hou, Ji, Zhang, and Stefanidis,
[*Kolmogorov--Arnold Networks: A Critical Assessment of Claims,
Performance, and Practical Viability*](https://arxiv.org/abs/2407.11075)
(arXiv:2407.11075v8, 2025), emphasize two evaluation principles that are
directly useful here:

1. compare models at a controlled parameter or computation budget; and
2. check whether the basis functions match the structure of the data.

The paper is a critical survey, not a source of a new hypercube estimate.  The
results below are elementary finite-cube consequences of applying those two
principles to binary inputs.  In particular, they do not rely on the survey's
compiled benchmark numbers.

The second suggested paper, Critch and Tsimerman,
[*A Taxonomy of Omnicidal Futures Involving Artificial
Intelligence*](https://arxiv.org/abs/2507.09369), supplies a qualitative
decision tree and narrative scenarios, not event probabilities or calibration
data.  Assigning probabilities to its five leaves would therefore be an
unsupported modeling choice.  It is not used below.

## 2. Binary collapse of a scalar one-layer KAN

Write the Boolean cube as \(\{-1,1\}^d\).  A scalar-output, one-layer KAN has
the additive form

\[
 A(x)=c+\sum_{i=1}^d \phi_i(x_i).
 \tag{1}
\]

The edge functions may be splines, wavelets, polynomials, or arbitrary
functions.  On a binary coordinate, however,

\[
 \phi_i(x_i)
 =
 \frac{\phi_i(1)+\phi_i(-1)}2
 +
 \frac{\phi_i(1)-\phi_i(-1)}2\,x_i.
 \tag{2}
\]

Consequently, the complete function class in (1) is exactly

\[
 V_{\leq 1}
 =
 \operatorname{span}\{1,x_1,\ldots,x_d\}.
 \tag{3}
\]

Adding more spline knots to the edge functions does not enlarge this class.
If an edge uses \(p\) basis coefficients, its values on cube data are obtained
through a \(2\times p\) evaluation matrix.  Its rank is at most \(2\), so at
least \(p-2\) coefficient directions are unidentifiable from the binary
training data alone.  Regularization or off-cube evaluation can distinguish
those directions, but the observed cube function cannot.

This statement is only about a single additive layer with scalar output.
Later nonlinear KAN layers can create coordinate interactions.

## 3. The exact optimal error certificate

For \(S\subseteq[d]\), let

\[
 \chi_S(x)=\prod_{i\in S}x_i,
 \qquad
 f(x)=\sum_{S\subseteq[d]}\widehat f(S)\chi_S(x),
 \]

where expectations and Fourier coefficients use the uniform law on the cube.
Orthogonal projection gives the best possible scalar one-layer KAN:

\[
 A_\star(x)
 =
 \widehat f(\varnothing)
 +
 \sum_{i=1}^d \widehat f(\{i\})x_i.
 \tag{4}
\]

Its exact squared error is

\[
 \inf_A \mathbb E[(f(X)-A(X))^2]
 =
 \sum_{|S|\geq2}\widehat f(S)^2.
 \tag{5}
\]

Thus no optimizer, initialization, spline grid, or GPU benchmark is needed to
decide whether this architecture can fit a finite Boolean observable.  For a
balanced \(\{-1,1\}\)-valued observable, its captured variance is precisely
the level-one Walsh mass.

There is also an exact probability-estimation consequence.  Couple \(X,Y\)
coordinatewise so that

\[
 \mathbb E[X_iY_i]=\rho,\qquad 0\leq\rho\leq1.
\]

Then

\[
 \operatorname{Stab}_\rho(f)
 =
 \mathbb E[f(X)f(Y)]
 =
 \sum_{k=0}^d \rho^k W_k,
 \qquad
 W_k=\sum_{|S|=k}\widehat f(S)^2.
 \tag{6}
\]

The optimal one-layer estimate and its certified remainder are

\[
 W_0+\rho W_1
 \leq
 \operatorname{Stab}_\rho(f)
 \leq
 W_0+\rho W_1+\rho^2 E_1,
 \qquad
 E_1=\sum_{k\geq2}W_k.
 \tag{7}
\]

For a Boolean winner, the disagreement probability is
\((1-\operatorname{Stab}_\rho(f))/2\), so (7) immediately gives a two-sided
probability interval.

## 4. The cube-native upgrade: retain interactions

The natural controlled-complexity replacement is not a denser univariate
spline grid.  It is the degree-\(r\) Walsh projection

\[
 P_{\leq r}f
 =
 \sum_{|S|\leq r}\widehat f(S)\chi_S,
 \qquad
 E_r=\|f-P_{\leq r}f\|_2^2=\sum_{k>r}W_k.
 \tag{8}
\]

This is the best \(L^2\) approximation using interactions of order at most
\(r\).  Truncating (6) gives the rigorous estimate

\[
 0
 \leq
 \operatorname{Stab}_\rho(f)
 -
 \sum_{k=0}^r \rho^kW_k
 \leq
 \rho^{r+1}E_r.
 \tag{9}
\]

The same certificate applies to the coordinate-refresh walk from
[the mixing note](HypercubeWalkMatrixMixing.md).  If one coordinate is chosen
uniformly and replaced by a fresh fair sign at every step, then, for
\(0\leq r<d\),

\[
 0
 \leq
 \mathbb E[f(X_0)f(X_t)]
 -
 \sum_{k=0}^r W_k\left(1-\frac{k}{d}\right)^t
 \leq
 E_r\left(1-\frac{r+1}{d}\right)^t.
 \tag{10}
\]

Equations (9) and (10) improve monotonically with \(r\), expose the exact
cost of neglected interactions, and use the same Walsh masses already
computed for the Hex experiments.

## 5. Exact small-Hex diagnosis

Let \(G_n\in\{-1,1\}\) be the tie-free winner on a fully colored
\(n\times n\) Hex board, as defined in
[the Hex note](HexWinnerNoiseMixing.md).  It is balanced, so \(W_0=0\).
Exact enumeration gives:

| board | \(d\) | one-layer mass \(W_1\) | \(E_1\) | \(E_2\) | \(E_3\) |
|---|---:|---:|---:|---:|---:|
| \(1\times1\) | 1 | \(1\) | \(0\) | \(0\) | \(0\) |
| \(2\times2\) | 4 | \(5/8\) | \(3/8\) | \(1/8\) | \(0\) |
| \(3\times3\) | 9 | \(7765/16384\) | \(8619/16384\) | \(4727/16384\) | \(1503/16384\) |

So the scalar one-layer KAN is exact only for the trivial board.  On the
\(3\times3\) board it misses more than half the winner's variance:

\[
 E_1=\frac{8619}{16384}\approx0.52606.
\]

Allowing pair and triple Walsh interactions reduces the certified error to
\(4727/16384\) and \(1503/16384\), respectively.  This is a genuine
improvement of the estimates in (9)--(10), but it is an exact spectral
improvement rather than evidence that KAN training itself is superior.

For example, the complete \(2\times2\) product-noise curve is

\[
 \operatorname{Stab}_\rho(G_2)
 =
 \frac58\rho+\frac14\rho^2+\frac18\rho^3.
 \tag{11}
\]

The additive KAN baseline keeps only \(5\rho/8\); adding pair interactions
keeps the first two terms and leaves the certified remainder \(\rho^3/8\).

## 6. A case where a deeper KAN is exactly aligned

Suppose \(f\) is invariant under every coordinate permutation.  It then
depends only on

\[
 s(x)=x_1+\cdots+x_d\in\{-d,-d+2,\ldots,d\}.
\]

Choose a continuous piecewise-linear function \(\psi\) that interpolates the
\(d+1\) required values.  Then

\[
 f(x)=\psi\!\left(\sum_i x_i\right)
\tag{12}
\]

is an exact width-one, two-layer KAN representation.  This includes majority
and parity.  It also illustrates the paper's data--architecture alignment
principle: the model succeeds because the target has a one-dimensional
sufficient statistic.

Hex is not invariant under arbitrary cell permutations; adjacency and the
two distinguished crossing directions matter.  Hamming weight alone is
therefore not sufficient for the Hex winner.

## 7. What should be benchmarked next

For larger boards, a fair empirical comparison should report:

1. test \(L^2\) error and winner-disagreement error;
2. parameter count, forward FLOPs, elapsed time, and peak memory;
3. several fixed seeds with uncertainty intervals;
4. the exact \(W_1\) projection as the one-layer lower bound;
5. degree-\(r\) Walsh models at matched parameter budgets; and
6. residual noise-stability error against (9), rather than prediction error
   alone.

The exact small-board calculations already reject the hypothesis that extra
univariate spline resolution repairs missing Hex interactions.  A deeper KAN
may still be worth testing, but it must beat the cube-native Walsh baselines
under matched resource budgets.

## 8. Verification

Run

```text
python verification/verify_kan_hypercube_baselines.py
```

The checker uses exact integer and rational arithmetic.  It exhausts every
Boolean observable through dimension \(3\), verifies the projection and
noise/refresh-remainder identities, checks exact two-layer interpolation for
every symmetric Boolean observable through dimension \(6\), and independently
enumerates the displayed Hex spectra through \(3\times3\).
