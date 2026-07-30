# Hex winners as hypercube observables

## Status and source

A fully colored \(n\times n\) Hex board is a point of a Boolean hypercube.
The topological no-tie theorem makes its winner a total sign-valued
observable, and the hypercube walk diagonalizes its relaxation by Walsh
level.

The Hex topology input is Theorems 1.2.6 and 1.2.8 in Anna Karlin and Yuval
Peres, [*Game Theory, Alive*](https://bookstore.ams.org/MBK/101).  This note
uses the theorem rather than reproducing its path-of-arrows or reduction-to-Y
proof.

## 1. A tie-free Boolean winner

Give the \(N=n^2\) cells signs

```math
x_v=
\begin{cases}
+1,&v\text{ is blue},\\
-1,&v\text{ is yellow}.
\end{cases}
```

Let \(B_n(x)\) indicate a blue left-to-right crossing, and let \(Y_n(x)\)
indicate a yellow top-to-bottom crossing.  Hex topology says

```math
\boxed{B_n(x)+Y_n(x)=1}
\tag{1}
```

for every completed standard board.  Define the signed winner

```math
G_n(x)=B_n(x)-Y_n(x)=2B_n(x)-1.
\tag{2}
```

Then

```math
\boxed{G_n(x)\in\{-1,1\},\qquad G_n(x)^2=1.}
\tag{3}
```

The game-theoretic strategy-stealing theorem concerns optimal alternating
play.  Equations (1)--(3) concern uniformly random completed boards; neither
statement implies the other.

## 2. The exact symmetry is twisted color reversal

On a symmetric rhombus, let \(J\) transpose the two board coordinates and
swap the two colors.  Then

```math
\boxed{G_n(Jx)=-G_n(x).}
\tag{4}
```

Because \(J\) preserves the uniform cube law, (4) gives

```math
\mathbb E G_n=0,\qquad
\Pr(G_n=1)=\Pr(G_n=-1)=\frac12.
\tag{5}
```

It is important that \(J\) includes the coordinate transpose.  In general,

```math
G_n(-x)\ne-G_n(x).
```

Thus the raw winner can have even Walsh levels.  They already occur for the
\(2\times2\) board.

## 3. Full winner spectrum

Use the unnormalized Walsh transform

```math
H_G(S)=\sum_{x\in\{-1,1\}^N}G_n(x)\chi_S(x)
```

and level masses

```math
W_k(G_n)=\sum_{|S|=k}H_G(S)^2.
\tag{6}
```

Parseval and (3) give the exact normalization

```math
\boxed{
\sum_{k=0}^N W_k(G_n)=4^N,\qquad W_0(G_n)=0.}
\tag{7}
```

Exact enumeration gives:

| Board | Nonzero full-winner masses |
| --- | --- |
| \(1\times1\) | \(W_1=4\) |
| \(2\times2\) | \(W_1=160,\ W_2=64,\ W_3=32\) |
| \(3\times3\) | \(W_1=124240,\ W_2=62272,\ W_3=51584,\ W_4=16320,\ W_5=5920,\ W_6=1216,\ W_7=512,\ W_8=64,\ W_9=16\) |

The \(W_2\) entries explicitly rule out the tempting but false inference
that the no-tie theorem makes the winner complement-odd.

## 4. Exact winner mixing

For the discrete coordinate-refresh walk of
[the walk note](HypercubeWalkMatrixMixing.md), stationarity and Walsh
diagonalization give

```math
\boxed{
\mathbb E_\pi[G_n(X_0)G_n(X_t)]
=
4^{-N}\sum_{k=1}^N
W_k(G_n)\left(1-\frac{k}{N}\right)^t.}
\tag{8}
```

Since \(G_n\) is sign-valued,

```math
\boxed{
\Pr_\pi(G_n(X_t)\ne G_n(X_0))
=
\frac12\left(
1-\mathbb E_\pi[G_n(X_0)G_n(X_t)]
\right).}
\tag{9}
```

For the Poissonized refresh chain with total update rate one, each coordinate
retains correlation

```math
\rho_t=e^{-t/N},
```

so

```math
\boxed{
\mathbb E_\pi[G_n(X_0)G_n(X_t)]
=4^{-N}\mathcal N_{G_n}(e^{-t/N}),}
\tag{10}
```

where

```math
\mathcal N_{G_n}(\rho)=\sum_kW_k(G_n)\rho^k.
```

Equations (8)--(10) turn the full spectral table into an exact
winner-mixing curve.  If a sequence of larger Hex boards moves its
normalized spectral mass toward high levels, the winner decorrelates before
the whole \(N\)-cube mixes.  Establishing that asymptotic movement is a
noise-sensitivity theorem, not a consequence of Hex topology alone.

The reflection-positive gluing argument in
[the picture-language note](ReflectionPositiveHypercubePictures.md) adds
Hankel positivity, complete monotonicity, and log-convexity constraints to
this curve.  In particular, the total influence alone gives certified
two-sided correlation and disagreement bounds even when the complete
spectrum is unavailable.

### Total influence and exact anchors

If \(X^{\oplus v}\) is obtained by flipping cell \(v\), the total influence
is

```math
\begin{aligned}
\operatorname{Inf}(G_n)
&=\sum_v\Pr\!\left(G_n(X)\ne G_n(X^{\oplus v})\right)\\
&=4^{-N}\sum_{k=1}^N kW_k(G_n).
\end{aligned}
\tag{10a}
```

It is also the left derivative at \(\rho=1\) of the product-noise
stability.  Exhaustive enumeration gives

| Board | Exact total influence |
| --- | ---: |
| \(1\times1\) | \(1\) |
| \(2\times2\) | \(3/2\) |
| \(3\times3\) | \(249/128\) |
| \(4\times4\) | \(2405/1024\) |

These values are useful anchors for checking a large-board estimator without
pretending that a finite simulation proves an asymptotic exponent.

### Seeded larger-board experiment

For a uniformly random cell \(V\), the estimator

```math
\widehat{\operatorname{Inf}}(G_n)
=N\,\frac1M\sum_{j=1}^M
\mathbf 1\!\left\{
G_n(X_j)\ne G_n(X_j^{\oplus V_j})
\right\}
\tag{10b}
```

is unbiased.  The experiment also couples \(Y\) to \(X\) by independently
flipping every cell with probability \((1-\rho)/2\) and estimates
\(\Pr(G_n(X)\ne G_n(Y))\).  It uses \(M=65{,}536\), seed `20260730`, board
sides \(5,7,9,11,15,21,31\), and ten values of \(\rho\) from \(0\) to \(1\).
The batched flood-fill kernel is exhaustively compared with an independent
scalar breadth-first search on every board through \(3\times3\).

The fixed run produced:

| \(n\) | \(\widehat{\operatorname{Inf}}(G_n)\) | one SE | disagreement at \(\rho=.9\) | disagreement at \(\rho=.99\) |
| ---: | ---: | ---: | ---: | ---: |
| 5 | 2.6684 | 0.0302 | 0.1192 | 0.0126 |
| 7 | 3.4648 | 0.0491 | 0.1407 | 0.0166 |
| 9 | 4.0811 | 0.0692 | 0.1583 | 0.0193 |
| 11 | 4.8410 | 0.0926 | 0.1756 | 0.0231 |
| 15 | 5.7129 | 0.1383 | 0.2030 | 0.0275 |
| 21 | 7.1127 | 0.2170 | 0.2347 | 0.0350 |
| 31 | 9.8687 | 0.3785 | 0.2719 | 0.0440 |

A least-squares fit of
\(\log\widehat{\operatorname{Inf}}(G_n)\) on \(\log n\) over only these
seven sizes has slope \(0.6971\).  This is a descriptive finite-size
diagnostic.  The monotone increase in fixed-\(\rho\) winner disagreement is
evidence of increasing noise sensitivity over this range, but neither the
fit nor the table establishes a limit law.

The source and fixed output are

```text
experiments/simulate_hex_noise_gpu.py
data/hex_noise_seed_20260730.csv
```

The script automatically uses CUDA when PyTorch has CUDA support.  The fixed
dataset records the actual execution device and library versions; the
machine used for the committed run had a visible GPU but only a CPU PyTorch
build, so the validated CPU fallback was used.  Reported standard errors are
marginal binomial errors.  Because the same base boards are reused across
the \(\rho\) values, differences between curve points are correlated.

## 5. The complement-odd projection

Define

```math
A_n(x)
=\frac{G_n(x)-G_n(-x)}2
=B_n(x)-B_n(-x).
\tag{11}
```

Then

```math
A_n(-x)=-A_n(x),
\tag{12}
```

so every even Walsh coefficient vanishes.  The odd masses of \(A_n\) are
exactly the odd masses displayed for \(G_n\), while all even masses are
deleted.

This is the signed crossing contrast used in the supercongruence
repository's
[planar-noise Walsh note](https://github.com/rbajaj5/a183068-supercongruence/blob/main/related-results/BlackNoiseWalshCongruences.md).
It is the correct input for the odd-chaos adjacent-scale congruence.  The raw
winner \(G_n\) is the correct input for tie-free winner mixing.  They answer
different questions.

## 6. The Y reduction is a majority renormalization

The reduction-to-Y proof in the cited source has an additional probabilistic
interpretation.  It maps a colored Y board of side \(n\) to one of side
\(n-1\) by replacing each target cell with the majority color of a
three-cell triangle.  The proof shows that this deterministic
coarse-graining preserves the existence and color of the unique Y.

In sign coordinates, the local gate is

```math
\boxed{
\operatorname{Maj}_3(x,y,z)
=\frac{x+y+z-xyz}{2}.}
\tag{13}
```

It has Walsh mass \(3/4\) at level one and \(1/4\) at level three.  If two
input triples are coupled coordinatewise with correlation \(\rho\), the
output correlation is exactly

```math
\boxed{
\operatorname{Stab}_\rho(\operatorname{Maj}_3)
=\frac{3\rho+\rho^3}{4}.}
\tag{14}
```

For \(0<\rho<1\), this is smaller than \(\rho\): a local majority gate
amplifies noise under this coupling.

Unlike fixed-boundary Hex, the two colors in Y have the same three target
sides.  The signed Y winner therefore satisfies raw complement oddness

```math
G_n^Y(-x)=-G_n^Y(x).
\tag{15}
```

The majority reduction commutes with color reversal, so it preserves this
oddness at every scale.  This supplies a direct topological-renormalization
source for an odd-only Walsh spectrum.

Neighboring triangles in the actual reduction overlap.  Their majority
outputs are consequently correlated even when the original cells are
independent.  Equation (14) is an exact local noise law, but iterating it as
though the entire reduced board remained a product measure would be
incorrect.  The correct global object is the pushforward measure under the
overlapping majority circuit.

## 7. What the Hex and Y excerpts add

- The winner is defined on every completed board and has no zero state.
- Its total Walsh mass is exactly \(4^N\).
- Symmetry makes the two winners equiprobable under uniform coloring.
- The full winner-mixing curve follows from its entire spectrum.
- Antisymmetrization, not topology by itself, removes even chaos.
- The Y reduction provides a winner-preserving majority
  coarse-graining and a genuinely complement-odd winner.

## Verification

Run

```text
python verification/verify_hex_winner_noise.py
python verification/verify_hex_noise_simulation.py
```

The checker exhaustively verifies existence and uniqueness through
\(4\times4\), the twisted color-reversal symmetry, the \(1\times1\) through
\(3\times3\) spectral tables, exact total influence through \(4\times4\),
the odd projection, and the discrete refresh-walk correlation formula.  It
also verifies the majority polynomial and its exact noise-stability map.
The separate data checker recomputes every probability, standard error,
correlation, and influence entry from the integer counts and audits the
symmetry and endpoint controls.
