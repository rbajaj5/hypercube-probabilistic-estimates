# Knot-volume phases as projector-hypercube mixtures

## 1. Making the proposal precise

A knot, viewed as a one-dimensional curve, has zero ordinary
three-dimensional volume. For a hyperbolic knot \(K\), however, the
complement \(S^3\setminus K\) has a canonical hyperbolic volume after the
curvature is normalized to \(-1\). Denote it by

```math
V(K)=\operatorname{Vol}(S^3\setminus K).
```

This note tests the following reproducible construction:

1. choose a knot \(K_j\) uniformly from a fixed finite sample;
2. form the circle phase

   ```math
   u_j=\{V(K_j)\}\in[0,1);
   ```

3. set

   ```math
   q_j=\cos(2\pi u_j);
   ```

4. use \(q_j\) as the angle parameter of the
   [two-basis projector hypercube](ProjectorMatrixHypercubeLaw.md); and
5. average the resulting exact conditional laws over the sampled knots.

The fractional parts themselves do not sum to one and are not probability
weights. The canonical probability object is their empirical circle measure

```math
\mu_N=\frac1N\sum_{j=1}^N\delta_{u_j}.
\tag{1}
```

The projector construction pushes (1), together with a uniform sign vector,
to a second probability measure on operator norms.

## 2. Reproducible knot sample

The experiment samples \(N=16\) indices uniformly without replacement from
[SnapPy 3.3.2](https://snappy.computop.org/)'s
[`CensusKnots`](https://snappy.computop.org/censuses.html) table of 3,116
hyperbolic knot exteriors, using Python's `random.Random` with seed
`20260730`. The fixed output is stored in
[`data/census_knot_volumes_seed_20260730.csv`](../data/census_knot_volumes_seed_20260730.csv).

The generator is
[`experiments/generate_knot_volume_sample.py`](../experiments/generate_knot_volume_sample.py).
It requires the optional SnapPy package; all downstream verification uses the
fixed CSV and the Python standard library only.

| Census index | Knot | \(V(K)\) | \(u=\{V\}\) | \(q=\cos(2\pi u)\) | Angle atoms |
|---:|---|---:|---:|---:|---:|
| 138 | `K7_67` | 5.489070 | 0.489070 | -0.997643 | 0.068659, 1.998821 |
| 243 | `K8_43` | 4.885012 | 0.885012 | 0.750159 | 1.870914, 0.706881 |
| 691 | `K9_190` | 5.982569 | 0.982569 | 0.994008 | 1.997002, 0.109469 |
| 705 | `K9_204` | 6.010158 | 0.010158 | 0.997964 | 1.998982, 0.063812 |
| 716 | `K9_215` | 6.053722 | 0.053722 | 0.943570 | 1.971583, 0.335946 |
| 1041 | `K9_540` | 7.606770 | 0.606770 | -0.783289 | 0.658348, 1.888539 |
| 1351 | `K10_85` | 5.087757 | 0.087757 | 0.851794 | 1.924471, 0.544437 |
| 1634 | `K10_368` | 6.528916 | 0.528916 | -0.983541 | 0.181433, 1.991753 |
| 2186 | `K10_920` | 7.952890 | 0.952890 | 0.956511 | 1.978136, 0.294921 |
| 2209 | `K10_943` | 7.989150 | 0.989150 | 0.997677 | 1.998838, 0.068162 |
| 2303 | `K10_1037` | 8.157653 | 0.157653 | 0.548219 | 1.759670, 0.950558 |
| 2304 | `K10_1038` | 8.158449 | 0.158449 | 0.544030 | 1.757288, 0.954955 |
| 2427 | `K10_1161` | 8.349217 | 0.349217 | -0.583798 | 0.912362, 1.779774 |
| 2449 | `K10_1183` | 8.371814 | 0.371814 | -0.692811 | 0.783823, 1.840006 |
| 2867 | `K10_1601` | 8.961972 | 0.961972 | 0.971590 | 1.985744, 0.238369 |
| 3005 | `K10_1739` | 9.251127 | 0.251127 | -0.007080 | 1.409199, 1.419211 |

The last column lists

```math
\sqrt{2+2q_j},\qquad \sqrt{2-2q_j}.
```

SnapPy's [`Manifold.volume`](https://snappy.computop.org/manifold.html)
reported ten digits of volume accuracy and positively oriented tetrahedra for
every sampled exterior. These are numerical hyperbolic-volume computations,
not interval-certified volumes.

## 3. The modulo-one probability distribution

Dividing the circle into eight half-open intervals of width \(1/8\) gives
the empirical categorical law

| Bin | Interval | Count | Probability |
|---:|---|---:|---:|
| 0 | \([0,1/8)\) | 3 | \(3/16\) |
| 1 | \([1/8,2/8)\) | 2 | \(2/16\) |
| 2 | \([2/8,3/8)\) | 3 | \(3/16\) |
| 3 | \([3/8,4/8)\) | 1 | \(1/16\) |
| 4 | \([4/8,5/8)\) | 2 | \(2/16\) |
| 5 | \([5/8,6/8)\) | 0 | \(0\) |
| 6 | \([6/8,7/8)\) | 0 | \(0\) |
| 7 | \([7/8,1)\) | 5 | \(5/16\) |

Its Shannon entropy is approximately \(2.43004\) bits, compared with
\(3\) bits for the uniform eight-bin law. The first four empirical circle
Fourier magnitudes are

```math
0.339322,\quad 0.376910,\quad 0.189610,\quad 0.078726.
```

This small sample is visibly nonuniform, but it is not evidence against
equidistribution of knot volumes in any larger family. The census ordering,
the finite sampling measure, the sample size, and the choice of multiplier
\(1\) before reduction modulo one all matter.

## 4. Sixteen specific matrix hypercubes

For each phase \(u_j\), define

```math
Z=
\begin{pmatrix}1&0\\0&-1\end{pmatrix},
\qquad
Q_j=
\begin{pmatrix}
\cos(2\pi u_j)&\sin(2\pi u_j)\\
\sin(2\pi u_j)&-\cos(2\pi u_j)
\end{pmatrix},
```

and

```math
A=\frac{I_2+Z}{2},\quad
B=\frac{I_2-Z}{2},\quad
C_j=\frac{I_2-Q_j}{2},\quad
D_j=\frac{I_2+Q_j}{2}.
```

The \(j\)-th knot labels the 16-vertex matrix hypercube

```math
\mathcal H_j
=
\left\{
\varepsilon_AA+\varepsilon_BB
+\varepsilon_CC_j+\varepsilon_DD_j:
\varepsilon\in\{-1,1\}^4
\right\}.
\tag{2}
```

Choose \(J\) uniformly from the 16 knots and choose the sign vector
independently and uniformly. If \(W\) is the operator norm of the selected
vertex, then the induced probability law is

```math
\boxed{
\nu_{16}
=
\frac18\delta_0+\frac58\delta_2
+\frac1{128}\sum_{j=1}^{16}
\left(
\delta_{\sqrt{2+2q_j}}
+\delta_{\sqrt{2-2q_j}}
\right).}
\tag{3}
```

The coefficients in (3) sum to one exactly. Thus the volume phases have
produced a proper probability distribution without treating arbitrary
fractional parts as unnormalized weights. Its second moment remains

```math
\mathbb E_{\nu_{16}}W^2=3,
\tag{4}
```

while its first moment is

```math
\mathbb E_{\nu_{16}}W
\approx1.55032865537.
\tag{5}
```

## 5. The induced black-swan coefficients

Let \(R\) be an independent standard Pareto variable with exponent
\(\alpha>0\). The
[heavy-tail projector theorem](ProjectorHeavyTailMixtures.md)
applied conditionally to each knot gives, for \(t\geq2\),

```math
\mathbb P(RW>t)
=
\overline C_{\alpha,16}t^{-\alpha},
\tag{6}
```

where

```math
\overline C_{\alpha,16}
=
\frac58\,2^\alpha
+\frac1{128}\sum_{j=1}^{16}
\left(
(2+2q_j)^{\alpha/2}
+(2-2q_j)^{\alpha/2}
\right).
\tag{7}
```

For this sample:

| Pareto exponent \(\alpha\) | \(\overline C_{\alpha,16}\) | Zhang \(q=0\) | Coincident \(|q|=1\) |
|---:|---:|---:|---:|
| 0.5 | 1.137495 | 1.181185 | 1.060660 |
| 1.0 | 1.550329 | 1.603553 | 1.500000 |
| 1.5 | 2.148813 | 2.188215 | 2.121320 |
| 2.0 | 3.000000 | 3.000000 | 3.000000 |
| 3.0 | 5.904117 | 5.707107 | 6.000000 |
| 4.0 | 11.688151 | 11.000000 | 12.000000 |

The sample respects the general phase transition: its coefficient lies
between the mutually unbiased and coincident extremes, agrees with both at
\(\alpha=2\), and reverses ordering above \(2\).

## 6. Interpretation and scope

This construction validates the proposed pipeline in a precise form:

```text
random census knot
    -> hyperbolic complement volume
    -> fractional circle phase
    -> angle-labeled 4-dimensional sign cube
    -> exact operator-norm probability law
    -> exact Pareto tail coefficient.
```

What is not justified is declaring the raw values
\(\{V(K_j)\}\) themselves to be probabilities. They do not sum to one, and
normalizing them would introduce an arbitrary weighting rule. The empirical
measure (1) and pushforward (3) are canonical once the sampling measure and
phase map have been chosen.

The calculation is a reproducible numerical experiment. It does not claim:

- that hyperbolic knot volumes are equidistributed modulo one;
- that SnapPy's census order defines a canonical random-knot model;
- that the displayed decimal volumes are rigorous interval certificates; or
- that knot geometry changes the Pareto exponent. It changes the exact
  angular coefficient through (7).

## Verification

Run

```text
python verification/verify_knot_volume_hypercubes.py
```

The checker is dependency-free. It verifies the seeded census indices, CSV
metadata, modulo-one bins, Fourier diagnostics, normalization of the
256-vertex mixture, its exact second moment, and the displayed tail
coefficients.
