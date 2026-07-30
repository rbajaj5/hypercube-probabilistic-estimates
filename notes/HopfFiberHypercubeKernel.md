# Hopf-fiber phase and a projective hypercube kernel

## Status and source boundary

This note extracts two standard geometric ingredients emphasized in Jennifer
Lorraine Nielsen's evolving preprint,
[*The Complex Hopf Fibration as the Canonical Space for Gauge-Gravity
Unification*](https://philpapers.org/versions/NIETTU):

1. Fourier orthogonality along the \(U(1)\) fiber of
   \(S^1\to S^{2d-1}\to\mathbb{CP}^{d-1}\);
2. the Fubini--Study transition kernel
   \(|\langle\psi,\phi\rangle|^2\) on projective space.

The archive had 289 versions when accessed on 2026-07-30.  The current
PhilArchive file was inspected as version 289; the stable formula references
used below are also visible in the open, non-peer-reviewed
[version 3](https://doi.org/10.20944/preprints202604.0315.v3), posted
2026-06-22.  This note does **not** import or endorse the preprint's
gauge-gravity, particle-spectrum, or uniqueness claims.  It uses only the
displayed fiber Fourier and projective-overlap identities, then proves the
Boolean-cube consequences independently.

## 1. Projectivizing a Boolean vertex

For \(x=(x_1,\ldots,x_d)\in\{-1,1\}^d\), set

```math
\psi_x=\frac1{\sqrt d}(x_1,\ldots,x_d)\in S^{d-1}
\subset S^{2d-1}\subset\mathbb C^d
```

and let \([\psi_x]\in\mathbb{CP}^{d-1}\) be its projective class.  The
Fubini--Study transition kernel becomes

```math
\boxed{
K_d(x,y)
=|\langle\psi_x,\psi_y\rangle|^2
=\left(\frac1d\sum_{i=1}^d x_i y_i\right)^2
=\left(1-\frac{2h(x,y)}d\right)^2,}
\tag{1}
```

where \(h(x,y)\) is Hamming distance.

This is a genuine projective quantity: multiplying either lift by a global
fiber phase does nothing,

```math
K_d(e^{i\theta}\psi_x,e^{i\varphi}\psi_y)=K_d(x,y).
\tag{2}
```

Consequently, attaching an independent global \(U(1)\) phase to a Hex board
cannot by itself improve a winner estimate.  Relative phases or a new
observable would be required.

## 2. The kernel is exactly degree two

Expanding the square gives

```math
\boxed{
K_d(x,y)
=\frac1d+\frac{2}{d^2}
\sum_{1\le i<j\le d}
\chi_{\{i,j\}}(x)\chi_{\{i,j\}}(y).}
\tag{3}
```

Thus the projective kernel contains only the constant and Walsh level-two
sectors.  It is not a mysterious new entropy: on the embedded cube it is a
completely explicit low-degree hypercube observable.

If \(X,Y\) are coupled coordinatewise with
\(\mathbb E[X_iY_i]=\rho\), then

```math
\boxed{
\mathbb E K_d(X,Y)
=\frac1d+\frac{d-1}{d}\rho^2
=\rho^2+\frac{1-\rho^2}{d}.}
\tag{4}
```

Equation (4) is an exact finite-dimensional correction to the naive
large-\(d\) answer \(\rho^2\).

## 3. Exact coordinate-refresh relaxation

Let \(X_t\) be the coordinate-refresh walk started at \(X_0=x\): at each
step, choose one coordinate uniformly and replace it by an independent fair
sign.  Walsh level two has eigenvalue \(1-2/d\).  Equations (3) and the walk
spectrum therefore give

```math
\boxed{
\mathbb E_x K_d(X_0,X_t)
=\frac1d+\frac{d-1}{d}
\left(1-\frac2d\right)^t.}
\tag{5}
```

The projective overlap relaxes on the \(d/2\) scale of a degree-two
observable, not on the \(d\log d\) scale required to mix the entire cube.

## 4. Fiber Fourier selection

For a uniform fiber phase \(\Theta\),

```math
\mathbb E e^{i(k-\ell)\Theta}=\mathbf1_{\{k=\ell\}},
\tag{6}
```

and a trilinear fiber average obeys

```math
\mathbb E e^{i(k+\ell+m)\Theta}
=\mathbf1_{\{k+\ell+m=0\}}.
\tag{7}
```

These are the rigorous contents of the fiber-mode decoupling and
charge-selection rule that are reusable here.  They are circle Fourier
orthogonality, not a new probability theorem.

## 5. Projectively weighted Hex stability

Let \(G:\{-1,1\}^d\to\{-1,1\}\) be a signed Boolean observable, such as the
tie-free Hex winner, and define

```math
C_G(\rho)
=\mathbb E_\rho[
G(X)G(Y)K_d(X,Y)].
\tag{8}
```

Write

```math
\mu_k=4^{-d}W_k(G),
\qquad \sum_k\mu_k=1,
\tag{9}
```

for its normalized Walsh level masses.  Multiplication by
\(\chi_{\{i,j\}}\) shifts a Walsh set \(S\) to
\(S\mathbin{\triangle}\{i,j\}\).  Averaging over pairs yields the exact
formula

```math
\boxed{
\begin{aligned}
C_G(\rho)
=\sum_{k=0}^d\mu_k\Bigg[
&\frac{\rho^k}{d}\\
&+\frac{2}{d^2}\left\{
\binom{k}{2}\rho^{k-2}
+k(d-k)\rho^k
+\binom{d-k}{2}\rho^{k+2}
\right\}
\Bigg].
\end{aligned}}
\tag{10}
```

Terms with zero binomial coefficient are omitted, so (10) is well defined at
\(\rho=0\).  This gives a new diagnostic at no simulation cost: the exact
\(1\times1\), \(2\times2\), and \(3\times3\) Hex spectra already committed
in [the Hex note](HexWinnerNoiseMixing.md) determine their complete
projectively weighted noise curves.

Substitution of those tables gives

```math
C_{G_1}(\rho)=\rho,
\tag{10a}
```

```math
C_{G_2}(\rho)
=\frac{\rho^4+10\rho^3+6\rho^2+14\rho+1}{32},
\tag{10b}
```

and

```math
\begin{aligned}
C_{G_3}(\rho)=\frac1{1327104}\big(&
73\rho^9+556\rho^8+5696\rho^7+24044\rho^6\\
&+116194\rho^5+215724\rho^4+587320\rho^3\\
&+156244\rho^2+213469\rho+7784\big).
\end{aligned}
\tag{10c}
```

Each polynomial equals \(1\) at \(\rho=1\).  The positive constant terms for
\(n=2,3\) show that even independent boards retain correlation after
weighting by projective/Hamming proximity; this is a property of the
weighted diagnostic, not residual dependence of the unweighted winners.

Equation (10) does not improve the unweighted Hex crossing probability or
total-influence bound.  Its gain is more precise: it converts the paper's
projective geometry into a rigorously computable, Hamming-sensitive Hex
correlation and makes the global-phase limitation explicit.

## 6. Finite Gaussian determinant on the cube

The preprint also motivates sectorwise Gaussian determinants.  On the finite
cube no zeta regularization is needed.  For \(L=I-P\), where \(P\) is the
coordinate-refresh operator, the nonzero eigenvalues are \(k/d\) with
multiplicity \(\binom dk\).  Hence

```math
\boxed{
\det{}'L
=\prod_{k=1}^d\left(\frac{k}{d}\right)^{\binom dk},
\qquad
Z_d=(\det{}'L)^{-1/2}
=\prod_{k=1}^d\left(\frac{d}{k}\right)^{\binom dk/2}.}
\tag{11}
```

This is the exact finite analogue of a Gaussian determinant factor.  It is a
normalization identity for the hypercube Laplacian, not evidence for the
preprint's physical determinant interpretation.

## Verification

Run

```text
python verification/verify_hopf_fiber_hypercube.py
```

The checker exhaustively verifies the Hamming/projective identity and its
Walsh decomposition, the product-noise and refresh-walk laws, fiber Fourier
selection, the finite determinant multiplicities, and formula (10) against
direct coupled enumeration for the exact Hex winner through \(3\times3\).
