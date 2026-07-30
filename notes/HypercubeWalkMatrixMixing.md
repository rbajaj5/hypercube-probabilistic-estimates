# Random walks on hypercubes and almost-orthogonal projector frames

## Status and boundary

This note adds a Markov-chain layer to the repository's static hypercube
laws.  The cube spectrum and mixing estimates are classical.  The matrix
application is an exact Walsh-degree calculation combined with the beta
overlap law for random directions.

The walk is an auxiliary sampler on parity or sign vectors.  It is not the
Collatz map: under the Bernstein--Lagarias coding, Collatz acts as a shift on
an infinite parity sequence, whereas the chain below refreshes a randomly
chosen coordinate of a fixed finite cube.

The
[GOE determinant-factor repository](https://github.com/rbajaj5/goe-determinant-factor-density)
is used for its radial-angular decomposition
\(B\sim\operatorname{Beta}(1/2,m)\).  Its “orthogonality” note primarily
concerns orthogonal polynomials; it does not itself claim that GOE matrices
are nearly orthogonal.

## 1. Choose an aperiodic cube walk

On

```math
\Omega_d=\{-1,1\}^d,
```

the ordinary walk that flips one uniformly selected coordinate has Walsh
eigenvalues

```math
1-\frac{2k}{d},\qquad 0\leq k\leq d.
\tag{1}
```

The top-degree character has eigenvalue \(-1\), so this walk has period two
and does not converge to the uniform law.

Use the coordinate-refresh walk instead:

1. choose \(I\) uniformly from \(\{1,\ldots,d\}\);
2. replace \(\epsilon_I\) by an independent fair sign; and
3. leave the other coordinates fixed.

This is the same kernel as the half-lazy coordinate-flip walk.  Its
stationary law \(\pi\) is uniform on \(\Omega_d\).

## 2. Exact spectrum and a mixing bound

For \(S\subseteq[d]\), let

```math
\chi_S(\epsilon)=\prod_{j\in S}\epsilon_j.
```

### Theorem 1 (Walsh diagonalization)

The refresh kernel \(P\) satisfies

```math
\boxed{
P\chi_S=
\left(1-\frac{|S|}{d}\right)\chi_S.}
\tag{2}
```

Indeed, refreshing a coordinate outside \(S\) preserves the character,
whereas refreshing one inside \(S\) gives conditional expectation zero.

Starting from any vertex \(x\), Parseval gives the exact chi-squared
distance

```math
\boxed{
\chi^2(P^t(x,\cdot)\Vert\pi)
=
\sum_{k=1}^d
\binom dk\left(1-\frac{k}{d}\right)^{2t}.}
\tag{3}
```

Consequently,

```math
\begin{aligned}
\|P^t(x,\cdot)-\pi\|_{\mathrm{TV}}
&\leq
\frac12
\left[
\sum_{k=1}^d
\binom dk\left(1-\frac{k}{d}\right)^{2t}
\right]^{1/2}\\
&\leq
\frac12
\left[
\exp\!\left(de^{-2t/d}\right)-1
\right]^{1/2}.
\end{aligned}
\tag{4}
```

At

```math
t=\frac d2(\log d+c),
```

the last bound becomes

```math
\frac12\sqrt{\exp(e^{-c})-1}.
\tag{5}
```

Thus the full cube mixes on the classical \(d\log d\) scale, with the
half-lazy normalization fixing the leading upper-bound scale at
\((d/2)\log d\).

## 3. Pushforwards can only mix faster

Let \(F:\Omega_d\to\mathcal Y\) be any observable: an affine Collatz
correction, an operator norm, a determinant, a knot-labeled matrix, or an
entropy bin.  Total variation contracts under a measurable map:

```math
\boxed{
\|P^t(x,\cdot)\circ F^{-1}-\pi\circ F^{-1}\|_{\mathrm{TV}}
\leq
\|P^t(x,\cdot)-\pi\|_{\mathrm{TV}}.}
\tag{6}
```

The bound can be very loose because \(F\) may identify many cube vertices.
Walsh degree gives a sharper observable-specific rate.  If

```math
F=\sum_S\widehat F(S)\chi_S,
```

then

```math
\mathbb E_xF(\epsilon_t)-\mathbb E_\pi F
=
\sum_{S\ne\varnothing}
\widehat F(S)\chi_S(x)
\left(1-\frac{|S|}{d}\right)^t.
\tag{7}
```

An observable supported at Walsh level \(k\) relaxes at the exact rate
\((1-k/d)^t\), often much faster than the whole chain.

## 4. Matrix-valued sign cubes

Let \(u_1,\ldots,u_d\) be unit vectors in \(\mathbb R^D\), put

```math
P_i=u_iu_i^\mathsf T,
\qquad
S(\epsilon)=\sum_{i=1}^d\epsilon_iP_i,
\tag{8}
```

and define the squared overlaps

```math
B_{ij}
=\operatorname{tr}(P_iP_j)
=|\langle u_i,u_j\rangle|^2.
\tag{9}
```

The Frobenius energy has the exact Boolean expansion

```math
\boxed{
\operatorname{tr}S(\epsilon)^2
=d+2\sum_{i<j}B_{ij}\epsilon_i\epsilon_j.}
\tag{10}
```

Its centered part is pure Walsh degree two.  Therefore

```math
\boxed{
\mathbb E_x[\operatorname{tr}S(\epsilon_t)^2-d]
=
\left(1-\frac2d\right)^t
\bigl(\operatorname{tr}S(x)^2-d\bigr),}
\tag{11}
```

while at stationarity

```math
\boxed{
\operatorname{Var}_\pi(\operatorname{tr}S^2)
=4\sum_{i<j}B_{ij}^2.}
\tag{12}
```

This is the useful matrix-mixing refinement: the quadratic energy relaxes
on the \(d\) scale even though the entire sign vector needs \(d\log d\)
steps.

## 5. What “almost orthogonal matrices” means here

Rank-one projectors are not orthogonal matrices.  The precise statement is
that they can be almost orthogonal as vectors in the Hilbert--Schmidt space.
Their Gram matrix is

```math
H_{ij}=\operatorname{tr}(P_iP_j)=B_{ij}.
\tag{13}
```

If

```math
\max_{i\ne j}|\langle u_i,u_j\rangle|\leq\eta,
```

then Gershgorin's theorem gives

```math
\boxed{
\|H-I_d\|_{\mathrm{op}}\leq(d-1)\eta^2.}
\tag{14}
```

Thus \((d-1)\eta^2\ll1\) is one checkable almost-orthogonality regime.

For the repository's two-basis \(2\times2\) frame, the normalized
Hilbert--Schmidt Gram matrix of the traceless involutions \(Z,Q\) is

```math
\begin{pmatrix}1&q\\q&1\end{pmatrix}.
```

Its orthogonality defect is exactly \(|q|\), and its condition number is

```math
\frac{1+|q|}{1-|q|}
\tag{15}
```

when \(|q|<1\).  The mutually unbiased case \(q=0\) is exactly orthogonal in
this matrix geometry.

There is also an exact vertex statement.  For \(|q|<1\), six of the sixteen
vertices of the two-basis projector cube are nonzero scalar multiples of
orthogonal matrices:

- two are \(\pm2I_2\);
- four satisfy \(S^2=(2\pm2q)I_2\).

The remaining eight nonzero vertices are singular, so a claim that every
vertex is “almost orthogonal” would be false.

The same matrix-valued sign cube also admits an
[exact discrete resolvent loop
equation](DiscreteLoopEquationsProjectorCubes.md).  Whereas (10)--(12)
measure only the degree-two energy, that identity retains every repeated
projector overlap through rank-one cavity resolvents.

## 6. The GOE beta law supplies random-frame estimates

The GOE factor repository writes a standard Gaussian direction
\(\omega\in S^{2m}\) as a radial-angular pair and proves

```math
B=\omega_0^2\sim\operatorname{Beta}\!\left(\frac12,m\right).
\tag{16}
```

Equivalently, if \(u,v\) are independent Haar unit vectors in dimension

```math
D=2m+1,
```

then

```math
B=|\langle u,v\rangle|^2
```

has the same law.  Its first two moments are

```math
\boxed{
\mathbb EB=\frac1D,\qquad
\mathbb EB^2=\frac{3}{D(D+2)}.}
\tag{17}
```

For \(d\) independent Haar directions, averaging (12) over the random frame
therefore gives

```math
\boxed{
\mathbb E_{\mathrm{frame}}
\operatorname{Var}_\pi(\operatorname{tr}S^2)
=
\frac{6d(d-1)}{D(D+2)}.}
\tag{18}
```

This is small when \(d\ll D\).  Chebyshev yields the joint frame/sign bound

```math
\Pr\!\left(
|\operatorname{tr}S^2-d|\geq a
\right)
\leq
\frac{6d(d-1)}{D(D+2)a^2}.
\tag{19}
```

For pairwise coherence, the elementary beta-moment bound

```math
\Pr(|\langle u,v\rangle|\geq\eta)
\leq\frac1{D\eta^2}
\tag{20}
```

and a union bound give

```math
\Pr\!\left(
\max_{i<j}|\langle u_i,u_j\rangle|\geq\eta
\right)
\leq
\frac{\binom d2}{D\eta^2}.
\tag{21}
```

The exact beta tail can replace (20) when sharper constants matter.

## 7. Collatz and finite-field mixing remain separate chains

Three notions of mixing now coexist:

1. the coordinate-refresh sampler on the static parity cube, with exact
   spectrum (2);
2. the Bernoulli shift that codes Haar-random \(2\)-adic Collatz dynamics;
3. the time-inhomogeneous affine correction
   \(C_{j+1}=3^{\epsilon_j}C_j+\epsilon_j2^j\bmod p\).

The third chain converges exponentially to uniform for \(p\geq5\), after
periodic blocking, while \(p=3\) has the exact nonuniform limiting profile
proved in [the Collatz parity note](CollatzParityHypercube.md).  These chains
can be compared through observables, but their mixing times are not
interchangeable.

A nonuniform heat-bath variant is developed in
[the novelty-refresh note](NoveltyRefreshHypercube.md).  Centering each
coordinate at its Bernoulli novelty probability preserves the exact
subset-spectrum formula while changing the invariant entropy and the
escape time from the all-zero state.

## Verification

Run

```text
python verification/verify_hypercube_walk_matrix_mixing.py
```

The checker verifies the Walsh spectrum, exact chi-squared formula,
pushforward contraction, degree-two matrix relaxation, conformal-vertex
count, and GOE beta moments using only the Python standard library.
