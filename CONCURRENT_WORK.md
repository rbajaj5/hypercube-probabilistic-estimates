# Concurrent work from July 2026

Several short manuscripts appeared within hours of one another after the
Dirichlet-calibration result of
[Vlassis and Thomas](https://arxiv.org/abs/2607.08415). This note separates
their overlapping formulas, their different scopes, and the proof dependency
that they share.

Public timestamps this close do not establish discovery priority. The table
is only a chronology of public versions.

| Work | Public timestamp | Result used here |
|---|---:|---|
| [Fu--Han--Wang--Yan--Zhang--Zhou](https://arxiv.org/abs/2607.23980) | arXiv submission: 27 July 2026, 04:06 UTC | Bound for every \(n,\delta\), sharp for \(\delta\geq1\) |
| [Mark Stander](https://doi.org/10.5281/zenodo.21622951) | Zenodo record: 27 July 2026, 10:35 UTC; manuscript dated 26 July | Conditional derivation of the shared small-slack bound and the sharp unit-slack constant |
| [Zipei Nie--Jiaye Wei](https://arxiv.org/abs/2607.24528) | arXiv submission: 27 July 2026, 15:07 UTC | Shared small-slack bound, plus the exact \(n=2\) result from \(K_2^{\mathrm{ad}}\) |

## 1. The common small-slack formula

For \(0<\delta\leq1\), all three works obtain

```math
B_{n,\delta}
:=
\delta\left(\frac{n}{n+\delta}\right)^n.
```

Their common proof architecture is:

1. use the Vlassis--Thomas theorem to treat the Dirichlet statistic
   \(K_n\) as a valid merger;
2. bound \(K_n\) on a simplex halfspace using the centroid or off-centroid
   inequalities of Grunbaum and Letwin--Yaskin; and
3. translate the pointwise geometric bound into a probability bound.

Stander labels the conclusion *conditional* because the Vlassis--Thomas
theorem was then a new, unrefereed preprint and is not reproved in his note.
The other two papers cite the same theorem as an input. This is a difference
in presentation, not a different mathematical dependency.

## 2. Where the papers differ

For \(\delta>1\), Fu--Han--Wang--Yan--Zhang--Zhou prove the sharper
\(\delta\)-dependent estimate

```math
\left(1-\frac1{n+\delta}\right)^n,
```

whereas Stander's stated every-slack corollary retains the unit-slack value

```math
\left(\frac{n}{n+1}\right)^n.
```

The former is strictly larger when \(\delta>1\), since

```math
\frac{n+\delta-1}{n+\delta}
>
\frac{n}{n+1}
\quad\Longleftrightarrow\quad
\delta>1.
```

Nie--Wei focus on \(0<\delta\leq1\), and at \(n=2\) use the admissible merger
\(K_2^{\mathrm{ad}}\) to obtain the exact all-slack value

```math
\min\left\{
\frac{\delta}{1+\delta},
\left(\frac{1+\delta}{2+\delta}\right)^2
\right\}.
```

They explicitly note that this bivariate theorem is already known.

## 3. What the low-dimensional refinement adds

For \(1\leq n\leq4\), the classical Samuels theorem and the
finite-dimensional inequality in Paulin's proof give the exact value

```math
C_{n,\delta}
:=
\min\left\{
\frac{\delta}{1+\delta},
\left(1-\frac1{n+\delta}\right)^n
\right\}.
```

Consequently:

- for \(n=2\), this is exactly the bivariate theorem recorded by Nie--Wei;
- for \(n=3,4\) and \(0<\delta<1\), it is strictly larger than the common
  July-2026 bound \(B_{n,\delta}\);
- for \(n\leq4\), its proof does **not** depend on the
  Vlassis--Thomas calibration theorem; and
- it makes no claim about the unresolved general \(n\geq5\) small-slack
  problem.

Thus the value of this small-deviation module is a precise low-dimensional
comparison and an independent classical proof path, not a claim of a new
extremal theorem.
