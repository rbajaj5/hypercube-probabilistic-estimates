# Noisy reciprocity on a two-bit hypercube

## 1. Source and scope

The Veritasium video
[*What The Prisoner's Dilemma Reveals About Life, The Universe, and
Everything*](https://www.veritasium.com/videos/2024/1/15/what-the-prisoners-dilemma-reveals-about-life-the-universe-and-everything)
surveys Axelrod's tournaments, Tit-for-Tat, cooperation, and the damage caused
by implementation noise.  Its source list includes Wu and Axelrod,
[*How to Cope with Noise in the Iterated Prisoner's
Dilemma*](https://deepblue.lib.umich.edu/items/544c6210-b867-4e9d-8268-56fd4cd12a6d),
which compares generosity, contrition, and Pavlov under noisy play.

The video is motivation, not a theorem saying that one strategy is optimal in
every population.  Tournament rankings depend on the opponents, payoff
matrix, horizon, mutation rule, and noise model.

This note extracts an exact hypercube calculation.  The outcome of one round
is one of

\[
 CC,\ CD,\ DC,\ DD,
\]

which is the two-dimensional cube \(\{-1,1\}^2\), with cooperation encoded by
\(+1\).  A pair of memory-one strategies induces a four-state Markov chain.
Implementation errors make its stationary probabilities, payoffs, and mixing
rates exactly computable.

## 2. General memory-one transition law

For player \(X\), let

\[
 p=(p_{CC},p_{CD},p_{DC},p_{DD})
\tag{1}
\]

be the intended probability of cooperation after the previous outcome,
written from \(X\)'s viewpoint.  Define \(q\) similarly for player \(Y\).
Independent implementation error \(\varepsilon\) flips each intended action,
so the executed cooperation probability is

\[
 \widetilde p_s
 =
 \varepsilon+(1-2\varepsilon)p_s.
 \tag{2}
\]

From \(Y\)'s viewpoint the middle two states are reversed.  Thus, in the
state order \(CC,CD,DC,DD\), put

\[
 \widetilde q^{\,X}
 =
 (\widetilde q_{CC},\widetilde q_{DC},
   \widetilde q_{CD},\widetilde q_{DD}).
\]

The transition row from state \(s\) is

\[
 M_s=
 \left(
 \widetilde p_s\widetilde q^{\,X}_s,\quad
 \widetilde p_s(1-\widetilde q^{\,X}_s),\
 (1-\widetilde p_s)\widetilde q^{\,X}_s,\
 (1-\widetilde p_s)(1-\widetilde q^{\,X}_s)
 \right).
 \tag{3}
\]

When \(0<\varepsilon<1/2\), every entry is positive, so \(M\) has a unique
stationary law \(\pi\).  For payoff vectors

\[
 r_X=(R,S,T,P),\qquad r_Y=(R,T,S,P),
\]

the exact long-run payoffs are

\[
 u_X=\pi r_X,\qquad u_Y=\pi r_Y.
 \tag{4}
\]

For a finite horizon \(h\), initial law \(\nu\), and no discounting, the
expected average payoff is

\[
 \frac1h\,
 \nu\left(I+M+\cdots+M^{h-1}\right)r_X.
 \tag{5}
\]

Equations (3)--(5) replace a noisy tournament estimate by exact rational
arithmetic whenever the strategy and error probabilities are rational.

## 3. Tit-for-Tat and generosity

A symmetric reactive strategy cooperates with executed probability \(a\)
after the opponent cooperated and \(b\) after the opponent defected.  In spin
coordinates \(x,y\in\{-1,1\}\), define

\[
 \alpha=a+b-1,\qquad \beta=a-b.
\]

The Markov operator satisfies

\[
 Mx=\alpha+\beta y,\qquad
 My=\alpha+\beta x,
 \tag{6}
\]

and

\[
 M(xy)=\alpha^2+\alpha\beta(x+y)+\beta^2xy.
 \tag{7}
\]

In the ordered basis \(1,x+y,xy,x-y\), the operator is triangular.  Its
eigenvalues are

\[
 \boxed{1,\ \beta,\ \beta^2,\ -\beta.}
 \tag{8}
\]

The negative eigenvalue is the exact retaliatory oscillation: after one
player defects, the asymmetric \(CD/DC\) mode changes sign on the next round.

For Generous Tit-for-Tat with generosity \(g\), the intended response is

\[
 (p_C,p_D)=(1,g).
\]

Writing \(\delta=1-2\varepsilon\), its executed parameters are

\[
 a=1-\varepsilon,\qquad
 b=\varepsilon+\delta g,
\]

so

\[
 \alpha=\delta g,\qquad
 \boxed{\beta=\delta(1-g).}
 \tag{9}
\]

The unique stationary law is a product law.  Each player cooperates with
probability

\[
 \boxed{
 c=
 \frac{\varepsilon+\delta g}
 {2\varepsilon+\delta g}.}
 \tag{10}
\]

The centered symmetric mode decays as \(\beta^t\), while the retaliation mode
decays as \((-\beta)^t\).  Generosity therefore has two exact effects in this
model:

1. it raises the stationary cooperation rate; and
2. it shortens the error cascade from \(\delta^t\) to
   \([\delta(1-g)]^t\).

For ordinary Tit-for-Tat, \(g=0\).  For every
\(0<\varepsilon<1/2\),

\[
 c=\frac12,\qquad
 \pi=(1/4,1/4,1/4,1/4).
 \tag{11}
\]

Thus arbitrarily rare errors eventually destroy mutual cooperation in the
infinite-horizon stationary limit, and the slow mode approaches \(1\) as
\(\varepsilon\downarrow0\).  This does not contradict successful finite
tournaments: the order of the small-noise and long-time limits matters.

## 4. Relation to reflection positivity

Let \(A=x-y\) be the antisymmetric retaliation observable.  In stationary
symmetric reactive play,

\[
 \frac{\mathbb E[A_0A_t]}{\mathbb E[A_0^2]}
 =
 (-\beta)^t.
 \tag{12}
\]

The negative odd correlations show why the one-step strategy chain does not
satisfy the positive-spectrum hypothesis in
[the reflection-positive note](ReflectionPositiveHypercubePictures.md).
However \(M^2\) has eigenvalues

\[
 1,\ \beta^2,\ \beta^2,\ \beta^4\geq0.
\]

Consequently the even-lag sequence is reflection positive and obeys all the
Hankel, complete-monotonicity, and moment constraints proved there.

This supplies a useful perturbation test: negative Hankel or odd-lag
correlations need not indicate bad data; they may identify an authentic
period-two strategic mode.  Subsampling every two rounds removes that sign
oscillation.

## 5. Exact payoff example

Use the common payoff convention

\[
 (R,S,T,P)=(3,0,5,1)
\]

and implementation error \(\varepsilon=1/20\).  Exact self-play values are:

| strategy | stationary cooperation | stationary payoff |
|---|---:|---:|
| Always Cooperate | \(19/20\) | \(1179/400\) |
| Always Defect | \(1/20\) | \(459/400\) |
| Tit-for-Tat | \(1/2\) | \(9/4\) |
| Generous Tit-for-Tat, \(g=1/3\) | \(7/8\) | \(183/64\) |
| Pavlov | \(1729/2000\) | \(5553/2000\) |

For symmetric reactive self-play, (10) gives the product stationary law and

\[
 u(c)
 =
 Rc^2+(S+T)c(1-c)+P(1-c)^2
 =
 1+3c-c^2.
 \tag{13}
\]

Self-play is not a tournament ranking.  At the same error rate, Always Defect
scores \(1881/400\) against Always Cooperate, while Always Cooperate scores
only \(81/400\) in the reverse role.  Against Generous Tit-for-Tat with
\(g=1/3\), Always Defect scores \(2451/1000\), while the generous player
scores \(801/1000\).  These exact values expose the cooperation/exploitation
tradeoff that a self-play table hides.

## 6. What scales to larger cubes

With \(d\) simultaneous players, one round of cooperate/defect actions lies
in \(\{-1,1\}^d\).  The two-player calculation does **not** imply that an
arbitrary network game is Walsh diagonal: game-dependent responses couple
coordinates and can introduce negative or complex modes.

The reusable procedure is:

1. specify the strategy state, payoff rule, update schedule, and error model;
2. construct the finite Markov operator;
3. compute its stationary law and signed spectrum;
4. separate positive relaxation modes from oscillatory modes; and
5. compare simulations with exact low-dimensional anchors.

Spatial evolutionary games on the Hex adjacency graph are a plausible next
experiment, but their conclusions will depend on the reproduction/update
rule.  The Veritasium video alone does not select that rule or prove an Ising,
percolation, or evolutionary scaling law.

## 7. Verification

Run

```text
python verification/verify_noisy_reciprocity_hypercube.py
```

The checker uses exact rational arithmetic.  It verifies the general
memory-one transition and stationary equations, the complete reactive
spectrum, stationary cooperation and retaliation correlations, the displayed
payoff table and asymmetric examples, finite-horizon matrix sums, and
reflection positivity of the two-step chain.
