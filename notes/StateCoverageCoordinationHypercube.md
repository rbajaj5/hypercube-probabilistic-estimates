# State coverage, supervision, and coordination on a hypercube

## Scope and source boundary

[OvercookedV2](https://arxiv.org/abs/2503.17821) argues empirically that poor
state coverage explains much of the zero-shot-coordination gap in the original
Overcooked benchmark, then introduces tasks in which partial information and
test-time protocol formation remain essential.  Its state-augmentation
procedure samples starting states from a cross-play rollout buffer.

[AlphaZero in Sparsely Rewarded Games](https://arxiv.org/abs/2607.08984)
contrasts strong empirical play with exact oracle consistency in Connect Four
and Chomp.  It reports that an oracle-derived auxiliary policy loss improves
consistency, while high aggregate match rates still need not yield a single
perfect trace.

This note extracts exact finite probability statements suggested by that
distinction.  It does not reproduce the paper's reinforcement-learning
experiments, prove convergence of its training algorithm, or claim an
improvement to either paper.  The results below concern arbitrary probability
laws and finite action sets on a Boolean cube.

## The sharp state-coverage certificate

Let \(\Omega=\{0,1\}^d\), let \(P\) be a training-state law, and let \(Q\) be a
test-state law.  For every statistic \(f:\Omega\to[0,1]\),

\[
\left|\mathbb E_Q f-\mathbb E_P f\right|
\leq
\operatorname{TV}(Q,P)
:=
\frac12\sum_{x\in\Omega}|Q(x)-P(x)|.
\tag{1}
\]

The constant is exact:

\[
\sup_{0\leq f\leq1}
\left|\mathbb E_Q f-\mathbb E_P f\right|
=\operatorname{TV}(Q,P).
\tag{2}
\]

Indeed, (1) follows by separating the positive and negative parts of \(Q-P\),
and equality is attained by the indicator of
\(\{x:Q(x)>P(x)\}\).  For \(f\in[a,b]\), rescaling gives the sharp bound
\((b-a)\operatorname{TV}(Q,P)\).

This is a coverage certificate, not a policy-learning theorem: it controls
the shift of a *fixed bounded state statistic*.

## What distribution augmentation guarantees

Model an augmented training law by

\[
R_\alpha=(1-\alpha)P+\alpha A,\qquad 0\leq\alpha\leq1,
\tag{3}
\]

where \(A\) is the law represented by an augmentation buffer.  Convexity of
the \(\ell_1\) norm gives

\[
\operatorname{TV}(Q,R_\alpha)
\leq
(1-\alpha)\operatorname{TV}(Q,P)
+\alpha\operatorname{TV}(Q,A).
\tag{4}
\]

If the augmentation law is exactly the target law, \(A=Q\), then

\[
\operatorname{TV}(Q,R_\alpha)
=(1-\alpha)\operatorname{TV}(Q,P).
\tag{5}
\]

Thus ideal target-state augmentation contracts every bound in (1) by the
exact factor \(1-\alpha\).  An approximate buffer still helps whenever the
right side of (4) is smaller than \(\operatorname{TV}(Q,P)\).  Equation (4)
does not say that an optimizer will learn a good action at the newly covered
states.

## Exact finite-buffer missed mass

Let \(X_1,\ldots,X_N\) be independent samples from \(A\).  The test mass not
represented in the buffer is

\[
U_N
=
\sum_{x\in\Omega}
Q(x)\mathbf 1\{x\notin\{X_1,\ldots,X_N\}\}.
\tag{6}
\]

Its expectation is exactly

\[
\mathbb E U_N
=
\sum_{x\in\Omega}Q(x)(1-A(x))^N
\leq
\sum_{x\in\Omega}Q(x)e^{-NA(x)}.
\tag{7}
\]

Consequently,

\[
\Pr(U_N\geq\varepsilon)
\leq
\frac{\mathbb E U_N}{\varepsilon}.
\tag{8}
\]

For \(A=Q\) uniform on \(M=2^d\) states,

\[
\mathbb E U_N
=
(1-M^{-1})^N
\leq e^{-N/M}.
\tag{9}
\]

The convenient sufficient condition
\(N\geq M\log(1/\varepsilon)\) therefore makes the expected missed mass at
most \(\varepsilon\).  Covering *every* state is more expensive:

\[
\Pr(\text{some state is unseen})
\leq
M(1-M^{-1})^N
\leq Me^{-N/M}.
\tag{10}
\]

Hence

\[
N\geq
2^d\bigl(d\log2+\log(1/\delta)\bigr)
\tag{11}
\]

is sufficient for full coverage with probability at least \(1-\delta\).
The exponential dependence is not removed by calling the samples a replay
buffer.

There is a similarly transparent partner-imbalance cost.  Suppose partner
type \(h\) has a disjoint support of \(s_h\) states and is sampled with weight
\(w_h\), uniformly within that support.  Its expected missed test mass is

\[
\left(1-\frac{w_h}{s_h}\right)^N.
\tag{12}
\]

Rare partner types therefore cost on the scale \(s_h/w_h\), even when the
combined buffer is large.

## Sparse supervision has an exponential blind spot

The same calculation is a lower bound for learning from sparse feedback.
Suppose a uniformly sampled cube state is informative only when it belongs to
\(G\subseteq\Omega\), where \(|G|=k\).  After \(N\) independent samples,

\[
\Pr(\text{no informative state is observed})
=
\left(1-\frac{k}{2^d}\right)^N.
\tag{13}
\]

Consider two equally likely label or reward rules that agree outside \(G\)
but require opposite answers on \(G\).  Conditional on observing no point of
\(G\), their data distributions are identical, so any learner has conditional
error at least \(1/2\) in distinguishing them.  Its minimax identification
error is therefore at least

\[
\frac12\left(1-\frac{k}{2^d}\right)^N.
\tag{14}
\]

For a single informative vertex, the natural sample scale is \(2^d\).  Search,
structure, or auxiliary labels may change the sampling law or reveal labels
away from \(G\), but sparse terminal feedback alone does not remove (14).

## Bit-flip augmentation is a Walsh low-pass filter

Another natural augmentation starts from \(X\sim P\) and independently flips
each bit with probability \(\eta\in[0,1/2]\).  Write \(A=P*T_\eta\) for the
law of the perturbed state and

\[
\chi_S(x)=(-1)^{\sum_{i\in S}x_i},
\qquad
\widehat P(S)=\mathbb E_P\chi_S(X).
\]

Independence of the noise bits gives the exact multiplier

\[
\widehat A(S)
=(1-2\eta)^{|S|}\widehat P(S).
\tag{15}
\]

If \(V\) is uniform on the cube, Parseval's identity yields

\[
\chi^2(A\|V)
=
\sum_{\varnothing\ne S\subseteq[d]}
(1-2\eta)^{2|S|}\widehat P(S)^2,
\tag{16}
\]

and Cauchy--Schwarz gives

\[
\operatorname{TV}(A,V)
\leq \frac12\sqrt{\chi^2(A\|V)}.
\tag{17}
\]

Noise augmentation suppresses high-degree state structure particularly
quickly.  At \(\eta=1/2\), it produces the uniform law exactly.  This is a
coverage-smoothing statement; it need not move \(P\) toward a nonuniform
partner law \(Q\).

## What an oracle auxiliary target certifies

Let a state have \(m\) legal actions, of which the set \(B\) of size \(b\) is
oracle-optimal.  The smoothed target used in the AlphaZero auxiliary-loss
paper has weights proportional to

\[
\bigl(\mathbf 1\{a\in B\}+\varepsilon\bigr)
\mathbf 1\{a\text{ is legal}\}.
\]

After normalization,

\[
q_\varepsilon(B)
=
\frac{b(1+\varepsilon)}{b+m\varepsilon}
=1-\beta,
\qquad
\beta=\frac{(m-b)\varepsilon}{b+m\varepsilon}.
\tag{18}
\]

Thus smoothing itself assigns the exact mass \(\beta\) to nonoptimal actions.
If a learned policy \(p\) satisfies
\(\operatorname{TV}(p,q_\varepsilon)\leq\tau\), then

\[
p(B)\geq1-\beta-\tau.
\tag{19}
\]

Equivalently, if the excess cross-entropy is
\(\operatorname{KL}(q_\varepsilon\|p)\leq\Delta\), Pinsker's inequality gives
\(\tau\leq\sqrt{\Delta/2}\).  This certifies the success probability when an
action is *sampled* from \(p\); it is not automatically a certificate for the
argmax action.

There is a sharp elementary certificate for the unsmoothed target.  If \(q\)
is uniform on \(B\), its auxiliary cross-entropy is

\[
L_B(p)=-\frac1b\sum_{a\in B}\log p(a).
\]

The arithmetic--geometric mean inequality gives

\[
p(B)\geq b e^{-L_B(p)}
=e^{-(L_B(p)-\log b)}.
\tag{20}
\]

Equality holds when the probabilities of the optimal actions are equal.  A
simple sufficient condition for the greedy action to be optimal is

\[
p(B)>\frac{b}{b+1},
\tag{21}
\]

because then some optimal action has mass greater than \(1-p(B)\), the total
mass of all nonoptimal actions.

## Match rate does not certify a perfect trace

Let \(F_t\) be the event of an oracle mistake at labeled ply \(t\), and put
\(e_t=\Pr(F_t)\).  Without any independence assumption,

\[
\Pr(\text{the length-\(T\) trace is perfect})
\geq
1-\sum_{t=1}^T e_t.
\tag{22}
\]

The bound is sharp when the mistake events are disjoint.  More starkly, a
dataset in which every \(T\)-ply trace contains exactly one mistake has pooled
oracle-match rate \(1-1/T\) but perfect-trace rate zero.  Under the additional
assumption of independent mistakes with common match probability \(r\), the
perfect-trace probability is \(r^T\).

These identities explain why move-level match, first failure, longest
consistent chain, and perfect-trace rate measure genuinely different things.
They do not assert that mistakes made by a learned game-playing agent are
independent.

## Perfect coverage need not produce coordination

Coverage and protocol information can be separated exactly.  Let the observed
state \(X\) have any full-support distribution on \(\Omega\).  Independently,
draw a hidden protocol

\[
\Theta\sim\operatorname{Unif}(\{0,1\}^b).
\]

The correct action is \(\Theta\).  Any policy whose action is based only on
\(X\), deterministic or randomized, has success probability

\[
\Pr(\text{action}=\Theta)=2^{-b}.
\tag{23}
\]

Here train and test state laws may be identical, so their total variation is
zero and every observed state may have been visited.  Nevertheless, the
coordination failure probability is \(1-2^{-b}\).  If the agent instead
receives a grounded signal \(Y=\Theta\), the decoder \(Y\mapsto Y\) succeeds
with probability one.

The button-game example in the OvercookedV2 paper illustrates the same design
principle.  If a lamp index is encoded as

\[
L=2a+\theta,
\tag{24}
\]

then the decoder \(L\bmod2=\theta\) is invariant to the otherwise arbitrary
button index \(a\).  The protocol is grounded in a shared feature, so no
convention about \(a\) must be learned.

Equations (1)--(17) quantify state coverage and sparse feedback.
Equations (18)--(22) separate move-level supervision from trajectory-level
perfection, while (23)--(24) show why coverage alone cannot certify zero-shot
coordination.

## What is verified

`verification/verify_state_coverage_coordination.py` uses exact rational
arithmetic and exhaustive finite enumeration to check:

1. sharpness of the total-variation observable bound;
2. augmentation convexity and the exact target-mixture contraction;
3. the missed-mass formula by enumerating all small sample buffers;
4. the uniform coupon bound and the disjoint partner-mixture formula;
5. the sparse-signal blind-spot probability;
6. every Walsh multiplier, the chi-squared Parseval identity, and (17);
7. target-smoothing mass, the exact polynomial form of (20), and the
   trajectory counterexample;
8. the hidden-protocol success law and the parity decoder.

The checker proves the finite identities it implements.  It does not evaluate
learned policies or establish that a particular augmentation distribution
matches an empirical deployment distribution.

## Relation to the other modules

- [Latent environments](LatentEnvironmentHypercube.md) quantify dependence
  caused by hidden variables.  The protocol example here instead makes the
  hidden variable decision-relevant and proves a success lower ceiling.
- [Hypercube walks and matrix mixing](HypercubeWalkMatrixMixing.md) analyze a
  dynamical refresh chain.  Equation (13) is the one-step product-noise
  analogue.
- [Noisy reciprocity](NoisyReciprocityHypercube.md) studies coordination when
  actions affect future states.  The present results deliberately isolate
  static coverage before any strategic dynamics are introduced.
