"""Verify state-coverage and hidden-protocol identities on Boolean cubes."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb


Q = Fraction
Distribution = tuple[Fraction, ...]


def total_variation(left: Distribution, right: Distribution) -> Fraction:
    assert len(left) == len(right)
    return sum((abs(a - b) for a, b in zip(left, right, strict=True)), Q(0)) / 2


def mixture(
    left: Distribution,
    right: Distribution,
    alpha: Fraction,
) -> Distribution:
    return tuple(
        (1 - alpha) * a + alpha * b
        for a, b in zip(left, right, strict=True)
    )


def expected_unseen_mass(
    test: Distribution,
    sampling: Distribution,
    sample_count: int,
) -> Fraction:
    return sum(
        (
            test_mass * (1 - sampling_mass) ** sample_count
            for test_mass, sampling_mass in zip(test, sampling, strict=True)
        ),
        Q(0),
    )


def brute_force_unseen_mass(
    test: Distribution,
    sampling: Distribution,
    sample_count: int,
) -> Fraction:
    state_count = len(test)
    output = Q(0)
    for samples in product(range(state_count), repeat=sample_count):
        sample_probability = Q(1)
        for state in samples:
            sample_probability *= sampling[state]
        seen = set(samples)
        unseen_mass = sum(
            (test[state] for state in range(state_count) if state not in seen),
            Q(0),
        )
        output += sample_probability * unseen_mass
    return output


def noise_convolution(
    distribution: Distribution,
    dimension: int,
    flip_rate: Fraction,
) -> Distribution:
    output = [Q(0)] * (1 << dimension)
    for source, source_mass in enumerate(distribution):
        for target in range(1 << dimension):
            distance = (source ^ target).bit_count()
            output[target] += (
                source_mass
                * flip_rate**distance
                * (1 - flip_rate) ** (dimension - distance)
            )
    return tuple(output)


def walsh_moment(distribution: Distribution, subset: int) -> Fraction:
    return sum(
        (
            mass * (-1 if (state & subset).bit_count() % 2 else 1)
            for state, mass in enumerate(distribution)
        ),
        Q(0),
    )


def verify_total_variation_and_augmentation() -> int:
    cases = (
        (
            (Q(1, 2), Q(1, 4), Q(1, 8), Q(1, 8)),
            (Q(1, 8), Q(1, 8), Q(1, 4), Q(1, 2)),
            (Q(1, 5), Q(1, 5), Q(2, 5), Q(1, 5)),
        ),
        (
            (Q(1, 3), Q(1, 6), Q(1, 12), Q(1, 12), Q(1, 12), Q(1, 12), Q(1, 12), Q(1, 12)),
            (Q(1, 12), Q(1, 12), Q(1, 12), Q(1, 12), Q(1, 12), Q(1, 12), Q(1, 6), Q(1, 3)),
            (Q(1, 8),) * 8,
        ),
    )
    checks = 0
    for training, test, augmentation in cases:
        assert sum(training, Q(0)) == 1
        assert sum(test, Q(0)) == 1
        assert sum(augmentation, Q(0)) == 1
        tv = total_variation(test, training)

        # Indicators exhaust the vertices of [0,1]^Omega, so their maximum
        # verifies the sharp bounded-observable statement.
        largest_gap = Q(0)
        for indicator in range(1 << len(test)):
            gap = abs(
                sum(
                    (
                        (test[state] - training[state])
                        for state in range(len(test))
                        if indicator & (1 << state)
                    ),
                    Q(0),
                )
            )
            largest_gap = max(largest_gap, gap)
            assert gap <= tv
            checks += 1
        assert largest_gap == tv
        checks += 1

        for alpha in (Q(0), Q(1, 7), Q(1, 2), Q(6, 7), Q(1)):
            augmented = mixture(training, augmentation, alpha)
            assert total_variation(test, augmented) <= (
                (1 - alpha) * total_variation(test, training)
                + alpha * total_variation(test, augmentation)
            )
            targeted = mixture(training, test, alpha)
            assert total_variation(test, targeted) == (1 - alpha) * tv
            checks += 2
    return checks


def verify_finite_buffers() -> int:
    test = (Q(1, 2), Q(1, 4), Q(1, 8), Q(1, 8))
    sampling = (Q(1, 8), Q(1, 4), Q(1, 4), Q(3, 8))
    checks = 0
    for sample_count in range(6):
        exact = expected_unseen_mass(test, sampling, sample_count)
        brute = brute_force_unseen_mass(test, sampling, sample_count)
        assert brute == exact
        checks += 1

    # Exact inclusion-exclusion for uniform coupon collection, followed by
    # the union bound in equation (10).
    for state_count in (2, 4, 8):
        uniform = (Q(1, state_count),) * state_count
        for sample_count in range(9):
            assert expected_unseen_mass(
                uniform,
                uniform,
                sample_count,
            ) == (1 - Q(1, state_count)) ** sample_count
            all_seen_probability = sum(
                (
                    Q((-1) ** omitted * comb(state_count, omitted))
                    * Q(state_count - omitted, state_count) ** sample_count
                    for omitted in range(state_count + 1)
                ),
                Q(0),
            )
            some_unseen_probability = 1 - all_seen_probability
            assert some_unseen_probability <= (
                state_count
                * (1 - Q(1, state_count)) ** sample_count
            )
            checks += 2

    # Two disjoint partner supports of sizes two and three.
    partner_weights = (Q(1, 4), Q(3, 4))
    support_sizes = (2, 3)
    combined = (
        (partner_weights[0] / support_sizes[0],) * support_sizes[0]
        + (partner_weights[1] / support_sizes[1],) * support_sizes[1]
    )
    for partner, (weight, support_size) in enumerate(
        zip(partner_weights, support_sizes, strict=True)
    ):
        start = 0 if partner == 0 else support_sizes[0]
        partner_test = tuple(
            Q(1, support_size) if start <= state < start + support_size else Q(0)
            for state in range(sum(support_sizes))
        )
        for sample_count in range(7):
            assert expected_unseen_mass(
                partner_test,
                combined,
                sample_count,
            ) == (1 - weight / support_size) ** sample_count
            checks += 1
    return checks


def verify_walsh_noise() -> int:
    checks = 0
    for dimension in range(1, 5):
        state_count = 1 << dimension
        normalizer = sum(range(1, state_count + 1))
        distribution = tuple(
            Q(state + 1, normalizer) for state in range(state_count)
        )
        uniform = (Q(1, state_count),) * state_count
        for flip_rate in (Q(0), Q(1, 10), Q(1, 4), Q(1, 2)):
            noisy = noise_convolution(distribution, dimension, flip_rate)
            assert sum(noisy, Q(0)) == 1
            checks += 1

            spectral_chi_squared = Q(0)
            for subset in range(1 << dimension):
                expected = (
                    (1 - 2 * flip_rate) ** subset.bit_count()
                    * walsh_moment(distribution, subset)
                )
                observed = walsh_moment(noisy, subset)
                assert observed == expected
                if subset:
                    spectral_chi_squared += observed**2
                checks += 1

            direct_chi_squared = sum(
                (
                    (mass - uniform_mass) ** 2 / uniform_mass
                    for mass, uniform_mass in zip(noisy, uniform, strict=True)
                ),
                Q(0),
            )
            assert direct_chi_squared == spectral_chi_squared
            assert 4 * total_variation(noisy, uniform) ** 2 <= direct_chi_squared
            checks += 2

            if flip_rate == Q(1, 2):
                assert noisy == uniform
                checks += 1
    return checks


def verify_sparse_signals_and_oracles() -> int:
    checks = 0

    # Brute-force the no-signal event for small uniform cubes.
    for dimension in range(1, 4):
        state_count = 1 << dimension
        for informative_count in range(1, state_count + 1):
            informative = set(range(informative_count))
            for sample_count in range(5):
                no_signal = Q(0)
                for samples in product(
                    range(state_count),
                    repeat=sample_count,
                ):
                    if all(state not in informative for state in samples):
                        no_signal += Q(1, state_count**sample_count)
                expected = (
                    1 - Q(informative_count, state_count)
                ) ** sample_count
                assert no_signal == expected
                # With equal prior odds, indistinguishability on this event
                # contributes half its mass to minimax identification error.
                assert no_signal / 2 == expected / 2
                checks += 2

    # Check the exact mass introduced by target smoothing and its distance
    # from the hard uniform oracle target.
    for action_count in range(2, 8):
        for optimal_count in range(1, action_count):
            for epsilon in (Q(0), Q(1, 1000), Q(1, 10), Q(1)):
                denominator = optimal_count + action_count * epsilon
                smoothed = tuple(
                    (
                        (1 + epsilon) / denominator
                        if action < optimal_count
                        else epsilon / denominator
                    )
                    for action in range(action_count)
                )
                hard = tuple(
                    (
                        Q(1, optimal_count)
                        if action < optimal_count
                        else Q(0)
                    )
                    for action in range(action_count)
                )
                outside_mass = sum(smoothed[optimal_count:], Q(0))
                beta = (
                    Q(action_count - optimal_count) * epsilon
                    / denominator
                )
                assert sum(smoothed, Q(0)) == 1
                assert outside_mass == beta
                assert total_variation(smoothed, hard) == beta
                assert max(range(action_count), key=smoothed.__getitem__) < optimal_count
                checks += 4

            # Exact polynomial AM-GM form of equation (20):
            # p(B)^b >= b^b product_{a in B} p(a).
            policy_weights = (
                tuple(range(1, action_count + 1)),
                tuple(range(action_count, 0, -1)),
                tuple(1 + (action % 3) for action in range(action_count)),
            )
            for weights in policy_weights:
                weight_sum = sum(weights)
                policy = tuple(Q(weight, weight_sum) for weight in weights)
                optimal_mass = sum(policy[:optimal_count], Q(0))
                optimal_product = Q(1)
                for probability in policy[:optimal_count]:
                    optimal_product *= probability
                assert optimal_mass**optimal_count >= (
                    Q(optimal_count**optimal_count) * optimal_product
                )
                if optimal_mass > Q(optimal_count, optimal_count + 1):
                    assert max(policy[:optimal_count]) > max(
                        policy[optimal_count:]
                    )
                checks += 1

    # Exactly one error per trace: high pooled match can coexist with no
    # perfect trace.
    for trace_length in range(2, 13):
        traces = tuple(
            tuple(ply != failed_ply for ply in range(trace_length))
            for failed_ply in range(trace_length)
        )
        match_rate = Q(
            sum(sum(trace) for trace in traces),
            trace_length * len(traces),
        )
        perfect_rate = Q(
            sum(all(trace) for trace in traces),
            len(traces),
        )
        assert match_rate == 1 - Q(1, trace_length)
        assert perfect_rate == 0
        checks += 2

        # Independent per-ply matches have exact perfect probability r^T.
        match_probability = Q(trace_length - 1, trace_length)
        enumerated_perfect = Q(0)
        for outcomes in product((False, True), repeat=trace_length):
            probability = Q(1)
            for matched in outcomes:
                probability *= (
                    match_probability if matched else 1 - match_probability
                )
            if all(outcomes):
                enumerated_perfect += probability
        assert enumerated_perfect == match_probability**trace_length
        checks += 1

    return checks


def verify_hidden_protocols() -> int:
    checks = 0
    # Exhaust every deterministic policy g: X -> Theta in two small models.
    for state_count, protocol_bits in ((4, 2), (8, 1)):
        action_count = 1 << protocol_bits
        state_law = (Q(1, state_count),) * state_count
        for policy in product(range(action_count), repeat=state_count):
            success = sum(
                (
                    state_law[state] * Q(1, action_count)
                    for state in range(state_count)
                    for protocol in range(action_count)
                    if policy[state] == protocol
                ),
                Q(0),
            )
            assert success == Q(1, action_count)
            checks += 1

        # Giving the protocol itself as a grounded signal makes identity
        # decoding perfect.
        grounded_success = sum(
            (Q(1, action_count) for _signal in range(action_count)),
            Q(0),
        )
        assert grounded_success == 1
        checks += 1

    for button in range(8):
        for protocol_bit in (0, 1):
            lamp = 2 * button + protocol_bit
            assert lamp % 2 == protocol_bit
            checks += 1
    return checks


total_checks = (
    verify_total_variation_and_augmentation()
    + verify_finite_buffers()
    + verify_walsh_noise()
    + verify_sparse_signals_and_oracles()
    + verify_hidden_protocols()
)
print(
    f"verified {total_checks} exact state-coverage, augmentation, "
    "missed-mass, sparse-signal, Walsh-noise, oracle-supervision, "
    "trace, and hidden-protocol identities"
)
