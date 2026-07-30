"""Verify exact latent-environment identities on finite Boolean cubes."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import comb


Q = Fraction
State = tuple[int, ...]


def falling(value: int, order: int) -> int:
    output = 1
    for offset in range(order):
        output *= value - offset
    return output


def mixture_probability(
    state: State,
    alpha: Fraction,
    high_rate: Fraction,
    low_rate: Fraction,
) -> Fraction:
    weight = sum(state)
    dimension = len(state)
    return (
        alpha
        * high_rate**weight
        * (1 - high_rate) ** (dimension - weight)
        + (1 - alpha)
        * low_rate**weight
        * (1 - low_rate) ** (dimension - weight)
    )


def verify_case(
    dimension: int,
    alpha: Fraction,
    high_rate: Fraction,
    low_rate: Fraction,
) -> int:
    assert dimension >= 2
    assert 0 < alpha < 1
    assert 0 < low_rate < high_rate < 1

    states = tuple(product((0, 1), repeat=dimension))
    probabilities = [
        mixture_probability(state, alpha, high_rate, low_rate)
        for state in states
    ]
    mean_rate = alpha * high_rate + (1 - alpha) * low_rate
    environment_variance = (
        alpha * (1 - alpha) * (high_rate - low_rate) ** 2
    )
    checks = 0

    assert sum(probabilities, Q(0)) == 1
    checks += 1

    # Exchangeability: the mass of a vertex depends only on its Hamming weight.
    for weight in range(dimension + 1):
        masses = {
            probability
            for state, probability in zip(
                states,
                probabilities,
                strict=True,
            )
            if sum(state) == weight
        }
        assert len(masses) == 1
        checks += 1

    for coordinate in range(dimension):
        observed_mean = sum(
            (
                probability * state[coordinate]
                for state, probability in zip(
                    states,
                    probabilities,
                    strict=True,
                )
            ),
            Q(0),
        )
        assert observed_mean == mean_rate
        checks += 1

    for left, right in combinations(range(dimension), 2):
        observed_product = sum(
            (
                probability * state[left] * state[right]
                for state, probability in zip(
                    states,
                    probabilities,
                    strict=True,
                )
            ),
            Q(0),
        )
        assert (
            observed_product - mean_rate**2
            == environment_variance
        )
        checks += 1

    # Every centered interaction is the corresponding centered moment of the
    # latent rate.
    for subset_bits in product((0, 1), repeat=dimension):
        subset = tuple(
            coordinate
            for coordinate, included in enumerate(subset_bits)
            if included
        )
        order = len(subset)
        if order == 0:
            continue
        observed = sum(
            (
                probability
                * product_value(
                    tuple(
                        Q(state[coordinate]) - mean_rate
                        for coordinate in subset
                    )
                )
                for state, probability in zip(
                    states,
                    probabilities,
                    strict=True,
                )
            ),
            Q(0),
        )
        expected = (high_rate - low_rate) ** order * (
            alpha * (1 - alpha) ** order
            + (1 - alpha) * (-alpha) ** order
        )
        assert observed == expected
        checks += 1

    count_law = [Q(0)] * (dimension + 1)
    for state, probability in zip(states, probabilities, strict=True):
        count_law[sum(state)] += probability

    expected_count_law = [
        Q(comb(dimension, weight))
        * (
            alpha
            * high_rate**weight
            * (1 - high_rate) ** (dimension - weight)
            + (1 - alpha)
            * low_rate**weight
            * (1 - low_rate) ** (dimension - weight)
        )
        for weight in range(dimension + 1)
    ]
    assert count_law == expected_count_law
    checks += dimension + 1

    observed_mean_count = sum(
        (
            Q(weight) * probability
            for weight, probability in enumerate(count_law)
        ),
        Q(0),
    )
    observed_second_moment = sum(
        (
            Q(weight * weight) * probability
            for weight, probability in enumerate(count_law)
        ),
        Q(0),
    )
    expected_variance = (
        dimension * mean_rate * (1 - mean_rate)
        + dimension * (dimension - 1) * environment_variance
    )
    assert observed_mean_count == dimension * mean_rate
    assert (
        observed_second_moment - observed_mean_count**2
        == expected_variance
    )
    checks += 2

    for order in range(dimension + 1):
        observed_factorial_moment = sum(
            (
                Q(falling(weight, order)) * probability
                for weight, probability in enumerate(count_law)
            ),
            Q(0),
        )
        expected_factorial_moment = Q(falling(dimension, order)) * (
            alpha * high_rate**order
            + (1 - alpha) * low_rate**order
        )
        assert observed_factorial_moment == expected_factorial_moment
        checks += 1

    # Bayes posterior and the constant adjacent-count odds multiplier.
    posterior_odds = []
    for weight in range(dimension + 1):
        high_mass = (
            alpha
            * high_rate**weight
            * (1 - high_rate) ** (dimension - weight)
        )
        low_mass = (
            (1 - alpha)
            * low_rate**weight
            * (1 - low_rate) ** (dimension - weight)
        )
        posterior = high_mass / (high_mass + low_mass)
        assert posterior / (1 - posterior) == high_mass / low_mass
        posterior_odds.append(high_mass / low_mass)
        checks += 1

    odds_multiplier = (
        high_rate * (1 - low_rate)
        / (low_rate * (1 - high_rate))
    )
    for weight in range(dimension):
        assert (
            posterior_odds[weight + 1]
            == posterior_odds[weight] * odds_multiplier
        )
        assert posterior_odds[weight + 1] > posterior_odds[weight]
        checks += 2

    return checks


def product_value(values: tuple[Fraction, ...]) -> Fraction:
    output = Q(1)
    for value in values:
        output *= value
    return output


CASES = (
    (2, Q(1, 2), Q(3, 4), Q(1, 4)),
    (3, Q(1, 3), Q(2, 3), Q(1, 5)),
    (4, Q(2, 5), Q(5, 6), Q(1, 6)),
    (5, Q(3, 7), Q(7, 8), Q(2, 9)),
    (6, Q(4, 9), Q(9, 10), Q(1, 10)),
    (7, Q(5, 11), Q(10, 11), Q(3, 20)),
)


total_checks = sum(
    verify_case(dimension, alpha, high_rate, low_rate)
    for dimension, alpha, high_rate, low_rate in CASES
)
print(
    f"verified {total_checks} exact latent-environment normalization, "
    "exchangeability, covariance, Fourier-moment, count, overdispersion, "
    f"factorial-moment, and posterior identities across {len(CASES)} cubes"
)
