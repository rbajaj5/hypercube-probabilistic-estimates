"""Verify novelty-refresh identities on finite Boolean hypercubes."""

from __future__ import annotations

from fractions import Fraction
from itertools import product as cartesian_product
from math import prod


Q = Fraction
State = tuple[int, ...]
Matrix = list[list[Fraction]]


def state_space(dimension: int) -> tuple[State, ...]:
    return tuple(cartesian_product((0, 1), repeat=dimension))


def stationary_probability(
    state: State,
    novelty: tuple[Fraction, ...],
) -> Fraction:
    output = Q(1)
    for bit, probability in zip(state, novelty, strict=True):
        output *= probability if bit else 1 - probability
    return output


def transition_matrix(
    selection: tuple[Fraction, ...],
    novelty: tuple[Fraction, ...],
) -> tuple[tuple[State, ...], Matrix]:
    states = state_space(len(selection))
    index = {state: position for position, state in enumerate(states)}
    matrix = [
        [Q(0) for _ in states]
        for _ in states
    ]
    for source_position, source in enumerate(states):
        for coordinate, selection_probability in enumerate(selection):
            for new_bit, reset_probability in (
                (0, 1 - novelty[coordinate]),
                (1, novelty[coordinate]),
            ):
                target = list(source)
                target[coordinate] = new_bit
                matrix[source_position][index[tuple(target)]] += (
                    selection_probability * reset_probability
                )
    return states, matrix


def apply_transition(
    distribution: list[Fraction],
    matrix: Matrix,
) -> list[Fraction]:
    size = len(distribution)
    return [
        sum(
            (
                distribution[source] * matrix[source][target]
                for source in range(size)
            ),
            Q(0),
        )
        for target in range(size)
    ]


def operator_on_function(
    matrix: Matrix,
    values: list[Fraction],
) -> list[Fraction]:
    return [
        sum(
            (
                matrix[source][target] * values[target]
                for target in range(len(values))
            ),
            Q(0),
        )
        for source in range(len(values))
    ]


def total_variation(
    left: list[Fraction],
    right: list[Fraction],
) -> Fraction:
    return sum(
        (abs(a - b) for a, b in zip(left, right, strict=True)),
        Q(0),
    ) / 2


def chi_squared(
    distribution: list[Fraction],
    stationary: list[Fraction],
) -> Fraction:
    return sum(
        (
            (observed - expected) ** 2 / expected
            for observed, expected in zip(
                distribution,
                stationary,
                strict=True,
            )
        ),
        Q(0),
    )


def poisson_binomial(novelty: tuple[Fraction, ...]) -> list[Fraction]:
    coefficients = [Q(1)]
    for probability in novelty:
        updated = [Q(0)] * (len(coefficients) + 1)
        for count, coefficient in enumerate(coefficients):
            updated[count] += coefficient * (1 - probability)
            updated[count + 1] += coefficient * probability
        coefficients = updated
    return coefficients


def verify_chain(
    selection: tuple[Fraction, ...],
    novelty: tuple[Fraction, ...],
) -> int:
    dimension = len(selection)
    assert len(novelty) == dimension
    assert sum(selection, Q(0)) == 1
    assert all(value > 0 for value in selection)
    assert all(0 < value < 1 for value in novelty)

    states, matrix = transition_matrix(selection, novelty)
    stationary = [
        stationary_probability(state, novelty)
        for state in states
    ]
    zero_position = states.index((0,) * dimension)
    checks = 0

    assert sum(stationary, Q(0)) == 1
    assert all(sum(row, Q(0)) == 1 for row in matrix)
    assert apply_transition(stationary, matrix) == stationary
    checks += len(states) + 2

    for left in range(len(states)):
        for right in range(len(states)):
            assert (
                stationary[left] * matrix[left][right]
                == stationary[right] * matrix[right][left]
            )
            checks += 1

    # Every centered product is an exact eigenfunction.  Square-root
    # normalization is unnecessary for this algebraic check.
    for subset_bits in cartesian_product((0, 1), repeat=dimension):
        subset = tuple(
            index
            for index, included in enumerate(subset_bits)
            if included
        )
        values = [
            prod(
                (
                    Q(state[index]) - novelty[index]
                    for index in subset
                ),
                start=Q(1),
            )
            for state in states
        ]
        eigenvalue = 1 - sum(
            (selection[index] for index in subset),
            Q(0),
        )
        assert operator_on_function(matrix, values) == [
            eigenvalue * value for value in values
        ]
        checks += len(states)

    # Stationary novelty count and exact entropy coefficients.
    stationary_count_law = [Q(0)] * (dimension + 1)
    for state, probability in zip(states, stationary, strict=True):
        stationary_count_law[sum(state)] += probability
    assert stationary_count_law == poisson_binomial(novelty)
    assert sum(
        (
            count * probability
            for count, probability in enumerate(stationary_count_law)
        ),
        Q(0),
    ) == sum(novelty, Q(0))
    assert 1 - stationary_count_law[0] == (
        1
        - prod(
            (1 - probability for probability in novelty),
            start=Q(1),
        )
    )
    for coordinate in range(dimension):
        coefficient_log_nu = sum(
            (
                probability * state[coordinate]
                for state, probability in zip(
                    states,
                    stationary,
                    strict=True,
                )
            ),
            Q(0),
        )
        coefficient_log_one_minus_nu = sum(
            (
                probability * (1 - state[coordinate])
                for state, probability in zip(
                    states,
                    stationary,
                    strict=True,
                )
            ),
            Q(0),
        )
        assert coefficient_log_nu == novelty[coordinate]
        assert coefficient_log_one_minus_nu == 1 - novelty[coordinate]
        checks += 2
    checks += 3

    escape_probability = sum(
        (
            selection[index] * novelty[index]
            for index in range(dimension)
        ),
        Q(0),
    )
    assert 1 - matrix[zero_position][zero_position] == escape_probability
    for time in range(9):
        assert (
            matrix[zero_position][zero_position] ** time
            == (1 - escape_probability) ** time
        )
        checks += 1

    # Finite-time mean, chi-squared identity, and two TV upper bounds.
    distribution = [Q(0)] * len(states)
    distribution[zero_position] = 1
    for time in range(9):
        observed_mean = sum(
            (
                probability * sum(state)
                for state, probability in zip(
                    states,
                    distribution,
                    strict=True,
                )
            ),
            Q(0),
        )
        expected_mean = sum(
            (
                novelty[index]
                * (1 - (1 - selection[index]) ** time)
                for index in range(dimension)
            ),
            Q(0),
        )
        assert observed_mean == expected_mean

        observed_chi_squared = chi_squared(distribution, stationary)
        spectral_chi_squared = Q(0)
        for subset_bits in cartesian_product((0, 1), repeat=dimension):
            subset = tuple(
                index
                for index, included in enumerate(subset_bits)
                if included
            )
            if not subset:
                continue
            initial_coefficient_squared = prod(
                (
                    novelty[index] / (1 - novelty[index])
                    for index in subset
                ),
                start=Q(1),
            )
            eigenvalue = 1 - sum(
                (selection[index] for index in subset),
                Q(0),
            )
            spectral_chi_squared += (
                initial_coefficient_squared * eigenvalue ** (2 * time)
            )
        assert observed_chi_squared == spectral_chi_squared

        observed_tv = total_variation(distribution, stationary)
        assert 4 * observed_tv * observed_tv <= observed_chi_squared
        coupon_bound = sum(
            (
                (1 - selection[index]) ** time
                for index in range(dimension)
            ),
            Q(0),
        )
        assert observed_tv <= coupon_bound
        checks += 4
        distribution = apply_transition(distribution, matrix)

    return checks


CHAINS = (
    ((Q(1),), (Q(1, 3),)),
    (
        (Q(1, 2), Q(1, 2)),
        (Q(1, 5), Q(2, 5)),
    ),
    (
        (Q(1, 2), Q(1, 3), Q(1, 6)),
        (Q(1, 7), Q(2, 7), Q(3, 7)),
    ),
    (
        (Q(1, 4),) * 4,
        (Q(1, 10), Q(1, 5), Q(3, 10), Q(2, 5)),
    ),
    (
        (Q(1, 15), Q(2, 15), Q(3, 15), Q(4, 15), Q(5, 15)),
        (Q(1, 11), Q(2, 11), Q(3, 11), Q(4, 11), Q(5, 11)),
    ),
)


total_checks = sum(
    verify_chain(selection, novelty)
    for selection, novelty in CHAINS
)
print(
    f"verified {total_checks} exact novelty-refresh transition, "
    "stationary, entropy-coefficient, spectral, escape, and mixing "
    f"identities across {len(CHAINS)} nonuniform hypercubes"
)
