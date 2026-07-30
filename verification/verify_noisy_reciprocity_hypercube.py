"""Verify noisy memory-one Prisoner's Dilemma hypercube identities."""

from __future__ import annotations

from fractions import Fraction
from itertools import product


F = Fraction
STATES = ((1, 1), (1, -1), (-1, 1), (-1, -1))
PAYOFF_X = (F(3), F(0), F(5), F(1))
PAYOFF_Y = (F(3), F(5), F(0), F(1))


def solve(matrix: list[list[Fraction]], vector: list[Fraction]) -> list[Fraction]:
    size = len(vector)
    augmented = [matrix[row][:] + [vector[row]] for row in range(size)]
    for column in range(size):
        pivot = next(
            row
            for row in range(column, size)
            if augmented[row][column]
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        for entry in range(column, size + 1):
            augmented[column][entry] /= pivot_value
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                for entry in range(column, size + 1):
                    augmented[row][entry] -= factor * augmented[column][entry]
    return [augmented[row][-1] for row in range(size)]


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    size = len(work)
    output = F(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return F(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            output = -output
        pivot_value = work[column][column]
        output *= pivot_value
        for entry in range(column, size):
            work[column][entry] /= pivot_value
        for row in range(column + 1, size):
            factor = work[row][column]
            if factor:
                for entry in range(column, size):
                    work[row][entry] -= factor * work[column][entry]
    return output


def transition(
    strategy_x: tuple[Fraction, ...],
    strategy_y: tuple[Fraction, ...],
    error: Fraction,
) -> list[list[Fraction]]:
    executed_x = [error + (1 - 2 * error) * value for value in strategy_x]
    y_from_x_view = (
        strategy_y[0],
        strategy_y[2],
        strategy_y[1],
        strategy_y[3],
    )
    executed_y = [
        error + (1 - 2 * error) * value
        for value in y_from_x_view
    ]
    output = []
    for px, py in zip(executed_x, executed_y, strict=True):
        output.append(
            [
                px * py,
                px * (1 - py),
                (1 - px) * py,
                (1 - px) * (1 - py),
            ]
        )
    return output


def stationary(matrix: list[list[Fraction]]) -> list[Fraction]:
    equations = [
        [
            matrix[column][row] - (1 if row == column else 0)
            for column in range(4)
        ]
        for row in range(3)
    ]
    equations.append([F(1)] * 4)
    return solve(equations, [F(0), F(0), F(0), F(1)])


def row_times_matrix(
    row: list[Fraction], matrix: list[list[Fraction]]
) -> list[Fraction]:
    return [
        sum(row[index] * matrix[index][column] for index in range(4))
        for column in range(4)
    ]


def operator_on_function(
    matrix: list[list[Fraction]], values: list[Fraction]
) -> list[Fraction]:
    return [
        sum(matrix[row][column] * values[column] for column in range(4))
        for row in range(4)
    ]


def inner(
    stationary_law: list[Fraction],
    left: list[Fraction],
    right: list[Fraction],
) -> Fraction:
    return sum(
        stationary_law[index] * left[index] * right[index]
        for index in range(4)
    )


def payoff(
    stationary_law: list[Fraction],
) -> tuple[Fraction, Fraction]:
    return (
        sum(probability * reward for probability, reward in zip(
            stationary_law, PAYOFF_X, strict=True
        )),
        sum(probability * reward for probability, reward in zip(
            stationary_law, PAYOFF_Y, strict=True
        )),
    )


def matrix_multiply(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [
            sum(left[row][index] * right[index][column] for index in range(4))
            for column in range(4)
        ]
        for row in range(4)
    ]


def identity_matrix() -> list[list[Fraction]]:
    return [
        [F(int(row == column)) for column in range(4)]
        for row in range(4)
    ]


def finite_horizon_payoff(
    initial: list[Fraction],
    matrix: list[list[Fraction]],
    horizon: int,
) -> Fraction:
    distribution = initial[:]
    total = F(0)
    for _ in range(horizon):
        total += sum(
            distribution[index] * PAYOFF_X[index]
            for index in range(4)
        )
        distribution = row_times_matrix(distribution, matrix)
    return total / horizon


ALL_C = (F(1), F(1), F(1), F(1))
ALL_D = (F(0), F(0), F(0), F(0))
TFT = (F(1), F(0), F(1), F(0))
GTFT = (F(1), F(1, 3), F(1), F(1, 3))
PAVLOV = (F(1), F(0), F(0), F(1))
STRATEGIES = (ALL_C, ALL_D, TFT, GTFT, PAVLOV)


# General transition and stationarity checks.
transition_checks = 0
stationary_checks = 0
for error in (F(1, 100), F(1, 20), F(1, 5)):
    for strategy_x, strategy_y in product(STRATEGIES, repeat=2):
        matrix = transition(strategy_x, strategy_y, error)
        assert all(sum(row) == 1 for row in matrix)
        assert all(0 < entry < 1 for row in matrix for entry in row)
        transition_checks += 20
        law = stationary(matrix)
        assert sum(law) == 1
        assert all(value > 0 for value in law)
        assert row_times_matrix(law, matrix) == law
        stationary_checks += 9


# Symmetric Generous Tit-for-Tat formulas and full characteristic polynomial.
reactive_checks = 0
x_values = [F(state[0]) for state in STATES]
y_values = [F(state[1]) for state in STATES]
antisymmetric = [
    x_values[index] - y_values[index]
    for index in range(4)
]
for error in (F(1, 100), F(1, 20), F(1, 5)):
    delta = 1 - 2 * error
    for generosity in (F(0), F(1, 10), F(1, 3), F(1)):
        strategy = (F(1), generosity, F(1), generosity)
        matrix = transition(strategy, strategy, error)
        alpha = delta * generosity
        beta = delta * (1 - generosity)
        cooperation = (error + delta * generosity) / (
            2 * error + delta * generosity
        )
        mean_spin = 2 * cooperation - 1
        expected_law = [
            cooperation * cooperation,
            cooperation * (1 - cooperation),
            (1 - cooperation) * cooperation,
            (1 - cooperation) * (1 - cooperation),
        ]
        assert row_times_matrix(expected_law, matrix) == expected_law
        assert stationary(matrix) == expected_law

        centered_symmetric = [
            x_values[index] + y_values[index] - 2 * mean_spin
            for index in range(4)
        ]
        assert operator_on_function(matrix, antisymmetric) == [
            -beta * value for value in antisymmetric
        ]
        assert operator_on_function(matrix, centered_symmetric) == [
            beta * value for value in centered_symmetric
        ]

        # Equality at five lambda values proves equality of the two monic
        # degree-four characteristic polynomials.
        for eigenvalue_parameter in (F(-2), F(-1), F(0), F(1), F(2)):
            characteristic_matrix = [
                [
                    eigenvalue_parameter * int(row == column)
                    - matrix[row][column]
                    for column in range(4)
                ]
                for row in range(4)
            ]
            expected_characteristic = (
                (eigenvalue_parameter - 1)
                * (eigenvalue_parameter - beta)
                * (eigenvalue_parameter + beta)
                * (eigenvalue_parameter - beta * beta)
            )
            assert determinant(characteristic_matrix) == expected_characteristic
            reactive_checks += 1

        # Exact retaliation correlation and two-step positivity.
        evolved = antisymmetric[:]
        variance = inner(expected_law, antisymmetric, antisymmetric)
        for time in range(9):
            if variance:
                assert (
                    inner(expected_law, antisymmetric, evolved) / variance
                    == (-beta) ** time
                )
            evolved = operator_on_function(matrix, evolved)
            reactive_checks += 1

        even_sequence = [beta ** (2 * time) for time in range(10)]
        current = even_sequence
        for order in range(6):
            sign = 1 if order % 2 == 0 else -1
            assert all(sign * value >= 0 for value in current)
            current = [
                current[index + 1] - current[index]
                for index in range(len(current) - 1)
            ]
            reactive_checks += len(current) + 1


# Displayed exact benchmark at five-percent implementation error.
benchmark_error = F(1, 20)
expected_self_play = {
    ALL_C: ([F(361, 400), F(19, 400), F(19, 400), F(1, 400)], F(1179, 400)),
    ALL_D: ([F(1, 400), F(19, 400), F(19, 400), F(361, 400)], F(459, 400)),
    TFT: ([F(1, 4)] * 4, F(9, 4)),
    GTFT: ([F(49, 64), F(7, 64), F(7, 64), F(1, 64)], F(183, 64)),
    PAVLOV: (
        [F(817, 1000), F(19, 400), F(19, 400), F(11, 125)],
        F(5553, 2000),
    ),
}
benchmark_checks = 0
for strategy, (expected_law, expected_payoff) in expected_self_play.items():
    law = stationary(transition(strategy, strategy, benchmark_error))
    assert law == expected_law
    assert payoff(law) == (expected_payoff, expected_payoff)
    benchmark_checks += 6

all_d_vs_all_c = payoff(stationary(transition(ALL_D, ALL_C, benchmark_error)))
assert all_d_vs_all_c == (F(1881, 400), F(81, 400))
all_d_vs_gtft = payoff(stationary(transition(ALL_D, GTFT, benchmark_error)))
assert all_d_vs_gtft == (F(2451, 1000), F(801, 1000))
benchmark_checks += 4


# Finite-horizon matrix-sum identity.  Compare iterative distributions with
# an explicitly accumulated matrix geometric sum.
finite_checks = 0
initial = [F(1), F(0), F(0), F(0)]
for strategy_x, strategy_y in product(STRATEGIES, repeat=2):
    matrix = transition(strategy_x, strategy_y, benchmark_error)
    power = identity_matrix()
    geometric = [[F(0) for _ in range(4)] for _ in range(4)]
    for horizon in range(1, 9):
        geometric = [
            [
                geometric[row][column] + power[row][column]
                for column in range(4)
            ]
            for row in range(4)
        ]
        exact_sum = sum(
            initial[row] * geometric[row][column] * PAYOFF_X[column]
            for row in range(4)
            for column in range(4)
        ) / horizon
        assert exact_sum == finite_horizon_payoff(
            initial, matrix, horizon
        )
        power = matrix_multiply(power, matrix)
        finite_checks += 1


print(
    f"verified {transition_checks} transition probabilities, "
    f"{stationary_checks} stationary-law identities, "
    f"{reactive_checks} reactive spectra/correlations, "
    f"{benchmark_checks} exact payoff benchmarks, and "
    f"{finite_checks} finite-horizon matrix sums"
)
