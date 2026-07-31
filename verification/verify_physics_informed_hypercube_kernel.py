"""Verify the exact physics-informed Walsh filter on finite cubes."""

from __future__ import annotations

from fractions import Fraction
from itertools import product


Q = Fraction


def character(state: int, subset: int) -> int:
    return -1 if (state & subset).bit_count() % 2 else 1


def walsh(values: tuple[Fraction, ...], subset: int) -> Fraction:
    return sum(
        (
            value * character(state, subset)
            for state, value in enumerate(values)
        ),
        Q(0),
    ) / len(values)


def inverse_walsh(
    coefficients: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    return tuple(
        sum(
            (
                coefficient * character(state, subset)
                for subset, coefficient in enumerate(coefficients)
            ),
            Q(0),
        )
        for state in range(len(coefficients))
    )


def kernel(state: int, target: int, dimension: int, rho: Fraction) -> Fraction:
    output = Q(1)
    for coordinate in range(dimension):
        same = ((state ^ target) & (1 << coordinate)) == 0
        output *= 1 + rho * (1 if same else -1)
    return output


def laplacian(values: tuple[Fraction, ...], dimension: int) -> tuple[Fraction, ...]:
    return tuple(
        sum(
            (
                values[state] - values[state ^ (1 << coordinate)]
                for coordinate in range(dimension)
            ),
            Q(0),
        )
        for state in range(1 << dimension)
    )


def verify_spectra() -> int:
    checks = 0
    for dimension in range(1, 6):
        state_count = 1 << dimension
        for rho in (Q(1, 5), Q(1, 2), Q(1)):
            for subset in range(state_count):
                kernel_action = tuple(
                    sum(
                        (
                            kernel(x, y, dimension, rho)
                            * character(y, subset)
                            for y in range(state_count)
                        ),
                        Q(0),
                    )
                    / state_count
                    for x in range(state_count)
                )
                eigenvalue = rho ** subset.bit_count()
                assert kernel_action == tuple(
                    eigenvalue * character(x, subset)
                    for x in range(state_count)
                )
                assert eigenvalue > 0
                checks += 2

        for subset in range(state_count):
            mode = tuple(
                Q(character(state, subset)) for state in range(state_count)
            )
            assert laplacian(mode, dimension) == tuple(
                Q(2 * subset.bit_count()) * value for value in mode
            )
            checks += 1
    return checks


def verify_filter() -> int:
    checks = 0
    for dimension in range(1, 6):
        state_count = 1 << dimension
        target = tuple(
            Q(((7 * state + 3) % 11) - 5, 7)
            for state in range(state_count)
        )
        target_hat = tuple(walsh(target, subset) for subset in range(state_count))
        assert inverse_walsh(target_hat) == target
        checks += 1

        for rho, gamma, ridge in (
            (Q(1, 2), Q(0), Q(1, 7)),
            (Q(1, 3), Q(2, 5), Q(1, 4)),
            (Q(3, 4), Q(5, 3), Q(2, 7)),
        ):
            estimate_hat = []
            expected_error = Q(0)
            expected_residual = Q(0)
            for subset, coefficient in enumerate(target_hat):
                degree = subset.bit_count()
                mu = Q(2 * degree)
                kappa = rho**degree
                denominator = 1 + gamma * mu**2 + ridge / kappa
                estimate = (1 + gamma * mu**2) * coefficient / denominator
                estimate_hat.append(estimate)

                # Coefficientwise first-order condition.
                assert (
                    (estimate - coefficient)
                    + gamma * mu * (mu * estimate - mu * coefficient)
                    + ridge * estimate / kappa
                    == 0
                )
                error = ridge / kappa * coefficient / denominator
                expected_error += error**2
                expected_residual += mu**2 * error**2
                checks += 1

            estimate = inverse_walsh(tuple(estimate_hat))
            difference = tuple(
                observed - truth
                for observed, truth in zip(estimate, target, strict=True)
            )
            observed_error = sum((value**2 for value in difference), Q(0)) / state_count
            observed_residual_values = laplacian(difference, dimension)
            observed_residual = (
                sum((value**2 for value in observed_residual_values), Q(0))
                / state_count
            )
            assert observed_error == expected_error
            assert observed_residual == expected_residual
            checks += 2
    return checks


def verify_noisy_risk() -> int:
    checks = 0
    for signal, mu, kappa, ridge, value_noise, operator_noise in (
        (Q(2, 3), Q(2), Q(1, 2), Q(1, 5), Q(1, 7), Q(1, 4)),
        (Q(-3, 5), Q(4), Q(1, 4), Q(2, 7), Q(2, 9), Q(1, 3)),
        (Q(5, 8), Q(6), Q(1, 8), Q(3, 10), Q(1, 6), Q(2, 5)),
    ):
        r = ridge / kappa
        value_variance = value_noise**2
        operator_variance = operator_noise**2
        optimum = (
            r**2 * signal**2 + value_variance
        ) / (operator_variance * (1 + r))

        for gamma in (Q(0), Q(1, 3), optimum, 2 * optimum):
            enumerated_risk = Q(0)
            for value_sign, operator_sign in product((-1, 1), repeat=2):
                observed_value = signal + value_sign * value_noise
                observed_operator = mu * signal + operator_sign * operator_noise
                estimate = (
                    observed_value + gamma * mu * observed_operator
                ) / (1 + gamma * mu**2 + r)
                enumerated_risk += (estimate - signal) ** 2 / 4

            formula = (
                r**2 * signal**2
                + value_variance
                + gamma**2 * mu**2 * operator_variance
            ) / (1 + gamma * mu**2 + r) ** 2
            assert enumerated_risk == formula
            checks += 1

        # The derivative numerator is linear in gamma and changes sign at
        # the stated global minimizer.
        base = r**2 * signal**2 + value_variance
        assert operator_variance * (1 + r) * optimum == base
        risk_at_optimum = (
            base + optimum**2 * mu**2 * operator_variance
        ) / (1 + optimum * mu**2 + r) ** 2
        for competitor in (Q(0), optimum / 2, 2 * optimum, 4 * optimum):
            competitor_risk = (
                base + competitor**2 * mu**2 * operator_variance
            ) / (1 + competitor * mu**2 + r) ** 2
            assert risk_at_optimum <= competitor_risk
            checks += 1
    return checks


total_checks = verify_spectra() + verify_filter() + verify_noisy_risk()
print(
    f"verified {total_checks} exact Hamming-kernel, cube-Laplacian, "
    "physics-informed filter, bias, residual, and noisy-risk identities"
)
