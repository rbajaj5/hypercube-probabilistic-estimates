"""Verify the KAN/Boolean-cube collapse and Walsh error certificates."""

from __future__ import annotations

from collections import defaultdict, deque
from fractions import Fraction
from itertools import product


def fwht(values: list[int | Fraction]) -> list[int | Fraction]:
    output = values[:]
    block = 1
    while block < len(output):
        for start in range(0, len(output), 2 * block):
            for offset in range(block):
                left = output[start + offset]
                right = output[start + offset + block]
                output[start + offset] = left + right
                output[start + offset + block] = left - right
        block *= 2
    return output


def character(mask: int, subset: int) -> int:
    return -1 if (mask & subset).bit_count() % 2 else 1


def normalized_coefficients(
    values: list[int | Fraction],
) -> list[Fraction]:
    size = len(values)
    return [Fraction(value, size) for value in fwht(values)]


def level_masses(coefficients: list[Fraction]) -> dict[int, Fraction]:
    masses: dict[int, Fraction] = defaultdict(Fraction)
    for subset, coefficient in enumerate(coefficients):
        masses[subset.bit_count()] += coefficient * coefficient
    return dict(masses)


# Any univariate edge function on {-1, 1} is affine on its observed domain.
edge_checks = 0
for minus, plus in product(range(-4, 5), repeat=2):
    intercept = Fraction(plus + minus, 2)
    slope = Fraction(plus - minus, 2)
    assert intercept - slope == minus
    assert intercept + slope == plus
    edge_checks += 2

# Polynomial coefficient directions (t^2 - 1)t^j are invisible at both
# binary input points.  This witnesses p - 2 null directions for p
# polynomial coefficients.
nullity_checks = 0
for parameter_count in range(2, 11):
    for degree in range(parameter_count - 2):
        for point in (-1, 1):
            assert (point * point - 1) * point**degree == 0
            nullity_checks += 1


# Exhaust every Boolean observable through dimension three.  Verify that its
# best additive approximation is the level-at-most-one Walsh projection, and
# verify every degree truncation/noise certificate directly.
projection_checks = 0
noise_checks = 0
refresh_checks = 0
rho_values = (
    Fraction(0),
    Fraction(1, 7),
    Fraction(1, 3),
    Fraction(1, 2),
    Fraction(4, 5),
    Fraction(1),
)
for dimension in range(1, 4):
    size = 1 << dimension
    for function_mask in range(1 << size):
        values = [
            1 if function_mask & (1 << state) else -1
            for state in range(size)
        ]
        coefficients = normalized_coefficients(values)
        masses = level_masses(coefficients)
        assert sum(masses.values()) == 1

        additive = [
            sum(
                coefficient * character(state, subset)
                for subset, coefficient in enumerate(coefficients)
                if subset.bit_count() <= 1
            )
            for state in range(size)
        ]
        direct_error = sum(
            (Fraction(value) - prediction) ** 2
            for value, prediction in zip(values, additive, strict=True)
        ) / size
        spectral_error = sum(
            mass for level, mass in masses.items() if level >= 2
        )
        assert direct_error == spectral_error
        projection_checks += 1

        for rho in rho_values:
            exact_stability = sum(
                rho**level * mass for level, mass in masses.items()
            )
            for degree in range(dimension + 1):
                estimate = sum(
                    rho**level * mass
                    for level, mass in masses.items()
                    if level <= degree
                )
                tail = sum(
                    mass
                    for level, mass in masses.items()
                    if level > degree
                )
                assert 0 <= exact_stability - estimate
                assert exact_stability - estimate <= rho ** (degree + 1) * tail
                noise_checks += 1

        for time in range(6):
            exact_correlation = sum(
                mass * Fraction(dimension - level, dimension) ** time
                for level, mass in masses.items()
            )
            for degree in range(dimension):
                estimate = sum(
                    mass * Fraction(dimension - level, dimension) ** time
                    for level, mass in masses.items()
                    if level <= degree
                )
                tail = sum(
                    mass
                    for level, mass in masses.items()
                    if level > degree
                )
                remainder = exact_correlation - estimate
                rate = Fraction(dimension - degree - 1, dimension)
                assert 0 <= remainder <= tail * rate**time
                refresh_checks += 1


# Every permutation-symmetric Boolean function factors exactly through the
# Hamming sum and a piecewise-linear univariate interpolant.
def linear_interpolant(nodes: list[int], values: list[int], point: int) -> Fraction:
    if point <= nodes[0]:
        return Fraction(values[0])
    if point >= nodes[-1]:
        return Fraction(values[-1])
    for index in range(len(nodes) - 1):
        left, right = nodes[index : index + 2]
        if left <= point <= right:
            weight = Fraction(point - left, right - left)
            return (1 - weight) * values[index] + weight * values[index + 1]
    raise AssertionError("point not bracketed")


symmetric_checks = 0
for dimension in range(1, 7):
    nodes = list(range(-dimension, dimension + 1, 2))
    for function_mask in range(1 << (dimension + 1)):
        outputs = [
            1 if function_mask & (1 << index) else -1
            for index in range(dimension + 1)
        ]
        for state in range(1 << dimension):
            hamming_sum = dimension - 2 * state.bit_count()
            node_index = (hamming_sum + dimension) // 2
            assert linear_interpolant(nodes, outputs, hamming_sum) == outputs[node_index]
            symmetric_checks += 1


# Independent small-Hex enumeration.
DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1))


def crossing(mask: int, side: int, color: int, axis: int) -> bool:
    def has_color(x: int, y: int) -> bool:
        return ((mask >> (y * side + x)) & 1) == color

    starts = []
    for coordinate in range(side):
        x, y = (0, coordinate) if axis == 0 else (coordinate, 0)
        if has_color(x, y):
            starts.append((x, y))
    queue = deque(starts)
    seen = set(starts)
    while queue:
        x, y = queue.popleft()
        if (x if axis == 0 else y) == side - 1:
            return True
        for dx, dy in DIRECTIONS:
            xx, yy = x + dx, y + dy
            if (
                0 <= xx < side
                and 0 <= yy < side
                and (xx, yy) not in seen
                and has_color(xx, yy)
            ):
                seen.add((xx, yy))
                queue.append((xx, yy))
    return False


def winner(mask: int, side: int) -> int:
    blue = crossing(mask, side, color=1, axis=0)
    yellow = crossing(mask, side, color=0, axis=1)
    assert blue != yellow
    return 1 if blue else -1


expected_unnormalized_masses = {
    1: {1: 4},
    2: {1: 160, 2: 64, 3: 32},
    3: {
        1: 124240,
        2: 62272,
        3: 51584,
        4: 16320,
        5: 5920,
        6: 1216,
        7: 512,
        8: 64,
        9: 16,
    },
}
expected_errors = {
    1: (Fraction(0), Fraction(0), Fraction(0)),
    2: (Fraction(3, 8), Fraction(1, 8), Fraction(0)),
    3: (
        Fraction(8619, 16384),
        Fraction(4727, 16384),
        Fraction(1503, 16384),
    ),
}
hex_checks = 0
for side, expected_masses in expected_unnormalized_masses.items():
    dimension = side * side
    size = 1 << dimension
    values = [winner(mask, side) for mask in range(size)]
    transform = fwht(values)
    integer_masses: dict[int, int] = defaultdict(int)
    for subset, coefficient in enumerate(transform):
        integer_masses[subset.bit_count()] += coefficient * coefficient
    integer_masses = {
        level: mass for level, mass in integer_masses.items() if mass
    }
    assert integer_masses == expected_masses
    assert sum(values) == 0

    masses = {
        level: Fraction(mass, 4**dimension)
        for level, mass in integer_masses.items()
    }
    errors = tuple(
        sum(mass for level, mass in masses.items() if level > degree)
        for degree in (1, 2, 3)
    )
    assert errors == expected_errors[side]

    for rho in rho_values:
        exact = sum(rho**level * mass for level, mass in masses.items())
        for degree, error in zip((1, 2, 3), errors, strict=True):
            estimate = sum(
                rho**level * mass
                for level, mass in masses.items()
                if level <= degree
            )
            assert 0 <= exact - estimate <= rho ** (degree + 1) * error
            hex_checks += 1


print(
    f"verified {edge_checks} binary edge evaluations, "
    f"{nullity_checks} invisible coefficient directions, "
    f"{projection_checks} exhaustive additive projections, "
    f"{noise_checks} Walsh noise certificates, "
    f"{refresh_checks} Walsh refresh certificates, "
    f"{symmetric_checks} exact symmetric two-layer values, and "
    f"{hex_checks} exact small-Hex KAN/Walsh bounds"
)
