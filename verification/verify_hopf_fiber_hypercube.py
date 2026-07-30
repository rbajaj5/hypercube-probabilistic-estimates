"""Verify the Hopf-fiber/projective hypercube identities."""

from __future__ import annotations

import cmath
import math
from collections import defaultdict, deque
from fractions import Fraction
from itertools import combinations


DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1))


def signs(mask: int, dimension: int) -> tuple[int, ...]:
    return tuple(1 if mask & (1 << i) else -1 for i in range(dimension))


def projective_kernel(left: tuple[int, ...], right: tuple[int, ...]) -> Fraction:
    dimension = len(left)
    overlap = sum(x * y for x, y in zip(left, right, strict=True))
    return Fraction(overlap * overlap, dimension * dimension)


def kernel_walsh_form(
    left: tuple[int, ...], right: tuple[int, ...]
) -> Fraction:
    dimension = len(left)
    level_two = sum(
        left[i] * left[j] * right[i] * right[j]
        for i, j in combinations(range(dimension), 2)
    )
    return Fraction(1, dimension) + Fraction(2 * level_two, dimension**2)


kernel_checks = 0
for dimension in range(1, 7):
    for left_mask in range(1 << dimension):
        left = signs(left_mask, dimension)
        for right_mask in range(1 << dimension):
            right = signs(right_mask, dimension)
            hamming = (left_mask ^ right_mask).bit_count()
            expected = Fraction((dimension - 2 * hamming) ** 2, dimension**2)
            assert projective_kernel(left, right) == expected
            assert kernel_walsh_form(left, right) == expected
            kernel_checks += 2


noise_checks = 0
for dimension in range(1, 8):
    fixed = (1,) * dimension
    for rho in (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)):
        direct = Fraction(0)
        for flip_mask in range(1 << dimension):
            hamming = flip_mask.bit_count()
            probability = (
                ((1 - rho) / 2) ** hamming
                * ((1 + rho) / 2) ** (dimension - hamming)
            )
            direct += probability * projective_kernel(
                fixed, signs(((1 << dimension) - 1) ^ flip_mask, dimension)
            )
        expected = rho**2 + (1 - rho**2) / dimension
        assert direct == expected
        noise_checks += 1


walk_checks = 0
for dimension in range(1, 7):
    distribution = {0: Fraction(1)}
    for time in range(11):
        direct = sum(
            probability
            * Fraction((dimension - 2 * state.bit_count()) ** 2, dimension**2)
            for state, probability in distribution.items()
        )
        expected = Fraction(1, dimension)
        expected += Fraction(dimension - 1, dimension) * (
            Fraction(dimension - 2, dimension) ** time
        )
        assert direct == expected
        walk_checks += 1

        updated = defaultdict(Fraction)
        for state, probability in distribution.items():
            updated[state] += probability / 2
            for coordinate in range(dimension):
                updated[state ^ (1 << coordinate)] += probability / (2 * dimension)
        distribution = dict(updated)


fiber_checks = 0
root_count = 97
for frequency in range(-12, 13):
    average = sum(
        cmath.exp(2j * math.pi * frequency * index / root_count)
        for index in range(root_count)
    ) / root_count
    expected = 1 if frequency == 0 else 0
    assert abs(average - expected) < 2e-14
    fiber_checks += 1
for first in range(-4, 5):
    for second in range(-4, 5):
        for third in range(-4, 5):
            average = sum(
                cmath.exp(
                    2j
                    * math.pi
                    * (first + second + third)
                    * index
                    / root_count
                )
                for index in range(root_count)
            ) / root_count
            expected = 1 if first + second + third == 0 else 0
            assert abs(average - expected) < 2e-14
            fiber_checks += 1


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


def hex_winner(mask: int, side: int) -> int:
    blue = crossing(mask, side, color=1, axis=0)
    yellow = crossing(mask, side, color=0, axis=1)
    assert blue != yellow
    return 1 if blue else -1


def fwht(values: list[int]) -> list[int]:
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


def level_masses(values: list[int]) -> dict[int, int]:
    masses = defaultdict(int)
    for subset, coefficient in enumerate(fwht(values)):
        masses[subset.bit_count()] += coefficient * coefficient
    return dict(masses)


def projective_weighted_formula(
    masses: dict[int, int], dimension: int, rho: Fraction
) -> Fraction:
    total = Fraction(0)
    normalization = 4**dimension
    for level, mass in masses.items():
        bracket = Fraction(rho**level, dimension)
        pair_sum = Fraction(level * (dimension - level)) * rho**level
        if level >= 2:
            pair_sum += Fraction(math.comb(level, 2)) * rho ** (level - 2)
        if dimension - level >= 2:
            pair_sum += (
                Fraction(math.comb(dimension - level, 2))
                * rho ** (level + 2)
            )
        bracket += Fraction(2, dimension**2) * pair_sum
        total += Fraction(mass, normalization) * bracket
    return total


def projective_weighted_coefficients(
    masses: dict[int, int], dimension: int
) -> dict[int, Fraction]:
    coefficients = defaultdict(Fraction)
    normalization = 4**dimension
    for level, mass in masses.items():
        weight = Fraction(mass, normalization)
        coefficients[level] += weight / dimension
        coefficients[level] += (
            weight
            * Fraction(2 * level * (dimension - level), dimension**2)
        )
        if level >= 2:
            coefficients[level - 2] += (
                weight * Fraction(2 * math.comb(level, 2), dimension**2)
            )
        if dimension - level >= 2:
            coefficients[level + 2] += (
                weight
                * Fraction(
                    2 * math.comb(dimension - level, 2),
                    dimension**2,
                )
            )
    return {degree: value for degree, value in coefficients.items() if value}


expected_hex_coefficients = {
    1: {1: Fraction(1)},
    2: {
        0: Fraction(1, 32),
        1: Fraction(7, 16),
        2: Fraction(3, 16),
        3: Fraction(5, 16),
        4: Fraction(1, 32),
    },
    3: {
        0: Fraction(7784, 1327104),
        1: Fraction(213469, 1327104),
        2: Fraction(156244, 1327104),
        3: Fraction(587320, 1327104),
        4: Fraction(215724, 1327104),
        5: Fraction(116194, 1327104),
        6: Fraction(24044, 1327104),
        7: Fraction(5696, 1327104),
        8: Fraction(556, 1327104),
        9: Fraction(73, 1327104),
    },
}


hex_checks = 0
for side in (1, 2, 3):
    dimension = side * side
    state_count = 1 << dimension
    values = [hex_winner(mask, side) for mask in range(state_count)]
    masses = level_masses(values)
    assert (
        projective_weighted_coefficients(masses, dimension)
        == expected_hex_coefficients[side]
    )
    hex_checks += 1
    for rho in (Fraction(0), Fraction(1, 2), Fraction(3, 4), Fraction(1)):
        direct = Fraction(0)
        for state, value in enumerate(values):
            for flip_mask in range(state_count):
                hamming = flip_mask.bit_count()
                probability = Fraction(1, state_count)
                probability *= ((1 - rho) / 2) ** hamming
                probability *= ((1 + rho) / 2) ** (dimension - hamming)
                kernel = Fraction((dimension - 2 * hamming) ** 2, dimension**2)
                direct += (
                    probability
                    * value
                    * values[state ^ flip_mask]
                    * kernel
                )
        assert direct == projective_weighted_formula(masses, dimension, rho)
        hex_checks += 1


determinant_checks = 0
for dimension in range(1, 13):
    multiplicities = defaultdict(int)
    for subset in range(1 << dimension):
        multiplicities[subset.bit_count()] += 1
    assert all(
        multiplicities[level] == math.comb(dimension, level)
        for level in range(dimension + 1)
    )
    log_determinant = sum(
        math.comb(dimension, level) * math.log(level / dimension)
        for level in range(1, dimension + 1)
    )
    log_partition = sum(
        Fraction(math.comb(dimension, level), 2)
        * math.log(dimension / level)
        for level in range(1, dimension + 1)
    )
    assert abs(log_partition + log_determinant / 2) < 1e-12
    determinant_checks += dimension + 2


print(
    f"verified {kernel_checks} projective/Hamming identities, "
    f"{noise_checks} product-noise laws, {walk_checks} refresh laws, "
    f"{fiber_checks} fiber Fourier selections, "
    f"{hex_checks} projectively weighted Hex correlations, and "
    f"{determinant_checks} finite determinant identities"
)
