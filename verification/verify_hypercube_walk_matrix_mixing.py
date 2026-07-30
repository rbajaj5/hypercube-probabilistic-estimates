"""Verify hypercube-walk mixing and almost-orthogonal frame formulas."""

from __future__ import annotations

import itertools
import math
from collections import defaultdict
from fractions import Fraction


def character(state: int, subset: int) -> int:
    return -1 if (state & subset).bit_count() % 2 else 1


def refresh_step(
    distribution: dict[int, Fraction], dimension: int
) -> dict[int, Fraction]:
    updated: dict[int, Fraction] = defaultdict(Fraction)
    weight = Fraction(1, 2 * dimension)
    for state, probability in distribution.items():
        for coordinate in range(dimension):
            updated[state] += probability * weight
            updated[state ^ (1 << coordinate)] += probability * weight
    return dict(updated)


# Exact Walsh eigenvalues.
eigenvalue_checks = 0
for dimension in range(1, 9):
    for subset in range(1 << dimension):
        level = subset.bit_count()
        expected_ratio = Fraction(dimension - level, dimension)
        for state in range(1 << dimension):
            averaged = Fraction(0)
            for coordinate in range(dimension):
                averaged += Fraction(character(state, subset), 2 * dimension)
                averaged += Fraction(
                    character(state ^ (1 << coordinate), subset),
                    2 * dimension,
                )
            assert averaged == expected_ratio * character(state, subset)
            eigenvalue_checks += 1


# Exact chi-squared formula from a point mass.
chi_squared_checks = 0
for dimension in range(1, 8):
    distribution = {0: Fraction(1)}
    uniform = Fraction(1, 2**dimension)
    for time in range(9):
        direct = sum(
            (distribution.get(state, Fraction(0)) - uniform) ** 2 / uniform
            for state in range(1 << dimension)
        )
        spectral = sum(
            math.comb(dimension, level)
            * Fraction(dimension - level, dimension) ** (2 * time)
            for level in range(1, dimension + 1)
        )
        assert direct == spectral
        chi_squared_checks += 1
        distribution = refresh_step(distribution, dimension)


# Non-lazy single-coordinate flip remains in one parity class at each time.
periodicity_checks = 0
for dimension in range(1, 9):
    states = {0}
    for time in range(8):
        parities = {state.bit_count() % 2 for state in states}
        assert parities == {time % 2}
        states = {
            state ^ (1 << coordinate)
            for state in states
            for coordinate in range(dimension)
        }
        periodicity_checks += 1


def projector_overlap(
    first: tuple[Fraction, ...], second: tuple[Fraction, ...]
) -> Fraction:
    inner = sum(a * b for a, b in zip(first, second, strict=True))
    return inner * inner


vectors = (
    (Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(1)),
    (Fraction(3, 5), Fraction(4, 5)),
    (Fraction(4, 5), Fraction(-3, 5)),
)
dimension = len(vectors)
overlaps = {
    (i, j): projector_overlap(vectors[i], vectors[j])
    for i in range(dimension)
    for j in range(i + 1, dimension)
}


def centered_energy(state: int) -> Fraction:
    return 2 * sum(
        overlap
        * (1 if ((state >> i) & 1) == ((state >> j) & 1) else -1)
        for (i, j), overlap in overlaps.items()
    )


# The matrix energy is pure Walsh degree two.
stationary_variance = sum(
    centered_energy(state) ** 2 for state in range(1 << dimension)
) / (2**dimension)
expected_variance = 4 * sum(value * value for value in overlaps.values())
assert stationary_variance == expected_variance

energy_distribution = {0: Fraction(1)}
initial_energy = centered_energy(0)
energy_mixing_checks = 0
for time in range(11):
    actual = sum(
        probability * centered_energy(state)
        for state, probability in energy_distribution.items()
    )
    expected = Fraction(dimension - 2, dimension) ** time * initial_energy
    assert actual == expected
    energy_mixing_checks += 1
    energy_distribution = refresh_step(energy_distribution, dimension)


# Exact norm-category pushforward for the two-basis projector cube.
q = Fraction(3, 5)


def norm_square_label(state: int) -> Fraction:
    signs = tuple(1 if state & (1 << j) else -1 for j in range(4))
    a, b, c, d = signs
    first_aligned = a == b
    second_aligned = c == d
    if first_aligned and second_aligned:
        return Fraction(4) if a == c else Fraction(0)
    if first_aligned != second_aligned:
        return Fraction(4)
    x = Fraction(a - b, 2)
    y = Fraction(d - c, 2)
    return 2 + 2 * x * y * q


stationary_categories: dict[Fraction, Fraction] = defaultdict(Fraction)
conformal_nonzero = 0
for state in range(16):
    stationary_categories[norm_square_label(state)] += Fraction(1, 16)
    signs = tuple(1 if state & (1 << j) else -1 for j in range(4))
    a, b, c, d = signs
    if a == b and c == d and a == c:
        conformal_nonzero += 1
    elif a != b and c != d:
        conformal_nonzero += 1
assert conformal_nonzero == 6
assert stationary_categories == {
    Fraction(0): Fraction(1, 8),
    Fraction(4): Fraction(5, 8),
    Fraction(16, 5): Fraction(1, 8),
    Fraction(4, 5): Fraction(1, 8),
}

pushforward_distribution = {0: Fraction(1)}
pushforward_checks = 0
for _time in range(13):
    uniform = Fraction(1, 16)
    cube_tv = sum(
        abs(pushforward_distribution.get(state, Fraction(0)) - uniform)
        for state in range(16)
    ) / 2
    categories: dict[Fraction, Fraction] = defaultdict(Fraction)
    for state, probability in pushforward_distribution.items():
        categories[norm_square_label(state)] += probability
    category_tv = sum(
        abs(categories.get(label, Fraction(0)) - stationary_probability)
        for label, stationary_probability in stationary_categories.items()
    ) / 2
    assert category_tv <= cube_tv
    pushforward_checks += 1
    pushforward_distribution = refresh_step(pushforward_distribution, 4)


# GOE angular beta moments and the averaged frame-variance formula.
beta_checks = 0
for m in range(1, 101):
    ambient_dimension = 2 * m + 1
    mean = Fraction(1, 2) / Fraction(2 * m + 1, 2)
    second_moment = (
        Fraction(1, 2)
        * Fraction(3, 2)
        / (
            Fraction(2 * m + 1, 2)
            * Fraction(2 * m + 3, 2)
        )
    )
    assert mean == Fraction(1, ambient_dimension)
    assert second_moment == Fraction(
        3, ambient_dimension * (ambient_dimension + 2)
    )
    for frame_size in (2, 3, 5, 8):
        averaged_variance = (
            4 * math.comb(frame_size, 2) * second_moment
        )
        assert averaged_variance == Fraction(
            6 * frame_size * (frame_size - 1),
            ambient_dimension * (ambient_dimension + 2),
        )
        beta_checks += 1


print(
    f"verified {eigenvalue_checks} Walsh eigenfunction values, "
    f"{chi_squared_checks} chi-squared identities, "
    f"{periodicity_checks} periodic flip supports, "
    f"{energy_mixing_checks} degree-two matrix relaxations, "
    f"{pushforward_checks} norm-pushforward contractions, "
    "6 conformal nonzero projector vertices, and "
    f"{beta_checks} GOE beta/frame moment identities"
)
