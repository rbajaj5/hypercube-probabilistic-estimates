"""Verify Hex topology consequences and winner-mixing identities."""

from __future__ import annotations

import math
from collections import defaultdict, deque
from fractions import Fraction


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


def transpose_complement(mask: int, side: int) -> int:
    output = 0
    for y in range(side):
        for x in range(side):
            old_bit = (mask >> (x * side + y)) & 1
            new_bit = 1 - old_bit
            output |= new_bit << (y * side + x)
    return output


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


def level_masses(values: list[int | Fraction]) -> dict[int, int | Fraction]:
    transform = fwht(values)
    masses: dict[int, int | Fraction] = defaultdict(int)
    for subset, coefficient in enumerate(transform):
        masses[subset.bit_count()] += coefficient * coefficient
    return {level: mass for level, mass in masses.items() if mass}


# Exhaustive no-tie topology and twisted antisymmetry.
topology_checks = 0
winner_tables = {}
for side in range(1, 5):
    dimension = side * side
    blue_wins = 0
    values = []
    for mask in range(1 << dimension):
        value = winner(mask, side)
        values.append(value)
        blue_wins += value == 1
        assert winner(transpose_complement(mask, side), side) == -value
        topology_checks += 2
    assert blue_wins == 2 ** (dimension - 1)
    winner_tables[side] = values


# The derivative of noise stability at rho=1 is total influence.  These
# exact small-board values anchor the larger seeded Monte Carlo experiment.
expected_influences = {
    1: Fraction(1),
    2: Fraction(3, 2),
    3: Fraction(249, 128),
    4: Fraction(2405, 1024),
}
influence_checks = 0
exact_influences = {}
for side, values in winner_tables.items():
    dimension = side * side
    directed_pivots = sum(
        values[mask] != values[mask ^ (1 << coordinate)]
        for mask in range(1 << dimension)
        for coordinate in range(dimension)
    )
    influence = Fraction(directed_pivots, 1 << dimension)
    exact_influences[side] = influence
    if side in expected_influences:
        assert influence == expected_influences[side]
    influence_checks += 1


expected_winner_masses = {
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
expected_odd_masses = {
    1: {1: 4},
    2: {1: 160, 3: 32},
    3: {1: 124240, 3: 51584, 5: 5920, 7: 512, 9: 16},
}

spectral_checks = 0
for side, expected in expected_winner_masses.items():
    dimension = side * side
    size = 1 << dimension
    full_mask = size - 1
    values = winner_tables[side]
    masses = level_masses(values)
    assert masses == expected
    assert sum(masses.values()) == 4**dimension
    assert masses.get(0, 0) == 0

    odd_projection = [
        Fraction(values[mask] - values[full_mask ^ mask], 2)
        for mask in range(size)
    ]
    odd_masses = level_masses(odd_projection)
    assert odd_masses == expected_odd_masses[side]
    assert all(level % 2 for level in odd_masses)
    spectral_checks += len(masses) + len(odd_masses)


def refresh_function(values: list[Fraction], dimension: int) -> list[Fraction]:
    output = []
    for state in range(1 << dimension):
        total = Fraction(values[state], 2)
        total += sum(
            values[state ^ (1 << coordinate)]
            for coordinate in range(dimension)
        ) / (2 * dimension)
        output.append(total)
    return output


# Directly verify the stationary autocorrelation formula for small boards.
mixing_checks = 0
for side in (1, 2, 3):
    dimension = side * side
    size = 1 << dimension
    winner_values = [winner(mask, side) for mask in range(size)]
    masses = expected_winner_masses[side]
    evolved = [Fraction(value) for value in winner_values]
    for time in range(8):
        direct = sum(
            Fraction(winner_values[state]) * evolved[state]
            for state in range(size)
        ) / size
        spectral = sum(
            Fraction(mass, 4**dimension)
            * Fraction(dimension - level, dimension) ** time
            for level, mass in masses.items()
        )
        assert direct == spectral
        disagreement = (1 - direct) / 2
        assert 0 <= disagreement <= 1
        evolved = refresh_function(evolved, dimension)
        mixing_checks += 1


def majority(first: int, second: int, third: int) -> int:
    return 1 if first + second + third > 0 else -1


# Verify the exact majority polynomial and its noise-stability map by direct
# enumeration of every coordinatewise correlated sign coupling.
majority_checks = 0
for first in (-1, 1):
    for second in (-1, 1):
        for third in (-1, 1):
            assert majority(first, second, third) == (
                first + second + third - first * second * third
            ) // 2
            majority_checks += 1

for rho in (
    Fraction(-1),
    Fraction(-3, 4),
    Fraction(-1, 3),
    Fraction(0),
    Fraction(1, 5),
    Fraction(2, 3),
    Fraction(1),
):
    direct = Fraction(0)
    for x_mask in range(8):
        x = tuple(1 if x_mask & (1 << j) else -1 for j in range(3))
        for y_mask in range(8):
            y = tuple(1 if y_mask & (1 << j) else -1 for j in range(3))
            probability = Fraction(1, 8)
            for xx, yy in zip(x, y, strict=True):
                probability *= (1 + rho * xx * yy) / 2
            direct += probability * majority(*x) * majority(*y)
    assert direct == (3 * rho + rho**3) / 4
    majority_checks += 1


print(
    f"verified {topology_checks} Hex crossing/symmetry facts, "
    f"{influence_checks} exact influence values, "
    f"{spectral_checks} full/odd spectral masses, and "
    f"{mixing_checks} exact winner autocorrelations, plus "
    f"{majority_checks} majority-renormalization identities"
)
