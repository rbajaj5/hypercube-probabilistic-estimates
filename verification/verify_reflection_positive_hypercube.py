"""Verify reflection-positive hypercube moment and Hex bounds."""

from __future__ import annotations

from collections import defaultdict, deque
from fractions import Fraction
from itertools import combinations


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


def walsh_masses(values: list[int]) -> dict[int, Fraction]:
    size = len(values)
    masses: dict[int, Fraction] = defaultdict(Fraction)
    for subset, coefficient in enumerate(fwht(values)):
        masses[subset.bit_count()] += Fraction(coefficient * coefficient, size * size)
    return dict(masses)


def correlations(
    masses: dict[int, Fraction], dimension: int, count: int
) -> list[Fraction]:
    return [
        sum(
            mass * Fraction(dimension - level, dimension) ** time
            for level, mass in masses.items()
        )
        for time in range(count)
    ]


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    size = len(work)
    output = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
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


def matrix_rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        for entry in range(column, columns):
            work[pivot_row][entry] /= pivot_value
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = work[row][column]
            if factor:
                for entry in range(column, columns):
                    work[row][entry] -= factor * work[pivot_row][entry]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def principal_minors_nonnegative(matrix: list[list[Fraction]]) -> int:
    size = len(matrix)
    checks = 0
    for order in range(1, size + 1):
        for indices in combinations(range(size), order):
            minor = [[matrix[i][j] for j in indices] for i in indices]
            assert determinant(minor) >= 0
            checks += 1
    return checks


def solve(matrix: list[list[Fraction]], vector: list[Fraction]) -> list[Fraction]:
    size = len(vector)
    augmented = [
        matrix[row][:] + [vector[row]]
        for row in range(size)
    ]
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


def verify_moment_law(
    masses: dict[int, Fraction], dimension: int
) -> tuple[int, int, int, int]:
    sequence = correlations(masses, dimension, 18)
    difference_checks = 0
    current = sequence
    for order in range(8):
        sign = 1 if order % 2 == 0 else -1
        assert all(sign * value >= 0 for value in current)
        difference_checks += len(current)
        current = [
            current[index + 1] - current[index]
            for index in range(len(current) - 1)
        ]

    hankel_checks = 0
    rank_checks = 0
    active_levels = sum(mass > 0 for mass in masses.values())
    for size in range(1, 6):
        hankel = [
            [sequence[i + j] for j in range(size)]
            for i in range(size)
        ]
        shifted = [
            [sequence[i + j + 1] for j in range(size)]
            for i in range(size)
        ]
        difference = [
            [hankel[i][j] - shifted[i][j] for j in range(size)]
            for i in range(size)
        ]
        hankel_checks += principal_minors_nonnegative(hankel)
        hankel_checks += principal_minors_nonnegative(difference)
        assert matrix_rank(hankel) == min(size, active_levels)
        rank_checks += 1

    log_convex_checks = 0
    for time in range(len(sequence) - 2):
        assert sequence[time] * sequence[time + 2] >= sequence[time + 1] ** 2
        log_convex_checks += 1

    return difference_checks, hankel_checks, rank_checks, log_convex_checks


# Exhaust all Boolean observables through dimension three.
exhaustive_functions = 0
difference_checks = 0
hankel_checks = 0
rank_checks = 0
log_convex_checks = 0
bound_checks = 0
for dimension in range(1, 4):
    size = 1 << dimension
    for function_mask in range(1 << size):
        values = [
            1 if function_mask & (1 << state) else -1
            for state in range(size)
        ]
        masses = walsh_masses(values)
        assert sum(masses.values()) == 1
        counts = verify_moment_law(masses, dimension)
        difference_checks += counts[0]
        hankel_checks += counts[1]
        rank_checks += counts[2]
        log_convex_checks += counts[3]

        mean_square = masses.get(0, Fraction(0))
        variance = 1 - mean_square
        if variance:
            influence = sum(level * mass for level, mass in masses.items())
            first_active = min(level for level, mass in masses.items() if level and mass)
            centered = {
                level: mass for level, mass in masses.items() if level
            }
            centered_sequence = correlations(centered, dimension, 12)
            first = variance - influence / dimension
            assert centered_sequence[0] == variance
            assert centered_sequence[1] == first
            mean_node = first / variance
            upper_node = Fraction(dimension - first_active, dimension)
            for time in range(1, 12):
                lower = variance * mean_node**time
                upper = first * upper_node ** (time - 1)
                assert lower <= centered_sequence[time] <= upper
                bound_checks += 1
        exhaustive_functions += 1


# Small Hex boards.
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


expected_influences = {
    1: Fraction(1),
    2: Fraction(3, 2),
    3: Fraction(249, 128),
    4: Fraction(2405, 1024),
}
expected_second_correlations = {
    1: Fraction(0),
    2: Fraction(27, 64),
    3: Fraction(2177, 3456),
    4: Fraction(1544691, 2097152),
}
expected_bounds_at_two = {
    2: (Fraction(25, 64), Fraction(15, 32)),
    3: (Fraction(90601, 147456), Fraction(301, 432)),
    4: (Fraction(195412441, 268435456), Fraction(209685, 262144)),
}

hex_checks = 0
reconstruction_checks = 0
for side in range(1, 5):
    dimension = side * side
    values = [winner(mask, side) for mask in range(1 << dimension)]
    assert sum(values) == 0
    masses = walsh_masses(values)
    sequence = correlations(masses, dimension, max(dimension + 1, 3))
    influence = sum(level * mass for level, mass in masses.items())
    assert influence == expected_influences[side]
    assert sequence[1] == 1 - influence / dimension
    assert sequence[2] == expected_second_correlations[side]
    hex_checks += 4

    counts = verify_moment_law(masses, dimension)
    difference_checks += counts[0]
    hankel_checks += counts[1]
    rank_checks += counts[2]
    log_convex_checks += counts[3]

    if side in expected_bounds_at_two:
        lower, upper = expected_bounds_at_two[side]
        assert lower <= sequence[2] <= upper
        hex_checks += 1

    # Recover every level mass from C_0,...,C_d by the known-node
    # Vandermonde system.  Keep this exact solve to sides at most three;
    # side four is already independently checked above.
    if side <= 3:
        nodes = [
            Fraction(dimension - level, dimension)
            for level in range(dimension + 1)
        ]
        vandermonde = [
            [node**time for node in nodes]
            for time in range(dimension + 1)
        ]
        recovered = solve(vandermonde, sequence[: dimension + 1])
        for level, mass in enumerate(recovered):
            assert mass == masses.get(level, 0)
            reconstruction_checks += 1

    # Consecutive-lag interpolation using the exact spectral ceiling.
    ceiling = Fraction(dimension - 1, dimension)
    for time in range(1, min(8, dimension + 1)):
        if sequence[time - 1] == 0:
            continue
        ratio = sequence[time] / sequence[time - 1]
        for extra in range(5):
            target = correlations(masses, dimension, time + extra + 1)[-1]
            assert sequence[time] * ratio**extra <= target
            assert target <= sequence[time] * ceiling**extra
            hex_checks += 2


# Product-noise semigroup sampled at tau multiples with rho = 3/4.  Its
# nodes are rho^k, so the same reflected-gluing Hankel law must hold.
product_checks = 0
base_rho = Fraction(3, 4)
sample_masses = {
    0: Fraction(1, 10),
    1: Fraction(2, 10),
    3: Fraction(3, 10),
    6: Fraction(4, 10),
}
product_sequence = [
    sum(mass * base_rho ** (level * time) for level, mass in sample_masses.items())
    for time in range(12)
]
for size in range(1, 6):
    matrix = [
        [product_sequence[i + j] for j in range(size)]
        for i in range(size)
    ]
    product_checks += principal_minors_nonnegative(matrix)


print(
    f"verified {exhaustive_functions} exhaustive Boolean observables, "
    f"{difference_checks} complete-monotonicity values, "
    f"{hankel_checks} reflection-positive principal minors, "
    f"{rank_checks} exact Hankel ranks, "
    f"{log_convex_checks} log-convexity inequalities, "
    f"{bound_checks} influence-only bounds, "
    f"{hex_checks} Hex correlation/interpolation facts, "
    f"{reconstruction_checks} recovered Hex Walsh masses, and "
    f"{product_checks} product-noise principal minors"
)
