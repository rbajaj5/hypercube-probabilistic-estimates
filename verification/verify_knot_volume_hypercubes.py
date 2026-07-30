"""Verify the fixed knot-volume sample and its induced hypercube law."""

from __future__ import annotations

import cmath
import csv
import itertools
import math
import random
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "census_knot_volumes_seed_20260730.csv"

SEED = 20260730
CENSUS_SIZE = 3116
SAMPLE_SIZE = 16
EXPECTED_BINS = [3, 2, 3, 1, 2, 0, 0, 5]
EXPECTED_ENTROPY = 2.4300365325772657
EXPECTED_FOURIER_MAGNITUDES = [
    0.33932203154673535,
    0.3769104576367034,
    0.18961048131400243,
    0.07872620347365585,
]
EXPECTED_TAIL_COEFFICIENTS = {
    0.5: 1.1374945282870066,
    1.0: 1.5503286553694133,
    1.5: 2.1488130561381107,
    2.0: 3.0,
    3.0: 5.904117239592693,
    4.0: 11.688151048065665,
}


def close(actual: float, expected: float, tolerance: float = 2e-11) -> None:
    assert math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance), (
        actual,
        expected,
    )


with DATA.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))

assert len(rows) == SAMPLE_SIZE
assert {int(row["seed"]) for row in rows} == {SEED}
assert {int(row["census_size"]) for row in rows} == {CENSUS_SIZE}
assert {int(row["sample_size"]) for row in rows} == {SAMPLE_SIZE}
assert len({row["name"] for row in rows}) == SAMPLE_SIZE

indices = [int(row["index"]) for row in rows]
expected_indices = sorted(
    random.Random(SEED).sample(range(CENSUS_SIZE), SAMPLE_SIZE)
)
assert indices == expected_indices
assert all(int(row["volume_accuracy"]) >= 10 for row in rows)
assert all(
    row["solution_type"] == "all tetrahedra positively oriented"
    for row in rows
)

phases = [float(row["volume"]) % 1.0 for row in rows]
angles = [math.cos(2.0 * math.pi * phase) for phase in phases]
angle_atoms = [
    (math.sqrt(2.0 + 2.0 * q), math.sqrt(2.0 - 2.0 * q))
    for q in angles
]


def symmetric_norm(matrix: tuple[tuple[float, float], tuple[float, float]]) -> float:
    a, b = matrix[0]
    _, d = matrix[1]
    midpoint = (a + d) / 2.0
    radius = math.hypot((a - d) / 2.0, b)
    return max(abs(midpoint + radius), abs(midpoint - radius))


for phase, (first, second) in zip(phases, angle_atoms, strict=True):
    theta = 2.0 * math.pi * phase
    q = math.cos(theta)
    sine = math.sin(theta)
    matrices = (
        ((1.0, 0.0), (0.0, 0.0)),
        ((0.0, 0.0), (0.0, 1.0)),
        (((1.0 - q) / 2.0, -sine / 2.0), (-sine / 2.0, (1.0 + q) / 2.0)),
        (((1.0 + q) / 2.0, sine / 2.0), (sine / 2.0, (1.0 - q) / 2.0)),
    )
    actual_norms = []
    for signs in itertools.product((-1.0, 1.0), repeat=4):
        matrix = tuple(
            tuple(
                sum(sign * source[row][column] for sign, source in zip(signs, matrices))
                for column in range(2)
            )
            for row in range(2)
        )
        actual_norms.append(symmetric_norm(matrix))
    expected_norms = [0.0] * 2 + [first] * 2 + [second] * 2 + [2.0] * 10
    for actual, expected in zip(
        sorted(actual_norms), sorted(expected_norms), strict=True
    ):
        close(actual, expected)

bins = [0] * 8
for phase in phases:
    bins[min(int(8 * phase), 7)] += 1
assert bins == EXPECTED_BINS

probabilities = [Fraction(count, SAMPLE_SIZE) for count in bins]
assert sum(probabilities) == 1
entropy = -sum(
    float(probability) * math.log2(float(probability))
    for probability in probabilities
    if probability
)
close(entropy, EXPECTED_ENTROPY)

fourier_magnitudes = [
    abs(
        sum(
            cmath.exp(2j * math.pi * frequency * phase)
            for phase in phases
        )
        / SAMPLE_SIZE
    )
    for frequency in range(1, 5)
]
for actual, expected in zip(
    fourier_magnitudes, EXPECTED_FOURIER_MAGNITUDES, strict=True
):
    close(actual, expected)

# The law has fixed masses 1/8 at zero and 5/8 at two. Each of the
# 2 * SAMPLE_SIZE angle atoms has mass 1/(8 * SAMPLE_SIZE).
assert Fraction(1, 8) + Fraction(5, 8) + (
    2 * SAMPLE_SIZE * Fraction(1, 8 * SAMPLE_SIZE)
) == 1
assert 16 * SAMPLE_SIZE == 256

for first, second in angle_atoms:
    close(first * first + second * second, 4.0)

mean = (5.0 / 8.0) * 2.0 + sum(
    first + second for first, second in angle_atoms
) / (8.0 * SAMPLE_SIZE)
close(mean, EXPECTED_TAIL_COEFFICIENTS[1.0])

second_moment = (5.0 / 8.0) * 4.0 + sum(
    first * first + second * second for first, second in angle_atoms
) / (8.0 * SAMPLE_SIZE)
close(second_moment, 3.0)


def tail_coefficient(alpha: float) -> float:
    return (5.0 / 8.0) * (2.0**alpha) + sum(
        first**alpha + second**alpha for first, second in angle_atoms
    ) / (8.0 * SAMPLE_SIZE)


for alpha, expected in EXPECTED_TAIL_COEFFICIENTS.items():
    close(tail_coefficient(alpha), expected)

print(
    "verified 16 seeded knot volumes, 8-bin/Fourier diagnostics, "
    "256 equally likely knot-sign configurations, and 6 tail coefficients"
)
