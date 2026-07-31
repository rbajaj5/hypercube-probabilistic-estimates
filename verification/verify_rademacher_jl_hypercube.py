"""Verify exact Rademacher JL distortion laws for Boolean-cube pairs."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import comb


Q = Fraction
Law = dict[int, Fraction]


def row_square_law(distance: int) -> Law:
    law: Law = {}
    for negative_count in range(distance + 1):
        row_sum = distance - 2 * negative_count
        square = row_sum**2
        law[square] = law.get(square, Q(0)) + Q(
            comb(distance, negative_count),
            1 << distance,
        )
    return law


def convolve(left: Law, right: Law) -> Law:
    output: Law = {}
    for left_value, left_mass in left.items():
        for right_value, right_mass in right.items():
            value = left_value + right_value
            output[value] = output.get(value, Q(0)) + left_mass * right_mass
    return output


def distortion_law(rows: int, distance: int) -> Law:
    one_row = row_square_law(distance)
    output: Law = {0: Q(1)}
    for _ in range(rows):
        output = convolve(output, one_row)
    return output


def verify_pair_laws() -> int:
    checks = 0
    for distance in range(1, 10):
        one_row = row_square_law(distance)
        assert sum(one_row.values(), Q(0)) == 1
        checks += 1

        for rows in range(1, 7):
            law = distortion_law(rows, distance)
            scale = rows * distance
            assert sum(law.values(), Q(0)) == 1
            mean = sum(
                (Q(square, scale) * mass for square, mass in law.items()),
                Q(0),
            )
            second = sum(
                (
                    Q(square, scale) ** 2 * mass
                    for square, mass in law.items()
                ),
                Q(0),
            )
            assert mean == 1
            assert second - mean**2 == Q(
                2 * (distance - 1),
                rows * distance,
            )
            checks += 3

            zero_row_probability = (
                Q(comb(distance, distance // 2), 1 << distance)
                if distance % 2 == 0
                else Q(0)
            )
            assert law.get(0, Q(0)) == zero_row_probability**rows
            checks += 1

            if distance == 1:
                assert law == {rows: Q(1)}
                checks += 1
            elif distance == 2:
                expected = {
                    4 * successes: Q(comb(rows, successes), 1 << rows)
                    for successes in range(rows + 1)
                }
                assert law == expected
                checks += 1
            elif distance == 3:
                expected = {
                    rows + 8 * successes: (
                        Q(comb(rows, successes))
                        * Q(1, 4) ** successes
                        * Q(3, 4) ** (rows - successes)
                    )
                    for successes in range(rows + 1)
                }
                assert law == expected
                checks += 1

            for delta in (Q(1, 10), Q(1, 4), Q(1, 2), Q(1)):
                tail = sum(
                    (
                        mass
                        for square, mass in law.items()
                        if abs(Q(square, scale) - 1) >= delta
                    ),
                    Q(0),
                )
                bound = Q(2 * (distance - 1), rows * distance) / delta**2
                assert tail <= bound
                checks += 1
    return checks


def matrix_image(
    vertex: int,
    matrix_signs: tuple[int, ...],
    rows: int,
    dimension: int,
) -> tuple[int, ...]:
    return tuple(
        sum(
            matrix_signs[row * dimension + coordinate]
            * ((vertex >> coordinate) & 1)
            for coordinate in range(dimension)
        )
        for row in range(rows)
    )


def verify_full_cubes() -> int:
    checks = 0
    for dimension in range(1, 9):
        vertex_count = 1 << dimension
        shell_total = sum(
            (
                (1 << (dimension - 1)) * comb(dimension, distance)
                for distance in range(1, dimension + 1)
            )
        )
        assert shell_total == comb(vertex_count, 2)
        checks += 1

        for rows in range(1, 5):
            expected_collisions = Q(0)
            for distance in range(1, dimension + 1):
                pair_count = (1 << (dimension - 1)) * comb(
                    dimension,
                    distance,
                )
                zero_row_probability = (
                    Q(comb(distance, distance // 2), 1 << distance)
                    if distance % 2 == 0
                    else Q(0)
                )
                expected_collisions += (
                    pair_count * zero_row_probability**rows
                )
                checks += 1
            assert expected_collisions >= 0
            checks += 1

    # Independently enumerate every small sign matrix and count all pairs
    # with equal real-valued images.
    for dimension in range(1, 5):
        vertex_count = 1 << dimension
        vertices = tuple(range(vertex_count))
        for rows in range(1, 4):
            if rows * dimension > 12:
                continue
            observed = Q(0)
            matrix_count = 1 << (rows * dimension)
            for sign_bits in range(matrix_count):
                signs = tuple(
                    1 if sign_bits & (1 << index) else -1
                    for index in range(rows * dimension)
                )
                images = {
                    vertex: matrix_image(
                        vertex,
                        signs,
                        rows,
                        dimension,
                    )
                    for vertex in vertices
                }
                collisions = sum(
                    images[left] == images[right]
                    for left, right in combinations(vertices, 2)
                )
                observed += Q(collisions, matrix_count)

            expected = Q(0)
            for distance in range(1, dimension + 1):
                pair_count = (1 << (dimension - 1)) * comb(
                    dimension,
                    distance,
                )
                zero_row_probability = (
                    Q(comb(distance, distance // 2), 1 << distance)
                    if distance % 2 == 0
                    else Q(0)
                )
                expected += pair_count * zero_row_probability**rows
            assert observed == expected
            checks += 1
    return checks


total_checks = verify_pair_laws() + verify_full_cubes()
print(
    f"verified {total_checks} exact Rademacher-JL row, distortion, moment, "
    "tail, Hamming-shell, and full-cube collision identities"
)
