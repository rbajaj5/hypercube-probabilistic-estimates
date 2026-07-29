"""Exact checks for notes/AffineSpectrumHashing.md.

Run with:
    python verification/verify_affine_spectrum_hashing.py

All probabilities use fractions; there is no floating-point tolerance.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


def parity(x: int) -> int:
    return x.bit_count() & 1


def apply_map(rows: tuple[int, ...], x: int) -> int:
    """Apply the binary matrix represented by its row masks."""
    out = 0
    for i, row in enumerate(rows):
        out |= parity(row & x) << i
    return out


def all_maps(n: int, t: int):
    """Enumerate all linear maps F_2^n -> F_2^t as tuples of row masks."""
    total = 1 << (n * t)
    mask = (1 << n) - 1
    for code in range(total):
        yield tuple((code >> (i * n)) & mask for i in range(t))


def is_injective_on(rows: tuple[int, ...], points: tuple[int, ...]) -> bool:
    images = {apply_map(rows, x) for x in points}
    return len(images) == len(points)


def rank_probability_formula(r: int, t: int) -> Fraction:
    if t < r:
        return Fraction(0)
    ans = Fraction(1)
    for j in range(r):
        ans *= Fraction((1 << t) - (1 << j), 1 << t)
    return ans


def exact_injection_probability(points: tuple[int, ...], n: int, t: int) -> Fraction:
    good = 0
    total = 0
    for rows in all_maps(n, t):
        total += 1
        good += is_injective_on(rows, points)
    return Fraction(good, total)


def difference_set(points: tuple[int, ...]) -> set[int]:
    return {a ^ b for a in points for b in points}


def binary_rank(vectors: tuple[int, ...]) -> int:
    """Rank of binary vectors represented by integers."""
    basis: dict[int, int] = {}
    for vector in vectors:
        x = vector
        while x:
            pivot = x.bit_length() - 1
            if pivot in basis:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                break
    return len(basis)


def characteristic_polynomial_value(elements: tuple[int, ...], q: int) -> int:
    """Evaluate the vector-matroid characteristic polynomial at q."""
    rho = binary_rank(elements)
    total = 0
    for mask in range(1 << len(elements)):
        subset = tuple(
            elements[i] for i in range(len(elements)) if (mask >> i) & 1
        )
        term = q ** (rho - binary_rank(subset))
        total += -term if mask.bit_count() & 1 else term
    return total


def matroid_injection_probability(points: tuple[int, ...], t: int) -> Fraction:
    elements = tuple(sorted(difference_set(points) - {0}))
    rho = binary_rank(elements)
    q = 1 << t
    return Fraction(characteristic_polynomial_value(elements, q), q**rho)


def check_affine_formula() -> int:
    checks = 0
    # It suffices to use the standard r-dimensional subspace. Translation and
    # ambient zero coordinates cannot change collision probabilities.
    for r in range(0, 5):
        points = tuple(range(1 << r))
        for t in range(0, 6):
            exact = exact_injection_probability(points, r, t)
            expected = rank_probability_formula(r, t)
            matroid = matroid_injection_probability(points, t)
            assert exact == expected, (r, t, exact, expected)
            assert exact == matroid, (r, t, exact, matroid)
            if t >= r:
                failure = 1 - exact
                assert failure <= Fraction((1 << r) - 1, 1 << t)
                assert failure < Fraction(1 << r, 1 << t)
            checks += 1
    return checks


def check_difference_set_bound() -> int:
    checks = 0
    # Exhaust every nonempty subset of F_2^n for n <= 3 and every t <= 3.
    for n in range(0, 4):
        universe = tuple(range(1 << n))
        for size in range(1, len(universe) + 1):
            for points in combinations(universe, size):
                distinct_differences = len(difference_set(points)) - 1
                for t in range(0, 4):
                    exact_success = exact_injection_probability(points, n, t)
                    matroid_success = matroid_injection_probability(points, t)
                    assert exact_success == matroid_success, (
                        n,
                        points,
                        t,
                        exact_success,
                        matroid_success,
                    )
                    exact_failure = 1 - exact_success
                    bound = Fraction(distinct_differences, 1 << t)
                    assert exact_failure <= bound, (
                        n,
                        points,
                        t,
                        exact_failure,
                        bound,
                    )
                    checks += 1
    return checks


def walsh_support_matching_phase(k: int, delta: int) -> set[int]:
    """Compute the Walsh support by exact summation.

    Coordinates are x_0,y_0,...,x_{k-1},y_{k-1},z.
    """
    n = 2 * k + 1
    support: set[int] = set()
    for frequency in range(1 << n):
        coefficient_numerator = 0
        for point in range(1 << n):
            quadratic = 0
            for j in range(k):
                x = (point >> (2 * j)) & 1
                y = (point >> (2 * j + 1)) & 1
                quadratic ^= x & y
            z = (point >> (2 * k)) & 1
            exponent = quadratic ^ (delta & z) ^ parity(frequency & point)
            coefficient_numerator += -1 if exponent else 1
        if coefficient_numerator:
            support.add(frequency)
    return support


def expected_matching_support(k: int, delta: int) -> set[int]:
    offset = delta << (2 * k)
    return {offset | low for low in range(1 << (2 * k))}


def check_matching_spectra() -> int:
    checks = 0
    for k in range(0, 4):
        for delta in (0, 1):
            actual = walsh_support_matching_phase(k, delta)
            expected = expected_matching_support(k, delta)
            assert actual == expected
            assert len(actual) == 4**k
            checks += 1
    return checks


def main() -> None:
    affine_checks = check_affine_formula()
    difference_checks = check_difference_set_bound()
    spectrum_checks = check_matching_spectra()
    total = affine_checks + difference_checks + spectrum_checks
    print(f"affine formula checks: {affine_checks}")
    print(f"difference-set and matroid-law checks: {difference_checks}")
    print(f"matching-spectrum checks: {spectrum_checks}")
    print(f"all {total} exact checks passed")


if __name__ == "__main__":
    main()
