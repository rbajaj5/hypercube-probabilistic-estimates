"""Verify the exact finite claims in the Collatz parity-hypercube note."""

from __future__ import annotations

import cmath
import itertools
import math
from collections import Counter, deque


def collatz_step(n: int) -> int:
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def parity_vector(n: int, horizon: int) -> tuple[int, ...]:
    bits = []
    for _ in range(horizon):
        bits.append(n % 2)
        n = collatz_step(n)
    return tuple(bits)


def affine_data(bits: tuple[int, ...]) -> tuple[int, int]:
    ones = 0
    correction = 0
    for j, bit in enumerate(bits):
        correction = (3**bit) * correction + bit * (2**j)
        ones += bit
    return ones, correction


# Parity vectors are exactly the Boolean cube, and the affine identity holds.
for horizon in range(1, 11):
    seen = set()
    for seed in range(2**horizon):
        bits = parity_vector(seed, horizon)
        seen.add(bits)
        ones, correction = affine_data(bits)
        iterate = seed
        for _ in range(horizon):
            iterate = collatz_step(iterate)
        assert (2**horizon) * iterate == (3**ones) * seed + correction
    assert len(seen) == 2**horizon


# Exhaustively check Kolmogorov's maximal inequality for these finite cubes.
maximal_checks = 0
for horizon in range(1, 13):
    paths = []
    for bits in itertools.product((0, 1), repeat=horizon):
        walk = 0
        maximum = 0
        for bit in bits:
            walk += 2 * bit - 1
            maximum = max(maximum, abs(walk))
        paths.append(maximum)
    for threshold in range(1, horizon + 1):
        exceptional = sum(value >= threshold for value in paths)
        assert exceptional * threshold**2 <= (2**horizon) * horizon
        maximal_checks += 1


def normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def rademacher_kolmogorov_distance(horizon: int) -> float:
    cumulative = 0.0
    distance = 0.0
    for ones in range(horizon + 1):
        location = (2 * ones - horizon) / math.sqrt(horizon)
        distance = max(distance, abs(cumulative - normal_cdf(location)))
        cumulative += math.comb(horizon, ones) / (2**horizon)
        distance = max(distance, abs(cumulative - normal_cdf(location)))
    return distance


# Numerical audit of the exact binomial CDF.  This is not used as a proof of
# the Berry--Esseen theorem; it records its 1/sqrt(m) scale in this model.
berry_esseen_checks = 0
for horizon in range(1, 501):
    distance = rademacher_kolmogorov_distance(horizon)
    assert distance <= 0.4 / math.sqrt(horizon)
    berry_esseen_checks += 1


# Check the exact fixed-frequency circle characteristic function.
base = 10.0
alpha = math.log(3.0, base)
beta = math.log(2.0, base)
circle_checks = 0
for horizon in range(1, 11):
    for frequency in range(1, 6):
        direct = sum(
            cmath.exp(
                2j
                * math.pi
                * frequency
                * (sum(bits) * alpha - horizon * beta)
            )
            for bits in itertools.product((0, 1), repeat=horizon)
        ) / (2**horizon)
        formula = cmath.exp(
            -2j * math.pi * frequency * horizon * beta
        ) * (
            (1.0 + cmath.exp(2j * math.pi * frequency * alpha)) / 2.0
        ) ** horizon
        assert abs(direct - formula) < 2e-13
        circle_checks += 1


def residue_counts(horizon: int, prime: int) -> list[int]:
    counts = [0] * prime
    counts[0] = 1
    for j in range(horizon):
        updated = [0] * prime
        shift = pow(2, j, prime)
        for value, count in enumerate(counts):
            updated[value] += count
            updated[(3 * value + shift) % prime] += count
        counts = updated
    return counts


# Prime 3 has the exact last-selected-index law in equations (21)--(22).
prime_three_checks = 0
for horizon in range(1, 21):
    counts = residue_counts(horizon, 3)
    if horizon % 2 == 0:
        expected = [
            1,
            (2**horizon - 1) // 3,
            2 * (2**horizon - 1) // 3,
        ]
    else:
        expected = [
            1,
            (2 ** (horizon + 1) - 1) // 3,
            (2**horizon - 2) // 3,
        ]
    assert counts == expected
    prime_three_checks += 1


def multiplicative_order(value: int, prime: int) -> int:
    state = 1
    for exponent in range(1, prime):
        state = state * value % prime
        if state == 1:
            return exponent
    raise AssertionError((value, prime))


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, math.isqrt(value) + 1))


# Verify transitivity of the block generators for a representative exact
# range of odd nonsingular primes.
mixing_checks = 0
for prime in range(5, 100):
    if not is_prime(prime):
        continue
    period = multiplicative_order(2, prime)
    generators = [
        lambda x, shift=pow(2, j, prime): (3 * x + shift) % prime
        for j in range(period)
    ]
    reached = {0}
    queue = deque([0])
    while queue:
        value = queue.popleft()
        for generator in generators:
            image = generator(value)
            if image not in reached:
                reached.add(image)
                queue.append(image)
    assert len(reached) == prime
    mixing_checks += 1


# Verify the finite-field Fourier recursion against direct residue counts.
fourier_checks = 0
for prime in (3, 5, 7, 11, 13):
    transform = [1.0 + 0.0j] * prime
    for j in range(8):
        updated = [
            (
                transform[frequency]
                + cmath.exp(
                    2j * math.pi * frequency * pow(2, j, prime) / prime
                )
                * transform[(3 * frequency) % prime]
            )
            / 2.0
            for frequency in range(prime)
        ]
        transform = updated
        counts = residue_counts(j + 1, prime)
        direct = [
            sum(
                count * cmath.exp(2j * math.pi * frequency * value / prime)
                for value, count in enumerate(counts)
            )
            / (2 ** (j + 1))
            for frequency in range(prime)
        ]
        assert max(
            abs(actual - expected)
            for actual, expected in zip(transform, direct, strict=True)
        ) < 3e-13
        fourier_checks += prime


print(
    "verified parity bijections and affine iterates through horizon 10; "
    f"{maximal_checks} maximal bounds; "
    f"{berry_esseen_checks} binomial-normal audits; "
    f"{circle_checks} circle Fourier identities; "
    f"{prime_three_checks} prime-3 laws; "
    f"{mixing_checks} nonsingular-prime transitivity checks; and "
    f"{fourier_checks} finite-field Fourier coefficients"
)
