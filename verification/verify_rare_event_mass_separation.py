"""Verify exact rare-event mass separation on finite Boolean cubes."""

from __future__ import annotations

from fractions import Fraction
from math import comb


Q = Fraction


def character(state: int, subset: int) -> int:
    return -1 if (state & subset).bit_count() % 2 else 1


def in_event(state: int, constrained: int, pattern: int) -> bool:
    return (state & constrained) == (pattern & constrained)


def verify_case(
    dimension: int,
    constrained: int,
    pattern: int,
) -> int:
    state_count = 1 << dimension
    codimension = constrained.bit_count()
    uniform_mass = Q(1, state_count)
    event = tuple(
        in_event(state, constrained, pattern)
        for state in range(state_count)
    )
    density = tuple(Q(1 << codimension) if hit else Q(0) for hit in event)
    checks = 0

    support_size = sum(event)
    assert support_size == 1 << (dimension - codimension)
    assert Q(support_size, state_count) == Q(1, 1 << codimension)
    assert sum((uniform_mass * value for value in density), Q(0)) == 1
    checks += 3

    for order in range(1, 6):
        observed = sum(
            (uniform_mass * value**order for value in density),
            Q(0),
        )
        assert observed == Q(1 << (codimension * (order - 1)))
        checks += 1

    mean = Q(1)
    variance = sum(
        (uniform_mass * (value - mean) ** 2 for value in density),
        Q(0),
    )
    assert variance == (1 << codimension) - 1
    checks += 1

    tilted = tuple(uniform_mass * value for value in density)
    total_variation = sum(
        (abs(mass - uniform_mass) for mass in tilted),
        Q(0),
    ) / 2
    assert total_variation == 1 - Q(1, 1 << codimension)
    assert sum(mass > 0 for mass in tilted) == support_size
    checks += 2

    total_influence = Q(0)
    for coordinate in range(dimension):
        influence = sum(
            (
                uniform_mass
                for state in range(state_count)
                if event[state] != event[state ^ (1 << coordinate)]
            ),
            Q(0),
        )
        expected = (
            Q(1, 1 << (codimension - 1))
            if constrained & (1 << coordinate)
            else Q(0)
        )
        assert influence == expected
        total_influence += influence
        checks += 1
    assert total_influence == Q(codimension, 1 << (codimension - 1))
    checks += 1

    squared_fourier_mass = Q(0)
    for subset in range(state_count):
        coefficient = sum(
            (
                uniform_mass * density[state] * character(state, subset)
                for state in range(state_count)
            ),
            Q(0),
        )
        if subset & ~constrained:
            expected = Q(0)
        else:
            expected = Q(character(pattern, subset))
        assert coefficient == expected
        squared_fourier_mass += coefficient**2
        checks += 1
    assert squared_fourier_mass == sum(
        (uniform_mass * value**2 for value in density),
        Q(0),
    )
    checks += 1

    for flip_rate in (Q(0), Q(1, 10), Q(1, 4), Q(1, 2)):
        joint = Q(0)
        density_correlation = Q(0)
        for source in range(state_count):
            for target in range(state_count):
                distance = (source ^ target).bit_count()
                transition = (
                    flip_rate**distance
                    * (1 - flip_rate) ** (dimension - distance)
                )
                pair_mass = uniform_mass * transition
                if event[source] and event[target]:
                    joint += pair_mass
                density_correlation += (
                    pair_mass * density[source] * density[target]
                )
        assert joint == Q(1 - flip_rate, 2) ** codimension
        assert density_correlation == (2 * (1 - flip_rate)) ** codimension
        assert joint / Q(1, 1 << codimension) == (
            1 - flip_rate
        ) ** codimension
        checks += 3

    for sample_count in range(7):
        hit_probability = Q(1, 1 << codimension)
        hit_count_law = tuple(
            Q(comb(sample_count, hits))
            * hit_probability**hits
            * (1 - hit_probability) ** (sample_count - hits)
            for hits in range(sample_count + 1)
        )
        miss_probability = (
            1 - hit_probability
        ) ** sample_count
        assert sum(hit_count_law, Q(0)) == 1
        assert hit_count_law[0] == miss_probability
        checks += 2

    if codimension % 2 == 0:
        stage = codimension // 2
        assert 1 << codimension == 4**stage
        assert Q(1, 1 << codimension) == Q(1, 4**stage)
        assert Q(1, 1 << stage) ** 2 == Q(1, 4**stage)
        checks += 3

    return checks


total_checks = 0
case_count = 0
for dimension in range(1, 8):
    for codimension in range(1, dimension + 1):
        constrained = (1 << codimension) - 1
        pattern = sum(
            1 << coordinate
            for coordinate in range(codimension)
            if coordinate % 2
        )
        total_checks += verify_case(dimension, constrained, pattern)
        case_count += 1

print(
    f"verified {total_checks} exact normalization, rare-event, influence, "
    f"Fourier, noise, and sampling identities across {case_count} subcubes"
)
