"""Verify fixed-drawing plane-saturation laws on Boolean edge cubes."""

from __future__ import annotations

from fractions import Fraction
from functools import cache
from itertools import combinations


Q = Fraction
Polynomial = tuple[int, ...]


def add_edge(adjacency: list[int], left: int, right: int) -> None:
    adjacency[left] |= 1 << right
    adjacency[right] |= 1 << left


def conflict_graph(size: int, conflicts: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    adjacency = [0] * size
    for left, right in conflicts:
        assert 0 <= left < right < size
        add_edge(adjacency, left, right)
    return tuple(adjacency)


def chords_cross(
    first: tuple[int, int],
    second: tuple[int, int],
    vertex_count: int,
) -> bool:
    a, b = sorted(first)
    c, d = sorted(second)
    if len({a, b, c, d}) < 4:
        return False

    def between(start: int, point: int, end: int) -> bool:
        return (point - start) % vertex_count < (end - start) % vertex_count

    return between(a, c, b) != between(a, d, b)


def convex_complete_conflict_graph(vertex_count: int) -> tuple[int, ...]:
    chords = tuple(combinations(range(vertex_count), 2))
    conflicts = tuple(
        (left, right)
        for left, right in combinations(range(len(chords)), 2)
        if chords_cross(chords[left], chords[right], vertex_count)
    )
    return conflict_graph(len(chords), conflicts)


def is_independent(mask: int, adjacency: tuple[int, ...]) -> bool:
    remaining = mask
    while remaining:
        bit = remaining & -remaining
        vertex = bit.bit_length() - 1
        if adjacency[vertex] & mask:
            return False
        remaining ^= bit
    return True


def is_maximal_independent(mask: int, adjacency: tuple[int, ...]) -> bool:
    if not is_independent(mask, adjacency):
        return False
    full = (1 << len(adjacency)) - 1
    missing = full ^ mask
    while missing:
        bit = missing & -missing
        vertex = bit.bit_length() - 1
        if not adjacency[vertex] & mask:
            return False
        missing ^= bit
    return True


def planarity_indicator(mask: int, adjacency: tuple[int, ...]) -> int:
    output = 1
    for left in range(len(adjacency)):
        for right in range(left + 1, len(adjacency)):
            if adjacency[left] & (1 << right):
                output *= 1 - (
                    ((mask >> left) & 1) * ((mask >> right) & 1)
                )
    return output


def saturation_indicator(mask: int, adjacency: tuple[int, ...]) -> int:
    output = planarity_indicator(mask, adjacency)
    for vertex, neighbors in enumerate(adjacency):
        selected = (mask >> vertex) & 1
        none_of_neighbors = 1
        remaining = neighbors
        while remaining:
            bit = remaining & -remaining
            neighbor = bit.bit_length() - 1
            none_of_neighbors *= 1 - ((mask >> neighbor) & 1)
            remaining ^= bit
        output *= selected + (1 - selected) * (1 - none_of_neighbors)
    return output


def count_polynomials(
    adjacency: tuple[int, ...],
) -> tuple[Polynomial, Polynomial]:
    size = len(adjacency)
    independent = [0] * (size + 1)
    maximal = [0] * (size + 1)
    for mask in range(1 << size):
        weight = mask.bit_count()
        if is_independent(mask, adjacency):
            independent[weight] += 1
        if is_maximal_independent(mask, adjacency):
            maximal[weight] += 1
    return tuple(independent), tuple(maximal)


def evaluate(coefficients: Polynomial, value: Fraction) -> Fraction:
    output = Q(0)
    power = Q(1)
    for coefficient in coefficients:
        output += coefficient * power
        power *= value
    return output


def derivative_evaluate(
    coefficients: Polynomial,
    value: Fraction,
) -> Fraction:
    return sum(
        (
            Q(degree * coefficient) * value ** (degree - 1)
            for degree, coefficient in enumerate(coefficients)
            if degree
        ),
        Q(0),
    )


def convolve(left: Polynomial, right: Polynomial) -> Polynomial:
    output = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            output[left_degree + right_degree] += left_value * right_value
    return tuple(output)


def greedy_distribution(adjacency: tuple[int, ...]) -> tuple[Fraction, ...]:
    size = len(adjacency)

    @cache
    def recurse(vertices_mask: int) -> tuple[Fraction, ...]:
        if vertices_mask == 0:
            return (Q(1),)
        vertex_count = vertices_mask.bit_count()
        output = [Q(0)] * (vertex_count + 1)
        remaining = vertices_mask
        while remaining:
            bit = remaining & -remaining
            vertex = bit.bit_length() - 1
            next_mask = vertices_mask & ~(adjacency[vertex] | bit)
            child = recurse(next_mask)
            for child_size, probability in enumerate(child):
                output[child_size + 1] += probability / vertex_count
            remaining ^= bit
        return tuple(output)

    return recurse((1 << size) - 1)


def components(adjacency: tuple[int, ...]) -> tuple[int, ...]:
    unseen = (1 << len(adjacency)) - 1
    output = []
    while unseen:
        seed = unseen & -unseen
        component = 0
        frontier = seed
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            component |= bit
            frontier |= adjacency[vertex] & unseen & ~component
        unseen &= ~component
        output.append(component)
    return tuple(output)


def compressed_component(
    component_mask: int,
    adjacency: tuple[int, ...],
) -> tuple[int, ...]:
    vertices = [
        vertex
        for vertex in range(len(adjacency))
        if component_mask & (1 << vertex)
    ]
    positions = {vertex: index for index, vertex in enumerate(vertices)}
    output = [0] * len(vertices)
    for vertex in vertices:
        for neighbor in vertices:
            if adjacency[vertex] & (1 << neighbor):
                output[positions[vertex]] |= 1 << positions[neighbor]
    return tuple(output)


def verify_graph(name: str, adjacency: tuple[int, ...]) -> int:
    size = len(adjacency)
    independent_polynomial, maximal_polynomial = count_polynomials(adjacency)
    checks = 0

    for mask in range(1 << size):
        assert planarity_indicator(mask, adjacency) == int(
            is_independent(mask, adjacency)
        )
        assert saturation_indicator(mask, adjacency) == int(
            is_maximal_independent(mask, adjacency)
        )
        checks += 2

    for probability in (Q(1, 4), Q(1, 3), Q(1, 2), Q(2, 3)):
        activity = probability / (1 - probability)
        direct_plane = Q(0)
        direct_saturated = Q(0)
        direct_weighted_size = Q(0)
        for mask in range(1 << size):
            mass = (
                probability ** mask.bit_count()
                * (1 - probability) ** (size - mask.bit_count())
            )
            if is_independent(mask, adjacency):
                direct_plane += mass
                direct_weighted_size += mask.bit_count() * mass
            if is_maximal_independent(mask, adjacency):
                direct_saturated += mass

        partition = evaluate(independent_polynomial, activity)
        maximal_partition = evaluate(maximal_polynomial, activity)
        assert direct_plane == (1 - probability) ** size * partition
        assert (
            direct_saturated
            == (1 - probability) ** size * maximal_partition
        )
        assert direct_saturated / direct_plane == (
            maximal_partition / partition
        )
        expected_size = (
            activity
            * derivative_evaluate(independent_polynomial, activity)
            / partition
        )
        assert direct_weighted_size / direct_plane == expected_size
        checks += 4

    conflict_edges = sum(mask.bit_count() for mask in adjacency) // 2
    maximum_degree = max((mask.bit_count() for mask in adjacency), default=0)
    every_maximal_size = [
        mask.bit_count()
        for mask in range(1 << size)
        if is_maximal_independent(mask, adjacency)
    ]
    if size:
        assert every_maximal_size
        assert all(
            Q(value) >= Q(size, maximum_degree + 1)
            for value in every_maximal_size
        )
        checks += len(every_maximal_size)

    greedy = greedy_distribution(adjacency)
    assert sum(greedy, Q(0)) == 1
    greedy_expected = sum(
        (
            Q(output_size) * probability
            for output_size, probability in enumerate(greedy)
        ),
        Q(0),
    )
    degree_bound = sum(
        (Q(1, neighbors.bit_count() + 1) for neighbors in adjacency),
        Q(0),
    )
    density_bound = (
        Q(size * size, size + 2 * conflict_edges)
        if size
        else Q(0)
    )
    maximum_degree_bound = (
        Q(size, maximum_degree + 1)
        if size
        else Q(0)
    )
    assert greedy_expected >= degree_bound
    assert degree_bound >= density_bound
    assert density_bound >= maximum_degree_bound
    assert all(
        probability == 0 or output_size in every_maximal_size
        for output_size, probability in enumerate(greedy)
    )
    checks += len(greedy) + 4

    component_masks = components(adjacency)
    component_independent = (1,)
    component_maximal = (1,)
    component_greedy: tuple[Fraction, ...] = (Q(1),)
    for component_mask in component_masks:
        component_graph = compressed_component(component_mask, adjacency)
        local_independent, local_maximal = count_polynomials(component_graph)
        component_independent = convolve(
            component_independent,
            local_independent,
        )
        component_maximal = convolve(component_maximal, local_maximal)
        local_greedy = greedy_distribution(component_graph)
        next_greedy = [Q(0)] * (
            len(component_greedy) + len(local_greedy) - 1
        )
        for left_degree, left_probability in enumerate(component_greedy):
            for right_degree, right_probability in enumerate(local_greedy):
                next_greedy[left_degree + right_degree] += (
                    left_probability * right_probability
                )
        component_greedy = tuple(next_greedy)
    assert component_independent == independent_polynomial
    assert component_maximal == maximal_polynomial
    assert component_greedy == greedy
    checks += 3

    print(
        f"{name}: {size} candidate edges, {conflict_edges} crossings, "
        f"{sum(independent_polynomial)} plane states, "
        f"{sum(maximal_polynomial)} saturated states"
    )
    return checks


def matching_with_isolates(
    pair_count: int,
    isolate_count: int,
) -> tuple[int, ...]:
    size = 2 * pair_count + isolate_count
    conflicts = tuple(
        (2 * pair, 2 * pair + 1)
        for pair in range(pair_count)
    )
    adjacency = conflict_graph(size, conflicts)
    independent, maximal = count_polynomials(adjacency)

    expected_independent = (1,)
    expected_maximal = (1,)
    for _ in range(pair_count):
        expected_independent = convolve(expected_independent, (1, 2))
        expected_maximal = convolve(expected_maximal, (0, 2))
    for _ in range(isolate_count):
        expected_independent = convolve(expected_independent, (1, 1))
        expected_maximal = convolve(expected_maximal, (0, 1))
    expected_independent += (0,) * (size + 1 - len(expected_independent))
    expected_maximal += (0,) * (size + 1 - len(expected_maximal))
    assert independent == expected_independent
    assert maximal == expected_maximal
    return adjacency


GRAPHS = (
    ("single prescribed edge", conflict_graph(1, ())),
    ("three mutually crossing choices", conflict_graph(
        3,
        ((0, 1), (0, 2), (1, 2)),
    )),
    ("two conflict pairs plus isolates", matching_with_isolates(2, 2)),
    ("five-edge conflict path", conflict_graph(
        5,
        ((0, 1), (1, 2), (2, 3), (3, 4)),
    )),
    ("six-edge conflict star", conflict_graph(
        6,
        ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5)),
    )),
    ("convex K4", convex_complete_conflict_graph(4)),
    ("convex K5", convex_complete_conflict_graph(5)),
    ("convex K6", convex_complete_conflict_graph(6)),
)


total_checks = sum(
    verify_graph(name, adjacency)
    for name, adjacency in GRAPHS
)

# Explicit geometric identifications used in the note.
k4 = convex_complete_conflict_graph(4)
assert sorted(mask.bit_count() for mask in k4) == [0, 0, 0, 0, 1, 1]
k5 = convex_complete_conflict_graph(5)
assert sorted(mask.bit_count() for mask in k5) == [0] * 5 + [2] * 5

print(
    f"verified {total_checks} exact fixed-drawing planarity, saturation, "
    "Bernoulli, hard-core, Boolean-indicator, greedy, bound, and component "
    f"identities across {len(GRAPHS)} conflict hypercubes"
)
