"""Verify hard-core heat-bath identities on finite conflict hypercubes."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


Q = Fraction
SparseRow = dict[int, Fraction]


def add_edge(adjacency: list[int], left: int, right: int) -> None:
    adjacency[left] |= 1 << right
    adjacency[right] |= 1 << left


def conflict_graph(
    size: int,
    conflicts: tuple[tuple[int, int], ...],
) -> tuple[int, ...]:
    adjacency = [0] * size
    for left, right in conflicts:
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
    missing = ((1 << len(adjacency)) - 1) ^ mask
    while missing:
        bit = missing & -missing
        vertex = bit.bit_length() - 1
        if adjacency[vertex] & mask == 0:
            return False
        missing ^= bit
    return True


def independent_states(adjacency: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        mask
        for mask in range(1 << len(adjacency))
        if is_independent(mask, adjacency)
    )


def update_probability(
    state: int,
    vertex: int,
    adjacency: tuple[int, ...],
    occupancy_probability: Fraction,
) -> Fraction:
    if adjacency[vertex] & state:
        return Q(0)
    return occupancy_probability


def transition_rows(
    states: tuple[int, ...],
    adjacency: tuple[int, ...],
    occupancy_probability: Fraction,
) -> tuple[SparseRow, ...]:
    size = len(adjacency)
    positions = {state: index for index, state in enumerate(states)}
    rows = []
    for state in states:
        row: SparseRow = {}
        for vertex in range(size):
            bit = 1 << vertex
            probability_on = update_probability(
                state,
                vertex,
                adjacency,
                occupancy_probability,
            )
            off_state = state & ~bit
            on_state = state | bit
            row[positions[off_state]] = (
                row.get(positions[off_state], Q(0))
                + (1 - probability_on) / size
            )
            if probability_on:
                row[positions[on_state]] = (
                    row.get(positions[on_state], Q(0))
                    + probability_on / size
                )
        rows.append(row)
    return tuple(rows)


def apply_transition(
    distribution: list[Fraction],
    rows: tuple[SparseRow, ...],
) -> list[Fraction]:
    output = [Q(0)] * len(distribution)
    for source, probability in enumerate(distribution):
        for target, transition_probability in rows[source].items():
            output[target] += probability * transition_probability
    return output


def total_variation(
    left: list[Fraction],
    right: list[Fraction],
) -> Fraction:
    return sum(
        (abs(a - b) for a, b in zip(left, right, strict=True)),
        Q(0),
    ) / 2


def saturation_partition(
    states: tuple[int, ...],
    adjacency: tuple[int, ...],
    activity: Fraction,
) -> Fraction:
    return sum(
        (
            activity ** state.bit_count()
            for state in states
            if is_maximal_independent(state, adjacency)
        ),
        Q(0),
    )


def verify_chain(
    name: str,
    adjacency: tuple[int, ...],
    occupancy_probability: Fraction,
) -> int:
    size = len(adjacency)
    assert size >= 1
    assert 0 < occupancy_probability < 1
    activity = occupancy_probability / (1 - occupancy_probability)
    states = independent_states(adjacency)
    positions = {state: index for index, state in enumerate(states)}
    rows = transition_rows(states, adjacency, occupancy_probability)
    partition = sum(
        (activity ** state.bit_count() for state in states),
        Q(0),
    )
    stationary = [
        activity ** state.bit_count() / partition
        for state in states
    ]
    checks = 0

    assert all(sum(row.values(), Q(0)) == 1 for row in rows)
    assert all(
        target in row
        for source, row in enumerate(rows)
        for target in (source,)
    )
    checks += 2 * len(states)

    observed_stationary = apply_transition(stationary, rows)
    assert observed_stationary == stationary
    checks += len(states)

    for source, row in enumerate(rows):
        for target, probability in row.items():
            reverse = rows[target].get(source, Q(0))
            assert (
                stationary[source] * probability
                == stationary[target] * reverse
            )
            checks += 1

    # Connectivity follows constructively: delete to zero, then add the
    # target independent set.  Verify every required transition is positive.
    zero_position = positions[0]
    for state in states:
        current = state
        while current:
            bit = current & -current
            next_state = current ^ bit
            assert rows[positions[current]].get(
                positions[next_state],
                Q(0),
            ) > 0
            current = next_state
            checks += 1
        current = 0
        remaining = state
        while remaining:
            bit = remaining & -remaining
            next_state = current | bit
            assert rows[positions[current]].get(
                positions[next_state],
                Q(0),
            ) > 0
            current = next_state
            remaining ^= bit
            checks += 1

    stationary_mean = Q(0)
    stationary_eligible = Q(0)
    for source, state in enumerate(states):
        eligible = sum(
            1
            for vertex in range(size)
            if adjacency[vertex] & state == 0
        )
        direct_drift = sum(
            (
                probability
                * (states[target].bit_count() - state.bit_count())
                for target, probability in rows[source].items()
            ),
            Q(0),
        )
        expected_drift = (
            occupancy_probability * eligible - state.bit_count()
        ) / size
        assert direct_drift == expected_drift
        stationary_mean += stationary[source] * state.bit_count()
        stationary_eligible += stationary[source] * eligible
        checks += 1
    assert stationary_mean == occupancy_probability * stationary_eligible
    checks += 1

    # Empty-state first escape.
    assert (
        rows[zero_position][zero_position]
        == 1 - occupancy_probability
    )
    for time in range(9):
        assert (
            rows[zero_position][zero_position] ** time
            == (1 - occupancy_probability) ** time
        )
        checks += 1

    maximum_degree = max(neighbors.bit_count() for neighbors in adjacency)
    contraction = 1 - (
        1 - maximum_degree * occupancy_probability
    ) / size

    # Exact common-uniform coupling for every adjacent pair.
    for left in states:
        for differing_vertex in range(size):
            bit = 1 << differing_vertex
            if left & bit:
                continue
            right = left | bit
            if right not in positions:
                continue

            expected_distance = Q(0)
            for selected_vertex in range(size):
                left_probability = update_probability(
                    left,
                    selected_vertex,
                    adjacency,
                    occupancy_probability,
                )
                right_probability = update_probability(
                    right,
                    selected_vertex,
                    adjacency,
                    occupancy_probability,
                )
                distance_without_selected = (
                    (left & ~(1 << selected_vertex))
                    ^ (right & ~(1 << selected_vertex))
                ).bit_count()
                expected_distance += (
                    distance_without_selected
                    + abs(left_probability - right_probability)
                ) / size
            assert expected_distance <= contraction
            checks += 1

    # Exact finite-time TV checks in the strict contraction regime.
    if maximum_degree * occupancy_probability < 1:
        starts = [0, states[-1]]
        if len(states) > 2:
            starts.append(states[len(states) // 2])
        for start in dict.fromkeys(starts):
            distribution = [Q(0)] * len(states)
            distribution[positions[start]] = 1
            for time in range(11):
                observed_tv = total_variation(distribution, stationary)
                assert observed_tv <= size * contraction**time
                distribution = apply_transition(distribution, rows)
                checks += 1

    maximal_partition = saturation_partition(
        states,
        adjacency,
        activity,
    )
    stationary_saturation = sum(
        (
            stationary[index]
            for index, state in enumerate(states)
            if is_maximal_independent(state, adjacency)
        ),
        Q(0),
    )
    assert stationary_saturation == maximal_partition / partition
    checks += 1

    print(
        f"{name}: {size} coordinates, {len(states)} plane states, "
        f"Delta={maximum_degree}, rho={occupancy_probability}, "
        f"contraction={contraction}"
    )
    return checks


CHAINS = (
    (
        "four nonconflicting edges",
        conflict_graph(4, ()),
        Q(2, 3),
    ),
    (
        "one conflict pair plus isolates",
        conflict_graph(5, ((0, 1),)),
        Q(1, 2),
    ),
    (
        "five-edge conflict path",
        conflict_graph(
            5,
            ((0, 1), (1, 2), (2, 3), (3, 4)),
        ),
        Q(1, 4),
    ),
    (
        "convex K4",
        convex_complete_conflict_graph(4),
        Q(1, 2),
    ),
    (
        "convex K5",
        convex_complete_conflict_graph(5),
        Q(1, 3),
    ),
    (
        "convex K6",
        convex_complete_conflict_graph(6),
        Q(1, 6),
    ),
)


total_checks = sum(
    verify_chain(name, adjacency, occupancy_probability)
    for name, adjacency, occupancy_probability in CHAINS
)
print(
    f"verified {total_checks} exact hard-core transition, reversibility, "
    "connectivity, drift, escape, path-coupling, finite-time mixing, and "
    f"saturation identities across {len(CHAINS)} conflict hypercubes"
)
