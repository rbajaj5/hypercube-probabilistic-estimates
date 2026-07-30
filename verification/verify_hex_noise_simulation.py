"""Audit the fixed seeded Hex noise/pivotality simulation dataset."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


DATA = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "hex_noise_seed_20260730.csv"
)
EXPECTED_SIDES = (5, 7, 9, 11, 15, 21, 31)
EXPECTED_RHOS = (0.0, 0.25, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0)
EXPECTED_SEED = 20_260_730
EXPECTED_SAMPLES = 65_536


def close(left: float, right: float, tolerance: float = 5e-12) -> bool:
    return abs(left - right) <= tolerance


with DATA.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))

assert len(rows) == len(EXPECTED_SIDES) * len(EXPECTED_RHOS)
groups = defaultdict(list)
arithmetic_checks = 0
for row in rows:
    side = int(row["board_side"])
    cells = int(row["cells"])
    samples = int(row["samples"])
    rho = float(row["rho"])
    blue_wins = int(row["blue_wins"])
    pivots = int(row["pivotal_samples"])
    disagreements = int(row["noise_disagreements"])

    assert int(row["seed"]) == EXPECTED_SEED
    assert samples == EXPECTED_SAMPLES
    assert cells == side * side
    assert row["device_type"] in {"cpu", "cuda"}
    assert row["torch_version"]
    assert 0 <= blue_wins <= samples
    assert 0 <= pivots <= samples
    assert 0 <= disagreements <= samples

    blue_probability = blue_wins / samples
    blue_se = math.sqrt(blue_probability * (1 - blue_probability) / samples)
    pivotal_probability = pivots / samples
    pivotal_se = math.sqrt(
        pivotal_probability * (1 - pivotal_probability) / samples
    )
    disagreement = disagreements / samples
    disagreement_se = math.sqrt(disagreement * (1 - disagreement) / samples)

    assert close(float(row["blue_win_probability"]), blue_probability)
    assert close(float(row["blue_win_se"]), blue_se)
    assert close(
        float(row["random_coordinate_pivotal_probability"]),
        pivotal_probability,
    )
    assert close(
        float(row["total_influence_estimate"]),
        cells * pivotal_probability,
    )
    assert close(float(row["total_influence_se"]), cells * pivotal_se)
    assert close(float(row["noise_disagreement_probability"]), disagreement)
    assert close(float(row["noise_disagreement_se"]), disagreement_se)
    assert close(float(row["winner_correlation"]), 1 - 2 * disagreement)
    assert close(float(row["winner_correlation_se"]), 2 * disagreement_se)
    arithmetic_checks += 8
    groups[side].append(row)


statistical_checks = 0
influences = []
for side in EXPECTED_SIDES:
    group = sorted(groups[side], key=lambda item: float(item["rho"]))
    assert tuple(float(row["rho"]) for row in group) == EXPECTED_RHOS
    assert len({row["blue_wins"] for row in group}) == 1
    assert len({row["pivotal_samples"] for row in group}) == 1
    assert len({row["side_runtime_seconds"] for row in group}) == 1

    blue_probability = float(group[0]["blue_win_probability"])
    blue_se = float(group[0]["blue_win_se"])
    assert abs(blue_probability - 0.5) <= 6 * blue_se

    independent = group[0]
    independent_disagreement = float(
        independent["noise_disagreement_probability"]
    )
    independent_se = float(independent["noise_disagreement_se"])
    assert abs(independent_disagreement - 0.5) <= 6 * independent_se

    identical = group[-1]
    assert int(identical["noise_disagreements"]) == 0
    assert float(identical["noise_disagreement_probability"]) == 0
    assert float(identical["winner_correlation"]) == 1

    # Sampling noise can very slightly disturb monotonicity, so use a
    # four-standard-error tolerance for adjacent estimates.
    for left, right in zip(group[:-1], group[1:], strict=True):
        left_p = float(left["noise_disagreement_probability"])
        right_p = float(right["noise_disagreement_probability"])
        combined_se = math.hypot(
            float(left["noise_disagreement_se"]),
            float(right["noise_disagreement_se"]),
        )
        assert right_p <= left_p + 4 * combined_se

    influence = float(group[0]["total_influence_estimate"])
    influence_se = float(group[0]["total_influence_se"])
    near_one = next(row for row in group if float(row["rho"]) == 0.99)
    derivative_proxy = (
        2 * float(near_one["noise_disagreement_probability"]) / 0.01
    )
    assert derivative_proxy > 0
    influences.append((side, influence, influence_se, derivative_proxy))
    statistical_checks += 5 + len(group) - 1


log_sides = [math.log(side) for side, *_ in influences]
log_influences = [math.log(influence) for _, influence, *_ in influences]
mean_x = sum(log_sides) / len(log_sides)
mean_y = sum(log_influences) / len(log_influences)
slope = sum(
    (x - mean_x) * (y - mean_y)
    for x, y in zip(log_sides, log_influences, strict=True)
) / sum((x - mean_x) ** 2 for x in log_sides)
assert 0 < slope < 2

print(
    f"verified {len(rows)} seeded Hex simulation rows, "
    f"{arithmetic_checks} arithmetic identities, and "
    f"{statistical_checks} symmetry/noise diagnostics; "
    f"log-log influence slope {slope:.4f}"
)
