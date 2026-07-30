"""Seeded GPU/CPU Monte Carlo for Hex winner noise and pivotality.

PyTorch is an optional experiment dependency.  The committed exact checker
does not import it.  Typical use:

    python experiments/simulate_hex_noise_gpu.py --self-test
    python experiments/simulate_hex_noise_gpu.py --output data/hex_noise_seed_20260730.csv

The Boolean board tensor has shape ``(batch, side, side)`` and is true on
blue cells.  Blue wins by connecting left to right on the six-neighbour
rhombus.  Hex's no-tie theorem then determines the signed winner.
"""

from __future__ import annotations

import argparse
import csv
import math
import platform
import time
from collections import deque
from pathlib import Path

try:
    import torch
except ImportError as exc:  # pragma: no cover - depends on optional package
    raise SystemExit(
        "PyTorch is required for this experiment; install a CUDA build to "
        "use the GPU, or a CPU build for the fallback."
    ) from exc


DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1))
DEFAULT_RHOS = (0.0, 0.25, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sides",
        type=int,
        nargs="+",
        default=(5, 7, 9, 11, 15, 21, 31),
    )
    parser.add_argument("--samples", type=int, default=65_536)
    parser.add_argument("--batch-size", type=int, default=2_048)
    parser.add_argument("--seed", type=int, default=20_260_730)
    parser.add_argument("--rhos", type=float, nargs="+", default=DEFAULT_RHOS)
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="exhaustively compare the batched kernel with scalar BFS through 3x3",
    )
    return parser.parse_args()


def resolve_device(requested: str) -> torch.device:
    if requested == "cuda" and not torch.cuda.is_available():
        raise SystemExit("--device cuda requested, but CUDA is unavailable")
    if requested == "auto":
        requested = "cuda" if torch.cuda.is_available() else "cpu"
    return torch.device(requested)


def batched_crossing(boards: torch.Tensor, axis: int) -> torch.Tensor:
    """Return crossing indicators for Boolean boards.

    ``axis=0`` means a left-right connection and ``axis=1`` a top-bottom
    connection.  Each iteration expands the full frontier by one graph step.
    """

    if boards.dtype != torch.bool or boards.ndim != 3:
        raise ValueError("boards must be a Boolean tensor of shape (batch,n,n)")
    if boards.shape[1] != boards.shape[2]:
        raise ValueError("boards must be square")
    if axis not in (0, 1):
        raise ValueError("axis must be 0 or 1")

    visited = torch.zeros_like(boards)
    if axis == 0:
        visited[:, :, 0] = boards[:, :, 0]
    else:
        visited[:, 0, :] = boards[:, 0, :]

    # A simple path has at most n^2 vertices, so this cap is rigorous.
    for _ in range(boards.shape[1] ** 2):
        neighbours = torch.zeros_like(visited)
        neighbours[:, :, 1:] |= visited[:, :, :-1]
        neighbours[:, :, :-1] |= visited[:, :, 1:]
        neighbours[:, 1:, :] |= visited[:, :-1, :]
        neighbours[:, :-1, :] |= visited[:, 1:, :]
        neighbours[:, 1:, :-1] |= visited[:, :-1, 1:]
        neighbours[:, :-1, 1:] |= visited[:, 1:, :-1]
        expanded = boards & (visited | neighbours)
        if torch.equal(expanded, visited):
            break
        visited = expanded

    if axis == 0:
        return visited[:, :, -1].any(dim=1)
    return visited[:, -1, :].any(dim=1)


def batched_winner(boards: torch.Tensor) -> torch.Tensor:
    """Return the signed Hex winner (+1 blue, -1 yellow)."""

    blue = batched_crossing(boards, axis=0)
    return blue.to(torch.int8).mul_(2).sub_(1)


def scalar_crossing(board: list[list[bool]], color: bool, axis: int) -> bool:
    side = len(board)
    starts = [
        ((0, coordinate) if axis == 0 else (coordinate, 0))
        for coordinate in range(side)
    ]
    starts = [(x, y) for x, y in starts if board[y][x] is color]
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
                and board[yy][xx] is color
            ):
                seen.add((xx, yy))
                queue.append((xx, yy))
    return False


def exhaustive_self_test(device: torch.device) -> int:
    checks = 0
    for side in (1, 2, 3):
        cells = side * side
        states = torch.arange(1 << cells, dtype=torch.int64)
        bit_positions = torch.arange(cells, dtype=torch.int64)
        boards_cpu = (
            (states[:, None] >> bit_positions[None, :]) & 1
        ).bool().reshape(-1, side, side)
        boards = boards_cpu.to(device)
        blue_gpu = batched_crossing(boards, axis=0).cpu()
        yellow_gpu = batched_crossing(~boards, axis=1).cpu()
        if not torch.all(blue_gpu ^ yellow_gpu):
            raise AssertionError(f"Hex no-tie kernel check failed at side {side}")

        for state, board_tensor in enumerate(boards_cpu):
            board = board_tensor.tolist()
            blue = scalar_crossing(board, color=True, axis=0)
            yellow = scalar_crossing(board, color=False, axis=1)
            if blue != bool(blue_gpu[state]) or yellow != bool(yellow_gpu[state]):
                raise AssertionError(
                    f"batched/scalar mismatch at side {side}, state {state}"
                )
            checks += 2

        transformed = (~boards.transpose(1, 2)).contiguous()
        transformed_winner = batched_winner(transformed)
        if not torch.equal(transformed_winner, -batched_winner(boards)):
            raise AssertionError(f"twisted symmetry failed at side {side}")
        checks += 1 << cells
    return checks


def binomial_se(successes: int, trials: int) -> float:
    proportion = successes / trials
    return math.sqrt(proportion * (1.0 - proportion) / trials)


def device_name(device: torch.device) -> str:
    if device.type == "cuda":
        return torch.cuda.get_device_name(device)
    return platform.processor() or "CPU"


def simulate_side(
    side: int,
    samples: int,
    batch_size: int,
    rhos: tuple[float, ...],
    generator: torch.Generator,
    device: torch.device,
) -> list[dict[str, int | float | str]]:
    if side < 1:
        raise ValueError("board sides must be positive")
    if samples < 1 or batch_size < 1:
        raise ValueError("samples and batch size must be positive")
    if any(not 0.0 <= rho <= 1.0 for rho in rhos):
        raise ValueError("noise correlations must lie in [0,1]")

    blue_wins = 0
    pivotal = 0
    disagreements = [0] * len(rhos)
    completed = 0
    cells = side * side
    flip_probabilities = (
        (1.0 - torch.tensor(rhos, device=device, dtype=torch.float32)) / 2.0
    )[:, None, None, None]
    started = time.perf_counter()

    while completed < samples:
        current = min(batch_size, samples - completed)
        boards = torch.rand(
            (current, side, side), device=device, generator=generator
        ) < 0.5
        winners = batched_winner(boards)
        blue_wins += int((winners == 1).sum().item())

        coordinates = torch.randint(
            cells, (current,), device=device, generator=generator
        )
        flipped = boards.reshape(current, cells).clone()
        rows = torch.arange(current, device=device)
        flipped[rows, coordinates] = ~flipped[rows, coordinates]
        flipped_winners = batched_winner(flipped.reshape_as(boards))
        pivotal += int((flipped_winners != winners).sum().item())

        uniforms = torch.rand(
            (len(rhos), current, side, side),
            device=device,
            generator=generator,
        )
        coupled = boards.unsqueeze(0) ^ (uniforms < flip_probabilities)
        coupled_winners = batched_winner(
            coupled.reshape(len(rhos) * current, side, side)
        ).reshape(len(rhos), current)
        batch_disagreements = (coupled_winners != winners.unsqueeze(0)).sum(dim=1)
        for index, count in enumerate(batch_disagreements.tolist()):
            disagreements[index] += count
        completed += current

    if device.type == "cuda":
        torch.cuda.synchronize(device)
    elapsed = time.perf_counter() - started

    blue_probability = blue_wins / samples
    blue_se = binomial_se(blue_wins, samples)
    pivotal_probability = pivotal / samples
    pivotal_se = binomial_se(pivotal, samples)
    influence = cells * pivotal_probability
    influence_se = cells * pivotal_se

    rows = []
    for rho, count in zip(rhos, disagreements):
        disagreement = count / samples
        disagreement_se = binomial_se(count, samples)
        rows.append(
            {
                "board_side": side,
                "cells": cells,
                "samples": samples,
                "rho": rho,
                "blue_wins": blue_wins,
                "blue_win_probability": blue_probability,
                "blue_win_se": blue_se,
                "pivotal_samples": pivotal,
                "random_coordinate_pivotal_probability": pivotal_probability,
                "total_influence_estimate": influence,
                "total_influence_se": influence_se,
                "noise_disagreements": count,
                "noise_disagreement_probability": disagreement,
                "noise_disagreement_se": disagreement_se,
                "winner_correlation": 1.0 - 2.0 * disagreement,
                "winner_correlation_se": 2.0 * disagreement_se,
                "side_runtime_seconds": elapsed,
            }
        )
    return rows


def write_csv(
    destination: Path,
    rows: list[dict[str, int | float | str]],
    metadata: dict[str, int | str],
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "seed",
        "device_type",
        "device_name",
        "torch_version",
        "cuda_version",
        "batch_size",
        *rows[0].keys(),
    ]
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({**metadata, **row})


def main() -> None:
    args = parse_args()
    device = resolve_device(args.device)
    rhos = tuple(args.rhos)
    if len(set(rhos)) != len(rhos):
        raise SystemExit("--rhos values must be distinct")

    checks = exhaustive_self_test(device)
    print(
        f"self-test: {checks} scalar/kernel and topology checks passed",
        flush=True,
    )
    if args.self_test and args.output is None:
        return

    generator = torch.Generator(device=device)
    generator.manual_seed(args.seed)
    all_rows: list[dict[str, int | float | str]] = []
    print(
        f"device: {device.type} ({device_name(device)}); "
        f"samples/side: {args.samples:,}",
        flush=True,
    )
    for side in args.sides:
        rows = simulate_side(
            side,
            args.samples,
            args.batch_size,
            rhos,
            generator,
            device,
        )
        all_rows.extend(rows)
        first = rows[0]
        print(
            f"{side:>2}x{side:<2}: influence "
            f"{first['total_influence_estimate']:.4f} +/- "
            f"{first['total_influence_se']:.4f}; "
            f"blue {first['blue_win_probability']:.4f}; "
            f"{first['side_runtime_seconds']:.2f}s",
            flush=True,
        )

    if args.output is None:
        raise SystemExit("--output is required unless only --self-test is used")
    metadata: dict[str, int | str] = {
        "seed": args.seed,
        "device_type": device.type,
        "device_name": device_name(device),
        "torch_version": torch.__version__,
        "cuda_version": torch.version.cuda or "",
        "batch_size": args.batch_size,
    }
    write_csv(args.output, all_rows, metadata)
    print(f"wrote {len(all_rows)} rows to {args.output}", flush=True)


if __name__ == "__main__":
    main()
