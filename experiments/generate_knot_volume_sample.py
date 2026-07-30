"""Generate a reproducible SnapPy CensusKnots volume sample as CSV.

SnapPy is an optional experiment dependency:
    python -m pip install snappy

The committed verification suite reads the fixed CSV output and does not
require SnapPy.
"""

from __future__ import annotations

import argparse
import csv
import random
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument("--sample-size", type=int, default=16)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        import snappy
    except ImportError as exc:
        raise SystemExit(
            "SnapPy is required for generation; install it with "
            "`python -m pip install snappy`."
        ) from exc

    census_size = len(snappy.CensusKnots)
    if not 0 < args.sample_size <= census_size:
        raise SystemExit(
            f"sample size must lie between 1 and {census_size}"
        )

    rng = random.Random(args.seed)
    indices = sorted(rng.sample(range(census_size), args.sample_size))

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(
        (
            "seed",
            "census_size",
            "sample_size",
            "index",
            "name",
            "volume",
            "volume_accuracy",
            "solution_type",
        )
    )
    for index in indices:
        manifold = snappy.CensusKnots[index]
        volume, accuracy = manifold.volume(accuracy=True)
        writer.writerow(
            (
                args.seed,
                census_size,
                args.sample_size,
                index,
                manifold.name(),
                f"{float(volume):.12f}",
                accuracy,
                manifold.solution_type(),
            )
        )


if __name__ == "__main__":
    main()
