"""Exact checks for notes/LowDimensionalSmallDeviations.md.

Run with:
    python verification/verify_low_dimensional_bounds.py

The checks use fractions only.  They audit the displayed formulas and
finite grids; the cited extremal theorems remain the proof of the universal
probability statements.
"""

from __future__ import annotations

from fractions import Fraction


def samuels_candidate(k: int, delta: Fraction) -> Fraction:
    """The k-active-coordinate lower-tail probability."""
    return (Fraction(k - 1, 1) + delta) ** k / (Fraction(k, 1) + delta) ** k


def endpoint_bound(n: int, delta: Fraction) -> Fraction:
    """The sharp endpoint formula recorded for n <= 4."""
    return min(samuels_candidate(1, delta), samuels_candidate(n, delta))


def paper_bound(n: int, delta: Fraction) -> Fraction:
    """The dimension-dependent bound in arXiv:2607.23980."""
    if delta < 1:
        return delta * Fraction(n, 1) ** n / (Fraction(n, 1) + delta) ** n
    return samuels_candidate(n, delta)


def shared_small_slack_bound(n: int, delta: Fraction) -> Fraction:
    """The 0 < delta <= 1 formula in the three July 2026 manuscripts."""
    assert 0 < delta <= 1
    return delta * Fraction(n, 1) ** n / (Fraction(n, 1) + delta) ** n


def stander_every_slack_bound(n: int, delta: Fraction) -> Fraction:
    """The piecewise every-slack corollary in the Stander manuscript."""
    if delta <= 1:
        return shared_small_slack_bound(n, delta)
    return Fraction(n, n + 1) ** n


def all_positive_grid(max_denominator: int, max_value: int) -> list[Fraction]:
    values = {
        Fraction(numerator, denominator)
        for denominator in range(1, max_denominator + 1)
        for numerator in range(1, max_value * denominator + 1)
    }
    return sorted(values)


def check_samuels_endpoint_reduction() -> int:
    """Check min_k q_k = min(q_1, q_n) on a large exact grid."""
    checks = 0
    for n in range(1, 5):
        for delta in all_positive_grid(max_denominator=35, max_value=6):
            candidates = [samuels_candidate(k, delta) for k in range(1, n + 1)]
            assert min(candidates) == endpoint_bound(n, delta), (
                n,
                delta,
                candidates,
            )
            checks += 1
    return checks


def check_strict_small_slack_improvement() -> int:
    """Check c_{n,delta} > b_{n,delta} for n=2,3,4 and 0<delta<1."""
    checks = 0
    deltas = sorted(
        {
            Fraction(numerator, denominator)
            for denominator in range(2, 101)
            for numerator in range(1, denominator)
        }
    )
    for n in range(2, 5):
        for delta in deltas:
            sharp = endpoint_bound(n, delta)
            baseline = paper_bound(n, delta)
            assert sharp > baseline, (n, delta, sharp, baseline)
            checks += 1
    return checks


def check_large_slack_agreement() -> int:
    """Check that the endpoint formula selects the paper's sharp branch."""
    checks = 0
    deltas = all_positive_grid(max_denominator=30, max_value=8)
    for n in range(1, 5):
        for delta in deltas:
            if delta < 1:
                continue
            assert endpoint_bound(n, delta) == paper_bound(n, delta), (n, delta)
            checks += 1
    return checks


def check_concurrent_work_comparisons() -> int:
    """Audit the exact overlap and scope differences of the July 2026 work."""
    checks = 0
    deltas = all_positive_grid(max_denominator=40, max_value=8)
    for n in range(1, 13):
        for delta in deltas:
            fu_han_et_al = paper_bound(n, delta)
            stander = stander_every_slack_bound(n, delta)
            if delta <= 1:
                nie_wei = shared_small_slack_bound(n, delta)
                assert fu_han_et_al == nie_wei == stander, (n, delta)
            else:
                assert fu_han_et_al > stander, (
                    n,
                    delta,
                    fu_han_et_al,
                    stander,
                )
            checks += 1
    return checks


def check_sharpness_families() -> int:
    """Audit the exact means and event probabilities of both extremizers."""
    checks = 0
    for n in range(1, 5):
        for delta in all_positive_grid(max_denominator=30, max_value=6):
            # One active coordinate:
            high_one = 1 + delta
            high_probability_one = 1 / high_one
            low_probability_one = 1 - high_probability_one
            assert high_one * high_probability_one == 1
            assert low_probability_one == samuels_candidate(1, delta)

            # All active coordinates:
            high_all = n + delta
            high_probability_all = 1 / high_all
            zero_probability_all = 1 - high_probability_all
            assert high_all * high_probability_all == 1
            assert zero_probability_all**n == samuels_candidate(n, delta)

            assert min(low_probability_one, zero_probability_all**n) == (
                endpoint_bound(n, delta)
            )
            checks += 1
    return checks


def check_bivariate_formula() -> int:
    """Check the exact n=2 formula recorded by Nie--Wei."""
    checks = 0
    for delta in all_positive_grid(max_denominator=50, max_value=6):
        t = 2 + delta
        merger_level_at_zero = Fraction(2, 1) / t - 1 / t**2
        merger_level_at_one = 1 / (t - 1)
        complement = 1 - max(merger_level_at_zero, merger_level_at_one)
        expected = min(
            delta / (1 + delta),
            ((1 + delta) / (2 + delta)) ** 2,
        )
        assert complement == expected
        assert expected == endpoint_bound(2, delta)
        checks += 1
    return checks


def main() -> None:
    endpoint_checks = check_samuels_endpoint_reduction()
    improvement_checks = check_strict_small_slack_improvement()
    agreement_checks = check_large_slack_agreement()
    concurrent_checks = check_concurrent_work_comparisons()
    sharpness_checks = check_sharpness_families()
    bivariate_checks = check_bivariate_formula()
    total = (
        endpoint_checks
        + improvement_checks
        + agreement_checks
        + concurrent_checks
        + sharpness_checks
        + bivariate_checks
    )
    print(f"Samuels endpoint-grid checks: {endpoint_checks}")
    print(f"strict small-slack comparisons: {improvement_checks}")
    print(f"large-slack agreement checks: {agreement_checks}")
    print(f"concurrent-work comparisons: {concurrent_checks}")
    print(f"sharpness-family checks: {sharpness_checks}")
    print(f"Nie--Wei bivariate-formula checks: {bivariate_checks}")
    print(f"all {total} exact checks passed")


if __name__ == "__main__":
    main()
