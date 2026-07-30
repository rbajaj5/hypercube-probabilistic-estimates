"""Exact checks for notes/ProjectorHeavyTailMixtures.md."""

from fractions import Fraction
from itertools import product

from verify_projector_matrix_hypercube import (
    I,
    Q,
    ZERO,
    expected_norm_squared,
    madd,
    matrix,
    mscale,
    phase_expression,
    projector_frame,
    trace,
    determinant,
)


def moment_from_cube(q: Fraction, alpha: int) -> Fraction:
    """Return E[W_q^alpha] exactly for positive even alpha."""
    assert alpha > 0 and alpha % 2 == 0
    exponent = alpha // 2
    total = ZERO
    for signs in product((-1, 1), repeat=4):
        total += expected_norm_squared(q, signs) ** exponent
    return total / 16


def formula_moment(q: Fraction, alpha: int) -> Fraction:
    """The exact Pareto-mixture coefficient for positive even alpha."""
    assert alpha > 0 and alpha % 2 == 0
    exponent = alpha // 2
    return (
        Q(5, 8) * 2**alpha
        + Q(1, 8) * (2 + 2 * q) ** exponent
        + Q(1, 8) * (2 - 2 * q) ** exponent
    )


def polynomial_derivative(coefficients):
    return tuple(
        (index + 1) * coefficients[index + 1]
        for index in range(len(coefficients) - 1)
    )


def polynomial_add(left, right):
    degree = max(len(left), len(right))
    return tuple(
        (left[i] if i < len(left) else ZERO)
        + (right[i] if i < len(right) else ZERO)
        for i in range(degree)
    )


def polynomial_scale(scalar, coefficients):
    return tuple(scalar * coefficient for coefficient in coefficients)


def scaled_frame(frame, scales):
    return tuple(
        mscale(scale, projector)
        for scale, projector in zip(scales, frame)
    )


def main() -> None:
    q_values = (Q(-1), Q(-5, 13), Q(0), Q(3, 5), Q(1))
    moment_checks = 0
    for q in q_values:
        for alpha in (2, 4, 6, 8):
            direct = moment_from_cube(q, alpha)
            formula = formula_moment(q, alpha)
            assert direct == formula
            moment_checks += 1

        # The second-moment boundary is angle-independent.
        assert formula_moment(q, 2) == 3
        # Above the boundary, q=0 is the minimum.
        assert formula_moment(q, 4) == 11 + q * q
        assert formula_moment(q, 6) == 42 + 6 * q * q

    assert formula_moment(Q(0), 4) < formula_moment(Q(1), 4)
    assert formula_moment(Q(0), 6) < formula_moment(Q(1), 6)
    # At alpha=1, the Zhang-minus-coincident coefficient is
    # (sqrt(2)-1)/4 > 0; squaring certifies sqrt(2)>1.
    assert Q(2) > Q(1) ** 2

    zhang_frame, _ = projector_frame(Q(0), Q(1))
    base_quartic = phase_expression(*zhang_frame)
    assert base_quartic == matrix(
        Q(-1, 2), Q(-1, 4), Q(-1, 4), ZERO
    )

    common_scale_checks = 0
    for scale in (Q(1), Q(3, 2), Q(2), Q(7, 3)):
        scaled = scaled_frame(zhang_frame, (scale,) * 4)
        assert phase_expression(*scaled) == mscale(
            scale**4, base_quartic
        )
        assert madd(*scaled) == mscale(2 * scale, I)
        benchmark = Q(1, 64) * (2 * scale) ** 4
        assert benchmark == scale**4 / 4
        common_scale_checks += 1

    independent_scale_checks = 0
    scale_vectors = (
        (Q(1), Q(1), Q(1), Q(1)),
        (Q(2), Q(3), Q(5), Q(7)),
        (Q(3, 2), Q(5, 4), Q(7, 3), Q(11, 5)),
    )
    for scales in scale_vectors:
        scaled = scaled_frame(zhang_frame, scales)
        scale_product = scales[0] * scales[1] * scales[2] * scales[3]
        assert phase_expression(*scaled) == mscale(
            scale_product, base_quartic
        )
        independent_scale_checks += 1

    # The base eigenvalue discriminant gives kappa=(1+sqrt(2))/4.
    assert trace(base_quartic) == Q(-1, 2)
    assert determinant(base_quartic) == Q(-1, 16)
    assert trace(base_quartic) ** 2 - 4 * determinant(base_quartic) == Q(1, 2)

    # If P(z)=sum_{j=0}^3 (alpha*z)^j/j!, then
    # d[e^{-alpha z}P(z)]/dz =
    # -e^{-alpha z} alpha^4 z^3/3!, the Erlang(4) density identity.
    erlang_checks = 0
    for alpha in (Q(1, 2), Q(1), Q(3), Q(7, 2)):
        survival_polynomial = (
            Q(1),
            alpha,
            alpha**2 / 2,
            alpha**3 / 6,
        )
        derivative = polynomial_derivative(survival_polynomial)
        logarithmic_derivative = polynomial_add(
            derivative,
            polynomial_scale(-alpha, survival_polynomial),
        )
        assert logarithmic_derivative == (
            ZERO,
            ZERO,
            ZERO,
            -(alpha**4) / 6,
        )
        assert survival_polynomial[0] == 1
        erlang_checks += 1

    print("projector heavy-tail checks passed")
    print(f"  {moment_checks} exact Pareto-coefficient moment checks")
    print(f"  {common_scale_checks} common-shock quartic checks")
    print(f"  {independent_scale_checks} independent-shock quartic checks")
    print(f"  {erlang_checks} exact Erlang survival-polynomial checks")


if __name__ == "__main__":
    main()
