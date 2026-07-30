"""Exact checks for notes/ProjectorMatrixHypercubeLaw.md."""

from fractions import Fraction
from itertools import product


Q = Fraction
ZERO = Q(0)
ONE = Q(1)


def matrix(a, b, c, d):
    return ((a, b), (c, d))


I = matrix(ONE, ZERO, ZERO, ONE)
Z = matrix(ONE, ZERO, ZERO, -ONE)


def madd(*matrices):
    return tuple(
        tuple(sum((m[i][j] for m in matrices), ZERO) for j in range(2))
        for i in range(2)
    )


def mscale(scalar, m):
    return tuple(tuple(scalar * entry for entry in row) for row in m)


def mmul(left, right):
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(2)), ZERO)
            for j in range(2)
        )
        for i in range(2)
    )


def transpose(m):
    return tuple(tuple(m[j][i] for j in range(2)) for i in range(2))


def trace(m):
    return m[0][0] + m[1][1]


def determinant(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def is_scalar_matrix(m, scalar):
    return m == mscale(scalar, I)


def projector_frame(q, r):
    """Return A,B,C,D with Q=D-C and Q^2=I."""
    assert q * q + r * r == 1
    involution = matrix(q, r, r, -q)
    a = mscale(Q(1, 2), madd(I, Z))
    b = mscale(Q(1, 2), madd(I, mscale(-1, Z)))
    c = mscale(Q(1, 2), madd(I, mscale(-1, involution)))
    d = mscale(Q(1, 2), madd(I, involution))
    return (a, b, c, d), involution


def signed_sum(frame, signs):
    return madd(*(mscale(sign, projector) for sign, projector in zip(signs, frame)))


def expected_norm_squared(q, signs):
    ea, eb, ec, ed = signs
    ab_aligned = ea == eb
    cd_aligned = ec == ed

    if ab_aligned and cd_aligned:
        return Q(ea + ec) ** 2
    if ab_aligned != cd_aligned:
        return Q(4)

    x = Q(ea - eb, 2)
    y = Q(ed - ec, 2)
    return 2 + 2 * x * y * q


def verify_norm_squared(frame, involution, q, signs):
    s = signed_sum(frame, signs)
    ea, eb, ec, ed = signs
    ab_aligned = ea == eb
    cd_aligned = ec == ed
    expected = expected_norm_squared(q, signs)

    if ab_aligned and cd_aligned:
        assert is_scalar_matrix(s, Q(ea + ec))
    elif ab_aligned != cd_aligned:
        # S=mI+R with m=+/-1 and R a signed involution.
        m = Q(ea) if ab_aligned else Q(ec)
        remainder = madd(s, mscale(-m, I))
        assert mmul(remainder, remainder) == I
        assert trace(remainder) == 0
        assert determinant(s) == 0
        assert abs(trace(s)) == 2
    else:
        assert trace(s) == 0
        assert is_scalar_matrix(mmul(s, s), expected)

    return expected


def phase_expression(a, b, c, d):
    anticommutator = madd(mmul(b, c), mmul(c, b))
    return madd(
        mmul(mmul(a, anticommutator), d),
        mmul(mmul(d, anticommutator), a),
    )


def main() -> None:
    rational_frames = (
        (Q(1), Q(0)),
        (Q(0), Q(1)),
        (Q(3, 5), Q(4, 5)),
        (Q(-5, 13), Q(12, 13)),
    )
    total_sign_checks = 0

    for q, r in rational_frames:
        frame, involution = projector_frame(q, r)
        a, b, c, d = frame

        for projector in frame:
            assert transpose(projector) == projector
            assert mmul(projector, projector) == projector
            assert trace(projector) == 1
            assert determinant(projector) == 0

        assert madd(a, b) == I
        assert madd(c, d) == I
        assert madd(*frame) == mscale(2, I)
        assert mmul(involution, involution) == I
        assert madd(mmul(Z, involution), mmul(involution, Z)) == mscale(
            2 * q, I
        )
        assert madd(*(mmul(projector, projector) for projector in frame)) == (
            mscale(2, I)
        )

        observed = {}
        norm_square_sum = ZERO
        for signs in product((-1, 1), repeat=4):
            norm_squared = verify_norm_squared(
                frame, involution, q, signs
            )
            observed[norm_squared] = observed.get(norm_squared, 0) + 1
            norm_square_sum += norm_squared
            total_sign_checks += 1

        expected_atoms = {}
        for value, multiplicity in (
            (Q(0), 2),
            (Q(4), 10),
            (2 + 2 * q, 2),
            (2 - 2 * q, 2),
        ):
            expected_atoms[value] = expected_atoms.get(value, 0) + multiplicity

        assert observed == expected_atoms
        assert norm_square_sum / 16 == 3

    # The mutually unbiased q=0 frame is exactly the Zhang projector quartet.
    zhang_frame, _ = projector_frame(Q(0), Q(1))
    a, b, c, d = zhang_frame
    assert a == matrix(ONE, ZERO, ZERO, ZERO)
    assert b == matrix(ZERO, ZERO, ZERO, ONE)
    assert c == matrix(Q(1, 2), Q(-1, 2), Q(-1, 2), Q(1, 2))
    assert d == matrix(Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2))

    norm_square_counts = {}
    for signs in product((-1, 1), repeat=4):
        value = expected_norm_squared(Q(0), signs)
        norm_square_counts[value] = norm_square_counts.get(value, 0) + 1
    assert norm_square_counts == {Q(0): 2, Q(2): 4, Q(4): 10}

    # Recover the exact quartic counterexample that supplies the frame.
    quartic = phase_expression(a, b, c, d)
    assert quartic == matrix(Q(-1, 2), Q(-1, 4), Q(-1, 4), ZERO)
    assert trace(quartic) == Q(-1, 2)
    assert determinant(quartic) == Q(-1, 16)
    assert trace(quartic) ** 2 - 4 * determinant(quartic) == Q(1, 2)

    print("projector matrix-hypercube checks passed")
    print(f"  {total_sign_checks} exact two-basis sign configurations")
    print("  exact general atom multiplicities and second moment")
    print("  Zhang law: 2 zero, 4 norm-sqrt(2), and 10 norm-2 vertices")
    print("  exact recovery of the quartic counterexample matrix")


if __name__ == "__main__":
    main()
