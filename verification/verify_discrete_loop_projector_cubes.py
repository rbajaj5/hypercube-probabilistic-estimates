"""Verify exact discrete loop equations for rank-one projector cubes."""

from __future__ import annotations

from fractions import Fraction
from itertools import product


Q = Fraction
Matrix = tuple[tuple[Fraction, ...], ...]
Vector = tuple[Fraction, ...]


def zero(size: int) -> Matrix:
    return tuple(tuple(Q(0) for _ in range(size)) for _ in range(size))


def identity(size: int) -> Matrix:
    return tuple(
        tuple(Q(int(row == column)) for column in range(size))
        for row in range(size)
    )


def add(*matrices: Matrix) -> Matrix:
    size = len(matrices[0])
    return tuple(
        tuple(
            sum((matrix[row][column] for matrix in matrices), Q(0))
            for column in range(size)
        )
        for row in range(size)
    )


def scale(value: Fraction, matrix: Matrix) -> Matrix:
    return tuple(
        tuple(value * entry for entry in row)
        for row in matrix
    )


def multiply(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    return tuple(
        tuple(
            sum(
                (
                    left[row][index] * right[index][column]
                    for index in range(size)
                ),
                Q(0),
            )
            for column in range(size)
        )
        for row in range(size)
    )


def trace(matrix: Matrix) -> Fraction:
    return sum(
        (matrix[index][index] for index in range(len(matrix))),
        Q(0),
    )


def frobenius_inner(left: Matrix, right: Matrix) -> Fraction:
    size = len(left)
    return sum(
        (
            left[row][column] * right[row][column]
            for row in range(size)
            for column in range(size)
        ),
        Q(0),
    )


def frobenius_squared(matrix: Matrix) -> Fraction:
    return frobenius_inner(matrix, matrix)


def determinant(matrix: Matrix) -> Fraction:
    work = [list(row) for row in matrix]
    size = len(work)
    output = Q(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return Q(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            output = -output
        pivot_value = work[column][column]
        output *= pivot_value
        for entry in range(column, size):
            work[column][entry] /= pivot_value
        for row in range(column + 1, size):
            factor = work[row][column]
            if factor:
                for entry in range(column, size):
                    work[row][entry] -= factor * work[column][entry]
    return output


def inverse(matrix: Matrix) -> Matrix:
    size = len(matrix)
    unit = identity(size)
    work = [
        list(matrix[row]) + list(unit[row])
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if work[row][column]
        )
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        for entry in range(2 * size):
            work[column][entry] /= pivot_value
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                for entry in range(2 * size):
                    work[row][entry] -= factor * work[column][entry]
    return tuple(
        tuple(work[row][size:])
        for row in range(size)
    )


def outer(vector: Vector) -> Matrix:
    return tuple(
        tuple(left * right for right in vector)
        for left in vector
    )


def quadratic(vector: Vector, matrix: Matrix) -> Fraction:
    size = len(vector)
    return sum(
        (
            vector[row] * matrix[row][column] * vector[column]
            for row in range(size)
            for column in range(size)
        ),
        Q(0),
    )


def signed_sum(
    projectors: tuple[Matrix, ...],
    signs: tuple[int, ...],
) -> Matrix:
    return add(
        *(
            scale(Q(sign), projector)
            for sign, projector in zip(signs, projectors, strict=True)
        )
    )


def resolvent(matrix: Matrix, spectral_parameter: Fraction) -> Matrix:
    return inverse(
        add(
            matrix,
            scale(-spectral_parameter, identity(len(matrix))),
        )
    )


def average(matrices: list[Matrix]) -> Matrix:
    return scale(Q(1, len(matrices)), add(*matrices))


def scalar_product(values: tuple[Fraction, ...]) -> Fraction:
    output = Q(1)
    for value in values:
        output *= value
    return output


def moment(
    projectors: tuple[Matrix, ...],
    sign_vectors: tuple[tuple[int, ...], ...],
    exponent: int,
) -> Fraction:
    dimension = len(projectors[0])
    total = Q(0)
    for signs in sign_vectors:
        matrix = signed_sum(projectors, signs)
        power = identity(dimension)
        for _ in range(exponent):
            power = multiply(power, matrix)
        total += trace(power) / dimension
    return total / len(sign_vectors)


def verify_frame(vectors: tuple[Vector, ...]) -> int:
    dimension = len(vectors[0])
    projectors = tuple(outer(vector) for vector in vectors)
    sign_vectors = tuple(product((-1, 1), repeat=len(vectors)))
    checks = 0

    for vector, projector in zip(vectors, projectors, strict=True):
        assert quadratic(vector, identity(dimension)) == 1
        assert multiply(projector, projector) == projector
        assert trace(projector) == 1
        checks += 3

    second_moment = moment(projectors, sign_vectors, 2)
    fourth_moment = moment(projectors, sign_vectors, 4)
    expected_fourth = Q(len(vectors))
    for left in range(len(projectors)):
        for right in range(left + 1, len(projectors)):
            overlap = trace(
                multiply(projectors[left], projectors[right])
            )
            expected_fourth += 4 * overlap + 2 * overlap * overlap
    assert moment(projectors, sign_vectors, 0) == 1
    assert moment(projectors, sign_vectors, 1) == 0
    assert second_moment == Q(len(vectors), dimension)
    assert moment(projectors, sign_vectors, 3) == 0
    assert fourth_moment == expected_fourth / dimension
    checks += 5

    transforms: dict[Fraction, Fraction] = {}
    for spectral_parameter in (
        Q(len(vectors) + 2),
        Q(-len(vectors) - 2),
    ):
        full_resolvents = {
            signs: resolvent(
                signed_sum(projectors, signs),
                spectral_parameter,
            )
            for signs in sign_vectors
        }
        expected_resolvent = average(list(full_resolvents.values()))
        transforms[spectral_parameter] = (
            trace(expected_resolvent) / dimension
        )

        coordinate_rhs: list[Matrix] = []
        scalar_terms: list[Fraction] = []
        for index, projector in enumerate(projectors):
            signed_resolvents = [
                scale(Q(signs[index]), full_resolvents[signs])
                for signs in sign_vectors
            ]
            integration_by_parts_left = average(signed_resolvents)

            cavity_products: list[Matrix] = []
            scalar_cavity_terms: list[Fraction] = []
            other_indices = tuple(
                other
                for other in range(len(vectors))
                if other != index
            )
            for other_signs in product((-1, 1), repeat=len(other_indices)):
                signs_by_index = dict(
                    zip(other_indices, other_signs, strict=True)
                )
                cavity_matrix = add(
                    *(
                        scale(
                            Q(signs_by_index[other]),
                            projectors[other],
                        )
                        for other in other_indices
                    ),
                    zero(dimension),
                )
                plus_resolvent = resolvent(
                    add(cavity_matrix, projector),
                    spectral_parameter,
                )
                minus_resolvent = resolvent(
                    add(cavity_matrix, scale(Q(-1), projector)),
                    spectral_parameter,
                )

                # Full finite-difference resolvent identity.
                assert add(
                    plus_resolvent,
                    scale(Q(-1), minus_resolvent),
                ) == scale(
                    Q(-2),
                    multiply(
                        multiply(plus_resolvent, projector),
                        minus_resolvent,
                    ),
                )

                product_term = multiply(
                    multiply(plus_resolvent, projector),
                    minus_resolvent,
                )
                cavity_products.append(product_term)

                cavity_resolvent = resolvent(
                    cavity_matrix,
                    spectral_parameter,
                )
                h_value = quadratic(vectors[index], cavity_resolvent)
                plus_h = quadratic(vectors[index], plus_resolvent)
                minus_h = quadratic(vectors[index], minus_resolvent)
                assert plus_h == h_value / (1 + h_value)
                assert minus_h == h_value / (1 - h_value)
                assert trace(
                    multiply(
                        multiply(projector, plus_resolvent),
                        multiply(projector, minus_resolvent),
                    )
                ) == h_value * h_value / (1 - h_value * h_value)

                base = add(
                    cavity_matrix,
                    scale(-spectral_parameter, identity(dimension)),
                )
                assert determinant(
                    add(base, projector)
                ) == determinant(base) * (1 + h_value)
                assert determinant(
                    add(base, scale(Q(-1), projector))
                ) == determinant(base) * (1 - h_value)
                scalar_cavity_terms.append(
                    h_value * h_value / (1 - h_value * h_value)
                )
                checks += 8

            averaged_product = average(cavity_products)
            assert integration_by_parts_left == scale(
                Q(-1), averaged_product
            )
            coordinate_rhs.append(
                multiply(projector, averaged_product)
            )
            scalar_terms.append(
                sum(scalar_cavity_terms, Q(0))
                / len(scalar_cavity_terms)
            )
            checks += 1

        # Exact matrix loop equation.
        loop_left = add(
            identity(dimension),
            scale(spectral_parameter, expected_resolvent),
        )
        loop_right = scale(Q(-1), add(*coordinate_rhs))
        assert loop_left == loop_right

        # Its trace and rank-one cavity forms.
        normalized_stieltjes = trace(expected_resolvent) / dimension
        scalar_left = 1 + spectral_parameter * normalized_stieltjes
        scalar_right = -sum(scalar_terms, Q(0)) / dimension
        assert trace(loop_left) / dimension == scalar_left
        assert trace(loop_right) / dimension == scalar_right
        assert scalar_left == scalar_right
        checks += 4

        # Paired one-sample residual: exact mean zero and a deterministic
        # Frobenius bound using the explicit outside-support gap.
        residuals: list[Matrix] = []
        for signs in sign_vectors:
            matrix = signed_sum(projectors, signs)
            full_resolvent = full_resolvents[signs]
            residual = add(
                identity(dimension),
                scale(spectral_parameter, full_resolvent),
            )
            for index, projector in enumerate(projectors):
                cavity_matrix = add(
                    *(
                        scale(Q(signs[other]), projectors[other])
                        for other in range(len(vectors))
                        if other != index
                    ),
                    zero(dimension),
                )
                plus_resolvent = resolvent(
                    add(cavity_matrix, projector),
                    spectral_parameter,
                )
                minus_resolvent = resolvent(
                    add(cavity_matrix, scale(Q(-1), projector)),
                    spectral_parameter,
                )
                residual = add(
                    residual,
                    multiply(
                        multiply(projector, plus_resolvent),
                        multiply(projector, minus_resolvent),
                    ),
                )
            residuals.append(residual)

        assert average(residuals) == zero(dimension)
        support_gap = abs(spectral_parameter) - len(vectors)
        assert support_gap > 0
        squared_bound = (
            dimension
            * len(vectors) ** 2
            * (1 / support_gap + 1 / (support_gap * support_gap)) ** 2
        )
        assert all(
            frobenius_squared(residual) <= squared_bound
            for residual in residuals
        )

        one_sample_mean_square = (
            sum(
                (frobenius_squared(residual) for residual in residuals),
                Q(0),
            )
            / len(residuals)
        )
        two_sample_mean_square = sum(
            (
                frobenius_squared(scale(Q(1, 2), add(left, right)))
                for left in residuals
                for right in residuals
            ),
            Q(0),
        ) / (len(residuals) ** 2)
        assert two_sample_mean_square == one_sample_mean_square / 2
        checks += len(residuals) + 2

    positive_parameter = Q(len(vectors) + 2)
    assert transforms[-positive_parameter] == -transforms[positive_parameter]
    checks += 1

    # Exact hierarchy with a product of three auxiliary Stieltjes transforms.
    base_parameter = positive_parameter
    external_parameters = (
        Q(len(vectors) + 3),
        Q(-len(vectors) - 4),
        Q(len(vectors) + 5),
    )
    hierarchy_data: dict[
        tuple[int, ...],
        tuple[Matrix, tuple[Fraction, ...], Fraction],
    ] = {}
    for signs in sign_vectors:
        matrix = signed_sum(projectors, signs)
        base_resolvent = resolvent(matrix, base_parameter)
        external_transforms = tuple(
            trace(resolvent(matrix, parameter)) / dimension
            for parameter in external_parameters
        )
        hierarchy_data[signs] = (
            base_resolvent,
            external_transforms,
            scalar_product(external_transforms),
        )

    hierarchy_left = average(
        [
            scale(
                phi,
                add(
                    identity(dimension),
                    scale(base_parameter, base_resolvent),
                ),
            )
            for base_resolvent, _, phi in hierarchy_data.values()
        ]
    )
    hierarchy_coordinate_terms: list[Matrix] = []
    for index, projector in enumerate(projectors):
        coordinate_left = average(
            [
                scale(
                    Q(signs[index]) * hierarchy_data[signs][2],
                    hierarchy_data[signs][0],
                )
                for signs in sign_vectors
            ]
        )
        cavity_terms: list[Matrix] = []
        other_indices = tuple(
            other
            for other in range(len(vectors))
            if other != index
        )
        for other_signs in product((-1, 1), repeat=len(other_indices)):
            signs_by_index = dict(
                zip(other_indices, other_signs, strict=True)
            )
            plus_signs = tuple(
                1 if coordinate == index else signs_by_index[coordinate]
                for coordinate in range(len(vectors))
            )
            minus_signs = tuple(
                -1 if coordinate == index else signs_by_index[coordinate]
                for coordinate in range(len(vectors))
            )
            plus_resolvent, plus_transforms, plus_phi = hierarchy_data[
                plus_signs
            ]
            minus_resolvent, minus_transforms, minus_phi = hierarchy_data[
                minus_signs
            ]
            average_resolvent = scale(
                Q(1, 2),
                add(plus_resolvent, minus_resolvent),
            )
            difference_phi = (plus_phi - minus_phi) / 2
            average_phi = (plus_phi + minus_phi) / 2

            average_transforms = tuple(
                (plus + minus) / 2
                for plus, minus in zip(
                    plus_transforms,
                    minus_transforms,
                    strict=True,
                )
            )
            difference_transforms = tuple(
                (plus - minus) / 2
                for plus, minus in zip(
                    plus_transforms,
                    minus_transforms,
                    strict=True,
                )
            )
            odd_subset_expansion = Q(0)
            for choices in product((0, 1), repeat=len(external_parameters)):
                if sum(choices) % 2 == 0:
                    continue
                factors = tuple(
                    difference_transforms[position]
                    if choice
                    else average_transforms[position]
                    for position, choice in enumerate(choices)
                )
                odd_subset_expansion += scalar_product(factors)
            assert difference_phi == odd_subset_expansion

            cavity_terms.append(
                add(
                    scale(difference_phi, average_resolvent),
                    scale(
                        -average_phi,
                        multiply(
                            multiply(plus_resolvent, projector),
                            minus_resolvent,
                        ),
                    ),
                )
            )
            checks += 1

        coordinate_right = average(cavity_terms)
        assert coordinate_left == coordinate_right
        hierarchy_coordinate_terms.append(
            multiply(projector, coordinate_right)
        )
        checks += 1

    hierarchy_right = add(*hierarchy_coordinate_terms)
    assert hierarchy_left == hierarchy_right
    checks += 1

    return checks


FRAMES: tuple[tuple[Vector, ...], ...] = (
    ((Q(1),),),
    (
        (Q(1), Q(0)),
        (Q(0), Q(1)),
        (Q(3, 5), Q(4, 5)),
        (Q(-4, 5), Q(3, 5)),
    ),
    (
        (Q(1), Q(0)),
        (Q(0), Q(1)),
        (Q(5, 13), Q(12, 13)),
        (Q(12, 13), Q(-5, 13)),
        (Q(7, 25), Q(24, 25)),
    ),
    (
        (Q(1), Q(0), Q(0)),
        (Q(0), Q(1), Q(0)),
        (Q(0), Q(0), Q(1)),
        (Q(3, 5), Q(4, 5), Q(0)),
        (Q(0), Q(5, 13), Q(12, 13)),
    ),
)


total_checks = sum(verify_frame(frame) for frame in FRAMES)
print(
    f"verified {total_checks} exact discrete integration-by-parts, "
    "resolvent, determinant, moment, residual, and loop-hierarchy identities "
    f"across {len(FRAMES)} rational projector frames"
)
