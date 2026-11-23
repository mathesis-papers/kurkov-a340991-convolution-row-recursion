import argparse
from functools import lru_cache
from itertools import product
import json
from math import comb, factorial, isqrt
from pathlib import Path


def check_sequence(sequence, zero, one, add, multiply, max_row, max_column):
    def term(index):
        return sequence[index - 1] if 1 <= index <= len(sequence) else zero

    source = [zero] + [term(index) for index in range(1, max_row + 1)]
    powers = [[one] + [zero] * max_row]
    for exponent in range(1, max_row + 1):
        coefficients = []
        for degree in range(max_row + 1):
            total = zero
            for split in range(degree + 1):
                total = add(total, multiply(powers[-1][split], source[degree - split]))
            coefficients.append(total)
        powers.append(coefficients)

    @lru_cache(maxsize=None)
    def array(row_index, column_index):
        if row_index == 1:
            return (term(column_index),) + (zero,) * max_row
        shifted = array(row_index - 1, column_index + 1)
        first_column = array(row_index - 1, 1)
        return tuple(add(shifted[degree], multiply(first_column[degree - 1], term(column_index)))
                     if degree else shifted[0] for degree in range(max_row + 1))

    rows = [[one] + [zero] * max_row]
    for row_index in range(1, max_row + 1):
        row = [powers[exponent][row_index] for exponent in range(max_row + 1)]
        assert tuple(row) == (zero,) + array(row_index, 1)[:-1], ("row", row_index)
        rows.append(row)
        for exponent in range(1, max_row + 1):
            left_total = zero
            right_total = zero
            for part in range(1, row_index + 1):
                left_total = add(left_total, multiply(term(part), rows[row_index - part][exponent - 1]))
                right_total = add(right_total, multiply(rows[row_index - part][exponent - 1], term(part)))
            assert row[exponent] == left_total == right_total, ('two-sided row recurrence', row_index, exponent)
        for column_index in range(1, max_column + 1):
            expected = []
            for exponent in range(max_row + 1):
                total = zero
                for degree in range(row_index):
                    total = add(total, multiply(powers[exponent][degree],
                                                term(column_index + row_index - 1 - degree)))
                expected.append(total)
            assert tuple(expected) == array(row_index, column_index), ("tail", row_index, column_index)
    array.cache_clear()
    return rows


def primes(count):
    values = []
    candidate = 2
    while len(values) < count:
        if all(candidate % divisor for divisor in range(2, isqrt(candidate) + 1)):
            values.append(candidate)
        candidate += 1
    return values


def matrix_add(left, right):
    return tuple((left[index] + right[index]) % 5 for index in range(4))


def matrix_multiply(left, right):
    return tuple(sum(left[2 * row + inner] * right[2 * inner + column] for inner in range(2)) % 5
                 for row in range(2) for column in range(2))


def word_add(left, right):
    result = dict(left)
    for word, coefficient in right.items():
        result[word] = result.get(word, 0) + coefficient
    return result


def word_multiply(left, right):
    result = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            word = left_word + right_word
            result[word] = result.get(word, 0) + left_coefficient * right_coefficient
    return result


def check_factorial_fractions(rows, bound):
    unit = {(0, 0): 1}
    variable_z = {(1, 0): 1}
    variable_x = {(0, 1): 1}

    def add(left, right):
        result = dict(left)
        for monomial, value in right.items():
            result[monomial] = result.get(monomial, 0) + value
        return {monomial: value for monomial, value in result.items() if value}

    def negative(value):
        return {monomial: -coefficient for monomial, coefficient in value.items()}

    def multiply(left, right):
        result = {}
        for (left_z, left_x), left_value in left.items():
            for (right_z, right_x), right_value in right.items():
                degree_z = left_z + right_z
                degree_x = left_x + right_x
                if degree_z <= bound and degree_x <= bound:
                    monomial = (degree_z, degree_x)
                    result[monomial] = result.get(monomial, 0) + left_value * right_value
        return {monomial: value for monomial, value in result.items() if value}

    def inverse(value):
        assert value.get((0, 0)) == 1
        remainder = add(unit, negative(value))
        assert all(degree_z > 0 for degree_z, degree_x in remainder)
        total = unit
        power = unit
        for exponent in range(1, bound + 1):
            power = multiply(power, remainder)
            total = add(total, power)
        return total

    coefficients = [(index + 2) // 2 for index in range(2 * bound + 4)]
    tail = unit
    tails = {}
    for index in reversed(range(len(coefficients))):
        tail = inverse(add(unit, negative(multiply({(1, 0): coefficients[index]}, tail))))
        tails[index] = tail
    factorial_series = {(degree, 0): factorial(degree) for degree in range(bound + 1)}
    assert tail == factorial_series, 'Euler factorial continued fraction'
    jacobi = unit
    for index in reversed(range(bound + 2)):
        denominator = add(unit, {(1, 0): -(2 * index + 1)})
        denominator = add(denominator, negative(multiply({(2, 0): (index + 1) ** 2}, jacobi)))
        jacobi = inverse(denominator)
    assert jacobi == factorial_series, 'Jacobi contraction'
    shifted_tail = multiply(variable_z, tails[2])
    inner = inverse(add(add(unit, negative(multiply(variable_x, variable_z))), negative(shifted_tail)))
    outer = inverse(add(unit, negative(multiply(variable_z, inner))))
    bala = add(add(unit, negative(variable_x)), multiply(variable_x, outer))
    rational = inverse(add(add(unit, variable_x), negative(multiply(variable_x, factorial_series))))
    expected = {(degree_z, degree_x): coefficient
                for degree_z, row in enumerate(rows[:bound + 1])
                for degree_x, coefficient in enumerate(row) if coefficient}
    assert bala == rational == expected, 'Bala factorial generating functions'


def check_siblings(fixtures, add, multiply):
    bell = [1]
    for index in range(20):
        bell.append(sum(comb(index, part) * bell[part] for part in range(index + 1)))
    catalan = [comb(2 * index, index) // (index + 1) for index in range(21)]
    inputs = {
        'A357583': bell[1:],
        'A090238': [factorial(index) for index in range(1, 21)],
        'A128899': catalan[1:],
        'A205574': bell,
        'A084938': [factorial(index) for index in range(20)],
        'A059365': catalan,
    }
    rows = {}
    for identifier, sequence in inputs.items():
        rows[identifier] = check_sequence(sequence, 0, 1, add, multiply, 12, 4)
        flattened = [value for index in range(10) for value in rows[identifier][index][:index + 1]]
        expected = fixtures['siblings'][identifier]['terms'][:55]
        if identifier == 'A059365':
            assert expected[0] == 0 and flattened[0] == 1
            assert flattened[1:] == expected[1:]
        else:
            assert flattened == expected, identifier
    for index in range(11):
        assert sum(rows['A357583'][index]) == fixtures['siblings']['A129247']['terms'][index]
    for index in range(7):
        assert rows['A357583'][2 * index][index] == fixtures['siblings']['A357584']['terms'][index]
    check_factorial_fractions(rows['A090238'], 12)
    for ordinary, shifted in (('A357583', 'A205574'), ('A090238', 'A084938'), ('A128899', 'A059365')):
        for degree in range(1, 9):
            for power in range(1, 5):
                expected = sum(comb(power, column) * rows[ordinary][degree][column]
                               for column in range(min(power, degree) + 1))
                assert rows[shifted][degree + power][power] == expected, ('binomial shift', ordinary, degree, power)
    for degree in range(1, 13):
        for power in range(1, 13):
            factorial_row = rows['A090238']
            assert power * factorial_row[degree][power] == (
                power * factorial_row[degree - 1][power - 1]
                + (degree + power - 1) * factorial_row[degree - 1][power])
            catalan_row = rows['A128899']
            next_column = catalan_row[degree - 1][power + 1] if power < 12 else 0
            assert catalan_row[degree][power] == (catalan_row[degree - 1][power - 1]
                                                 + 2 * catalan_row[degree - 1][power] + next_column)
    constant_source = [1, 0, 1]
    altered_row = [0, constant_source[2],
                   sum(constant_source[index] * constant_source[2 - index] for index in range(3))]
    ordinary_rows = check_sequence([0, 1], 0, 1, add, multiply, 2, 1)
    assert altered_row == [0, 1, 2]
    assert ordinary_rows[2][:3] == [0, 1, 0]


def main():
    parser = argparse.ArgumentParser(description="Exact convolution-row and shifted-tail identities.")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    fixtures = json.loads(Path(__file__).with_name("fixtures.json").read_text())
    integer_add = lambda left, right: left + right
    integer_multiply = lambda left, right: left * right
    prime_rows = check_sequence(primes(20), 0, 1, integer_add, integer_multiply, 16, 4)
    for row_index, expected in enumerate(fixtures["A340991_rows"]):
        assert prime_rows[row_index][:row_index + 1] == expected
        assert all(value == 0 for value in prime_rows[row_index][row_index + 1:])
    for row_index in range(11):
        assert sum(prime_rows[row_index]) == fixtures["A030017_n1_onward"][row_index]
        assert sum((-1) ** degree * coefficient for degree, coefficient in enumerate(prime_rows[row_index])) == fixtures["A030018_n0_onward"][row_index]
        if row_index:
            assert sum(degree * coefficient for degree, coefficient in enumerate(prime_rows[row_index])) == fixtures["A030281_n1_onward"][row_index - 1]
        assert prime_rows[row_index][row_index] == 2 ** row_index
    for index in range(6):
        assert prime_rows[2 * index][index] == fixtures["A340990_n0_onward"][index]
    for column_index, identifier in ((2, "A014342"), (3, "A014343"), (4, "A014344")):
        expected = fixtures["prime_columns"][identifier]
        assert [prime_rows[index + column_index][column_index] for index in range(10)] == expected
    prime_set = set(primes(20))
    weighted_indicator = [index if index in prime_set else 1 for index in range(1, 14)]
    indicator_rows = check_sequence(weighted_indicator, 0, 1, integer_add, integer_multiply, 10, 4)
    flattened = [value for row_index in range(10) for value in indicator_rows[row_index][:row_index + 1]]
    assert flattened == fixtures["A357368_flat_n0_onward"]

    for sequence in product((-1, 0, 1), repeat=6):
        check_sequence(sequence, 0, 1, integer_add, integer_multiply, 6, 4)
    for sequence in product((0, 1), repeat=6):
        check_sequence(sequence, 0, 1, lambda left, right: (left + right) % 2,
                       lambda left, right: left * right % 2, 6, 4)
        check_sequence(sequence, 0, 1, lambda left, right: int(bool(left or right)),
                       lambda left, right: left * right, 6, 4)

    matrix_zero = (0, 0, 0, 0)
    matrix_one = (1, 0, 0, 1)
    upper = (0, 1, 0, 0)
    lower = (0, 0, 1, 0)
    assert matrix_multiply(upper, lower) != matrix_multiply(lower, upper)
    for sequence in product((matrix_one, upper, lower), repeat=3):
        check_sequence(sequence, matrix_zero, matrix_one, matrix_add, matrix_multiply, 8, 4)
    symbols = [{(index,): 1} for index in range(1, 12)]
    check_sequence(symbols, {}, {(): 1}, word_add, word_multiply, 8, 4)
    check_siblings(fixtures, integer_add, integer_multiply)

    result = {
        "status": "pass",
        "arithmetic": "exact",
        "prime_rows": {"max_n": 16, "max_array_column": 4, "OEIS_rows_matched": [0, 10]},
        "signed_integer_sequences": {"alphabet": [-1, 0, 1], "support_bound": 6, "sequences": 729, "max_n": 6, "max_array_column": 4},
        "mod_2_sequences": {"support_bound": 6, "sequences": 64, "max_n": 6, "max_array_column": 4},
        "boolean_sequences": {"support_bound": 6, "sequences": 64, "max_n": 6, "max_array_column": 4},
        "noncommutative_matrices": {"ring": "2 by 2 matrices over Z/5Z", "support_bound": 3, "sequences": 27, "max_n": 8, "max_array_column": 4},
        "free_noncommutative_symbols": {"symbols": 11, "max_n": 8, "max_array_column": 4},
        "A357368": {"input": "a_j=j for prime j, otherwise 1", "max_n": 10, "OEIS_rows_matched": [0, 9]},
        "sibling_triangles": {
            "OEIS_rows_matched": [0, 9],
            "array_identity_max_n": 12,
            "max_array_column": 4,
            "inputs": {"A357583": "Bell(j)", "A090238": "j!", "A128899": "Catalan(j)",
                       "A205574": "Bell(j-1)", "A084938": "(j-1)!", "A059365": "Catalan(j-1)"},
            "A059365_boundary": "Positive rows match; the OEIS row-zero entry is 0, whereas P_0=1.",
            "A129247_row_sums": [0, 10], "A357584_central_terms": [0, 6],
        },
        "factorial_fractions": {"max_z_degree": 12, "max_x_degree": 12,
                                "identities": ["Euler Stieltjes fraction", "Jacobi contraction", "Bala weighted rational series", "Bala continued fraction"]},
        "two_sided_recurrence": "Both displayed multiplication orders checked in every sequence family above.",
        "nonzero_constant_variant": {"source": "1+z^2", "n": 2, "row_polynomial": "x+2x^2", "shifted_recursion": "x"},
        "variant_identities": {"binomial_shift": {"degree_n": [1, 8], "power_m": [1, 4]},
                               "factorial_recurrence": {"n": [1, 12], "r": [1, 12]},
                               "Catalan_recurrence": {"n": [1, 12], "r": [1, 12]}},
        "identities": ["P_n(x) = x R(n,1), n >= 1", "[x^r]R(n,k) = [z^(n-1)]A(z)^r B_k(z)"],
        "OEIS_specializations": ["A030017(n+1), 0 <= n <= 10", "A030018(n), 0 <= n <= 10", "A030281(n), 1 <= n <= 10", "A340990(n), 0 <= n <= 5", "T(n,n)=2^n, 0 <= n <= 10", "A014342 through A014344: ten terms per column with source offsets", "A357368 rows 0 through 9"],
        "scope": "Finite exact identities; the all-index identities and arbitrary-semiring extension have the written proof in note.tex.",
    }
    serialized = json.dumps(result, indent=2) + "\n"
    if arguments.output:
        arguments.output.write_text(serialized)
    else:
        print(serialized, end="")


if __name__ == "__main__":
    main()
