from math import comb, factorial
import sympy as symbolic


def rook_coefficient(index, degree):
    if index == 0:
        return int(degree == 0)
    return sum(
        (-1) ** (shift - degree)
        * (comb(2 * index - shift, shift)
           + (comb(2 * index - shift - 1, shift - 1) if shift else 0))
        * factorial(index - shift) * comb(shift, degree)
        for shift in range(degree, index + 1)
    )


def recurrence_coefficients(index):
    return [
        2 * (-252307 * index + 1041077),
        504614 * index**2 - 3362985 * index + 5118150,
        1280831 * index**2 - 7397886 * index + 6461565,
        746598 * index**2 - 2913543 * index - 1336090,
        -405481 * index**2 + 6175011 * index - 15469320,
        -375862 * index**2 + 4098537 * index - 8846430,
        2 * (-187931 * index + 560630),
    ]


def require_zero(expression):
    if symbolic.cancel(expression) != 0:
        raise ArithmeticError(str(expression))


def main():
    rows = [[1], [-1, 2], [0, 0, 2], [1, 0, 3, 2],
            [2, 8, 4, 8, 2], [13, 30, 40, 20, 15, 2]]
    for index, row in enumerate(rows):
        assert [rook_coefficient(index, degree)
                for degree in range(index + 1)] == row
    values = [rook_coefficient(index, 3) for index in range(301)]
    assert values[3:9] == [2, 8, 20, 152, 994, 7888]
    print('PASS: defining polynomial rows 0..5 and six initial values')
    for index in range(6, 301):
        assert sum(coefficient * values[index - shift]
                   for shift, coefficient in enumerate(
                       recurrence_coefficients(index))) == 0
    print('PASS: recurrence 9..300 and zero-extension cases 6..8')
    variable, argument, degree, index = symbolic.symbols('x z D n')
    numerator = -4*variable**5 + 6*variable**4 + 10*variable**3 + 8*variable**2 + 6*variable + 4
    constant = 4*variable**5 - 2*variable**4 - 16*variable**3 - 14*variable**2 - 10*variable - 4
    denominator = 3 * (1 + variable)**2
    pullback = variable / (1 + variable)**2
    homogeneous = numerator / denominator
    particular = constant / denominator
    multiplier = (1 - variable**3) / (variable**2 * (1 + variable))
    forcing = (variable**2 - 1) / variable**2
    theta_terms = [
        (1/argument - 1, -1/argument),
        (1 - 3/argument + 1/argument**2, 2/argument - 1/argument**2),
        (-1 + 7/argument - 6/argument**2 + 1/argument**3,
         -3/argument + 5/argument**2 - 1/argument**3),
    ]
    theta_homogeneous, theta_particular = symbolic.Integer(1), symbolic.Integer(0)
    for expected_homogeneous, expected_particular in theta_terms:
        theta_homogeneous, theta_particular = (
            argument * symbolic.diff(theta_homogeneous, argument)
            + theta_homogeneous * (1/argument - 1),
            argument * symbolic.diff(theta_particular, argument)
            - theta_homogeneous / argument,
        )
        require_zero(theta_homogeneous - expected_homogeneous)
        require_zero(theta_particular - expected_particular)
    theta_homogeneous = sum(weight * theta_terms[position][0]
                            for position, weight in enumerate(
                                [20-2*variable, 18-6*variable, 4-4*variable])) + 6
    theta_particular = sum(weight * theta_terms[position][1]
                           for position, weight in enumerate(
                               [20-2*variable, 18-6*variable, 4-4*variable]))
    require_zero(variable**3 / (3*(1+variable)**4)
                 * theta_homogeneous.subs(argument, pullback) - homogeneous)
    require_zero(variable**3 / (3*(1+variable)**4)
                 * theta_particular.subs(argument, pullback) - particular)
    require_zero(symbolic.diff(pullback, variable)
                 * (1-pullback) / pullback**2 - multiplier)
    require_zero(-symbolic.diff(pullback, variable) / pullback**2 - forcing)
    print('PASS: three Euler identities, affine reduction, and chain rule')
    first_homogeneous = symbolic.diff(homogeneous, variable) + homogeneous*multiplier
    first_particular = symbolic.diff(particular, variable) + homogeneous*forcing
    second_homogeneous = symbolic.diff(first_homogeneous, variable) + first_homogeneous*multiplier
    second_particular = symbolic.diff(first_particular, variable) + first_homogeneous*forcing
    second_operator = (-375862*variable**7 - 405481*variable**6 + 746598*variable**5
                       + 1280831*variable**4 + 504614*variable**3)
    first_operator = (-375862*variable**7 - 35945*variable**6 + 2525682*variable**5
                      + 2312643*variable**4 - 993731*variable**3
                      - 1849143*variable**2 - 504614*variable)
    zeroth_operator = (-1133912*variable**6 + 2249705*variable**5 + 2743028*variable**4
                       - 3357337*variable**3 - 3210883*variable**2
                       + 2259779*variable + 2082154)
    polynomial_forcing = 1136624*variable**3 - 10348*variable**4 + 1516100*variable**5
    require_zero(second_operator*second_homogeneous + first_operator*first_homogeneous
                 + zeroth_operator*homogeneous)
    require_zero(second_operator*second_particular + first_operator*first_particular
                 + zeroth_operator*particular - polynomial_forcing)
    operator = symbolic.Poly(sum(variable**shift * coefficient.subs(index, degree+shift)
                                for shift, coefficient in enumerate(
                                    recurrence_coefficients(index))), degree)
    require_zero(variable**2 * operator.nth(2) - second_operator)
    require_zero(variable * (operator.nth(2)+operator.nth(1)) - first_operator)
    require_zero(operator.nth(0) - zeroth_operator)
    print('PASS: both rational differential identities and all operator coefficients')
    primes = [candidate for candidate in range(2, 98)
              if candidate != 3 and all(candidate % divisor
                   for divisor in range(2, int(candidate**0.5)+1))]
    for prime in primes:
        for index in range(3, 101):
            assert (values[index] + 2*values[index+prime]
                    + values[index+2*prime]) % prime == 0
        polynomial = symbolic.Poly(
            numerator * sum(factorial(shift) * variable**shift
                            * (1+variable)**(2*prime-2-2*shift)
                            for shift in range(prime))
            + constant * (1+variable)**(2*prime-2), variable, modulus=prime)
        assert polynomial.degree() <= 2*prime+2
    assert (values[4] + 2*values[7] + values[10]) == 697756
    assert 697756 % 3 == 1
    print('PASS: prime-shift congruence for 24 primes <=97 and 3<=n<=100')
    print('PASS: truncated numerator degree bound for those primes; p=3 fails at n=4')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
