from math import comb, factorial, isqrt
import sympy as symbolic
from models import RECURRENCES, add, apply_recurrence, derivative, hit_model, scale, theta, variable, index
import verify_core as core

def main():
    core.main()
    limit = 300
    rows = {hits:[core.rook_coefficient(position,hits) for position in range(limit+1)]
            for hits in range(8)}
    for name,(hits,coefficients) in RECURRENCES.items():
        homogeneous,particular = apply_recurrence(hit_model(hits),coefficients)
        assert homogeneous == 0
        expected = -3202597896*variable**5 if name == 'A000185-six' else 0
        assert symbolic.expand(particular-expected) == 0
        start = max(hits+len(coefficients)-1,6 if name == 'A000159-short' else 0)
        for position in range(start,limit+1):
            assert sum(int(coefficient.subs(index,position))*rows[hits][position-offset]
                       for offset,coefficient in enumerate(coefficients)) == 0
        print('PASS:',name,'exact rational identity and recurrence',str(start)+'..300',flush=True)
    straight = [sum((-1)**offset*comb(2*position-offset,offset)*factorial(position-offset)
                    for offset in range(position+1)) for position in range(limit)]
    for position in range(1,limit+1):
        assert rows[1][position] == 2*position*straight[position-1]
    assert rows[1][4] == 8 and 2*4*straight[4] == 24
    linear = hit_model(1)
    straight_model = (1/(1+variable),symbolic.Integer(0))
    assert add(linear,scale(add(theta(straight_model),straight_model),-2*variable)) == (0,0)
    print('PASS: corrected one-hit identity 1..300; n=4 gives 8 versus 24',flush=True)
    model = hit_model(3)
    identity = add(scale(derivative(model),1-variable),scale(model,-1-1/variable))
    other = hit_model(2)
    assert add(identity,scale(add(theta(theta(other)),scale(other,-1)),-symbolic.Rational(2,3))) == (0,0)
    for position in range(3,limit):
        assert 3*(position*rows[3][position+1]-(position+1)*rows[3][position]) == 2*(position**2-1)*rows[2][position]
    print('PASS: two-hit/three-hit relation, exact rational identity and 3..299',flush=True)
    primes = [candidate for candidate in range(2,98) if all(candidate%divisor for divisor in range(2,isqrt(candidate)+1))]
    count = 0
    for prime in primes:
        for hits in range(min(prime,8)):
            for position in range(int(hits==0),101):
                assert (rows[hits][position]+2*rows[hits][position+prime]+rows[hits][position+2*prime])%prime == 0
                count += 1
    print('PASS: uniform hit congruence,',count,'pairs, 25 primes <=97 and hit degrees <min(p,8)',flush=True)
    for position in range(101):
        assert (rows[2][position]+2*rows[2][position+2]+rows[2][position+4])%2 == 0
    assert all(value%2 == 0 for degree in [1,3,5,7] for value in rows[degree])
    for prime in primes:
        if prime != 2:
            for position in range(301-prime**2):
                assert (rows[3][position+prime**2]+rows[3][position])%prime == 0
    for position in range(292):
        assert (rows[3][position+9]+rows[3][position])%3 == 0
    print('PASS: two-hit p=2 extension, odd-hit parity, three-hit antiperiods including p=3',flush=True)
    quotient = [0]+[rows[2][position]//position for position in range(1,301)]
    assert all(rows[2][position]%position==0 for position in range(1,301))
    quotient_model = (variable*(2+3*variable+2*variable**2)/(1+variable)**2,-2*variable)
    assert add(theta(quotient_model),scale(hit_model(2),-1)) == (0,0)
    for prime in primes:
        for position in range(1,101):
            assert (quotient[position]+2*quotient[position+prime]+quotient[position+2*prime])%prime == 0
    print('PASS: A000426 quotient identity and 2500 prime-shift checks',flush=True)
    for prime in primes:
        for position in range(1,max(1,301-4*prime)):
            assert (rows[0][position+4*prime]-2*rows[0][position+2*prime]+rows[0][position])%prime == 0
    print('PASS: A000179 four-prime-shift relation on all supported indices <=300',flush=True)
    hit_variable = symbolic.Symbol('t')
    shifted = variable*(hit_variable-1)
    binary_series = sum(factorial(degree)*variable**degree*(1+shifted)/(1-shifted)**(2*degree+1)
                        for degree in range(2))
    binary_target = 1+variable/(1+variable**2+variable**2*hit_variable**2)
    numerator = symbolic.cancel(binary_series-binary_target).as_numer_denom()[0]
    assert symbolic.Poly(numerator,variable,hit_variable,modulus=2).is_zero
    ternary_series = binary_series+2*variable**2*(1+shifted)/(1-shifted)**5
    cubic_coefficient = symbolic.cancel(symbolic.diff(ternary_series,hit_variable,3).subs(hit_variable,0)/6)
    ternary_target = (variable**7-variable**5-variable**3)/(1+variable)**8
    numerator = symbolic.cancel(cubic_coefficient-ternary_target).as_numer_denom()[0]
    assert symbolic.Poly(numerator,variable,modulus=3).is_zero
    ternary_polynomial = variable**8+variable**7-variable**6-variable**5-variable**4-variable**3
    numerator = symbolic.cancel((1+variable**9)*ternary_target-ternary_polynomial).as_numer_denom()[0]
    assert symbolic.Poly(numerator,variable,modulus=3).is_zero
    print('PASS: exact characteristic-two identity and characteristic-three antiperiod identity',flush=True)
    print('ALL_SIBLING_CHECKS_OK')

if __name__ == '__main__':
    main()
