from math import factorial
import sympy as symbolic

variable, argument, index, shift_degree = symbolic.symbols('x z n k')
pullback = variable / (1 + variable)**2
multiplier = (1-variable**3)/(variable**2*(1+variable))
forcing = (variable**2-1)/variable**2

def derivative(model):
    return (symbolic.cancel(symbolic.diff(model[0], variable)+model[0]*multiplier),
            symbolic.cancel(symbolic.diff(model[1], variable)+model[0]*forcing))

def add(left, right):
    return tuple(symbolic.cancel(left[entry]+right[entry]) for entry in range(2))

def scale(model, factor):
    return tuple(symbolic.cancel(part*factor) for part in model)

def theta(model):
    return scale(derivative(model), variable)

def hit_model(hits):
    homogeneous, particular = [symbolic.Integer(1)], [symbolic.Integer(0)]
    for degree in range(1, hits+1):
        homogeneous.append(symbolic.cancel(argument*symbolic.diff(homogeneous[-1],argument)
                                            +(1/argument-1)*homogeneous[-1]))
        particular.append(symbolic.cancel(argument*symbolic.diff(particular[-1],argument)
                                           -homogeneous[-2]/argument))
    alpha = symbolic.prod(2*shift_degree+offset for offset in range(1,hits+1))/symbolic.Integer(factorial(hits))
    beta = (symbolic.prod(2*shift_degree+offset for offset in range(1,hits))/symbolic.Integer(factorial(hits-1))
            if hits else symbolic.Integer(0))
    weights = symbolic.Poly((1-variable)*alpha+(1+variable)*beta,shift_degree)
    prefactor = variable**hits/(1+variable)**(hits+1)
    return tuple(symbolic.factor((prefactor*sum(weights.nth(degree)*terms[degree]
                              for degree in range(hits+1))).subs(argument,pullback))
                 for terms in [homogeneous,particular])

RECURRENCES = {
    'A000425': (1, [
        (index-1)*(index-3),
        -index*(index-2)*(index-3),
        -index*(index-1)*(index-3),
        -index*(index-1)]),
    'A000033': (2, [
        (index-3)*(index-2)*(2*index-5)*(2*index-7),
        -(index-3)*(index-2)*index*(2*index-7)**2,
        -(index-4)*(index-3)*index*(2*index-3)**2,
        -(index-2)*index*(2*index-5)*(2*index-3)]),
    'A000159-short': (3, [
        (index-3)*(index-1)*(index-4)*(2*index-7)**2,
        -(index-4)*(2*index-5)*index*(2*(index-3)**3-(index-3)**2-6*(index-3)+6),
        -(2*index-7)*index*(index-1)*(2*(index-3)**3+(index-3)**2-6*(index-3)-6),
        -index*(index-1)*(index-2)*(2*index-5)**2]),
    'A000181': (4, [
        3*(2111*index-8303)*(index-4)*(2*index-9)**2,
        -(index-3)*(25332*index**4-377236*index**3+1898681*index**2-3320738*index+484000),
        -2*(index-4)*(12140*index**4-118152*index**3+337063*index**2-377436*index+225720),
        1052*index**5-40656*index**4+266063*index**3-549153*index**2+49850*index+655200,
        (263*index+640)*(index-3)*(2*index-7)**2]),
    'A000185-six': (5, [
        4*(210968408*index**2-1603518486*index+2343057493),
        -843873632*index**3-4039256254*index**2+144382575631*index-553812368850,
        10453330198*index**3-175111274403*index**2+798927275864*index-639098546595,
        10059264970*index**3-98879552663*index**2+170576803994*index-134993524720,
        470894110*index**3-5178116941*index**2+108179055193*index-215961878286,
        1708832970*index**3-29554327949*index**2+137453332457*index-152801514054,
        3*(569610990*index**2-3742686463*index+4740040723)]),
    'A000185-four': (5, [
        (241*index-1066)*(2*index-11)*(index-5)**2,
        -482*index**5+10099*index**4-79756*index**3+285961*index**2-426904*index+149292,
        -(2*index-9)*(index-3)*(248*index**3-2229*index**2+5065*index-7134),
        -14*index**5-49*index**4-619*index**3+13174*index**2-51690*index+61248,
        -(index-3)*(index-4)*(7*index-87)*(2*index-7)]),
}

def apply_recurrence(model, coefficients):
    answer = (symbolic.Integer(0), symbolic.Integer(0))
    powers = [model]
    for degree in range(max(symbolic.degree(coefficient,index) for coefficient in coefficients)):
        powers.append(theta(powers[-1]))
    temporary = symbolic.Symbol('v')
    for offset, coefficient in enumerate(coefficients):
        polynomial = symbolic.Poly(coefficient.subs(index,temporary+offset),temporary)
        term = (symbolic.Integer(0),symbolic.Integer(0))
        for degree in range(polynomial.degree()+1):
            term = add(term,scale(powers[degree],polynomial.nth(degree)))
        answer = add(answer,scale(term,variable**offset))
    return tuple(symbolic.factor(part) for part in answer)

if __name__ == '__main__':
    for hits in range(6):
        print('MODEL',hits,hit_model(hits),flush=True)
    for name,(hits,coefficients) in RECURRENCES.items():
        certificate = apply_recurrence(hit_model(hits),coefficients)
        print('CERTIFICATE',name,certificate,flush=True)
        assert certificate[0] == 0
        assert symbolic.denom(certificate[1]) == 1
    model = hit_model(3)
    identity = add(scale(derivative(model),1-variable),scale(model,-1-1/variable))
    other = hit_model(2)
    print('C2 RELATION',add(identity,scale(add(theta(theta(other)),scale(other,-1)),-symbolic.Rational(2,3))))
    linear = hit_model(1)
    straight = (1/(1+variable),symbolic.Integer(0))
    print('C1 SHIFT IDENTITY', add(linear,scale(add(theta(straight),straight),-2*variable)))
