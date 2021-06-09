from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
K.<t> = NumberField(x^3-x+1)
P.<d> = NumberField(7*x^3 + x^2 + x - 1)
f = NumberFieldEmbedding(K, CLF, 1)
g = NumberFieldEmbedding(P, CLF, -1+I)
KC = ComplexBallField(255)
tg = KC(f(t))
td = KC(g(d))
L = 2 * arccosh(tg / 2)
R = arccosh(td)
params = {}
params['sinhdx'] = sinh(R.real() / 2)
params['sinhdy'] = sinh(R.real() / 2)
params['sintx2'] = sin(L.imag() / 2)
params['sinty2'] = sin(L.imag() / 2)
params['cosf'] = cos(R.imag())
def four_cosh_mu_sym(tg, td):
    coshP = td
    sinhP = sqrt(coshP^2 - 1)
    cosh2reP = abs(coshP)^2 + abs(sinhP)^2
    sinh2reP = sqrt(cosh2reP^2 - 1)
    x1 = abs(tg)^2
    x2 = abs(tg^2 - 4)
    return x1 + sinh2reP * x2 * x2 / sqrt(2 * (cosh2reP - 1) * x2 * x2)
params['coshmu'] = four_cosh_mu_sym(tg, td).real() / 4
F = RealBallField(255)
scale = list(map(lambda x : 8 * pow(F(2), F(x) / 6), range(0,-6,-1)))
from box_codes import get_box_codes
get_box_codes(params, scale=scale)
