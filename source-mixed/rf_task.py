import numpy as np
import Primitives as P
from scipy.optimize import root, show_options

# show_options(solver='root', method='hybr', disp=True), exit()

def fun(x):
    '''x is a 21-dimensional array with x[:11]=(11θ), x[11:18]=(7ϕ), and x[18:]=(x,y,z).
    Returns [σ12, σ34, σ13*σ24 - σ14*σ23].'''
    s = P.sigmas(x[:18],x[18:])
    return np.concatenate(( np.array([s[0], s[5], s[1]*s[4] - s[2]*s[3]]), np.zeros(18) ))

nsets = 10 # no. (x,y,z,11θ,7ϕ) sets to find s.t. fun(.)~0.

# initialize data output
h = open('data.txt','w')
h.write('The following set of points satisfy:\n')
h.write('  σ12 ~ 0.0\n')
h.write('  σ34 ~ 0.0\n')
h.write('  H := σ13*σ24 - σ14*σ23 ~ 0.0\n\n')

for n in range(nsets):
    # solve
    x0 = np.random.uniform(0., 1., 21); x0[:18] *= np.pi; x0[18:] = 10*x0[18:] - 5.;
    sol = root(fun=fun, x0=x0, args=(), method='hybr',
        jac=None, tol=0.0, callback=None, options={'maxfev':5_000})

    # write data to file
    temp = '(σ12,σ34,H) = ('+', '.join('{:.5e}'.format(elt) for elt in sol.fun[:3])+')\n'
    h.write('-'*(len(temp)-1) + '\n')
    h.write(temp)
    h.write('-'*(len(temp)-1) + '\n')
    h.write('(x,y,z) = ('+', '.join('{:.5f}'.format(elt) for elt in sol.x[18:])+')\n')
    h.write('θ/π = ('+', '.join('{:.5f}'.format(elt/np.pi) for elt in sol.x[:11])+')\n')
    if (n == nsets-1): end = '\n'
    else: end = '\n\n'
    h.write('ϕ/π = ('+', '.join('{:.5f}'.format(elt/np.pi) for elt in sol.x[11:18])+')'+end)

h.close()
