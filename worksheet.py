import sympy as sp
from matrix_utils import submat, dh, rx, ry, rz, txyz

#print(sp.simplify(sp.sqrt(157600)))
#print(260*260)

H=txyz(40,0,0)*ry(sp.pi/2)*rz(sp.pi)
sp.pprint(H)