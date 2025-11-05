import sympy as sp
import numpy as np

theta2, theta3, theta4, d1 = sp.symbols('theta2, theta3, theta4, d1')

def dh_matrix(theta, d, a, alpha):
    ct = sp.cos(theta)
    st = sp.sin(theta)
    ca = sp.cos(alpha)
    sa = sp.sin(alpha)

    return sp.Matrix([
        [ct, -st*ca, st*sa, a*ct],
        [st, ct*ca, -ct*sa, a*st],
        [0, sa, ca, d],
        [0, 0, 0, 1]
    ])

H10 = dh_matrix(0, 215+d1, 0, 0)
H21 = dh_matrix(theta2, 0, 800, 90)
H32 = dh_matrix(theta3+90, 0, 0, 90)
H43 = dh_matrix(theta4, 400, 0, 0)

#print("H10:")
#sp.pprint(H10);
#print("H21:")
#sp.pprint(H21);
#print("H32:")
#sp.pprint(H32);
#print("H43:")
#sp.pprint(H43);

H0R = sp.Matrix([
    [0, -1, 0, 0],
    [1, 0, 0, -500],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
    ])

Ho4 =  sp.Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 100],
    [0, 0, 0, 1]
    ])

#print("H0R:")
#sp.pprint(H0R);
#print("Ho4:")
#sp.pprint(Ho4);

H40 = H10 * H21 * H32 * H43
H = H0R * H40 * Ho4

#print("H:")
#sp.pprint(H);

H_outil_R = sp.Matrix([
    [cos(ez)*cos(ey), cos(ez)*sin(ey)*sin(ex)-sin(ez)*cos(ex), cos(ez)*sin(ey)*cos(ex)+sin(ez)*sin(ex),x],
    [sin(ez)*cos(ey), sin(ez)*sin(ey)*sin(ex)+cos(ez)*cos(ex), sin(ez)*sin(ey)*cos(ex)-cos(ez)*sin(ex), y],
    [-sin(ey), cos(ey)*sin(ex), cos(ey)*cos(ex), z],
    [0, 0, 0, 1]])

a1 = h10^-1
