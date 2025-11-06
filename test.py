from sympy import *

theta2, theta3, theta4, d1, nx, ox, ax, px, ny, oy, ay, py, nz, oz, az, pz = symbols('theta2, theta3, theta4, d1, n_x, o_x, a_x, p_x, n_y, o_y, a_y, p_y, n_z, o_z, a_z, p_z')

# s,c act as aliases to sin,cos
s = sin
c = cos

# the intended short-hand functions used to replace sin/cos later on
s1 = Function('s')
c1 = Function('c')

def dh_matrix(theta, d, a, alpha):
    ct = cos(theta)
    st = sin(theta)
    ca = cos(alpha)
    sa = sin(alpha)

    return Matrix([
        [ct, -st*ca, st*sa, a*ct],
        [st, ct*ca, -ct*sa, a*st],
        [0, sa, ca, d],
        [0, 0, 0, 1]
    ])

H1 = dh_matrix(0, 215+d1, 0, 0)
H2 = dh_matrix(theta2, 0, 800, pi/2)
H3 = dh_matrix(theta3+(pi/2), 0, 0, pi/2)
H4 = dh_matrix(theta4, 400, 0, 0)

#print("H10:")
#sp.pprint(H10);
#print("H21:")
#sp.pprint(H21);
#print("H32:")
#sp.pprint(H32);
#print("H43:")
#sp.pprint(H43);

H0R = Matrix([
    [0, -1, 0, 0],
    [1, 0, 0, -500],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
    ])

Ho4 = Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 100],
    [0, 0, 0, 1]
    ])

#print("H0R:")
#sp.pprint(H0R);
#print("Ho4:")
#sp.pprint(Ho4);

H40 = H1 * H2 * H3 * H4
H = H0R * H40 * Ho4

#print("H:")
#sp.pprint(H);

#H_outil_R = sp.Matrix([
#    [cos(ez)*cos(ey), cos(ez)*sin(ey)*sin(ex)-sin(ez)*cos(ex), cos(ez)*sin(ey)*cos(ex)+sin(ez)*sin(ex),x],
#    [sin(ez)*cos(ey), sin(ez)*sin(ey)*sin(ex)+cos(ez)*cos(ex), sin(ez)*sin(ey)*cos(ex)-cos(ez)*sin(ex), y],
#    [-sin(ey), cos(ey)*sin(ex), cos(ey)*cos(ex), z],
#    [0, 0, 0, 1]])

P = Matrix([
    [nx, px, ax, px],
    [ny, oy, ay, py],
    [nz, oz, az, pz],
    [0, 0, 0, 1]])

eq1 = simplify(H4)
eq2 = simplify((H1*H2*H3).inv()*P)
eq = Eq(eq1, eq2)

eq_latex = str(latex(eq)).replace('\sin{', 's_{').replace('\cos{','c_{').replace('\sin^{2}{','s^2_{')
print(eq_latex)

#pprint(eq.replace(s,s1).replace(c,c1))
