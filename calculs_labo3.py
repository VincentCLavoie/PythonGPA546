import sympy as sp

# Define symbols (all needed, including d1, x, y, z, ex, ey, ez)
theta2, theta3, theta4 = sp.symbols('theta2 theta3 theta4')
d1 = sp.symbols('d1')
x, y, z = sp.symbols('x y z')
ex, ey, ez = sp.symbols('ex ey ez')
nx, ny, nz = sp.symbols('n_x n_y n_z')
ox, oy, oz = sp.symbols('o_x o_y o_z')
ax, ay, az = sp.symbols('a_x a_y a_z')
px, py, pz = sp.symbols('p_x p_y p_z')

# Standard symbolic DH matrix function
def dh(theta, d, a, alpha):
    ct = sp.cos(theta)
    st = sp.sin(theta)
    ca = sp.cos(alpha)
    sa = sp.sin(alpha)
    return sp.Matrix([
        [ct, -st * ca, st * sa, a * ct],
        [st, ct * ca, -ct * sa, a * st],
        [0, sa, ca, d],
        [0, 0, 0, 1]
    ])

# Symbolic transformation matrices
H0R = sp.Matrix([
    [0, -1, 0, 0],
    [1,  0, 0, -500],
    [0,  0, 1, 0],
    [0,  0, 0, 1]
])

Ho4 = sp.Matrix([
    [1, 0, 0,   0],
    [0, 1, 0,   0],
    [0, 0, 1, 100],
    [0, 0, 0,   1]
])

H10 = dh(0, 215 + d1, 0, 0)
H21 = dh(theta2, 0, 800, sp.pi/2)
H32 = dh(theta3 + sp.pi/2, 0, 0, sp.pi/2)
H43 = dh(theta4, 400, 0, 0)

H_outil_R = sp.Matrix([
    [sp.cos(ez)*sp.cos(ey), sp.cos(ez)*sp.sin(ey)*sp.sin(ex)-sp.sin(ez)*sp.cos(ex), sp.cos(ez)*sp.sin(ey)*sp.cos(ex)+sp.sin(ez)*sp.sin(ex), x],
    [sp.sin(ez)*sp.cos(ey), sp.sin(ez)*sp.sin(ey)*sp.sin(ex)+sp.cos(ez)*sp.cos(ex), sp.sin(ez)*sp.sin(ey)*sp.cos(ex)-sp.cos(ez)*sp.sin(ex), y],
    [-sp.sin(ey), sp.cos(ey)*sp.sin(ex), sp.cos(ey)*sp.cos(ex), z],
    [0, 0, 0, 1]
])

# Example: multiply symbolic matrices (all must be SymPy matrices!)
# Combine transformations (example only)
P = H0R.inv() * H_outil_R * Ho4.inv()
 
# Construct the matrix
P_sym = sp.Matrix([
    [nx, ox, ax, px],
    [ny, oy, ay, py],
    [nz, oz, az, pz],
    [ 0,  0,  0,  1]
])


P_sym_reduced = H10.inv() * P_sym * H43.inv()
sp.pprint(P_sym_reduced)

H31=H21*H32
sp.pprint(H31)

