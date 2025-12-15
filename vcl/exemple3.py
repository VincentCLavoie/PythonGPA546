from GPA546 import *
from jacob import *

n_joints = 4
robconfig = ("P", "P" ,"P","R")

# Liste des matrices DH
H10 = dh(0, d1, 0, -pi/2)
H21 = dh(-pi/2, d2, 0, -pi/2)
H32 = dh(0, d3, 0, 0)
H43 = dh(theta4, 0, 0, 0)

# Matrice de transformation H0, atelier
H0a = Matrix([
    [0, 0, 1, 0],
    [1, 0, 0, 0],
    [0, 1, 0, -2000],
    [0, 0, 0, 1]
    ])

# Matrice de transformation Houtil, H4
Ho4 = Matrix([
    [0, -1, 0, -240],
    [1, 0, 0, 0],
    [0, 0, 1, 180],
    [0, 0, 0, 1]
    ])

Hoa = simplify(H0a*H10*H21*H32*H43*Ho4)

J = compute_jacobian(H0a, Hoa, H10, H21, H32, H43, n_joints, robconfig)
print_matrix(J)
