from GPA546 import *
from jacob import *

n_joints = 4
robconfig = ("R", "P" ,"R","R")

# Liste des matrices DH
H10 = dh(theta1, 650, 0, -pi/2)
H21 = dh(-pi/2, d2, 400, pi)
H32 = dh(theta3, 0, 600, 0)
H43 = dh(theta4, 0, 420, 0)

# Matrice de transformation H0, atelier
H0R = Matrix([
    [0, -1, 0, 0],
    [1, 0, 0, -2000],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
    ])

# Matrice de transformation Houtil, H4
Ho4 = Matrix([
    [0, 0, 1, 60],
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1]
    ])

Hoa = simplify(H0R*H10*H21*H32*H43*Ho4)

J = compute_jacobian(H0R, Hoa, H10, H21, H32, H43, n_joints, robconfig)
print_matrix(J)
