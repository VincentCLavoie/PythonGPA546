from GPA546 import *

# Liste des matrices DH
H1 = dh(0, 215+d1, 0, 0)
H2 = dh(theta2, 0, 800, pi/2)
H3 = dh(theta3+(pi/2), 0, 0, pi/2)
H4 = dh(theta4, 400, 0, 0)

# Matrice de transformation H0, atelier
H0R = Matrix([
    [0, -1, 0, 0],
    [1, 0, 0, -500],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
    ])

# Matrice de transformation Houtil, H4
Ho4 = Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 100],
    [0, 0, 0, 1]
    ])

eq1 = H1*H2*H3*H4
eq2 = matP

eq = Eq(eq1, eq2)

print_matrix(eq)

print_matrix_eqs(eq2, eq1)
