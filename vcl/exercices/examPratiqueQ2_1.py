from sympy import *
import pyperclip

# Symboles
theta2, theta3, theta4, d1, d2, d3, nx, ox, ax, px, ny, oy, ay, py, nz, oz, az, pz = symbols('theta2, theta3, theta4, d1, d2, d3, n_x, o_x, a_x, p_x, n_y, o_y, a_y, p_y, n_z, o_z, a_z, p_z')

# Aliases pour réduire la grosseur des résultats
s = sin
c = cos
s1 = Function('s')
c1 = Function('c')

# Fonction qui retourne une matrice de transformation selon des paramètres DH
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

# Liste des matrices DH
H1 = dh_matrix(0, d1, 0, -(pi/2))
H2 = dh_matrix(-(pi/2), d2, 0, -(pi/2))
H3 = dh_matrix(0, d3, 0, 0)
H4 = dh_matrix(theta4, 0, 0, 0)

Houtil = Matrix([
    [0, -1, 0, -240],
    [1, 0, 0, 0],
    [0, 0, 1, 180],
    [0, 0, 0, 1]])

H0 = Matrix([
    [0, 0, 1, 0],
    [1, 0, 0, 0],
    [0, 1, 0, -2000],
    [0, 0, 0, 1]])

# Matrice P 
P = Matrix([
    [nx, ox, ax, px],
    [ny, oy, ay, py],
    [nz, oz, az, pz],
    [0, 0, 0, 1]])

############################################
# Configuration de l'égalité
# 
# Exemple: eq1 = H1*H2*H3, eq2 = P*H4.inv()

eq1 = P
eq2 = H1*H2*H3*H4
############################################

# Création de l'équalité
eq = Eq(simplify(eq1), simplify(eq2))

# Transforme la réponse en code LaTeX et la copy sur le presse-papier
eq_latex = str(latex(eq)).replace('\sin{', 's_{').replace('\cos{','c_{').replace('\sin^{2}{','s^2_{').replace('\cos^{2}{','c^2_{').replace('\left(','').replace('right)','dummytext').replace(' \dummytext','')
print("latex version printed to clipboard")
pyperclip.copy(eq_latex)

# Imprime le résultat
pprint(eq.replace(s,s1).replace(c,c1))
