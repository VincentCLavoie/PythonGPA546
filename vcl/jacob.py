from GPA546 import *
from sympy import *

def compute_jacobian(H0a, Hoa, H10, H21, H32, H43, n_joints, robconfig):
    """Compute Jacobian once and cache it."""
    #kinematics = compute_kinematics()
    #Hoa = kinematics['Hoa']
    
    e0 = submat(H0a, 0, 2, 2, 2)
    e1 = submat(H0a * H10, 0, 2, 2, 2)
    e2 = submat(H0a * H10 * H21, 0, 2, 2, 2)
    e3 = submat(H0a * H10 * H21 * H32, 0, 2, 2, 2)

    p0 = submat(Hoa, 0, 3, 2, 3) - submat(H0a, 0, 3, 2, 3)
    p1 = submat(Hoa, 0, 3, 2, 3) - submat(H0a * H10, 0, 3, 2, 3)
    p2 = submat(Hoa, 0, 3, 2, 3) - submat(H0a * H10 * H21, 0, 3, 2, 3)
    p3 = submat(Hoa, 0, 3, 2, 3) - submat(H0a * H10 * H21 * H32, 0, 3, 2, 3)

    axes = [e0, e1, e2, e3]
    ps = [p0, p1, p2, p3]

    jl_cols = []
    ja_cols = []

    for i in range(n_joints):
        if robconfig[i] == "R":
            jl = axes[i].cross(ps[i])
            ja = axes[i]
        elif robconfig[i] == "P":
            jl = axes[i]
            ja = Matrix([0, 0, 0])
        else:
            raise ValueError(f"Unknown joint type at index {i}: {robconfig[i]}")
        jl_cols.append(jl)
        ja_cols.append(ja)

    J_matrix_top = Matrix.hstack(*jl_cols)
    J_matrix_bottom = Matrix.hstack(*ja_cols)
    J = Matrix.vstack(J_matrix_top, J_matrix_bottom)

    return J
