import sympy as sp

def txyz(x, y, z):
    """Pure translation homogeneous transform."""
    return sp.Matrix([
        [1, 0, 0, x],
        [0, 1, 0, y],
        [0, 0, 1, z],
        [0, 0, 0, 1],
    ])

def rx(alpha):
    """Rotation about x (radians) as homogeneous transform."""
    ca = sp.cos(alpha)
    sa = sp.sin(alpha)
    return sp.Matrix([
        [1,  0,   0, 0],
        [0, ca, -sa, 0],
        [0, sa,  ca, 0],
        [0,  0,   0, 1],
    ])

def ry(beta):
    """Rotation about y (radians) as homogeneous transform."""
    cb = sp.cos(beta)
    sb = sp.sin(beta)
    return sp.Matrix([
        [ cb, 0, sb, 0],
        [  0, 1,  0, 0],
        [-sb, 0, cb, 0],
        [  0, 0,  0, 1],
    ])

def rz(gamma):
    """Rotation about z (radians) as homogeneous transform."""
    cg = sp.cos(gamma)
    sg = sp.sin(gamma)
    return sp.Matrix([
        [cg, -sg, 0, 0],
        [sg,  cg, 0, 0],
        [ 0,   0, 1, 0],
        [ 0,   0, 0, 1],
    ])

def dh(theta, d, a, alpha):
    ct = sp.cos(theta)
    st = sp.sin(theta)
    ca = sp.cos(alpha)
    sa = sp.sin(alpha)
    return sp.Matrix([
        [ct, -st * ca,  st * sa, a * ct],
        [st,  ct * ca, -ct * sa, a * st],
        [ 0,       sa,      ca,      d],
        [ 0,        0,       0,      1]
    ])

def submat(matrix, row_start, col_start, row_end, col_end):
    """
    Extracts a submatrix from (row_start, col_start) to (row_end, col_end), inclusive.
    Works for numpy arrays and SymPy matrices.

    Parameters:
        matrix: 2D numpy array or sympy Matrix
        row_start: starting row index (int)
        col_start: starting column index (int)
        row_end: ending row index (int, inclusive)
        col_end: ending column index (int, inclusive)

    Returns:
        The selected submatrix
    """
    # In Python, slices are [start:end), so add +1 for inclusive ending
    if isinstance(matrix, sp.Matrix):
        return matrix[row_start:row_end+1, col_start:col_end+1]
    else:
        return matrix[row_start:row_end+1, col_start:col_end+1]

"""
# Example usage with numpy matrix
mat = np.arange(16).reshape(4,4)
sm = submat(mat, 0, 2, 2, 2)  # Third column, rows 0-2
print(sm)

# Example usage with sympy Matrix
mat_sp = sp.Matrix(mat)
sm_sp = submat(mat_sp, 0, 2, 2, 2)
sp.pprint(sm_sp)
"""