from sympy import *
import pyperclip

# Symboles
theta1, theta2, theta3, theta4, theta5, d1, d2, d3, d4, d5, nx, ox, ax, px, ny, oy, ay, py, nz, oz, az, pz = symbols('theta1, theta2, theta3, theta4, theta5, d1, d2, d3, d4, d5, n_x, o_x, a_x, p_x, n_y, o_y, a_y, p_y, n_z, o_z, a_z, p_z')

matP = Matrix([
    [nx, ox, ax, px],
    [ny, oy, ay, py],
    [nz, oz, az, pz],
    [0, 0, 0, 1]
])

def replace_trig_with_sc(expr, theta_syms):
    # Create s_i and c_i symbols
    N = len(theta_syms)
    s_syms = symbols(' '.join(f's{i+1}' for i in range(N)))
    c_syms = symbols(' '.join(f'c{i+1}' for i in range(N)))

    # Replace sin(theta_i) and cos(theta_i)
    repl = {sin(th): s_syms[i] for i, th in enumerate(theta_syms)}
    repl.update({cos(th): c_syms[i] for i, th in enumerate(theta_syms)})
    expr = expr.xreplace(repl)

    # Maps for grouping
    s_set, c_set = set(s_syms), set(c_syms)
    s_index_map = {s_syms[i]: str(i + 1) for i in range(N)}
    c_index_map = {c_syms[i]: str(i + 1) for i in range(N)}

    def group_mul(m):
        s_counts, c_counts, others = [], [], []
        for arg in m.args:
            if arg in s_set:
                s_counts.append((arg, 1))
            elif arg in c_set:
                c_counts.append((arg, 1))
            elif arg.is_Pow and arg.base in s_set and arg.exp.is_Integer and arg.exp > 0:
                s_counts.append((arg.base, int(arg.exp)))
            elif arg.is_Pow and arg.base in c_set and arg.exp.is_Integer and arg.exp > 0:
                c_counts.append((arg.base, int(arg.exp)))
            else:
                others.append(arg)

        new_args = []
        if sum(cnt for _, cnt in s_counts) >= 2:
            indices = []
            for base, cnt in s_counts:
                indices.extend([s_index_map[base]] * cnt)
            indices.sort(key=int)
            new_args.append(symbols('s' + ''.join(indices)))
        else:
            new_args.extend([base**cnt if cnt > 1 else base for base, cnt in s_counts])

        if sum(cnt for _, cnt in c_counts) >= 2:
            indices = []
            for base, cnt in c_counts:
                indices.extend([c_index_map[base]] * cnt)
            indices.sort(key=int)
            new_args.append(symbols('c' + ''.join(indices)))
        else:
            new_args.extend([base**cnt if cnt > 1 else base for base, cnt in c_counts])

        new_args.extend(others)
        return Mul(*new_args) if len(new_args) > 1 else new_args[0]

    return expr.replace(lambda ex: ex.is_Mul, group_mul)

# Fonction qui retourne une matrice de transformation selon des paramètres DH
def dh(theta, d, a, alpha):
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

def txyz(x, y, z):
    return Matrix([
        [1, 0, 0, x],
        [0, 1, 0, y],
        [0, 0, 1, z],
        [0, 0, 0, 1],
    ])
    
def rx(alpha):
    ca = cos(alpha)
    sa = sin(alpha)
    return Matrix([
        [1,  0,   0, 0],
        [0, ca, -sa, 0],
        [0, sa,  ca, 0],
        [0,  0,   0, 1],
    ])

def ry(beta):
    cb = cos(beta)
    sb = sin(beta)
    return Matrix([
        [ cb, 0, sb, 0],
        [  0, 1,  0, 0],
        [-sb, 0, cb, 0],
        [  0, 0,  0, 1],
    ])

def rz(gamma):
    cg = cos(gamma)
    sg = sin(gamma)
    return Matrix([
        [cg, -sg, 0, 0],
        [sg,  cg, 0, 0],
        [ 0,   0, 1, 0],
        [ 0,   0, 0, 1],
    ])

def print_matrix(matrix, do_simplify=True):

    if do_simplify:
        matrix = simplify(matrix)

    # Imprime le résultat
    disp = replace_trig_with_sc(matrix, [theta1, theta2, theta3, theta4, theta5])
    disp = simplify_sum_trig(disp, [theta1, theta2, theta3, theta4, theta5])
    pprint(disp)

    # Transforme la réponse en code LaTeX et la copy sur le presse-papier
    eq_latex = str(latex(disp))
    print("latex version printed to clipboard")
    pyperclip.copy(eq_latex)
        
    return disp


def equalities_from_matrix_eq(M_left, M_right, do_simplify=True):
    if M_left.shape != M_right.shape:
        raise ValueError(f"Shape mismatch: {M_left.shape} vs {M_right.shape}")

    rows, cols = M_left.shape
    eqs = []
    for i in range(rows):
        for j in range(cols):
            lhs = M_left[i, j]
            rhs = M_right[i, j]
            if do_simplify:
                lhs = simplify(lhs)
                rhs = simplify(rhs)
            eqs.append(Eq(lhs, rhs))
    return eqs


def print_matrix_eqs(M_left, M_right, do_simplify=True):
    if do_simplify:
        M1 = simplify(M_left)
        M2 = simplify(M_right)

    M1 = replace_trig_with_sc(M1, [theta1, theta2, theta3, theta4])
    M2 = replace_trig_with_sc(M2, [theta1, theta2, theta3, theta4])
    M1 = simplify_sum_trig(M1, [theta1, theta2, theta3, theta4, theta5])
    M2 = simplify_sum_trig(M2, [theta1, theta2, theta3, theta4, theta5])

    eqs = equalities_from_matrix_eq(M1, M2, do_simplify=do_simplify)

    # Best compatibility with SymPy's pretty printer and matrices
    for k, eq in enumerate(eqs, start=1):
        print(f"[{k}] ", end="")
        pprint(eq)



def simplify_sum_trig(expr, theta_syms):
    """
    Replace sin(theta_i + theta_j + ...) -> s_ij... and cos(...) -> c_ij...
    where indices are sorted and i,j,... correspond to positions in theta_syms.

    Rules:
    - Only rewrites when the angle is an Add of any subset of the provided theta_syms,
      each with coefficient +1 (e.g., theta1 + theta3 + theta5).
    - Order-insensitive: sin(theta3 + theta1) -> s_13.
    - If the angle contains anything else (constants, other symbols, or non-unit
      coefficients), it is left unchanged.
    """
    # Fast lookup: symbol -> 1-based index string
    index_map = {th: str(i + 1) for i, th in enumerate(theta_syms)}
    theta_set = set(theta_syms)

    def parse_indices_from_sum(arg):
        """
        Return sorted list of index strings if `arg` is a sum of unit-coefficient
        theta_k terms only; otherwise return None.
        """
        if not isinstance(arg, Add):
            return None

        indices = []
        for term in arg.args:
            # Term must be either 'theta_k' or '1*theta_k'
            if isinstance(term, Symbol) and term in theta_set:
                indices.append(index_map[term])
            else:
                # Check for unit-coefficient multiple like 1*theta_k
                if term.is_Mul and term.as_coeff_Mul()[0] == 1:
                    coeff, rest = term.as_coeff_Mul()
                    # rest must be exactly a single theta symbol
                    if rest in theta_set:
                        indices.append(index_map[rest])
                    else:
                        return None
                else:
                    return None  # Any other structure disqualifies
        # Normalize canonical order
        indices.sort(key=lambda s: int(s))
        return indices if indices else None

    def rewrite_func(f):
        """
        f is a SymPy Function like sin(<arg>) or cos(<arg>).
        If its arg matches a sum of thetas, replace with s_... or c_...
        """
        if f.func not in (sin, cos):
            return f

        indices = parse_indices_from_sum(f.args[0])
        if indices is None:
            return f  # leave as-is

        name = ('s_' if f.func is sin else 'c_') + ''.join(indices)
        return symbols(name)

    # Walk and rewrite all sin/cos nodes
    return expr.replace(lambda ex: ex.is_Function and ex.func in (sin, cos), rewrite_func)

def submat(matrix, row_start, col_start, row_end, col_end):
    # In Python, slices are [start:end), so add +1 for inclusive ending
    if isinstance(matrix, Matrix):
        return matrix[row_start:row_end+1, col_start:col_end+1]
    else:
        return matrix[row_start:row_end+1, col_start:col_end+1]
