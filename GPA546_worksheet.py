import sympy as sp
import numpy as np
from matrix_utils import submat, dh, rx, ry, rz, txyz
import streamlit as st

# ==============================================================================
# PAGE CONFIGURATION
# ==============================================================================
st.set_page_config(layout="wide", initial_sidebar_state="expanded")

# ==============================================================================
# SYMBOL DEFINITIONS
# ==============================================================================
theta1, theta2, theta3, theta4 = sp.symbols('theta1 theta2 theta3 theta4')
d1, d2, d3, d4 = sp.symbols('d1 d2 d3 d4')
x, y, z = sp.symbols('x y z')
ex, ey, ez = sp.symbols('ex ey ez')
nx, ny, nz = sp.symbols('n_x n_y n_z')
ox, oy, oz = sp.symbols('o_x o_y o_z')
ax, ay, az = sp.symbols('a_x a_y a_z')
px, py, pz = sp.symbols('p_x p_y p_z')

# ==============================================================================
# ROBOTICS SHORTHAND SYMBOLS + LATEX NAMES
# ==============================================================================
s1, s2, s3, s4 = sp.symbols('s1 s2 s3 s4')
c1, c2, c3, c4 = sp.symbols('c1 c2 c3 c4')

symbol_names = {
    s1: 's_{1}',  c1: 'c_{1}',
    s2: 's_{2}',  c2: 'c_{2}',
    s3: 's_{3}',  c3: 'c_{3}',
    s4: 's_{4}',  c4: 'c_{4}',
    sp.Symbol('s12'): 's_{12}',   sp.Symbol('c12'): 'c_{12}',
    sp.Symbol('s13'): 's_{13}',   sp.Symbol('c13'): 'c_{13}',
    sp.Symbol('s14'): 's_{14}',   sp.Symbol('c14'): 'c_{14}',
    sp.Symbol('s23'): 's_{23}',   sp.Symbol('c23'): 'c_{23}',
    sp.Symbol('s24'): 's_{24}',   sp.Symbol('c24'): 'c_{24}',
    sp.Symbol('s34'): 's_{34}',   sp.Symbol('c34'): 'c_{34}',
    sp.Symbol('s123'): 's_{123}', sp.Symbol('c123'): 'c_{123}',
    sp.Symbol('s124'): 's_{124}', sp.Symbol('c124'): 'c_{124}',
    sp.Symbol('s134'): 's_{134}', sp.Symbol('c134'): 'c_{134}',
    sp.Symbol('s234'): 's_{234}', sp.Symbol('c234'): 'c_{234}',
    sp.Symbol('s1234'): 's_{1234}', sp.Symbol('c1234'): 'c_{1234}',
}

def robotic_trig_shorthand(expr):
    """Replace sin/cos(theta combinations) by s.. / c.. Symbols."""
    repl = {
        # singles
        sp.sin(theta1): s1, sp.cos(theta1): c1,
        sp.sin(theta2): s2, sp.cos(theta2): c2,
        sp.sin(theta3): s3, sp.cos(theta3): c3,
        sp.sin(theta4): s4, sp.cos(theta4): c4,
        # 2-angle sums
        sp.sin(theta1 + theta2): sp.Symbol('s12'), sp.cos(theta1 + theta2): sp.Symbol('c12'),
        sp.sin(theta1 + theta3): sp.Symbol('s13'), sp.cos(theta1 + theta3): sp.Symbol('c13'),
        sp.sin(theta1 + theta4): sp.Symbol('s14'), sp.cos(theta1 + theta4): sp.Symbol('c14'),
        sp.sin(theta2 + theta3): sp.Symbol('s23'), sp.cos(theta2 + theta3): sp.Symbol('c23'),
        sp.sin(theta2 + theta4): sp.Symbol('s24'), sp.cos(theta2 + theta4): sp.Symbol('c24'),
        sp.sin(theta3 + theta4): sp.Symbol('s34'), sp.cos(theta3 + theta4): sp.Symbol('c34'),
        # 3-angle sums
        sp.sin(theta1 + theta2 + theta3): sp.Symbol('s123'), sp.cos(theta1 + theta2 + theta3): sp.Symbol('c123'),
        sp.sin(theta1 + theta2 + theta4): sp.Symbol('s124'), sp.cos(theta1 + theta2 + theta4): sp.Symbol('c124'),
        sp.sin(theta1 + theta3 + theta4): sp.Symbol('s134'), sp.cos(theta1 + theta3 + theta4): sp.Symbol('c134'),
        sp.sin(theta2 + theta3 + theta4): sp.Symbol('s234'), sp.cos(theta2 + theta3 + theta4): sp.Symbol('c234'),
        # 4-angle sum
        sp.sin(theta1 + theta2 + theta3 + theta4): sp.Symbol('s1234'),
        sp.cos(theta1 + theta2 + theta3 + theta4): sp.Symbol('c1234'),
    }
    return expr.xreplace(repl)

def latex_robotic(expr):
    disp = robotic_trig_shorthand(expr)
    return sp.latex(disp, symbol_names=symbol_names)

# ==============================================================================
# ROBOT PARAMETERS
# ==============================================================================
n_joints = 4
robconfig = ("R", "P", "P","R")
H10 = dh(theta1, 0, 0, 0)
H21 = dh(0, d2+405, 148, -sp.pi/2)
H32 = dh(0, d3+279.5, 0, sp.pi/2)
H43 = dh(theta4, 0, 0, -sp.pi/2)


H0a=txyz(-600,0,0)*rz(-sp.pi/2)
Ho4=txyz(0,0,58+75)*rz(sp.pi/2)

"""
H0a = sp.Matrix([
    [0, -1, 0, 0],
    [1, 0, 0, 600],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

Ho4 = sp.Matrix([
    [0, 1, 0, 0],
    [-1, 0, 0, 0],
    [0, 0, 1, 90],
    [0, 0, 0, 1]
])

"""

H_outil_R = sp.Matrix([
    [sp.cos(ez)*sp.cos(ey), sp.cos(ez)*sp.sin(ey)*sp.sin(ex)-sp.sin(ez)*sp.cos(ex), sp.cos(ez)*sp.sin(ey)*sp.cos(ex)+sp.sin(ez)*sp.sin(ex), x],
    [sp.sin(ez)*sp.cos(ey), sp.sin(ez)*sp.sin(ey)*sp.sin(ex)+sp.cos(ez)*sp.cos(ex), sp.sin(ez)*sp.sin(ey)*sp.cos(ex)-sp.cos(ez)*sp.sin(ex), y],
    [-sp.sin(ey), sp.cos(ey)*sp.sin(ex), sp.cos(ey)*sp.cos(ex), z],
    [0, 0, 0, 1]
])

P = H0a.inv() * H_outil_R * Ho4.inv()

P_sym = sp.Matrix([
    [nx, ox, ax, px],
    [ny, oy, ay, py],
    [nz, oz, az, pz],
    [0, 0, 0, 1]
])

# ==============================================================================
# CACHED KINEMATICS COMPUTATION
# ==============================================================================

@st.cache_data
def compute_kinematics():
    """Compute all kinematics matrices once and cache them."""
    H40 = H10 * H21 * H32 * H43
    Hoa = H0a * H40 * Ho4
    H30 = H10 * H21 * H32
    H20 = H10 * H21
    H42 = H32 * H43
    H41 = H21 * H32 * H43
    H31 = H21 * H32

    return {
        'H40': H40,
        'Hoa': Hoa,
        'H30': H30,
        'H20': H20,
        'H42': H42,
        'H41': H41,
        'H31': H31,
    }

@st.cache_data
def compute_jacobian():
    """Compute Jacobian once and cache it."""
    kinematics = compute_kinematics()
    Hoa = kinematics['Hoa']
    
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
            ja = sp.Matrix([0, 0, 0])
        else:
            raise ValueError(f"Unknown joint type at index {i}: {robconfig[i]}")
        jl_cols.append(jl)
        ja_cols.append(ja)

    J_matrix_top = sp.Matrix.hstack(*jl_cols)
    J_matrix_bottom = sp.Matrix.hstack(*ja_cols)
    J = sp.Matrix.vstack(J_matrix_top, J_matrix_bottom)

    return J

# Call once to cache
kinematics = compute_kinematics()
H40 = kinematics['H40']
H30 = kinematics['H30']
H20 = kinematics['H20']
H42 = kinematics['H42']
H41 = kinematics['H41']
H31 = kinematics['H31']
Hoa = kinematics['Hoa']

J = compute_jacobian()

# ==============================================================================
# CACHED MATRIX EQUALITIES
# ==============================================================================

@st.cache_data
def build_matrix_equalities():
    """Build all matrix equalities once and cache them."""
    return {
        "H40 = P_sym":              (H40, P_sym),
        "H41 = H10.inv() * P_sym":  (H41, H10.inv() * P_sym),
        "H30 = P_sym * H43.inv()":  (H30, P_sym * H43.inv()),
        "H42 = H20.inv() * P_sym":  (H42, H20.inv() * P_sym),
        "H43 = H30.inv() * P_sym":  (H43, H30.inv() * P_sym),
        "H20 = P_sym * H42.inv()":  (H20, P_sym * H42.inv()),
        "H10 = P_sym * H41.inv()":  (H10, P_sym * H41.inv()),
        "H31 = H10.inv() * P_sym * H43.inv()": (H31, H10.inv() * P_sym * H43.inv()),
        "H21 = H10.inv() * P_sym * H42.inv()": (H21, H10.inv() * P_sym * H42.inv()),
        "H32 = H20.inv() * P_sym * H43.inv()": (H32, H20.inv() * P_sym * H43.inv()),
        "Hoa": (Hoa, Hoa),
        "Ho4": (Ho4, Ho4),
        "H0a": (H0a, H0a),
        "Jacobien": (J,J),


    }

matrix_equalities = build_matrix_equalities()

def break_matrix_equality(lhs, rhs):
    """Return list of scalar equations Eq(lhs[i,j], rhs[i,j]) for all 16 entries."""
    eqs = []
    rows, cols = lhs.shape
    assert rows == 4 and cols == 4
    assert rhs.shape == lhs.shape
    for i in range(rows):
        for j in range(cols):
            eqs.append(sp.Eq(lhs[i, j], rhs[i, j]))
    return eqs

# ==============================================================================
# STREAMLIT UI
# ==============================================================================

# Top solution box
solution_box = st.container()

st.title("Robot Kinematics Equations")

# 1) Choose which matrix equality
eq_name = st.selectbox("Select a matrix equality", list(matrix_equalities.keys()))
lhs, rhs = matrix_equalities[eq_name]

st.subheader("Matrix equality")
st.latex(latex_robotic(sp.trigsimp(lhs)) + " = " + latex_robotic(sp.trigsimp(rhs)))

# 2) Element-wise equalities (always shown) with up to two selectable equations
st.subheader("Element-wise equalities")
scalar_eqs = break_matrix_equality(lhs, rhs)

selected_eq_indices = []

for idx, eq in enumerate(scalar_eqs):
    i, j = divmod(idx, 4)
    col_checkbox, col_eq = st.columns([1, 9])
    with col_checkbox:
        checked = st.checkbox(
            f"({i+1},{j+1})",
            key=f"eq_checkbox_{eq_name}_{idx}",
        )
        if checked:
            selected_eq_indices.append(idx)
    with col_eq:
        if eq is True or eq == True:
            st.latex(r"0 = 0")
        else:
            st.latex(latex_robotic(sp.trigsimp(eq.lhs)) + " = " + latex_robotic(sp.trigsimp(eq.rhs)))

# Limit to at most two equations
if len(selected_eq_indices) > 2:
    selected_eq_indices = selected_eq_indices[:2]

# ==============================================================================
# SIDEBAR: SOLVE BUTTON + VARIABLE SELECTION
# ==============================================================================

with st.sidebar:
    st.header("Solve")
    solve_clicked = st.button("Solve", use_container_width=True)

    st.markdown("**Variables to solve for**")

    variable_options = [
        (d1, r"d_1"), (d2, r"d_2"), (d3, r"d_3"), (d4, r"d_4"),
        (sp.sin(theta1), r"\sin(\theta_1)"),
        (sp.cos(theta1), r"\cos(\theta_1)"),
        (sp.tan(theta1), r"\tan(\theta_1)"),
        (sp.sin(theta2), r"\sin(\theta_2)"),
        (sp.cos(theta2), r"\cos(\theta_2)"),
        (sp.tan(theta2), r"\tan(\theta_2)"),
        (sp.sin(theta3), r"\sin(\theta_3)"),
        (sp.cos(theta3), r"\cos(\theta_3)"),
        (sp.tan(theta3), r"\tan(\theta_3)"),
        (sp.sin(theta4), r"\sin(\theta_4)"),
        (sp.cos(theta4), r"\cos(\theta_4)"),
        (sp.tan(theta4), r"\tan(\theta_4)"),
    ]

    if "selected_vars" not in st.session_state:
        st.session_state["selected_vars"] = {label: False for _, label in variable_options}

    selected_vars = []
    for sym, label in variable_options:
        checked = st.checkbox(
            f"${label}$",
            value=st.session_state["selected_vars"][label],
            key=f"var_{label}",
        )
        st.session_state["selected_vars"][label] = checked
        if checked:
            selected_vars.append(sym)

# ==============================================================================
# SOLVE AND DISPLAY RESULTS AT THE TOP
# ==============================================================================

with solution_box:
    if solve_clicked:
        if not selected_eq_indices:
            st.warning("Select at least one element-wise equation to solve from.")
        elif not selected_vars:
            st.warning("Select at least one variable to solve for in the sidebar.")
        else:
            st.markdown("### Solutions")
            with st.container(border=True):
                # Collect non-trivial selected equations
                selected_eqs = []
                for eq_idx in selected_eq_indices:
                    eq = scalar_eqs[eq_idx]
                    i, j = divmod(eq_idx, 4)
                    if eq is True or eq == True:
                        st.write(f"Equation ({i+1},{j+1}) is identically true; skipping.")
                        continue
                    selected_eqs.append((eq_idx, eq))

                if not selected_eqs:
                    st.write("No non-trivial equations selected.")
                else:
                    # Use all selected equations together as a system
                    eq_list = [eq for _, eq in selected_eqs]

                    for target in selected_vars:
                        sol = sp.solve(eq_list, target, dict=True)
                        if sol:
                            st.latex(
                                sp.latex(target)
                                + " = "
                                + latex_robotic(sp.trigsimp(sol[0][target]))
                            )
                        else:
                            idx_str = ", ".join(
                                f"({i+1},{j+1})"
                                for eq_idx, _ in selected_eqs
                                for i, j in [divmod(eq_idx, 4)]
                            )
                            st.write(
                                f"No solution for {sp.latex(target)} "
                                f"using equations {idx_str}."
                            )
