import streamlit as st
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="⚡ ComplexLab",
    page_icon="⚡",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- STYLE CONFIGURATION ---
st.markdown("""
<style>
    .main {
        background-color: #1e1e1e;
        color: #eeeeee;
    }
    .stButton>button {
        background-color: #00adb5;
        color: white;
        border: none;
        padding: 10px 24px;
        border-radius: 4px;
        font-weight: bold;
    }
    .stButton>button:hover {
        background-color: #007d85;
    }
    .stTextInput>div>div>input {
        background-color: #2e2e2e;
        color: #eeeeee;
    }
    .stSelectbox>div>div>select {
        background-color: #2e2e2e;
        color: #eeeeee;
    }
</style>
""", unsafe_allow_html=True)

# --- ENGINE FUNCTIONS (Same as original) ---

def moteur_suites_complexes(u0_str, q_str, n_str):
    u0 = sp.sympify(u0_str, locals={'i': sp.I})
    q = sp.sympify(q_str, locals={'i': sp.I})
    n = sp.sympify(n_str)
    
    un = u0 * (q ** n)
    if q == 1:
        sn = u0 * (n + 1)
    else:
        sn = u0 * (1 - q ** (n + 1)) / (1 - q)
    
    return sp.simplify(un), sp.simplify(sn)

def moteur_transformation(type_trans, z_init_str, parametre_str, centre_str="0"):
    z = sp.sympify(z_init_str, locals={'i': sp.I})
    omega = sp.sympify(centre_str, locals={'i': sp.I})
    p = sp.sympify(parametre_str, locals={'i': sp.I})
    
    if type_trans == "Translation":
        res = z + p
    elif type_trans == "Rotation":
        res = sp.exp(sp.I * p) * (z - omega) + omega
    elif type_trans == "Homothétie":
        res = p * (z - omega) + omega
    
    return sp.simplify(res)

def moteur_expo_vers_algebre(expression_expo):
    z = sp.sympify(expression_expo, locals={'i': sp.I})
    forme_alg = sp.expand(z, complex=True)
    return z, sp.simplify(forme_alg)

def moteur_racines_niemes(complexe_str, n_str):
    z = sp.sympify(complexe_str, locals={'i': sp.I})
    n = int(n_str)
    
    r = sp.Abs(z)
    theta = sp.arg(z)
    
    roots = []
    for k in range(n):
        angle = (theta + 2 * k * sp.pi) / n
        root = sp.root(r, n) * (sp.cos(angle) + sp.I * sp.sin(angle))
        roots.append(sp.simplify(root))
    
    return z, roots

def moteur_lieu_geometrique(left_side, right_side):
    x, y = sp.symbols('x y', real=True)
    z_sym = sp.Symbol('z')
    z_cart = x + sp.I * y
    
    G = sp.sympify(left_side, locals={'i': sp.I}).subs(z_sym, z_cart)
    D = sp.sympify(right_side, locals={'i': sp.I}).subs(z_sym, z_cart)
    
    sq_G = sp.expand(sp.re(G) ** 2 + sp.im(G) ** 2)
    sq_D = sp.expand(sp.re(D) ** 2 + sp.im(D) ** 2)
    
    equation = sp.simplify(sq_G - sq_D)
    
    return equation

def moteur_expansion_et_linearisation(mode, main_input, exponent=None):
    x = sp.Symbol('x', real=True)
    base = sp.sympify(main_input, locals={'i': sp.I})
    
    if mode == "binome":
        n = sp.sympify(exponent)
        original_expr = base ** n
        result = sp.expand(original_expr)
    else:
        if exponent and str(exponent).strip():
            n = sp.sympify(exponent)
            original_expr = base ** n
        else:
            original_expr = base
        
        step1 = original_expr.rewrite(sp.exp)
        step2 = sp.expand(step1)
        step3 = step2.rewrite(sp.cos)
        result = sp.expand(step3)
    
    return original_expr, result

# --- STREAMLIT APP ---

# Initialize session state for menu
if 'menu' not in st.session_state:
    st.session_state.menu = "2nd Degree Solver"

# Logo and Title
st.markdown("<h1 style='text-align: center; color: #00adb5;'>⚡ ComplexLab</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: #00adb5;'>POLYTECH IA</p>", unsafe_allow_html=True)

# Sidebar Menu
st.session_state.menu = st.sidebar.radio(
    "Select Module",
    [
        "2nd Degree Solver",
        "Quotient Z/Z'",
        "Simple Analyzer",
        "Modulus & Type",
        "Expansion & Linearization",
        "Geometry & Locus",
        "n-th Roots",
        "Exp_to_Alg Conversion",
        "Transformation",
        "Complex Sequences"
    ],
    key="menu_radio"
)

# --- 2ND DEGREE SOLVER ---
if st.session_state.menu == "2nd Degree Solver":
    st.header("Equation: ax² + bx + c = 0 (a, b, c in complex ℂ)")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        a = st.text_input("a", "1", key="solver_a")
    with col2:
        b = st.text_input("b", "0", key="solver_b")
    with col3:
        c = st.text_input("c", "0", key="solver_c")
    
    if st.button("CALCULATE / SOLVE", key="btn_solver"):
        try:
            raw_a = a.replace('i', 'I')
            raw_b = b.replace('i', 'I')
            raw_c = c.replace('i', 'I')
            
            a_val = sp.sympify(raw_a)
            b_val = sp.sympify(raw_b)
            c_val = sp.sympify(raw_c)
            z = sp.symbols('z')
            
            equation = a_val * z ** 2 + b_val * z + c_val
            st.latex(f"Equation : {sp.latex(equation)} = 0")
            
            delta = sp.simplify(b_val ** 2 - 4 * a_val * c_val)
            st.latex(f"Discriminant : \\Delta = {sp.latex(delta)}")
            
            solutions = sp.solve(equation, z)
            
            st.markdown("### --- SOLUTIONS ---")
            for i, sol in enumerate(solutions, 1):
                sol_simp = sp.simplify(sol)
                st.latex(f"z_{{{i}}} = {sp.latex(sol_simp)}")
                
                r = sp.Abs(sol_simp).simplify()
                theta = sp.arg(sol_simp).simplify()
                st.latex(f"Polar~form : {sp.latex(r)} e^{{i {sp.latex(theta)}}}")
        except Exception as e:
            st.error(f"Error : {e}")

# --- QUOTIENT Z/Z' ---
elif st.session_state.menu == "Quotient Z/Z'":
    st.header("Quotient Calculation: Z1 / Z2")
    
    col1, col2 = st.columns(2)
    with col1:
        z1 = st.text_input("Z1 (ex: 1+I)", key="quotient_z1")
    with col2:
        z2 = st.text_input("Z2 (ex: 1-I)", key="quotient_z2")
    
    if st.button("CALCULATE QUOTIENT", key="btn_quotient"):
        try:
            raw_z1 = z1.replace('i', 'I')
            raw_z2 = z2.replace('i', 'I')
            
            z1_val = sp.sympify(raw_z1)
            z2_val = sp.sympify(raw_z2)
            
            quotient_expr = z1_val / z2_val
            st.latex(f"Operation : \\frac{{{sp.latex(z1_val)}}}{{{sp.latex(z2_val)}}}")
            
            res_alg = sp.simplify(quotient_expr)
            st.latex(f"= {sp.latex(res_alg)}")
            
            r = sp.Abs(res_alg).simplify()
            theta = sp.arg(res_alg).simplify()
            st.latex(f"Exponential~form : {sp.latex(r)} e^{{i {sp.latex(theta)}}}")
        except Exception as e:
            st.error(f"Error : {e}")

# --- SIMPLE ANALYZER ---
elif st.session_state.menu == "Simple Analyzer":
    st.header("Analyzer: Real, Imaginary, Conjugate")
    
    nombre = st.text_input("Enter a complex (ex: 3 + 2*I)", key="analyzer_input")
    
    if st.button("ANALYZE", key="btn_analyze"):
        try:
            if nombre:
                raw_val = nombre.replace('i', 'I')
                nombre_val = sp.sympify(raw_val)
                
                real = sp.re(nombre_val)
                imag = sp.im(nombre_val)
                conj = sp.conjugate(nombre_val)
                
                st.latex(f"Number~z = {sp.latex(nombre_val)}")
                st.markdown("### --- RESULTS ---")
                st.latex(f"Re(z) = {sp.latex(real)}")
                st.latex(f"Im(z) = {sp.latex(imag)}")
                st.latex(f"Conjugate~\\bar{{z}} = {sp.latex(conj)}")
        except Exception as e:
            st.error(f"Error : {e}")

# --- MODULUS & TYPE ---
elif st.session_state.menu == "Modulus & Type":
    st.header("Analysis: Modulus, Argument, Exponential Form & Nature")
    
    z_type = st.text_input("Enter z (ex: 3*I, 1+I, 5)", key="type_input")
    
    if st.button("FULL DIAGNOSTIC", key="btn_type"):
        try:
            if z_type:
                raw_val = z_type.replace('i', 'I')
                z = sp.sympify(raw_val)
                r = sp.re(z)
                i = sp.im(z)
                mod = sp.Abs(z)
                angle = sp.arg(z)
                
                st.latex(f"z = {sp.latex(z)}")
                st.latex(f"|z| = {sp.latex(mod)} \\quad ; \\quad Arg(z) = {sp.latex(angle)}")
                
                st.latex(f"Trigonometric~Form : {sp.latex(mod)}(\\cos({sp.latex(angle)}) + i\\sin({sp.latex(angle)}))")
                st.latex(f"Exponential~Form : {sp.latex(mod)}(e^{{i {sp.latex(angle)}}})")
                
                st.markdown("### --- DIAGNOSTIC ---")
                if r != 0 and i == 0:
                    st.success("✅ It is a PURE REAL number because Im(z) = 0.")
                elif r == 0 and i != 0:
                    st.success("✅ It is a PURE IMAGINARY number because Re(z) = 0.")
                elif r != 0 and i != 0:
                    st.success("✅ It is a MIXED COMPLEX number (Real and Imaginary).")
                else:
                    st.warning("⚠️ The number is zero (z = 0).")
        except Exception as e:
            st.error(f"Error : {e}")

# --- EXPANSION & LINEARIZATION ---
elif st.session_state.menu == "Expansion & Linearization":
    st.header("Expansion & Linearization Pro")
    
    mode = st.radio("Select Mode", ["Binomial (ax+b)^n", "Total Linearization"], key="expansion_mode")
    
    if mode == "Binomial (ax+b)^n":
        main_input = st.text_input("Base (ex: x + 2*I)", key="expansion_base")
        exponent = st.text_input("n (exponent)", key="expansion_n")
        mode_val = "binome"
    else:
        main_input = st.text_input("Expression (ex: cos(x)**3 + sin(x)**3)", key="expansion_expr")
        exponent = None
        mode_val = "linear"
    
    if st.button("RUN AI ANALYSIS", key="btn_expansion"):
        try:
            if main_input:
                expr_org, result = moteur_expansion_et_linearisation(mode_val, main_input, exponent)
                
                st.latex(f"Input : {sp.latex(expr_org)}")
                st.markdown("### LINEARIZED / EXPANDED RESULT:")
                st.latex(sp.latex(result))
        except Exception as e:
            st.error(f"Error : {e}")

# --- GEOMETRY & LOCUS ---
elif st.session_state.menu == "Geometry & Locus":
    st.header("Locus: |f(z)| = k or |f(z)| = |g(z)|")
    
    col1, col2 = st.columns(2)
    with col1:
        left = st.text_input("Left side (ex: z - 2*I)", key="geo_left")
    with col2:
        right = st.text_input("Right side (ex: 2 or z + 1)", key="geo_right")
    
    if st.button("IDENTIFY LOCUS", key="btn_geo"):
        try:
            if left and right:
                eq_cartesian = moteur_lieu_geometrique(left, right)
                
                st.latex(f"Cartesian~Equation : {sp.latex(eq_cartesian)} = 0")
                
                x, y = sp.symbols('x y', real=True)
                try:
                    P = sp.Poly(eq_cartesian, x, y)
                    coeff_x2 = P.coeff_monomial(x**2)
                    coeff_y2 = P.coeff_monomial(y**2)
                    
                    st.markdown("### --- ANALYSIS ---")
                    
                    if coeff_x2 == 0 and coeff_y2 == 0:
                        st.success("🎯 Locus : LINE (Perpendicular bisector)")
                    elif coeff_x2 == coeff_y2 and coeff_x2 != 0:
                        eq_norm = eq_cartesian / coeff_x2
                        a = eq_norm.coeff(x, 1)
                        b = eq_norm.coeff(y, 1)
                        c = eq_norm.as_coefficients_dict()[1]
                        
                        center_x = -a / 2
                        center_y = -b / 2
                        r_sq = (a / 2) ** 2 + (b / 2) ** 2 - c
                        
                        if r_sq >= 0:
                            radius = sp.sqrt(r_sq)
                            st.success(f"🎯 CIRCLE with Center ({center_x}, {center_y}) and Radius {sp.simplify(radius)}")
                            
                            # Plot circle
                            fig, ax = plt.subplots(figsize=(6, 6))
                            circle = plt.Circle((float(center_x), float(center_y)), float(radius), color='#00adb5', fill=False, lw=2)
                            ax.add_patch(circle)
                            limit = float(max(abs(center_x), abs(center_y)) + radius + 1)
                            ax.set_xlim(-limit, limit)
                            ax.set_ylim(-limit, limit)
                            ax.axhline(0, color='white', lw=1)
                            ax.axvline(0, color='white', lw=1)
                            ax.set_aspect('equal')
                            ax.grid(True, linestyle='--', alpha=0.3)
                            fig.patch.set_facecolor('#1e1e1e')
                            ax.set_facecolor('#1e1e1e')
                            ax.tick_params(colors='white')
                            st.pyplot(fig)
                        else:
                            st.warning("⚠️ Imaginary circle (impossible radius)")
                    else:
                        st.info("ℹ️ ELLIPSE (Out of program?)")
                except:
                    st.latex(f"Cartesian~Equation : {sp.latex(eq_cartesian)} = 0")
        except Exception as e:
            st.error(f"Error : {e}")

# --- N-TH ROOTS ---
elif st.session_state.menu == "n-th Roots":
    st.header("n-th Root Calculation")
    
    col1, col2 = st.columns(2)
    with col1:
        z_root = st.text_input("Number z (ex: 1 + I)", key="roots_z")
    with col2:
        n_root = st.text_input("Value of n (ex: 3)", key="roots_n")
    
    if st.button("FIND ROOTS", key="btn_roots"):
        try:
            if z_root and n_root:
                z, roots = moteur_racines_niemes(z_root, n_root)
                
                st.latex(f"Roots~{n_root}{{-th~of~}} {sp.latex(z)} :")
                for i, r in enumerate(roots):
                    st.latex(f"z_{{{i}}} = {sp.latex(r)}")
        except Exception as e:
            st.error(f"Error : {e}")

# --- EXP TO ALG CONVERSION ---
elif st.session_state.menu == "Exp_to_Alg Conversion":
    st.header("Exponential ➔ Algebraic Form")
    
    expo = st.text_input("Ex: 2 * exp(i * pi/3)", key="conversion_input")
    
    if st.button("CONVERT", key="btn_conversion"):
        try:
            if expo:
                z, res = moteur_expo_vers_algebre(expo)
                st.latex(f"Z_{{expo}} = {sp.latex(z)}")
                st.latex(f"Z_{{alg}} = {sp.latex(res)}")
        except Exception as e:
            st.error(f"Error : {e}")

# --- TRANSFORMATION ---
elif st.session_state.menu == "Transformation":
    st.header("Plane Transformations")
    
    trans_type = st.selectbox("Transformation Type", ["Translation", "Rotation", "Homothétie"], key="trans_type")
    
    z_trans = st.text_input("Affix z (starting point)", key="trans_z")
    param = st.text_input("Vector b / Angle θ / Ratio k", key="trans_param")
    center = st.text_input("Center Ω (if Rotation/Homothety)", "0", key="trans_center")
    
    if st.button("CALCULATE IMAGE z'", key="btn_transform"):
        try:
            if z_trans and param:
                res = moteur_transformation(trans_type, z_trans, param, center)
                st.latex(f"The~image~of~the~point~is~z' = {sp.latex(res)}")
        except Exception as e:
            st.error(f"Error : {e}")

# --- COMPLEX SEQUENCES ---
elif st.session_state.menu == "Complex Sequences":
    st.header("Complex Geometric Sequences")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        u0 = st.text_input("First term u0 (ex: 1+I)", key="suite_u0")
    with col2:
        q = st.text_input("Ratio q (ex: I/2)", key="suite_q")
    with col3:
        n = st.text_input("Index n (ex: 5)", key="suite_n")
    
    if st.button("CALCULATE u_n AND S_n", key="btn_suite"):
        try:
            if u0 and q and n:
                un, sn = moteur_suites_complexes(u0, q, n)
                st.latex(f"Term~u_{{n}} = {sp.latex(un)}")
                st.latex(f"Sum~S_{{n}} = {sp.latex(sn)}")
        except Exception as e:
            st.error(f"Error : {e}")

# Footer
st.markdown("<hr>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: #00adb5; font-size: 12px;'>⚡ ComplexLab v2.0 - Streamlit Edition</p>", unsafe_allow_html=True)