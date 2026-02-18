import customtkinter as ctk
import sympy as sp
import matplotlib.pyplot as plt
from sympy.simplify.fu import TR8
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# --- STYLE CONFIGURATION ---
ctk.set_appearance_mode("Dark")  # Dark cyberpunk mode
ctk.set_default_color_theme("dark-blue")

# Custom "Mecha" look colors
COLOR_BG = "#1e1e1e"      # Very dark background
COLOR_ACCENT = "#00adb5"  # Futuristic cyan
COLOR_TEXT = "#eeeeee"    # Off-white for text


def moteur_suites_complexes(u0_str, q_str, n_str):
    u0 = sp.sympify(u0_str, locals={'i': sp.I})
    q = sp.sympify(q_str, locals={'i': sp.I})
    n = sp.sympify(n_str)

    un = u0 * (q ** n)
    # Sum Sn = u0 * (1 - q^(n+1)) / (1 - q)
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
        # z' = z + b
        res = z + p
    elif type_trans == "Rotation":
        # z' = e^(i t)(z - w) + w
        res = sp.exp(sp.I * p) * (z - omega) + omega
    elif type_trans == "Homothétie":
        # z' = k(z - w) + w
        res = p * (z - omega) + omega

    return sp.simplify(res)

def moteur_expo_vers_algebre(expression_expo):
    # Use our trick for 'i'
    z = sp.sympify(expression_expo, locals={'i': sp.I})
    # sp.expand(complex=True) forces form a + b i
    forme_alg = sp.expand(z, complex=True)
    return z, sp.simplify(forme_alg)

def moteur_racines_niemes(complexe_str, n_str):
    z = sp.sympify(complexe_str, locals={'i': sp.I})
    n = int(n_str)

    r = sp.Abs(z)
    theta = sp.arg(z)

    roots = []
    for k in range(n):
        # Formula: r^(1/n) * exp(i * (theta + 2kpi)/n)
        angle = (theta + 2 * k * sp.pi) / n
        root = sp.root(r, n) * (sp.cos(angle) + sp.I * sp.sin(angle))
        roots.append(sp.simplify(root))

    return z, roots

def moteur_lieu_geometrique(left_side, right_side):
    # Define variables
    x, y = sp.symbols('x y', real=True)
    z_sym = sp.Symbol('z')  # The symbol 'z' user types
    z_cart = x + sp.I * y  # x + i*y form

    # 1. Parsing and substitution
    # Replace user's 'z' with (x + I*y)
    G = sp.sympify(left_side, locals={'i': sp.I}).subs(z_sym, z_cart)
    D = sp.sympify(right_side, locals={'i': sp.I}).subs(z_sym, z_cart)

    # 2. Square-modulus trick (|Z|^2 = Re(Z)^2 + Im(Z)^2)
    # This removes square roots immediately
    sq_G = sp.expand(sp.re(G) ** 2 + sp.im(G) ** 2)
    sq_D = sp.expand(sp.re(D) ** 2 + sp.im(D) ** 2)

    # 3. Final equation: Left - Right = 0
    equation = sp.simplify(sq_G - sq_D)

    return equation


def moteur_expansion_et_linearisation(mode, main_input, exponent=None):
    # IMPORTANT: x must be real
    x = sp.Symbol('x', real=True)
    base = sp.sympify(main_input, locals={'i': sp.I})

    if mode == "binome":
        n = sp.sympify(exponent)
        original_expr = base ** n
        result = sp.expand(original_expr)
    else:
        # LINEARIZATION MODE
        if exponent and str(exponent).strip():
            n = sp.sympify(exponent)
            original_expr = base ** n
        else:
            original_expr = base

        # MECHANICAL METHOD (Fail-safe)
        # 1. Rewrite everything as exponentials (Euler)
        step1 = original_expr.rewrite(sp.exp)
        # 2. Expand parentheses
        step2 = sp.expand(step1)
        # 3. Convert back to Cosine forms
        step3 = step2.rewrite(sp.cos)
        # 4. Expand final result to separate terms
        result = sp.expand(step3)

        # The result will look like: 3/8 - cos(2*x)/2 + cos(4*x)/8
        # This is the correct linearized form.

    return original_expr, result

class MathDisplay(ctk.CTkFrame):
    """This frame only displays mathematical formulas in LaTeX"""

    def __init__(self, master, **kwargs):
        super().__init__(master, fg_color="transparent", **kwargs)
        self.figure = None
        self.canvas = None

    def afficher_latex(self, latex_str, font_size=12):
        # 1. Clear previous display if any
        if self.canvas:
            self.canvas.get_tk_widget().pack_forget()
            plt.close(self.figure)

        # 2. Create a transparent Matplotlib "image"
        # figsize=(width, height) in inches. dpi=100 for quality.
        self.figure = plt.figure(figsize=(5, 1), dpi=100)
        self.figure.patch.set_facecolor(COLOR_BG)  # Background same as app

        # 3. Write the math text
        ax = self.figure.add_subplot(111)
        ax.axis('off')  # Hide axes

        # Draw the formula in white centered
        ax.text(0.5, 0.5, f"${latex_str}$",
                fontsize=font_size,
                ha='center', va='center',
                color=COLOR_TEXT)

        # 4. Embed image into CustomTkinter
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(expand=True, fill='both')


class CyberMathApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("⚡ ComplexLab")
        self.geometry("1100x620")

        # Grid configuration (2 columns: Menu | Content)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # --- 1. LEFT SIDEBAR MENU ---
        self.sidebar = ctk.CTkFrame(self, width=200, corner_radius=0)
        self.sidebar.grid(row=0, column=0, sticky="nsew")

        self.logo = ctk.CTkLabel(self.sidebar, text="POLYTECH IA", font=("Orbitron", 20, "bold"))
        self.logo.pack(pady=30)

        btn_solveur = ctk.CTkButton(self.sidebar, text="2nd Degree Solver",
                                    fg_color=COLOR_ACCENT, hover_color="#007d85",
                                    command=self.show_solveur_frame)
        btn_solveur.pack(pady=10, padx=20, fill="x")

        btn_quotient = ctk.CTkButton(self.sidebar, text="Quotient Z/Z'",
                                     fg_color="transparent", hover_color="#3E444B",
                                     command=self.show_quotient_frame)
        btn_quotient.pack(pady=10, padx=20, fill="x")

        self.btn_analyse = ctk.CTkButton(self.sidebar, text="Simple Analyzer",
                                         fg_color="transparent", hover_color="#3E444B",
                                         command=self.show_analyse_frame)
        self.btn_analyse.pack(pady=10, padx=20, fill="x")

        self.btn_type = ctk.CTkButton(self.sidebar, text="Modulus & Type",
                                      fg_color="transparent", hover_color="#3E444B",
                                      command=self.show_type_frame)
        self.btn_type.pack(pady=10, padx=20, fill="x")

        self.btn_expansion = ctk.CTkButton(self.sidebar, text="Expansion",
                                           fg_color="transparent", hover_color="#3E444B",
                                           command=self.show_expansion_frame)
        self.btn_expansion.pack(pady=10, padx=20, fill="x")

        self.btn_geo = ctk.CTkButton(self.sidebar, text="Geometry & Locus",
                                     fg_color="transparent", hover_color="#3E444B",
                                     command=self.show_geo_frame)
        self.btn_geo.pack(pady=10, padx=20, fill="x")

        self.btn_geo = ctk.CTkButton(self.sidebar, text="n-th Roots",
                                     fg_color="transparent", hover_color="#3E444B",
                                     command=self.show_racines_frame)
        self.btn_geo.pack(pady=10, padx=20, fill="x")

        self.btn_conversion = ctk.CTkButton(self.sidebar, text="Exp_to_Alg Conversion",
                                           fg_color="transparent", hover_color="#3E444B",
                                            command=self.show_conversion_frame)
        self.btn_conversion.pack(pady=10, padx=20, fill="x")

        self.btn_transform = ctk.CTkButton(self.sidebar, text="Transformation",
                                           fg_color="transparent", hover_color="#3E444B",
                                           command=self.show_transform_frame)
        self.btn_transform.pack(pady=10, padx=20, fill="x")

        self.btn_suites = ctk.CTkButton(self.sidebar, text="Complex Sequences",
                                        fg_color="transparent", hover_color="#3E444B",
                                        command=self.show_suites_frame)
        self.btn_suites.pack(pady=10, padx=20, fill="x")

        # --- 2. MAIN AREA (Right) ---
        self.main_area = ctk.CTkScrollableFrame(self, fg_color="transparent")
        self.main_area.grid(row=0, column=1, sticky="nsew", padx=20, pady=20)

        self.setup_solveur_ui()
        self.setup_quotient_ui()
        self.show_solveur_frame()
        self.setup_analyse_ui()
        self.setup_type_ui()
        self.setup_expansion_ui()
        self.setup_racines_ui()
        self.setup_conversion_ui()
        self.setup_transform_ui()
        self.setup_suites_ui()

    def setup_solveur_ui(self):
        """Draw the inputs for a, b, c and the button"""

        # Title
        title = ctk.CTkLabel(self.main_area, text="Equation: ax² + bx + c = 0 (a, b, c in complex C)", font=("Arial", 22, "bold"))
        title.pack(pady=(0, 20))

        # Container for entries (horizontal)
        input_frame = ctk.CTkFrame(self.main_area, fg_color="transparent")
        input_frame.pack(pady=10)

        self.entry_a = self.creer_input(input_frame, "a")
        self.entry_b = self.creer_input(input_frame, "b")
        self.entry_c = self.creer_input(input_frame, "c")

        # Calculate button
        btn_calc = ctk.CTkButton(self.main_area, text="CALCULATE / SOLVE",
                                 font=("Arial", 16, "bold"),
                                 height=40,
                                 fg_color=COLOR_ACCENT,
                                 command=self.lancer_calcul_solveur)
        btn_calc.pack(pady=20)

        # Result area (LaTeX images)
        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def creer_input(self, master, placeholder):
        """Helper to create small nice entry boxes quickly"""
        entry = ctk.CTkEntry(master, placeholder_text=placeholder, width=80, justify="center")
        entry.pack(side="left", padx=10)
        return entry

    # =========================================================================
    # THE CORE: YOUR SOLVER LOGIC ADAPTED
    # =========================================================================
    def lancer_calcul_solveur(self):
        # 1. Clear old results
        for widget in self.result_container.winfo_children():
            widget.destroy()

        try:
            # Retrieve values (replace input())
            raw_a = self.entry_a.get().replace('i', 'I')
            raw_b = self.entry_b.get().replace('i', 'I')
            raw_c = self.entry_c.get().replace('i', 'I')

            # SymPy conversion
            a = sp.sympify(raw_a)
            b = sp.sympify(raw_b)
            c = sp.sympify(raw_c)
            z = sp.symbols('z')

            # --- STEP 1: Show equation ---
            equation = a * z ** 2 + b * z + c
            # Use sp.latex() to convert math object to LaTeX
            self.afficher_resultat_math(f"Equation : {sp.latex(equation)} = 0")

            # --- STEP 2: Discriminant ---
            delta = sp.simplify(b ** 2 - 4 * a * c)
            self.afficher_resultat_math(f"Discriminant : \\Delta = {sp.latex(delta)}")

            # --- STEP 3: Solutions ---
            solutions = sp.solve(equation, z)

            # Simple title text
            lbl_sol = ctk.CTkLabel(self.result_container, text="--- SOLUTIONS ---", text_color=COLOR_ACCENT)
            lbl_sol.pack(pady=10)

            for i, sol in enumerate(solutions, 1):
                sol_simp = sp.simplify(sol)
                # Algebraic form
                self.afficher_resultat_math(f"z_{i} = {sp.latex(sol_simp)}")

                # Polar form (optional)
                r = sp.Abs(sol_simp).simplify()
                theta = sp.arg(sol_simp).simplify()
                self.afficher_resultat_math(f"Polar form : {sp.latex(r)} e^{{i {sp.latex(theta)}}}")

        except Exception as e:
            err_label = ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red")
            err_label.pack()

    def afficher_resultat_math(self, latex_code):
        """Create a new line with MathDisplay rendering"""
        frame = MathDisplay(self.result_container)
        frame.pack(fill="x", pady=5)
        frame.afficher_latex(latex_code)

    def nettoyer_zone_principale(self):
        """Clear the right-hand area before showing a new module"""
        for widget in self.main_area.winfo_children():
            widget.destroy()

    def show_solveur_frame(self):
        self.nettoyer_zone_principale()  # 1. Clear everything (Quotient or other)
        self.setup_solveur_ui()  # 2. Redraw the Solver

    def show_quotient_frame(self):
        self.nettoyer_zone_principale()
        self.setup_quotient_ui()

    def show_analyse_frame(self):
        self.nettoyer_zone_principale()
        self.setup_analyse_ui()

    def show_type_frame(self):
        self.nettoyer_zone_principale()
        self.setup_type_ui()

    def show_expansion_frame(self):
        self.nettoyer_zone_principale()
        self.setup_expansion_ui()

    def show_geo_frame(self):
        self.nettoyer_zone_principale()
        self.setup_geo_ui()

    def show_racines_frame(self):
        self.nettoyer_zone_principale()
        self.setup_racines_ui()

    def show_conversion_frame(self):
        self.nettoyer_zone_principale()
        self.setup_conversion_ui()

    def show_transform_frame(self):
        self.nettoyer_zone_principale()
        self.setup_transform_ui()

    def show_suites_frame(self):
        self.nettoyer_zone_principale()
        self.setup_suites_ui()

    def setup_quotient_ui(self):
        # Title
        title = ctk.CTkLabel(self.main_area, text="Quotient Calculation: Z1 / Z2", font=("Arial", 22, "bold"))
        title.pack(pady=(0, 20))

        # Entry zone
        input_frame = ctk.CTkFrame(self.main_area, fg_color="transparent")
        input_frame.pack(pady=10)

        # Keep entries on self to retrieve later
        self.entry_z1 = self.creer_input(input_frame, "Z1 (ex: 1+I)")
        self.entry_z2 = self.creer_input(input_frame, "Z2 (ex: 1-I)")

        # Calculate button
        btn_calc = ctk.CTkButton(self.main_area, text="CALCULATE QUOTIENT",
                                 fg_color=COLOR_ACCENT,
                                 command=self.lancer_calcul_quotient)
        btn_calc.pack(pady=20)

        # Result area
        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def lancer_calcul_quotient(self):
        # 1. Clear previous results
        for widget in self.result_container.winfo_children():
            widget.destroy()

        try:
            # Retrieval
            raw_z1 = self.entry_z1.get().replace('i', 'I')
            raw_z2 = self.entry_z2.get().replace('i', 'I')

            # SymPy calculations
            z1 = sp.sympify(raw_z1)
            z2 = sp.sympify(raw_z2)

            # Show the operation
            quotient_expr = z1 / z2
            self.afficher_resultat_math(f"Operation : \\frac{{{sp.latex(z1)}}}{{{sp.latex(z2)}}}")

            # Simplified result (algebraic form)
            res_alg = sp.simplify(quotient_expr)
            self.afficher_resultat_math(f"= {sp.latex(res_alg)}")

            # Exponential form (bonus)
            r = sp.Abs(res_alg).simplify()
            theta = sp.arg(res_alg).simplify()
            self.afficher_resultat_math(f"Exponential form : {sp.latex(r)} e^{{i {sp.latex(theta)}}}")

        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

    def setup_analyse_ui(self):
        # Title
        title = ctk.CTkLabel(self.main_area, text="Analyzer: Real, Imaginary, Conjugate",
                             font=("Arial", 22, "bold"))
        title.pack(pady=(0, 20))

        # Input box
        self.entry_nombre = ctk.CTkEntry(self.main_area,
                                         placeholder_text="Enter a complex (ex: 3 + 2*I)",
                                         width=350, height=40)
        self.entry_nombre.pack(pady=10)

        # Analyze button
        btn_calc = ctk.CTkButton(self.main_area, text="ANALYZE",
                                 fg_color=COLOR_ACCENT,
                                 command=self.lancer_analyse_complexe)
        btn_calc.pack(pady=20)

        # LaTeX results area
        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def lancer_analyse_complexe(self):
        # 1. Clear
        for widget in self.result_container.winfo_children():
            widget.destroy()

        try:
            # 2. Get input
            raw_val = self.entry_nombre.get().replace('i', 'I')
            if not raw_val: return

            nombre = sp.sympify(raw_val)

            # 3. SymPy calculations
            real = sp.re(nombre)
            imag = sp.im(nombre)
            conj = sp.conjugate(nombre)

            # 4. Styled LaTeX display
            # Show the initial number
            self.afficher_resultat_math(f"Number~z = {sp.latex(nombre)}")

            # Small visual separator
            lbl = ctk.CTkLabel(self.result_container, text="--- RESULTS ---", text_color=COLOR_ACCENT)
            lbl.pack(pady=10)

            # Real part
            self.afficher_resultat_math(f"Re(z) = {sp.latex(real)}")

            # Imaginary part
            self.afficher_resultat_math(f"Im(z) = {sp.latex(imag)}")

            # Conjugate
            self.afficher_resultat_math(f"Conjugate~\\bar{{z}} = {sp.latex(conj)}")

        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

    def setup_type_ui(self):
        title = ctk.CTkLabel(self.main_area, text="Analysis: Modulus, Argument, Exponential Form & Nature", font=("Arial", 22, "bold"))
        title.pack(pady=(0, 20))

        self.entry_z_type = ctk.CTkEntry(self.main_area, placeholder_text="Enter z (ex: 3*I, 1+I, 5)", width=350)
        self.entry_z_type.pack(pady=10)

        btn_calc = ctk.CTkButton(self.main_area, text="FULL DIAGNOSTIC", fg_color=COLOR_ACCENT,
                                 command=self.lancer_analyse_type)
        btn_calc.pack(pady=20)

        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def lancer_analyse_type(self):
        for widget in self.result_container.winfo_children():
            widget.destroy()

        try:
            raw_val = self.entry_z_type.get().replace('i', 'I')
            if not raw_val: return

            z = sp.sympify(raw_val)
            r = sp.re(z)
            i = sp.im(z)
            mod = sp.Abs(z)
            angle = sp.arg(z)

            # --- MATH DISPLAY ---
            self.afficher_resultat_math(f"z = {sp.latex(z)}")
            self.afficher_resultat_math(f"|z| = {sp.latex(mod)} \\quad ; \\quad Arg(z) = {sp.latex(angle)}")

            # Trigonometric form for style
            self.afficher_resultat_math(
                f"Trigonometric~Form : {sp.latex(mod)}(\\cos({sp.latex(angle)}) + i\\sin({sp.latex(angle)}))")

            self.afficher_resultat_math(
                f"Exponential~Form : {sp.latex(mod)}(e^{{i {sp.latex(angle)}}})")
            # --- DIAGNOSTIC LOGIC (If-statements) ---
            diagnostic = ""
            if r != 0 and i == 0:
                diagnostic = f"✅ It is a PURE REAL number because Im(z) = 0."
            elif r == 0 and i != 0:
                diagnostic = f"✅ It is a PURE IMAGINARY number because Re(z) = 0."
            elif r != 0 and i != 0:
                diagnostic = f"✅ It is a MIXED COMPLEX number (Real and Imaginary)."
            else:
                diagnostic = f"⚠️ The number is zero (z = 0)."

            # Show diagnostic in a styled label
            lbl_diag = ctk.CTkLabel(self.result_container, text=diagnostic,
                                    font=("Arial", 16, "italic"), text_color=COLOR_ACCENT)
            lbl_diag.pack(pady=20)

        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

    def setup_expansion_ui(self):
        self.nettoyer_zone_principale()

        title = ctk.CTkLabel(self.main_area, text="Expansion & Linearization Pro", font=("Arial", 22, "bold"))
        title.pack(pady=20)

        # --- MODE SELECTOR ---
        self.mode_var = ctk.StringVar(value="binome")
        radio_frame = ctk.CTkFrame(self.main_area, fg_color="transparent")
        radio_frame.pack(pady=10)

        # Add "command=self.update_inputs" to update display live
        ctk.CTkRadioButton(radio_frame, text="Binomial (ax+b)^n", variable=self.mode_var,
                           value="binome", command=self.update_inputs).pack(side="left", padx=20)
        ctk.CTkRadioButton(radio_frame, text="Total Linearization", variable=self.mode_var,
                           value="linear", command=self.update_inputs).pack(side="left", padx=20)

        # --- INPUT ZONE ---
        self.input_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.input_container.pack(pady=10)

        # Main entry (Base or Expression)
        self.entry_main = ctk.CTkEntry(self.input_container, width=400, height=40)
        self.entry_main.pack(pady=5)
        self.label_help = ctk.CTkLabel(self.input_container, text="Enter the base (ex: 2*x + I)", font=("Arial", 11))
        self.label_help.pack()

        # Exponent entry (can be hidden)
        self.entry_n = ctk.CTkEntry(self.input_container, placeholder_text="n", width=60)
        self.entry_n.pack(pady=10)

        # Action button
        btn_calc = ctk.CTkButton(self.main_area, text="RUN AI ANALYSIS",
                                 fg_color=COLOR_ACCENT, command=self.lancer_calcul_expansion)
        btn_calc.pack(pady=20)

        # Result container
        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

        # Initialize correct display
        self.update_inputs()

    def update_inputs(self):
        """Change help texts and hide exponent if needed"""
        if self.mode_var.get() == "binome":
            self.entry_main.configure(placeholder_text="Base (ex: x + 2*I)")
            self.label_help.configure(text="The AI will compute (Base)^n")
            self.entry_n.pack(pady=10)  # Show n entry
        else:
            self.entry_main.configure(placeholder_text="Expression (ex: cos(x)**3 + sin(x)**3)")
            self.label_help.configure(text="The AI will linearize the entire expression")
            self.entry_n.pack_forget()  # Hide n entry

    def lancer_calcul_expansion(self):
        for widget in self.result_container.winfo_children():
            widget.destroy()

        try:
            mode = self.mode_var.get()
            main_val = self.entry_main.get()
            val_n = self.entry_n.get() if mode == "binome" else None

            if not main_val: return

            # Call engine
            expr_org, result = moteur_expansion_et_linearisation(mode, main_val, val_n)

            # LaTeX rendering
            self.afficher_resultat_math(f"Input : {sp.latex(expr_org)}")

            lbl_res = ctk.CTkLabel(self.result_container, text="LINEARIZED / EXPANDED RESULT:",
                                   font=("Arial", 13, "bold"), text_color=COLOR_ACCENT)
            lbl_res.pack(pady=10)

            expr_org, result = moteur_expansion_et_linearisation(mode, main_val, val_n)

            # LaTeX rendering
            self.afficher_resultat_math(f"Input : {sp.latex(expr_org)}")
            self.afficher_resultat_math(sp.latex(result))

        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

    def setup_geo_ui(self):
        """Interface for geometric loci"""
        self.nettoyer_zone_principale()

        title = ctk.CTkLabel(self.main_area, text="Locus: |f(z)| = k or |f(z)| = |g(z)|",
                             font=("Arial", 20, "bold"))
        title.pack(pady=20)

        # Input container
        input_frame = ctk.CTkFrame(self.main_area, fg_color="transparent")
        input_frame.pack(pady=10)

        # Left side
        ctk.CTkLabel(input_frame, text="|").pack(side="left")
        self.entry_geo_g = ctk.CTkEntry(input_frame, width=150, placeholder_text="z - 2*I")
        self.entry_geo_g.pack(side="left", padx=5)
        ctk.CTkLabel(input_frame, text="|  =  |").pack(side="left")

        # Right side (number or expression with z)
        self.entry_geo_d = ctk.CTkEntry(input_frame, width=150, placeholder_text="2  or  z + 1")
        self.entry_geo_d.pack(side="left", padx=5)
        ctk.CTkLabel(input_frame, text="|").pack(side="left")

        # Button
        btn = ctk.CTkButton(self.main_area, text="IDENTIFY LOCUS",
                            fg_color=COLOR_ACCENT, command=self.lancer_analyse_geo)
        btn.pack(pady=20)

        # Results
        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def lancer_analyse_geo(self):
        for widget in self.result_container.winfo_children():
            widget.destroy()

        try:
            left = self.entry_geo_g.get()
            right = self.entry_geo_d.get()

            if not left or not right: return

            # Call engine
            eq_cartesian = moteur_lieu_geometrique(left, right)

            # Show obtained equation
            self.afficher_resultat_math(f"Cartesian~Equation :")
            self.afficher_resultat_math(f"{sp.latex(eq_cartesian)} = 0")

            # --- SMART ANALYSIS ---
            x, y = sp.symbols('x y', real=True)
            # Convert to polynomial to read coefficients of x^2 and y^2
            P = sp.Poly(eq_cartesian, x, y)

            # Get coefficient in front of x^2
            coeff_x2 = P.coeff_monomial(x**2)
            # Get coefficient in front of y^2
            coeff_y2 = P.coeff_monomial(y**2)

            # Diagnostic
            nature = "UNKNOWN"

            # If both coefficients are zero (no squares), it's a line
            if coeff_x2 == 0 and coeff_y2 == 0:
                nature = "LINE (Perpendicular bisector)"
                lbl = ctk.CTkLabel(self.result_container, text=f"Locus : {nature}", font=("Arial", 16, "bold"),
                                   text_color="#00ff00")
                lbl.pack(pady=10)
            # If x^2 and y^2 have same coefficient (e.g. 1x^2 + 1y^2) it's a circle
            elif coeff_x2 == coeff_y2 and coeff_x2 != 0:
                nature = "CIRCLE"

                # Normalize by coeff_x2 to get x^2 + y^2 + ...
                eq_norm = eq_cartesian / coeff_x2

                # Extract coefficients a, b, c
                a = eq_norm.coeff(x, 1)
                b = eq_norm.coeff(y, 1)
                # Constant term
                c = eq_norm.as_coefficients_dict()[1]

                # Compute center and radius
                center_x = -a / 2
                center_y = -b / 2
                r_sq = (a / 2) ** 2 + (b / 2) ** 2 - c

                if r_sq >= 0:
                    radius = sp.sqrt(r_sq)
                    lbl = ctk.CTkLabel(self.result_container,
                                       text=f"CIRCLE with Center ({center_x}, {center_y}) and Radius {sp.simplify(radius)}",
                                       font=("Arial", 14, "bold"), text_color="#00ff00")
                    lbl.pack(pady=10)

                    # Button to draw
                    btn_plot = ctk.CTkButton(self.result_container, text="VIEW GRAPH",
                                             command=lambda: self.dessiner_locus("CIRCLE", center_x, center_y, float(radius)))
                    btn_plot.pack(pady=5)
                else:
                    details = "Imaginary circle (impossible radius)"

            # Bonus: If x^2 and y^2 differ (e.g. 2x^2 + 3y^2), it's an ellipse
            elif coeff_x2 != coeff_y2 and coeff_x2 != 0 and coeff_y2 != 0:
                nature = "ELLIPSE (Out of program?)"

            lbl_res = ctk.CTkLabel(self.result_container, text=f"Identified form : {nature}",
                                   font=("Arial", 16, "bold"), text_color="#00ff00")
            lbl_res.pack(pady=10)
        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

    def dessiner_locus(self, nature, cx=0, cy=0, r=1):
        # Create a new window
        plot_window = ctk.CTkToplevel(self)
        plot_window.title("Geometric Visualization")
        plot_window.geometry("500x500")

        fig, ax = plt.subplots(figsize=(5, 5))

        if nature == "CIRCLE":
            circle_visu = plt.Circle((cx, cy), r, color=COLOR_ACCENT, fill=False, lw=2)
            ax.add_patch(circle_visu)
            # Adjust limits to see entire circle
            limit = float(max(abs(cx), abs(cy)) + r + 1)
            ax.set_xlim(-limit, limit)
            ax.set_ylim(-limit, limit)

        ax.axhline(0, color='white', lw=1)  # X axis
        ax.axvline(0, color='white', lw=1)  # Y axis
        ax.set_aspect('equal')  # Keep circle round
        ax.grid(True, linestyle='--', alpha=0.3)
        fig.patch.set_facecolor(COLOR_BG)
        ax.set_facecolor(COLOR_BG)

        canvas = FigureCanvasTkAgg(fig, master=plot_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    def setup_racines_ui(self):
        self.nettoyer_zone_principale()
        title = ctk.CTkLabel(self.main_area, text="n-th Root Calculation", font=("Arial", 22, "bold"))
        title.pack(pady=20)

        self.entry_z_root = ctk.CTkEntry(self.main_area, placeholder_text="Number z (ex: 1 + I)", width=300)
        self.entry_z_root.pack(pady=10)

        self.entry_n_root = ctk.CTkEntry(self.main_area, placeholder_text="Value of n (ex: 3)", width=150)
        self.entry_n_root.pack(pady=10)

        btn = ctk.CTkButton(self.main_area, text="FIND ROOTS", fg_color=COLOR_ACCENT,
                            command=self.lancer_racines)
        btn.pack(pady=20)

        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def lancer_racines(self):
        for widget in self.result_container.winfo_children():
            widget.destroy()
        try:
            z_str = self.entry_z_root.get()
            n_str = self.entry_n_root.get()
            z, roots = moteur_racines_niemes(z_str, n_str)

            self.afficher_resultat_math(f"Roots~{n_str}{{-th~of~}} {sp.latex(z)} :")

            for i, r in enumerate(roots):
                self.afficher_resultat_math(f"z_{{{i}}} = {sp.latex(r)}")
        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

    def setup_conversion_ui(self):
        self.nettoyer_zone_principale()
        title = ctk.CTkLabel(self.main_area, text="Exponential ➔ Algebraic Form", font=("Arial", 22, "bold"))
        title.pack(pady=20)

        self.entry_expo = ctk.CTkEntry(self.main_area, placeholder_text="Ex: 2 * exp(i * pi/3)", width=350)
        self.entry_expo.pack(pady=10)

        btn = ctk.CTkButton(self.main_area, text="CONVERT", fg_color=COLOR_ACCENT, command=self.lancer_conversion)
        btn.pack(pady=20)

        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def lancer_conversion(self):
        for widget in self.result_container.winfo_children():
            widget.destroy()
        try:
            entry = self.entry_expo.get()
            z, res = moteur_expo_vers_algebre(entry)
            self.afficher_resultat_math(f"Z_{{expo}} = {sp.latex(z)}")
            self.afficher_resultat_math(f"Z_{{alg}} = {sp.latex(res)}")
        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

    def setup_transform_ui(self):
        self.nettoyer_zone_principale()
        title = ctk.CTkLabel(self.main_area, text="Plane Transformations", font=("Arial", 22, "bold"))
        title.pack(pady=20)

        self.trans_type = ctk.CTkComboBox(self.main_area, values=["Translation", "Rotation", "Homothétie"], width=200)
        self.trans_type.pack(pady=10)

        self.entry_z_trans = ctk.CTkEntry(self.main_area, placeholder_text="Affix z (starting point)", width=300)
        self.entry_z_trans.pack(pady=5)

        self.entry_param = ctk.CTkEntry(self.main_area, placeholder_text="Vector b / Angle θ / Ratio k", width=300)
        self.entry_param.pack(pady=5)

        self.entry_center = ctk.CTkEntry(self.main_area, placeholder_text="Center Ω (if Rotation/Homothety)",
                                         width=300)
        self.entry_center.pack(pady=5)

        btn = ctk.CTkButton(self.main_area, text="CALCULATE IMAGE z'", fg_color=COLOR_ACCENT,
                            command=self.lancer_transform)
        btn.pack(pady=20)

        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def lancer_transform(self):
        for widget in self.result_container.winfo_children(): widget.destroy()
        try:
            res = moteur_transformation(self.trans_type.get(), self.entry_z_trans.get(), self.entry_param.get(),
                                        self.entry_center.get() or "0")
            self.afficher_resultat_math(f"The image~of~the~point~is~z' = {sp.latex(res)}")
        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

    def setup_suites_ui(self):
        self.nettoyer_zone_principale()
        title = ctk.CTkLabel(self.main_area, text="Complex Geometric Sequences", font=("Arial", 22, "bold"))
        title.pack(pady=20)

        # Entries
        self.entry_u0 = ctk.CTkEntry(self.main_area, placeholder_text="First term u0 (ex: 1+I)", width=300)
        self.entry_u0.pack(pady=5)

        self.entry_q = ctk.CTkEntry(self.main_area, placeholder_text="Ratio q (ex: I/2)", width=300)
        self.entry_q.pack(pady=5)

        self.entry_n_suite = ctk.CTkEntry(self.main_area, placeholder_text="Index n (ex: 5)", width=100)
        self.entry_n_suite.pack(pady=5)

        btn = ctk.CTkButton(self.main_area, text="CALCULATE u_n AND S_n", fg_color=COLOR_ACCENT,
                            command=self.lancer_suites)
        btn.pack(pady=20)

        self.result_container = ctk.CTkFrame(self.main_area, fg_color="transparent")
        self.result_container.pack(fill="both", expand=True)

    def lancer_suites(self):
        for widget in self.result_container.winfo_children(): widget.destroy()
        try:
            un, sn = moteur_suites_complexes(self.entry_u0.get(), self.entry_q.get(), self.entry_n_suite.get())
            self.afficher_resultat_math(f"Term~u_{{n}} = {sp.latex(un)}")
            self.afficher_resultat_math(f"Sum~S_{{n}} = {sp.latex(sn)}")
        except Exception as e:
            ctk.CTkLabel(self.result_container, text=f"Error : {e}", text_color="red").pack()

if __name__ == "__main__":
    app = CyberMathApp()
    app.mainloop()