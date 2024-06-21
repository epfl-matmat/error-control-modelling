### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ c8805bcd-8bf0-43dd-b5d5-99668028a438
begin
	
toc = 
# sidebar --- DO NOT TOUCH THIS LINE
Markdown.parse(read("sidebar.md", String))
# sidebar --- DO NOT TOUCH THIS LINE 
	
using HypertextLiteral
	
end

# ╔═╡ 7c7e8792-c9c7-11ee-3f11-a1681da40276
# latex macros --- DO NOT TOUCH THIS LINE
include("latex_macros.jl") 
# latex macros --- DO NOT TOUCH THIS LINE

# ╔═╡ 48662bbd-d763-4e6b-957b-f4d6a783afc5
md"# Nomenclature"

# ╔═╡ 506de2af-faea-42e4-8f14-166a9ad600ed
md"""
Latin   | alphabet
--- | :---
$A$ 			| Generic matrix 
$\opA$ 			| Generic operator
$\bloch$ 	    | Bloch-Floquet transform
$\mathbb B$ 	| Plane wave basis
$\mathscr B$ 	| Set of all bounded operators
$\contour$ 		| Contour in the complex plane
$D^\alpha$ 		| Weak derivative
$D(\opA)$ 		| Domain of $\opA$
$\eigenspace _A(\lambda)$ | Eigenspace of $A$ associated with eigenvalue $\lambda$.
$G(\opA)$ 		| Graph of $\opA$
$H$    			| Sobolev space (see function spaces below)
$\opH$ 			| Schrödinger operator / hamiltonian ($- \laplacian / 2 + V$)
$\opH_k$ 		| Bloch fiber
$\hilbert$ 		| Hilbert space
$I$ 			| Identity matrix
$\mathbb K$ 	| $k$-grid or $k$-point mesh
$\mathbb L$ 	| Lattice
$\mathbb L^*$   | Reciprocal lattice
$q_A(u)$ 		| Quadratic form ($\langle u, Au \rangle$)
$Q(\opA)$ 		| Form domain of $\opA$
$a_A(u,v)$ 		| Sesquilinear form ($\langle u, Av \rangle$)
$R_A(u)$ 		| Rayleigh quotient ($\langle u , Au \rangle / \langle u, u \rangle$)
$R_z(A)$ 		| Resolvent ($(A-z I)^{-1}$)
$\mathcal T_R$ 	| Translation operator by $R$
"""

# ╔═╡ fc3cda3f-8ba8-4911-a8a4-ec46cf1b4dfb
md"""
Here is a non-comprehensive list of the notation used in the couse, with an emphasis on course-specific concepts, possible sources of confusion, and notation differing from typical physics / engineering notation.
"""

# ╔═╡ 7f43fd45-f6a5-449b-8401-627dcb047c80
md"""
Greek | alphabet
--- | :---
$\Delta$ 			| Laplace operator ($\text{div grad}$)
$\resolvent$ 		| Resovlent set
$\spectralradius$ 	| Spectral radius
$\sigma$ 			| Spectrum
$\Sigma$ 		| Bottom of the essential spectrum
"""

# ╔═╡ 1bb4d24e-f93c-4f9a-aaa3-1f799b5e2953
md"""
Other | - 
--- | :---
$\bullet ^*$ 				| Adjoint (operators)
$\tilde \bullet$ 			| Approximation of $\bullet$
$\overline{\bullet}$ 		| Complex conjugate (scalars), closure (sets)
$\dot \cup$ 				| Disjoint union
$\varnothing$ 				| Empty set
$\indicator_\Omega$  		| Indicator function over set $\Omega$
$\langle \bullet,\bullet \rangle$ | Inner product
$(\bullet,\bullet)$ 		| Open interval
$[\bullet, \bullet ]$ 		| Closed interval
$\leq$ 						| Vector subspace (sets), less or equal to (scalars)
$\to$ 						| Strong convergence
$\rightharpoonup$ 			| Weak convergence
"""

# ╔═╡ a1fc69f5-ce2e-4910-af4c-9db644a4dad5
md"""
Function (and other) spaces

Space | Definition | Hilbert space ?
--- | :--- | ---
$V_0$ 			| Space $V$ restricted to functions with compact support.
$C^0(\Omega,Y)$ |  $\{ f : \Omega \to Y \mid f \text{ is continuous} \}$
$C^1(\Omega,Y)$ |  $\{ f : \Omega \to Y \mid f' \text{ is continuous} \}$
$C^k(\Omega,Y)$ |  $\{ f : \Omega \to Y \mid f \text{ is infinitely differentiable} \}$
$C^\infty(\Omega,Y)$ |  $\{ f : \Omega \to Y \mid f^{(k)} \text{ is continuous} \}$
$L^2(\Omega)$ 	| $\{f : \Omega \to \mathbb C \mid \int_{\Omega} \vert f (x) \vert ^2 dx < \infty \}$ | ✓
$L^p(\Omega)$ 	| $\{f : \Omega \to \mathbb C \mid \int_{\Omega} \vert f (x) \vert ^p dx < \infty \}$ 
$L^p_{loc}(\Omega)$ 	| $\left \{  f : \Omega \rightarrow \mathbb{C} \ \middle \vert  \   f\rvert_K \in L^{p}(K) \quad \forall K \in \Omega, K \text { compact} \right \}$
$L^2_{per} (\Omega)$ | $\{ f \in L^2_{loc} (\mathbb R^3) \vert f \text{ is } \mathbb L \text{ periodic and } \mathbb L \text{ has unit cell } \Omega \}$ | ✓
$L^2_{qp} (\Omega^*, L^2_{per}(\Omega))$ | $\{ \mathbb R^d \ni k \mapsto u_k \in L^2_{per}(\Omega)  \vert \int_{\Omega^*} \| u_k \|^2_{L^2_{per}(Ω)} \ dk < ∞ \ \text{ and } u_{k+G} = u_k e^{-i G ⋅ x} \}$ | ✓
$\mathscr L(V)$ | $\{ f : V \to V \mid f \text{ linear} \}$
$\ell^2(\mathbb C)$ | $\{z : \mathbb N \to \mathbb C  \mid \ \sum_{n = 0}^\infty \vert z_i \vert^2 < \infty \}$ | ✓
$\ell^p(\mathbb C)$ | $\{z : \mathbb N \to \mathbb C  \mid \ \sum_{n = 0}^\infty \vert z_i \vert^p < \infty \}$
$H^n(\Omega)$ | $\{ f \in L^2(\Omega) \mid D^\alpha f \in L^2(\Omega) \ \forall \alpha \text{ s.t. } \Vert \alpha \Vert _1 \leq n \}$  | ✓
$H^S_{per} (\Omega)$ | $\{ f \in L^2_{per} (\Omega)  \vert \sum_{G \in \mathbb L^*} (1 + \vert G \vert^2)^S  \vert \hat f_G \vert^2 < \infty \}$ | ✓

Note that $\Omega$ here is used in most cases to denote the set on which the function is defined. 
However, in the case of periodic function spaces ($L^2_{per}(\Omega), H^S_{per}(\Omega)$), it denotes the unit cell.
"""

# ╔═╡ 14976bf5-5f50-4920-b8d1-632c440be175
md"""
Inner product of Hilbert spaces (for their definitions, see the table above).
To obtain the associated norm, recall $\| f \| = \sqrt{\langle f,f \rangle}$.

Space  | Inner Product $\langle f,g \rangle$ 
--- |  :--- | 
$L^2(\Omega)$ 	| $\int_\Omega \overline{f(x)} g(x) \ dx$
$\ell^2(\mathbb R)$ | $\sum_{i=0}^\infty \overline{f_i} g_i$ 
$H^n(\Omega)$ | $\sum_{\Vert \alpha \Vert _1 \leq n} \langle D^\alpha f, D^\alpha g \rangle_{L^2}$ 
$L^2_{per} (\Omega)$ | $\int_{\Omega} \overline{f(x)} g(x) dx$
$L^2_{qp} (\Omega^*, H^1_{per}(\Omega))$ | $\frac1{\vert \Omega^*\vert} \int_{\Omega^*} \langle f_k, g_k \rangle_{L^2_{per}(\Omega)} dk$
"""

# ╔═╡ 7c87a6e1-dfe6-4bc1-8857-3b4a55aa6c0b
begin

Sidebar(elts...; location="upper right") = @htl("""
<aside class="sidebar" style='top: 75px;right: 17px;'>$elts</aside>

<style>
aside.sidebar {
	position: fixed;
	max-width: min(30%, 300px, calc(100vw - 750px));
	padding: 0.4rem;
	border-radius: 10px;
	max-height: calc(100vh - 140px);
	overflow: auto;
	z-index: 10;
	background-color: rgba(0, 0, 0, 0.02);
}

aside.aside-sticky table {
	margin: 0.2rem 0;
}
</style>
""")

Sidebar(toc) 

end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"

[compat]
HypertextLiteral = "~0.9.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "5b37abdf7398dc5da4cd347d0609990238d895bb"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"
"""

# ╔═╡ Cell order:
# ╟─7c7e8792-c9c7-11ee-3f11-a1681da40276
# ╟─48662bbd-d763-4e6b-957b-f4d6a783afc5
# ╟─506de2af-faea-42e4-8f14-166a9ad600ed
# ╟─fc3cda3f-8ba8-4911-a8a4-ec46cf1b4dfb
# ╟─7f43fd45-f6a5-449b-8401-627dcb047c80
# ╟─1bb4d24e-f93c-4f9a-aaa3-1f799b5e2953
# ╟─a1fc69f5-ce2e-4910-af4c-9db644a4dad5
# ╟─14976bf5-5f50-4920-b8d1-632c440be175
# ╟─c8805bcd-8bf0-43dd-b5d5-99668028a438
# ╟─7c87a6e1-dfe6-4bc1-8857-3b4a55aa6c0b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
