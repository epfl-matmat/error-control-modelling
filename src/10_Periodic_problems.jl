### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ a0407e4d-3edc-48ea-a335-b815b796e378
begin
	import TikzPictures.TikzPicture
	using LaTeXStrings
	using PlutoUI
	using PlutoTeachingTools
	using HypertextLiteral
	using LinearAlgebra
	using DFTK
	using Plots

	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ e747d030-7487-4598-99d6-93afefd425e1
md"""
# Periodic Problems
"""

# ╔═╡ c7351c0e-c5d5-43e8-b0db-24c7b1985ae1
TODO(md"""Needs more visualisation and rationalisation.
  - Better use the route that Eric uses about introducing $k$-points as a way to exploit the periodicity in supercell calculations
""")

# ╔═╡ 1441bdfc-715e-46da-a4d2-8b788c71853f
md"""
In the previous lecture we saw that the discrete eigenvalues below the
essential spectrum of an operator can be approximated. 
Indeed, if $(\lambda_k, \phi_k)$ denotes the $k$-th such eigenpair (below $\Sigma(\mathcal A)$), then we can choose a *discrete basis* $\{ \phi_i \}_i$ spanning a subspace $S \leq D (\mathcal A)$ and
find an approximation $(\tilde \lambda_k , \tilde \phi_k)$ as the $k$-th
eigenpair of the matrix 
```math
\begin{align}
    A_{ij} \equiv \langle \phi_i | \mathcal A \phi_j \rangle.
\end{align}
```

Due to Courant-Fisher the approximation is from above and due to
Kato-Temple the error is controlled.

As an example for understanding this discretization procedure we will
discuss this in chapter the special setting of periodic problems with
Hamiltonian $\opH = - \frac1{2} \laplacian + V$, where $V$ is periodic.
Albeit such Hamiltonians have *no eigenvalues*, the Bloch-Floquet
transform will still allow us to construct a matching basis and
approximate the spectra of these operators in a meaningful way.
"""

# ╔═╡ 7971ffe3-39f0-4cf8-9cee-fd55483a8ed9
md"""
## Periodic functions and lattices

A periodic problem is characterized by being invariant to certain translations. 

- For example, the sine function satisfies 
  ```math
  \begin{align}
      \sin(x) = \sin(x + 2 \pi m) && \forall m \in \mathbb Z
  \end{align}
  ```
  making it $2 \pi$-periodic. Alternatively, we could denote this as 
  ```math
  \begin{align}
   \sin(x) = \sin(x+R) && \forall R \in \mathbb L
  \end{align}
  ```
  where $\mathbb L$ is the *lattice*, $\mathbb L = \{ 2 \pi m \mid m \in \mathbb Z \} = 2 \pi \mathbb Z$.


- The importance of periodicity for numerical computations is as follows : once an $\mathbb L$-periodic function, like $\sin$, is known inside the *unit cell* $\Omega$, the open and bounded set 
  ```math
  \begin{align}
      \Omega = \{ x \in \mathbb R \mid |x-R| > |x| \quad \forall R \in \mathbb L \setminus \{ 0 \} \} 
  \end{align}
  ```
  it is actually known for everywhere in $\mathbb R$. 
  This suggests that solving periodic problems inside the unit cell should be sufficient to know the solution everywhere.
  In the case of $\sin$, the unit cell is $(-\pi,\pi)$.
"""

# ╔═╡ 8d9fe76f-dfb3-4683-83c2-5d7cd2aff712
TikzPicture(L"""
        \draw (-5,0) node[left]{$\mathbb L$} 
        -- (-4,0) node{$\bullet$} node[below]{$-4 \pi$} 
        -- (-2,0) node{$\bullet$} node[below]{$-2 \pi$} 
        -- (-1,0)  node[above]{$- \pi$} 
        -- (0,0)  node{$\bullet$} node[below]{$0$} 
        -- (1,0)  node[above]{$\pi$} 
        -- (2,0)  node{$\bullet$} node[below]{$2 \pi$} 
        -- (4,0)  node{$\bullet$} node[below]{$4 \pi$} 
        -- (5,0)  
        ;

		\draw (-6.5,0) ;

        \draw (-1,0.5) node{$ ]$} -- (0,0.5) node[above]{$\Omega$} --  (1,0.5) node{$[$} ;
""",width="23cm",options="scale=1",preamble=raw"\usepackage{amsfonts}")

# ╔═╡ df23c23e-b509-4f48-9ab9-dd2b59834d01
md"""
- These concepts naturally generalize to higher dimensions, e.g. a 3D lattice is 
  ```math
  \begin{align}
      \mathbb L = \{ \mathbb Z a_1 + \mathbb Z a_2 + \mathbb Z a_3 \} 
  \end{align}
  ```
  where $a_1,a_2,a_3 \in \mathbb R ^3$ are three linearly independent vectors. They define a unit cell 
  ```math
  \begin{align}
    \Omega = \{ x \in \mathbb R^3 \mid |x-R| > |x| \quad \forall R \in \mathbb L \setminus \{ 0 \} \}.
  \end{align}
  ```

  Frequently, the matrix 
  ```math
  \begin{align}
    A = \begin{pmatrix}
        a_1 & a_2 & a_3\\
        \downarrow & \downarrow & \downarrow
        \end{pmatrix}
        \in \mathbb R ^{3 \times 3}    
  \end{align}
  ```
  is considered.

"""

# ╔═╡ 9b17dde4-6652-4800-a568-1244cf519593
md"""
- A function $f$ is now called $\mathbb L$-periodic if 
  ```math
  \begin{align}
     f(x) = f(x + R) && \forall R \in \mathbb L
  \end{align}
  ```
  and almost every $x \in \mathbb R^3$.


- Based on this observation we introduce the translation operator
  ```math
  \begin{align}
      \mathcal T_R : C_0^\infty (\mathbb R^3) \to C_0^\infty (\mathbb R^3), f \mapsto f(\bullet - R) 
  \end{align}
  ```
  which by denseness of $C_0^\infty$ on $L^2$ can be extended to all $L^2$ functions.
"""

# ╔═╡ fbbe239e-c06e-41b1-b1e3-857585f994bd
md"""
- Further we introduce the function space 
  ```math
  \begin{align}
    L^2_\text{per} (\Omega) = \{ f \in L^2_\text{loc} (\mathbb R^3) \mid f \text{ is } \mathbb L \text{ periodic} \}
  \end{align}
  ```
  with inner product 
  ```math
  \begin{align}
    \langle f, g \rangle _{L_\text{per}^2 (\Omega)} = \int_{\Omega} \overline{f(x)} g(x) dx.
  \end{align}
  ```


- A convenient basis for $L^2_\text{per} (\Omega)$ is the plane wave or Fourier basis 
  ```math
  \begin{align}
      \mathbb B = \left \{ e_G(x) \equiv \frac1{\sqrt{|\Omega|}} e^{i G x} \middle \vert G \in \mathbb L ^* \right \} 
  \end{align}
  ```
  where $\mathbb L^*$ denotes the *reciprocal lattice*

  ```math
  \begin{align}
      \mathbb L^* = \{ \mathbb Z b_1 + \mathbb Z b_2 + \mathbb Z b_3 \}
  \end{align}
  ```
  such that 
  ```math
  \begin{align}
      a_i \cdot b_j = 2 \pi \delta_{ij} && \text{for } i,j \in \{1,2,3\}
  \end{align}
  ```
  where $\delta_{ij}$ is the Kronecker delta, i.e. 
  ```math
  	B = 2 \pi A^{-1} 
  ```
  with $B = (b_1,b_2,b_3)$.
  We also note, that often we write with abuse of notation $\Omega\equiv |\Omega|$ (i.e. we think of $\Omega$ sometimes as the volume of the unit cell and not the unit cell itself).

"""

# ╔═╡ ea689cad-2859-41b0-aef3-36651a8d166b
md"""
- Note that the construction of $\mathbb L^*$ is such to ensure that $e_G(x)$ is $\mathbb L$-periodic by enforcing 
  ```math
  \begin{align}
      e^{i R \cdot G} = 1 && \forall \ R \in \mathbb L, G \in \mathbb L^*
  \end{align}
  ```


- We also define $\Omega^*$ as the unit cell of $\mathbb L^*$, also called the *first Brillouin zone*.
"""

# ╔═╡ 6ba1ce9a-6ced-4a82-9b5a-6b07a25c534b
md"""
### Results from Fourier analysis

We revise some concepts from Fourier analysis.

- For a function $f \in L^2_\text{per} (\Omega)$ the following expressions are frequently encountered :
  ```math
  \begin{align}
      f = \sum_{G \in \mathbb{L}^*} \hat f_G e_G \tag{Fourier Series Expansion}
  \end{align}
  ```
  with 
  ```math
  \begin{align}
      \hat f_G \equiv \langle e_G, f \rangle = \int_{\Omega} f(x) \overline{e_G(x)} dx    \tag{Fourier Coefficients}
  \end{align}
  ```
  and 
```math
\begin{align}
      \| f \|_{L^2_\text{per} (\Omega)} = \int_\Omega | f |^2 dx = \sum_{G \in \mathbb L^*} |\hat f_G|^2 < \infty. \tag{Parseval's Identity}
\end{align}
```

"""

# ╔═╡ 6c5c1231-9e16-4a43-8dec-68df585d78b9
md"""
- Employing a Fourier basis, **periodic Sobolev spaces** can be easily characterized as 
  ```math
  \begin{align}
      H^s_\text{per} (\Omega) = \left \{ f \in L^2_\text{per} (\Omega)  \middle | \sum_{G \in   \mathbb L^*} (1 + |G|^2)^s 
    |\hat f_G|^2 < \infty \right  \} 
  \end{align}
  ```
  where $s>0$.
"""

# ╔═╡ 3563cb4b-1e0d-47fd-8546-7a3885d397d1
md"""
- The **regularity of a function** $u$ is linked to the **decay rate of the Fourier coefficients**. To see this assume that $u : \Omega \to \mathbb{R}$ is $\mathbb{L}$-periodic and $m$ times differentiable with $\frac{d^m}{d x^m} u \in L^2_\text{per}(\Omega)$.

  Using that $\widehat{\frac{d}{dx} u}_G = i G \, \hat{u}_G$ we first obtain the two Fourier series
  ```math
  \begin{aligned}
  u(x) &= \sum_{G\in \mathbb{L}^\ast} \hat{u}_G e_G(x)\\
  \frac{d^m u}{dx^m}(x) &= \sum_{G\in \mathbb{L}^\ast} (iG)^m \hat{u}_G e_G(x)
  \end{aligned}
  ```
  By Parseval's theorem
  ```math
  \left\|\frac{d^m u}{dx^m}\right\|_{L^2_\text{per}}^2 = \sum_{G\in \mathbb{L}^\ast} |G|^{2m} \, |\hat{u}_G|^2 \geq \underbrace{|\widetilde{G}|^{2m} \, |\hat{u}_\widetilde{G}|^2}_\text{for one $\widetilde{G} \in \mathbb{L}^\ast$}
  ```
  Therefore $\forall G \in \mathbb{L}^\ast$ we have
  ```math
  \left\|\frac{d^m u}{dx^m}\right\|_{L^2_\text{per}} \geq |G|^m \,  |\hat{u}_G|
  \quad \Rightarrow \quad |\hat{u}_G| \leq \frac{\left\|\frac{d^m u}{dx^m}\right\|_{L^2_\text{per}}}{|G|^m}
  ```
  which implies that the decay of $\hat{u}_G$ with increasing $|G|$ is at least faster than a polynomial in $|G|$ of degree $m$.
- In particular if $u$ is **analytic** (i.e. all $m$ are possible), than the decay of $\hat{u}_G$ is faster than any polynomial in $|G|$.

More formally this is summarised by the following proposition:
"""

# ╔═╡ b90fe635-247e-4842-acbc-2a07204c89c8
#=
Theorem 1 (Fourier transforms and derivatives).
If $f \in L^1(\mathbb{R}^n)$ and $\frac{\partial f}{\partial x_j} \in L^1(\mathbb{R}^n)$, then
```math
\widehat{\frac{\partial f}{\partial x_j} (\xi)} = i \xi_j \hat{f}(\xi).
```
=#

md"""
!!! info "Proposition 0: Fourier decay and regularity" 
    If $f \in C^p(\mathbb{R}^n)$ and $\partial^α f \in L^1(\mathbb{R}^n)$ for $|α| \leq p$ then there is a $C>0$ such that for all $G \in \mathbb{R}^n$,
    ```math
    \left| \hat{f}(G) \right| \leq C (1 + |G|)^{-p}
    ```
    Reciprocally if $x^α \in L^1(\mathbb{R}^n)$ for all multi-indices $|α| \leq p$ then $\hat{f}$ is iof class $C^p$.

    Or colloquially: Regularity in real space implies decay in Fourier and vice versa.
"""

# ╔═╡ 20029c12-e608-4235-96c6-19fab70f02f3
md"""
### Example: Computing eigenvalues using plane-wave discretisations

For this example we consider the operator
```math
   \mathcal{H} = - \frac12 Δ + V(x) \quad \text{with $V(x) = \cos(x)$}
```
on $L^2_\text{per}(\Omega)$ with unit cell $\Omega = (-π, π)$ and periodic boundary conditions.
The lattice is $\mathbb{L} = 2π \mathbb{Z}$ and the lattice constant $a = 2π$.
We remark that this example is also discussed in the DFTK documentation, see [Comparing-discretization-techniques](https://docs.dftk.org/stable/guide/discretisation/).

The $\cos(x)$ takes the role of a lattice-periodic potential.
As we will see [in this later discussion](#Fibers-of-periodic-Schrödinger-operators) that sufficiently regular potentials (such as $\cos$) make this operator $\mathcal{H}$ self-adjoint on $L^2_\text{per}(\Omega)$ with domain $D(\opH) = H^2_\text{per}(\Omega)$ and form domain $Q(\opH) = H^1_\text{per}(\Omega)$. 
Moreover it actually has a compact resolvent, i.e. its **spectrum consists entirely of eigenvalues**.

Of note **the boundary condition matters**:
The spectrum we will get is different from e.g. the spectrum of $H$ on $L^2(\mathbb{R})$ as we will see later.
"""

# ╔═╡ 9cf1f790-835a-438d-be09-2b4eea978351
md"""
Above we introduced the Fourier basis
```math
\begin{align}
  \mathbb B = \left \{ e_G(x) \equiv \frac1{\sqrt{|\Omega|}} e^{i G x} \middle \vert G \in \mathbb L ^* \right \}
\end{align}
```
for $L^2_\text{per} (\Omega)$. In our 1D setting $\mathbb{L}^\ast = \frac{a}{2π} \mathbb{Z} = \mathbb{Z}$.

A natural way to form a **discretisation basis** for $L^2_\text{per} (\Omega)$ is to truncate this basis at a finite value for $|G|$.
For agreement with the physics literature one usually employs the condition that $\frac12 |G|^2 \leq E_\text{cut}$ where $E_\text{cut}$ is called the kinetic energy cutoff.

For our 1D problem an appropriate finite-sized plane-wave basis is thus
```math
\mathbb{B}^{E_\text{cut}} = \left\{ e_G \, \middle| \, G \in \mathbb{Z} \ \text{ and } \frac12 |G|^2 \leq E_\text{cut} \right\}
= \left\{ e_G \, \middle| \, G \in \mathbb{Z} \ \text{ and } |G| \leq \sqrt{2E_\text{cut}} \right\}
```
"""

# ╔═╡ 7c6dcdc6-692d-4a47-85ac-b9baee3a1be5
md"""
By construction the span of this basis provides a subspace to $L^2_\text{per} (\Omega)$. But notably since
```math
\left\| \sum_{|G| \leq \sqrt{2E_\text{cut}}} \hat{u}_G e_G \right\|^2_{H^1_\text{per}}
= \sum_{|G| \leq \sqrt{2E_\text{cut}}} (1 + |G|^2) |\hat{u}_G|^2
< \infty
```
every element of the span of $\mathbb{B}^{E_\text{cut}}$ also lives in $H^1_\text{per}$, such that $\mathbb{B}^{E_\text{cut}}$ also forms the basis for a proper subspace of the form domain $Q(\mathcal{H})$, i.e. $\text{span}\ \mathbb{B}^{E_\text{cut}} \leq H^1_\text{per} = Q(\mathcal{H})$.
In line with Courant-Fisher we thus consider the sequences
```math 
\tag{$\ast$}
\mu_j^{E_\text{cut}} (\opH) \coloneqq \min_{\substack{W \subset \textcolor{red}{\text{span}\ \mathbb{B}^{E_\text{cut}}} \\ \dim (W) = j}} \max_{0 \neq \varphi \in W}
\frac{\langle \phi, \mathcal{H} \phi \rangle}{\langle \phi, \phi\rangle}
```
Since $\mathcal{H}$ has only eigenvalues, we expect that
```math
\lim_{E_\text{cut} \to \infty} \mu_j^{E_\text{cut}} (\opH) = \lambda_j (\opH),
```
i.e. the $j$-th eigenvalue of $\opH$.

As discussed in [Operators and their spectra](https://teaching.matmat.org/error-control/09_Operators_Spectra.html), the RHS of ($\ast$) can also be obtained as the
$j$-th eigenvalue of a matrix $H \in \mathbb{C}^{N\times N}$
where $N = |\mathbb{B}^{E_\text{cut}}|$ (size of the basis) and the matrix has the elements
```math
H_{GG'} = \langle e_G, \opH e_{G'}\rangle,
```
since the $e_G$ are orthonormal.
"""

# ╔═╡ 696d9342-d8fd-442f-82ed-81032b6fcba2
md"""To compute these elements we first notice
```math
\langle e_G, e_{G'}\rangle = ∫_0^{2π} e_G^\ast(x) e_{G'}(x) d x = ∫_0^{2π} e^{i(G'-G)x} d x = δ_{G, G'}
```
and with $V(x) = \cos(x)$ further that
```math
\langle e_G, \opH e_{G'}\rangle = \frac 1 2 \left(|G|^2 \delta_{G,G'} + \delta_{G, G'+1} + \delta_{G, G'-1}\right)
```
"""

# ╔═╡ fec13421-c63f-4e96-9d0d-4d4f43ab172d
Foldable("Derivation of ⟨eG, H, eG'⟩", md"""
Next for $V(x) = \cos(x)$ we obtain
```math
\langle e_G, \mathcal{H} e_{G'}\rangle = \frac 1 2 ∫_0^{2π} e_G^\ast(x) H e_{G'}(x) d x
```
We start by applying the Hamiltonian to a plane-wave:
```math
\mathcal{H} e_{G'}(x) = - \frac{1}{2} (-|G|^2) \frac{1}{\sqrt{2π}} e^{iG'x} + \cos(x) \frac{1}{\sqrt{2π}} e{iG'x}
```
Then, using the result of the first part of the exercise and the fact that
$cos(x) = \frac 1 2 \left(e{ix} + e{-ix}\right)$, we get:
```math
\begin{align*}
⟨ e_G, H e_{G'}⟩
&= \frac 1 2 G^2 δ_{G, G'} + \frac 1 {4π} \left(∫_0^{2π} e^{ix ⋅ (G'-G+1)} d x + ∫_0^{2π} e^{ix ⋅ (G'-G-1)} d x \right) \\
&= \frac 1 2 \left(|G|^2 \delta_{G,G'} + \delta_{G, G'+1} + \delta_{G, G'-1}\right)
\end{align*}
```
""")

# ╔═╡ 2386c0dd-0b29-4168-9f25-2187d43b1461
md"""
For given `Ecut` the lowest $5$ eigenvalues are therefore obtained by the function:
"""

# ╔═╡ 3ceea351-5110-4a6d-a525-fdf259c74675
# Plane waves Hamiltonian -½Δ + cos on [0, 2pi].
function build_plane_waves_matrix_cos(Ecut)
	# Determine maximal value for G
	Gmax = floor(Int, sqrt(2Ecut))
    # Plane wave approximation to -½Δ
    Gsq = [float(i)^2 for i in -Gmax:Gmax]
    # Hamiltonian as derived above
    1/2 * Tridiagonal(ones(2Gmax), Gsq, ones(2Gmax))
end

# ╔═╡ 2a765fbb-759b-4480-913c-e8f01ef321a0
build_plane_waves_matrix_cos(10)

# ╔═╡ 5aeaa6b8-50b9-4c5f-bcce-592965b68353
let
	Ecut = 10.0
	H = build_plane_waves_matrix_cos(Ecut)
	L = eigvals(H)
	L[1:5]
end

# ╔═╡ 700ac798-5f2f-456e-8be8-cd40c76e6c96
md"""Plotting these as $\text{Ecut}$ (therefore the basis size) increases we see the expected convergence (by Courant-Fisher)
"""

# ╔═╡ 15fb9074-edbf-4458-9a41-bd3d8707a501
let
	reference = eigvals(build_plane_waves_matrix_cos(200))

	error_gs  = Float64[]
	error_5th = Float64[]

	Ecut_range = 2:2:20
	for Ecut in Ecut_range
		L = eigvals(build_plane_waves_matrix_cos(Ecut))

		push!(error_gs, abs(L[1] - reference[1]))
		push!(error_5th, abs(L[5] - reference[5]))
	end
	
	plot(Ecut_range, error_gs, mark=:x, yaxis=:log, lw=2, label="Ground state")
	plot!(Ecut_range, error_5th, label="5th state", lw=2, mark=:x)
end

# ╔═╡ 32086662-d60d-47cd-acb6-0f1a4f18e351
md"""
## Periodic operators and periodic Schrödinger problems

!!! note "Definition (L-periodicity)"
	A linear operator $\mathcal A : L^2(\mathbb R^d) \to L^2(\mathbb R^d)$ is called
	$\mathbb L$-periodic if it commutes with translations
	$\mathcal T_R$ for all $R \in \mathbb L$, i.e. 
	```math
	\begin{align}
	        \mathcal A (\mathcal T_R (\varphi)) = \mathcal T_R(\mathcal A (\varphi)) &&
	        \forall \ R \in \mathbb L, \varphi \in L^2(\mathbb R^d).
	    
	\end{align}
	```



!!! note "Theorem 1"
	Let $V \in L^{3/2}_\text{per} (\Omega)$ for $\Omega \subset \mathbb R^3$. The
	operator 
	```math
		\opH = - \frac1{2} \laplacian + V
	```
	is self-adjoint on
	$\hilbert = L^2(\mathbb R^3)$
	with $D(\opH) = H^2(\mathbb R^3)$ and $Q(\opH) = H^1(\mathbb R^3)$.
	$\opH$ is also bounded from below.

Notice that Theorem 1 makes a statement about $\opH$ on $L^2(\mathbb R^3)$ and not on $L^2_\text{per}(\mathbb R^3)$ as we considered in the preceeding numerical examples. While the **potential of the operator is periodic** the **domain is not periodic**. As a result, while $\opH$ is $\mathbb L$-periodic,
not all eigenfunctions of $\opH$ are $\mathbb L$-periodic.

While this sounds surprising, this is the setting that is physically most relevant. In fact constraining eigenfunctions to be of equal periodicity as $\opH$ turns out too restrictive.
"""

# ╔═╡ 1c07a977-68b7-44e1-9adf-ad1af1b44ebc
md"""
### High-level example: Approximating periodic operators using Bloch theory and DFTK

Before diving into the mathematical details, let us first discuss how to deal with periodic Schrödinger operators on a high level. This is a reformulation of the [Periodic problems and plane-wave discretisations](https://docs.dftk.org/stable/guide/periodic_problems/) tutorial of the [DFTK documentation](https://docs.dftk.org).

Let's consider the case where we want to find the spectrum of the **free-electron Hamiltonian** on $L^2(\mathbb{R})$, i.e. 
```math
\opH = - \frac12 \Delta.
```
Compared to Theorem 1 the potential $V(x) = 0$ is in $L^{3/2}_\text{per}(\Omega)$ for actually any valid unit cell $\Omega$ which tiles $\mathbb{R}$. For simplicity we choose $\Omega = (0, 2π)$, thus a lattice $\mathbb{L} = 2π \mathbb{Z}$ with lattice constant $a = 2π$.
"""

# ╔═╡ e2b26548-3f7d-4950-b22c-e7106912df49
md"""
Physically, in a free-electron model (which gives rise to this Hamiltonian) electron motion is only
by their own kinetic energy. As this model features no potential which could make one point
in space more preferred than another, we would expect this model to be periodic.
In fact the above Definition easy easy to verfiy in our setting,
i.e. one can easily show exactly that
```math
\mathcal{T}_{ma} \mathcal{H} = \mathcal{H} \mathcal{T}_{ma} \quad  ∀ m ∈ \mathbb{Z}.
```
"""

# ╔═╡ d7644757-dcee-4414-9dda-e09476779ab9
md"""
Without stating the details at this stage, 
**Bloch's theorem** now tells us that for periodic operators,
the solutions to the eigenproblem
```math
    \mathcal{H} ψ_{kn} = ε_{kn} ψ_{kn}
```
satisfy a factorization
```math
    ψ_{kn}(x) = e^{i k⋅x} u_{kn}(x)
```
into a plane wave ``e^{i k⋅x}`` and a lattice-periodic function
```math
   T_{ma} u_{kn}(x) = u_{kn}(x - ma) = u_{kn}(x) \quad ∀ m ∈ \mathbb{Z}.
```
In this ``n`` is a labeling integer index and ``k`` is a real number,
whose details will be clarified in the next section.
The index ``n`` is sometimes also called the **band index** and
functions ``ψ_{kn}`` satisfying this factorization are also known as
**Bloch functions** or **Bloch states**.
"""

# ╔═╡ 0d419731-99e9-4380-a667-d2f935044f17
TODO(md"Rationalise Bloch theory by solving a simple atomic chain model for increasing supercell sizes and then basically manually computing the Bloch decomposition / the mapping between supercells and k-points.")

# ╔═╡ 0b85ae26-7e44-4506-a096-71d2a49def36
md"""
Consider the application of ``2\mathcal{H} = -Δ = - \frac{d^2}{d x^2}``
to such a Bloch wave. First we notice for any function ``f``
```math
   -i∇ \left( e^{i k⋅x} f \right)
   = -i\frac{d}{dx} \left( e^{i k⋅x} f \right)
   = k e^{i k⋅x} f + e^{i k⋅x} (-i∇) f = e^{i k⋅x} (-i∇ + k) f.
```
Using this result twice one shows that applying ``-Δ`` yields
```math
\begin{aligned}
   -\Delta \left(e^{i k⋅x} u_{kn}(x)\right)
   &= -i∇ ⋅ \left[-i∇ \left(u_{kn}(x) e^{i k⋅x} \right) \right] \\
   &= -i∇ ⋅ \left[e^{i k⋅x} (-i∇ + k) u_{kn}(x) \right] \\
   &= e^{i k⋅x} (-i∇ + k)^2 u_{kn}(x) \\
   &= e^{i k⋅x} 2\mathcal{H}_k u_{kn}(x),
\end{aligned}
```
where we defined
```math
    \mathcal{H}_k = \frac12 (-i∇ + k)^2.
```
The action of this operator on a function ``u_{kn}`` is given by
```math
    \mathcal{H}_k u_{kn} = e^{-i k⋅x} \mathcal{H} e^{i k⋅x} u_{kn},
```
which in particular implies that
```math
   \mathcal{H}_k u_{kn} = ε_{kn} u_{kn} \quad ⇔ \quad \mathcal{H} (e^{i k⋅x} u_{kn}) = ε_{kn} (e^{i k⋅x} u_{kn}).
```
To seek the eigenpairs of ``\mathcal{H}`` we may thus equivalently
find the eigenpairs of *all* ``\mathcal{H}_k``.
The point of this is that the eigenfunctions ``u_{kn}`` of ``\mathcal{H}_k``
are periodic (unlike the eigenfunctions ``ψ_{kn}`` of ``\mathcal{H}``).
In contrast to ``ψ_{kn}`` the functions ``u_{kn}`` can thus be fully
represented considering the eigenproblem only on the unit cell.
"""

# ╔═╡ da661541-2007-4681-baea-36164b43efc7
md"""
A detailed mathematical analysis shows that the transformation from ``\mathcal{H}``
to the set of all ``\mathcal{H}_k`` for a suitable set of values for
``k`` (details below)
is actually a unitary transformation, the so-called **Bloch transform**.
This transform brings the Hamiltonian into the symmetry-adapted basis for
translational symmetry, which are exactly the Bloch functions.
Similar to the case of choosing a symmetry-adapted basis for other kinds of symmetries
(like the point group symmetry in molecules), the Bloch transform also makes
the Hamiltonian ``\mathcal{H}`` block-diagonal[^1]:
```math
    \mathcal{H} \stackrel{\bloch}{\longrightarrow} \left( \begin{array}{cccc} \mathcal{H}_1&&&0 \\ &\mathcal{H}_2\\&&\mathcal{H}_3\\0&&&\ddots \end{array} \right)
```
with each block ``\mathcal{H}_k`` taking care of one value ``k``.
This block-diagonal structure under the basis of Bloch functions lets us
completely describe the spectrum of ``\mathcal{H}`` by looking only at the spectrum
of all ``\mathcal{H}_k`` blocks, or written explicitly
```math
  	\sigma(\opH) = \bigcup_{k} \sigma(\opH_k).
```


[^1]: Notice that block-diagonal is a bit an abuse of terms here, since the Hamiltonian
      is not a matrix but an operator and the number of blocks is infinite.
      The mathematically precise term is that the Bloch transform reveals the fibers
      of the Hamiltonian.
"""

# ╔═╡ 0309d400-92ed-46ad-afce-a516a7f1380b
md"""
#### The Brillouin zone

We now consider the parameter ``k`` of the Hamiltonian blocks in detail.

- As discussed ``k`` is a real number. It turns out, however, that some of
  these ``k∈\mathbb{R}`` give rise to operators related by unitary transformations
  (again due to translational symmetry).
- Since such operators have the same eigenspectrum, only one version needs to be considered.
- The smallest subset from which ``k`` is chosen is the **Brillouin zone** (BZ) $\Omega^\ast$.

- The BZ is the unit cell of the **reciprocal lattice**, which may be constructed from
  the **real-space lattice** by a Fourier transform.
- In our simple 1D case the reciprocal lattice is just
  ```
    ... |--------|--------|--------| ...
           2π/a     2π/a     2π/a
  ```
  i.e. like the real-space lattice, but just with a different lattice constant
  ``b = 2π / a = 1``.
- The BZ in our example is thus ``\Omega^\ast = (-π/a, π/a) = (-½, ½)``.
  The members of ``\Omega^\ast`` are typically called **``k``-points**.

"""

# ╔═╡ a790a782-2252-4650-933a-4c002bd3e943
md"""
#### Discretization and plane-wave basis sets

With what we discussed so far the strategy to find all eigenpairs of a periodic
Hamiltonian ``H`` thus reduces to finding the eigenpairs of all ``H_k`` with ``k ∈ \overline{\Omega^\ast}`` (i.e. $\Omega^\ast$ plus its boundary to make a closed set).
This requires *two* discretisations:

  - ``\overline{\Omega^\ast}`` is infinite (and not countable). To discretize we first only pick a finite number
    of ``k``-points, $\mathbb K \subset \overline{\Omega^*}$.
    Usually this **``k``-point sampling** is done by picking ``k``-points
    along a regular grid inside the BZ, the **``k``-grid**.
  - Each ``\mathcal{H}_k`` is still an infinite-dimensional operator.
    Following a standard Ritz-Galerkin ansatz we project the operator into a finite basis
    and diagonalize the resulting matrix.

For the second step multiple types of bases are used in practice (finite differences,
finite elements, Gaussians, ...). Sticking to plane-waves, just as before, one usually employs a **$k$-adapted plane wave basis**
```math
  \mathbb B^{E_\text{cut}}_{k} = \left  \{ e_G \middle \vert G \in \mathbb L^* , \tfrac1{2} |G+k|^2 \leq E_\text{cut} \right \}
```
Note that thus a **different basis set** for each $k \in \mathbb K$ is employed,
which turns out to provide more consistent approximation properties compared to using the same plane-wave basis for all $k$-points.
"""

# ╔═╡ 3de086c4-877d-4d46-825a-376a7555bc7c
md"""
#### Solving the free-electron Hamiltonian numerically

One typical approach to get physical insight into a Hamiltonian ``\mathcal{H}`` is to plot
a so-called **band structure**, that is the eigenvalues of ``\mathcal{H}_k`` versus ``k``.
In DFTK we achieve this using the following steps:
"""

# ╔═╡ c0894d52-b2e1-4015-a56b-1f447099dbda
md"""
**Step 1**: Build the 1D lattice. DFTK is mostly tailored for 3D problems.
Therefore quantities related to the problem space are have a fixed
dimension 3. The convention is that for 1D / 2D problems the
trailing entries are always zero and ignored in the computation.
For the lattice we therefore construct a 3x3 matrix with only one entry.
"""

# ╔═╡ 47d5c6a3-294f-4fd6-be95-72209103db60
begin
	lattice_free = zeros(3, 3)
	lattice_free[1, 1] = 20.
end;

# ╔═╡ c1f717e0-fa07-4cbe-9105-a8745ffd053c
md"""
**Step 2**: Select a model. In this case we choose a free-electron model,
which is the same as saying that there is only a Kinetic term
(and no potential) in the model.
"""

# ╔═╡ 53aed70a-bdd0-4c28-8755-31ae96daf413
model_free = Model(lattice_free; terms=[Kinetic()])

# ╔═╡ f5c11992-a773-4420-9abd-9681e7ca9eed
md"""
**Step 3:** Define a plane-wave basis using this model and a cutoff ``E_\text{cut}``
of 300 Hartree. The ``k``-point grid is given as a regular grid in the BZ
(a so-called **Monkhorst-Pack** grid). Here we select only one ``k``-point (1x1x1),
see the note below for some details on this choice.
"""

# ╔═╡ 4779dfeb-c180-4d98-b006-571ca7adadcb
Foldable("Note on the selection of k-point grid", md"""
You might wonder why we only selected a single ``k``-point (clearly a very crude
and inaccurate approximation). In this example the `kgrid` parameter specified
in the construction of the `PlaneWaveBasis`
is not actually used for plotting the bands. It is only used when solving more
involved models like density-functional theory (DFT) where the Hamiltonian is
non-linear, which we will discuss later.
In these cases before plotting the bands the self-consistent field
equations (SCF) need to be solved first. This is typically done on
a different ``k``-point grid than the grid used for the bands later on.
In our case we don't need this extra step and therefore the `kgrid` value passed
to `PlaneWaveBasis` is actually arbitrary.
""")

# ╔═╡ ed486e03-f003-4267-851c-f8ce9713dfb2
basis_free = PlaneWaveBasis(model_free; Ecut=300, kgrid=(1, 1, 1))

# ╔═╡ 9864b4d7-baef-4798-982d-5a75cb3a2793
md"""
**Step 4**: Plot the bands! Select a density of ``k``-points for the ``k``-grid to use
for the bandstructure calculation, discretize the problem and diagonalize it.
Afterwards plot the bands.
"""

# ╔═╡ 3e5db56a-9b1c-4e15-b5ad-5ddd2f8c8dd8
let
	bands = compute_bands(basis_free; n_bands=6, kline_density=100)
	p = plot_bandstructure(bands)
	ylims!(p, -0.1, 0.5)
end

# ╔═╡ 762968dd-1d23-4ada-8f90-c6a98667efff
md"""
This plot is a plot of the eigenvalues of $\mathcal{H}_k$ in $y$-direction versus the $k$-point value in $x$-direction, going from $0$ all the way to $½$. The $n$-th eigenvalue of each $\mathcal{H}_k$ is connected via a blue line. These blue lines are referred to as **bands**.

We notice that at each selected $k$-point the operator $\mathcal{H}_k$ has a discrete spectrum. Considering the union of all bands, i.e. the eigenspectrum of the original $\mathcal{H}$ it is plausible that
```math
\sigma(\opH) = \bigcup_{k} \sigma(\opH_k) = [0, \infty)
```
as more detailed analysis shows.

This is in agreement with our discussion of the [Laplacian in an unbounded domain](https://teaching.matmat.org/error-control/02_Laplace_error_sources.html), where we found that $-\Delta$ on $L^2(\mathbb{R})$ indeed has a spectrum consisting of all positive real numbers.
"""

# ╔═╡ d2f9823d-0b80-400b-8ddd-bc029a24b987
md"""
#### Adding potentials

So far so good. But free electrons are actually a little boring,
so let's add a potential interacting with the electrons.

The modified problem we will look at consists of diagonalizing
```math
\mathcal{H}_k = \frac12 (-i \nabla + k)^2 + V
```
for all ``k \in \overline{\Omega^\ast}`` with a periodic potential ``V`` interacting with the electrons.

We will use `ElementGaussian`, which is an analytic potential describing a Gaussian
interaction with the electrons to DFTK. A single potential looks like:
"""

# ╔═╡ 0c996cdc-b54d-4167-b9e9-714ee4911025
let
	nucleus = ElementGaussian(0.3, 10.0)
	plot(r -> DFTK.local_potential_real(nucleus, norm(r)), xlims=(-50, 50))
end

# ╔═╡ a3f148b9-911e-4fc7-8de3-2c4ffcf97756
md"""
With this element at hand we can easily construct a setting
where two potentials of this form are located at positions
``20`` and ``80`` inside the lattice ``(0, 100)``:
"""

# ╔═╡ 61fed22d-8bed-4694-a74f-68b025b45d50
begin
	# Define the 1D lattice [0, 100]
	lattice = diagm([100., 0, 0])
	
	# Place them at 20 and 80 in *fractional coordinates*,
	# that is 0.2 and 0.8, since the lattice is 100 wide.
	nucleus   = ElementGaussian(0.3, 10.0)
	atoms     = [nucleus, nucleus]
	positions = [[0.2, 0, 0], [0.8, 0, 0]]
	
	# Assemble the model, discretize and build the Hamiltonian
	model = Model(lattice, atoms, positions; terms=[Kinetic(), AtomicLocal()])
	basis = PlaneWaveBasis(model; Ecut=300, kgrid=(1, 1, 1));
	ham   = Hamiltonian(basis)
	
	# Extract the total potential term of the Hamiltonian and plot it
	potential = DFTK.total_local_potential(ham)[:, 1, 1]
	rvecs = collect(r_vectors_cart(basis))[:, 1, 1]  # slice along the x axis
	x = [r[1] for r in rvecs]                        # only keep the x coordinate
	plot(x, potential, label="", xlabel="x", ylabel="V(x)")
	
end

# ╔═╡ 053648f0-d079-48f3-bc74-264738693d1f
md"""
This potential is the sum of two "atomic" potentials (the two "Gaussian" elements).
Due to the periodic setting we are considering interactions naturally also occur
across the unit cell boundary (i.e. wrapping from `100` over to `0`).
The required periodization of the atomic potential is automatically taken care,
such that the potential is smooth across the cell boundary at `100`/`0`.

With this setup, let's look at the bands:
"""

# ╔═╡ 8308ea7e-0f61-4630-ac03-e7c71c94fbd7
let
	bands = compute_bands(basis; n_bands=6, kline_density=500)
	plot_bandstructure(bands)
end

# ╔═╡ f0f4d26f-610c-4658-91ae-59a3b589b159
md"""
The bands are noticeably different.
 - The bands no longer overlap, meaning that the spectrum of $\mathcal{H}$ is
   no longer all positive numbers, but it has gaps.

 - The two lowest bands are almost flat. This is because they represent
   two tightly bound and localized electrons inside the two Gaussians.

 - The higher the bands are in energy, the more free-electron-like they are.
   In other words the higher the kinetic energy of the electrons, the less they feel
   the effect of the two Gaussian potentials. As it turns out the curvature of the bands,
   (the degree to which they are free-electron-like) is highly related to the delocalization
   of electrons in these bands: The more curved the more delocalized. In some sense
   "free electrons" correspond to perfect delocalization.

 - Try playing with the parameters of the Gaussian potentials by setting
   ```julia
   nucleus   = ElementGaussian(α, L)
   ```
   with different $α$ and $L$. You should notice the influence on the bands. Pay particular
   attention to the relation between the depth of the potential and the shape of the bands.
"""

# ╔═╡ 6092a639-f4ba-4833-9f2a-b9a513e30bb1
md"""
## Deep-dive: Bloch-Floquet theory (optional)

This section provides more mathematical details to the high-level discussion we provided above.
"""

# ╔═╡ 6bd7bc71-b827-4f33-8302-83520ba7cd7b
md"""
We start by considering the relationship between the Laplacian and the Fourier transform.
Since $\widehat{\laplacian
 f} (p) = - |p|^2 \hat f (p)$, the Fourier
transform 
```math
\begin{align}
    \hat f(p) = \frac1{(2 \pi)^{\frac{3}{2}}} \int_{\mathbb R^3} f(x) e^{-i p \cdot x} dx
\end{align}
```
 provides exactly the unitary transformation required to
diagonalize the Laplacian $\laplacian$. Notably eigenfunctions are the
plane wave $e^{i p \cdot x}$ with $p \in \mathbb R^3$. 
Note that other Fourier transform conventions exist. 
Here we explicitly choose the one which is unitary (which is *not* the
one used in DFTK).

In our problem $\opH = - \frac1{2} \laplacian+ V$ with $V$
$\mathbb L$-periodic. We thus need an equivalent operation, that is
better adapted to the $\mathbb L$-periodic setting. This turns out to be
the *Bloch transform*.

"""

# ╔═╡ 083e1f7a-e09e-479a-9cd9-83203e8c23dc
md"""
- First, we rewrite the plane waves as 
  ```math
  	\{ e^{i (k+G) \cdot x} \} _{k \in \Omega^*, G \in \mathbb L^*}
  ```
  where we note that for all $p \in \mathbb R^3$ a decomposition $p = k+G$ with $G \in \mathbb L^*, k \in \overline{\Omega^*}$ is always possible.


- Then, we use this to perform the inverse transform in two steps
  ```math
  \begin{align}
      f(x) &= \frac1{(2 \pi)^{\frac{3}{2}}} \int_{\mathbb R^3} \hat f (p) e^{i p \cdot x}
    \\
    &= \frac1{(2 \pi)^{\frac{3}{2}}} \sum_{G \in \mathbb L^*} \int_{\Omega^*} \hat f(k + G) e^{i (k+G)\cdot x} dk
  \end{align}
  ```
  where we take an $f \in C_0^\infty (\mathbb R^3)$. 
  By density arguments our discussion can be easily extended to $L^2 (\mathbb R^3)$.

"""

# ╔═╡ f8c35c68-0e60-45d7-ab93-5aace2bd1589
md"""
- From this, we introduce the **Bloch transform** $\bloch$ 
  ```math
  \begin{align}
      (\bloch f) (x) &= \Big(k \mapsto f_k(x)\Big)
  \end{align}
  ```
  where
  ```math
  \begin{align}
   f_k(x) &= \frac1{\sqrt{|\Omega|}} \sum_{G \in \mathbb L^*} \hat f (k + G) e^{i G\cdot x}
     \\
    &= \sum_{G \in \mathbb L^*} \hat f (k+G) e_G(x).
  \end{align}
  ```
  In other words the Bloch transform maps $f(x)$ to another function $k \mapsto f_k(x)$.
  We arrive at the reconstruction formula by considering the second step of the inverse Fourier transform
  ```math
  \begin{align}
    f(x) = \frac1{\sqrt{|\Omega^*|}} \int_{\Omega^*} f_k(x) e^{i k\cdot x} dk
  \end{align}
  ```
  where we used that 
  ```math
  \begin{align}
    \frac1{\sqrt{|\Omega|}} = \frac{\sqrt{|\Omega^*|}}{(2 \pi)^{\frac{3}{2}}}.
  \end{align}
  ```
  Notice that this normalization convention is chosen such that the normalized plane waves appear in the Bloch transform.

"""

# ╔═╡ db338525-379e-44e5-b3eb-8d98e33d171b
md"""
- The Bloch transform satisfies

  ```math
  \begin{align}
  f_k(x+R) &= f_k(x) && \text{(Periodicity)}
  \\
  f_{k+G} (x) &= e^{- i G \cdot x} f_k(x). &&\text{(Born-von Karman condition)}
  \end{align}
  ```
  
- As a result we can always restrict $f_k(x)$ to $k \in \Omega^*$ and $x \in \Omega$.

"""

# ╔═╡ 9a6793e3-a46a-4736-9521-c5c6a1636aa3
md"""
- An explicit construction of $f_k$ from $f$ is available via **Poisson's formula** : 
  ```math
  \begin{align}
    f_k(x) &= \frac1{\sqrt{|\Omega|}} \sum_{G \in \mathbb L^*} \hat f (k + G) e^{i G \cdot x} 
    \\
    &= \frac1{\sqrt{|\Omega^*|}} \sum_{R \in \mathbb L} f(x+R) e^{- i k \cdot (x +R)}
  \end{align}
  ```
  This can be seen from 
  ```math
  \begin{align}
    \hat f(k+G) &= \frac1{(2 \pi)^{\frac{3}{2}}} \int_{\mathbb R^3} f(x) e^{-i (k + G) \cdot x} dx
    \\
    &= \frac1{(2 \pi)^{\frac{3}{2}}} \sum_{R \in \mathbb L} \int_\Omega f(R+y) e^{-i (k+G) \cdot (R+y)} dy
    \\
    &= \frac1{\sqrt{|\Omega^*|}} \int_\Omega \left ( \sum_{R \in \mathbb L} f(R+y) e^{-i k \cdot  (R+y)} \right ) \frac{e^{-i G \cdot y}}{\sqrt{|\Omega|}} dy
  \end{align}
  ```
  and the fact that 
  ```math
  \begin{align}
    \sum_{G \in \mathbb L^*} e^{i G \cdot (x-y)} = \delta (x-y).
  \end{align}
  ```
  From Poisson's formula it is clear that *Parceval's identity* 
  ```math
    \int_{\mathbb R^3} | f(x) |^2 dx = \int_{\Omega^*} dk \int_\Omega dy \ |f_k(y)|^2 
  ```
  holds, since 
  ```math
  \begin{align}
    \int_{\Omega^*} dk \int_\Omega dy \ |f_k(y)|^2 &= \frac1{|\Omega^*|} \sum_{R, R' \in \mathbb L} \left ( \int_{\Omega^*} e^{i k \cdot (R-R')} dk \right ) \left ( \int_\Omega \overline{f(y+R')} f(y+R) dy \right )
    \\
    &= \sum_{R \in \mathbb L} \int_\Omega |f(y+R)|^2 \ dy 
    \\
    &= \int_{\mathbb R^3} |f(x)|^2 \  dx
  \end{align}
  ```
"""

# ╔═╡ 17d99340-383a-4b29-a1a2-244024c71155
md"""
We summarize our findings in a Theorem :

!!! note "Theorem 2 (Bloch-Floquet transform)"
	The Bloch-Floquet transform $\bloch$ is the unitary mapping from
	$L^2(\mathbb R^d)$ to the vector space
	$L^2_\text{qp} (\Omega^*, L^2_\text{per}(\Omega))$ (defined below)
	with the property that each $f \in C_0^\infty (\mathbb R^d)$ is mapped to 
	```math
	\Omega^* \ni K \mapsto f_k = \sum_{R \in \mathbb L } u(\bullet + R) e^{- k \cdot (\bullet + R)} \in L^2_\text{per}(\Omega).
	```
	The extension to $L^2(\mathbb R^d)$ follows by density of $C_0^\infty$ in $L^2$.

	In this $L^2_\text{qp} (\Omega^*, L^2_\text{per}(\Omega))$ is the space of $\mathbb L^*$-*quasiperiodic* functions which are $L^2_\text{per}(\Omega)$-valued :
	```math
	\begin{align}
	        L^2_\text{qp} (\Omega^*, L^2_\text{per}(\Omega)) = 
	         \left \{ \mathbb R^d \ni k \mapsto u_k \in L^2_\text{per}(\Omega) \middle \vert \int_{\Omega^*} \| u_k \|^2_{L^2_\text{per}(\Omega)} \ dk < \infty \right .
	        \\
	        \text{where } u_{k+G} = u_k e^{-i G \cdot x}  \text{ for all } G \in\mathbb L^*
	\\
		\text{and almost every } k \in \mathbb R^d \bigg \}.
	    
	\end{align}
	```
	On this vector space we have the inner product 
	```math
	\langle f,g \rangle_{L^2_\text{qp} (\Omega^*, L^2_\text{per}(\Omega))} = \frac1{|\Omega^*|} \int_{\Omega^*} \langle f_k, g_k \rangle_{L^2_\text{per}(\Omega)} dk.
	```
"""

# ╔═╡ fcf275be-b4c2-4d06-b502-a0c52e74b533
md"""
!!! tip "Remark"
	By way of Parceval's identity we can also understand $L^2_\text{qp} (\Omega^*, L^2_\text{per}(\Omega)) \equiv L^2(\Omega^* \times \Omega)$, i.e. $\bloch$ maps single variable functions of $y \in \mathbb R^d$ to functions of two variables : $k \in \Omega^*$ and $x \in \Omega$.

	We can illustrate it as
	```math
	\begin{align}
		L^2(\mathbb R^d) && \stackrel{\bloch}{\longrightarrow} && L^2_\text{qp} (\Omega^*, L^2_\text{per}(\Omega)) = L^2(\Omega^* \times \Omega)
	\\
		\mathbb R^d \ni x \mapsto  f(x) \in \mathbb C &&&&
		\Omega^* \times \Omega \ni (k,y) \mapsto f_k(y) \in \mathbb C
	\end{align}
	```
"""

# ╔═╡ 9d4a8178-52f4-4846-94d6-74b8d68abbb4
# use y instead of x??
md"""
In physics texts it is frequently written that the Bloch transform admits the decomposition $\psi(x) = u_{k'}(x) e^{i k' \cdot x}$ for an eigenstate $\psi$, where $u_k$ is $\mathbb L$-periodic. We want to make the connection of this with the present discussion.
As we will see in the following, an eigenvalue can be restricted to a single value of $k$.

Thus, we define $ψ_k(x)$ as 
```math
k \mapsto \begin{cases} u_{k'}(x) e^{- i G \cdot x} & \text{if } k = k'+G \text{ for } G \in \mathbb L \\ 0 & \text{otherwise} \end{cases}
```
The reconstruction formula then gives
```math
\begin{align}
	ψ(x) &= \frac1{\sqrt{|\Omega^*|}} \int_{\Omega^*} ψ_k(x) e^{ik \cdot x} dk
	\\ &= \frac1{\sqrt{|\Omega^*|}} u_{k'}(x) e^{i k' \cdot x}
\end{align}
```
which is just a scaled form of the commonly quoted representation.
"""

# ╔═╡ dede1f95-42ac-4325-ad20-07f884b58e2e
md"-----"

# ╔═╡ 21523a64-b1a9-422d-8dc0-ebb324f35380
#check the formula
md"""
The next natural question is to consider the Bloch-Floquet transform of operators.
We first consider the Laplace operator.
Using Poisson's formula,
```math
	f_k(x) = \frac1{\sqrt{|\Omega^*|}} \sum_{R \in \mathbb L} f(x + R) e^{- i k \cdot (x+R)}
```
we obtain
```math
\begin{align}
	(-i \nabla_x + k) f_k(x) &= \frac1{\sqrt{|\Omega^*|}} \sum_{R \in \mathbb L} [- i \nabla_x f(x+R)] e^{- ik \cdot (x+R)} +
	\\ &\quad + f(x+R) e^{-i k \cdot (x+R)} [(-i)^2 k + k] 
	\\ &= \frac1{\sqrt{|\Omega^*|}} \sum_{R \in \mathbb L} [- i ∇_x f(x+R)] e^{-i k (x+R)}
	\\ &= [- i ∇_x f]_k (x).
\end{align}
```
This is clearly true if $f \in C_0^\infty (\mathbb R^3)$.
Again by density arguments this shows 
```math
\bloch H^1(\mathbb R^3) = L^2_\text{qp} (\Omega^*, H^1_\text{per}(\Omega)).
```
Similarly
```math
\bloch (- \laplacian _x) f = (- \laplacian_x f)_k = \lvert -i \nabla_k + k |^2 f_k = \lvert -i \nabla_k + k |^2 \bloch f
```
suggesting
```math
\bloch H^2(\mathbb R^3) = L^2_\text{qp} (\Omega^*, H^2_\text{per} (\Omega))
```
i.e. that regularity in $f$ induces equal regularity in $f_k$ for all $k \in \Omega^*$.

Notably, this also suggests the Bloch transform of the kinetic energy operator $- \laplacian_x$ as
```math
\bloch (- \laplacian_x) \bloch^{-1} = \lvert - i ∇_x + k|^2
```
interpreted in the sense of distributions and defined on the Hilbert space $L^2_\text{qp} (\Omega^*, L^2_\text{per} (\Omega))$ with domain $D( \bloch (- \laplacian_x) \bloch^{-1} ) = L^2_\text{qp} (\Omega^*, H^2_\text{per} (\Omega))$.

This result is formalized in the following theorem.
"""

# ╔═╡ 62e5c394-dad2-4a62-a7ae-d4999f741505
md"""
Denoting $\mathscr L (V) = \{ f : V \to V | f \text{ linear} \}$ the space of linear maps between a function space and itself, we have :

!!! note "Theorem 3 (Bloch fibers)"
	Any bounded $\mathbb L$-periodic operator $\opA$ on $L^2(\mathbb R^3)$ is decomposed by the Bloch-Floquet transform in the following sense : there exists a function $k \mapsto \opA_k$ in $L^\infty_\text{qp} (\mathbb R^3, \mathscr L(L^2_\text{per} (\Omega))$ such that for any $u \in L^2(\mathbb R^3)$, any $G \in \mathbb L^*$, and almost any $k \in \Omega^*$ is holds that 
	```math
	\begin{align}
	(\opA u)_k = \opA_k u_k && \opA_{k+G} = e^{-i G \cdot x} \opA_k e^{i G \cdot x}
	\end{align}
	```
	where the operators $(\opA_k)_{k \in \Omega^*}$ on $\hilbert = L^2_\text{per} (\Omega)$ are called the *Bloch fibers* of $\opA$.

	The Bloch decomposition can be extended to unbounded operators on $L^2(\mathbb R^3)$, in which case $\opA_k$ are unbounded as well.

!!! warning "Example 1 (Laplace operator)"
	For $\opA = - \laplacian$, the Bloch fibers are $\opA_k = \lvert - i \nabla_x + k|^2$.
"""

# ╔═╡ 6da749bc-0181-49c2-8a7a-9c180930dce6
md"""
## Fibers of periodic Schrödinger operators

We consider $\opH = - \frac1{2} \laplacian + V$ on $L^2(\mathbb R^3)$ with domain $D(\opH) = H^2(\mathbb R^3)$ and $V \in L^2_\text{per} (\mathbb R^3)$ such that $\opH$ is $\mathbb L$-periodic.
Recall this operator is self-adjoint and as *no* eigenvalues.

- The fibers of $\opH$ are $\opH_k = ( - i \nabla_x + k)^2 + V$, where we note
  -  $\opH_k$ is a self-adjoint and bounded from below operator on $L^2_\text{per}(\Omega)$ with domain $H^2_\text{per}(\Omega)$ and form domain $H^1_\text{per}(\Omega)$
  -  $\opH_k$ has a compact resolvent $\resolvent(\opH_k)$, such that each $\opH_k$ has a purely discrete spectrum with eigenvalues accumulating at $+ \infty$ and eigenfunctions forming an orthonormal basis for $L^2_\text{per}(\Omega)$.


- To understand the spectrum of $\opH$, it suffices completely to study the spectrum of all $\opH_k$. 
  Indeed
  ```math
  	\sigma(\opH) = \bigcup_{k \in \overline{\Omega^*}} \sigma(\opH_k).
  ```


- We thus seek an $L^2_\text{per} (\Omega)$ orthonormal basis $\{ (λ_{k,n}, \phi_{k,n} )\}_{n \in \mathbb N} \subset (\mathbb R \times L^2_\text{per}(\Omega))^{\mathbb N}$ of eigenstates of $\opH_k$ :
  ```math
  \begin{align}
  	\opH_k \phi_{k,n} &= λ_{k,n} \phi_{k,n}
  	\\
  	\langle \phi_{k,n} , \phi_{k,m} \rangle_{L^2_\text{per}(\Omega)} &= \delta_{mn} && \forall n,m \in \mathbb N.
  \end{align}
  ```

  By convention, we sort for each $k \in \overline{\Omega^*}$ : $λ_{k,1} ≤ λ_{k,2} ≤ λ_{k,3} ≤ …$, thus counting multiplicities.
  Moreover, for each $n \in \mathbb N$ we define the mapping
  ```math
  	λ_n : \overline{\Omega^*} \to \mathbb R : k ↦ λ_{k,n}
  ```
  which is a Lipshitz function called the $n$-th **energy band**.


- Each band has a well-defined minimum/maximum, forming the following interval
  ```math
  	 λ_n \left (\overline{\Omega^*} \right ) = \left [ \min_{k ∈ \overline{\Omega^*}} λ_n(k) , \max_{k ∈ \overline{\Omega^*}} λ_n(k)  \right ].
  ```
  Moreover,
  ```math
  σ(\opH) = \bigcup_{n ∈ \mathbb N} λ_n \left ( \overline{\Omega^*} \right )
  ```

---

An image to have in the head is 
"""

# ╔═╡ 08b20ea7-f3c3-4302-bfc3-a13bbd57930a
TikzPicture(L"""
	%connecting dashes
	\draw[dashed,blue] (220,230) -- (250,230) --  (250,200) -- (189.5,200);
	\draw[dashed,cyan] (220,190) -- (249.6,190) -- (249.6,160) -- (190.5,160);
	\draw[dashed,teal] (125.5,170) -- (250.4,170) -- (250.4,130) -- (220,130);
	\draw[dashed,olive] (165.5,100) -- (250,100) -- (250,80) -- (190,80);

	
	%1st band
	\draw[blue,thick]   
	(100,230) .. controls (108.74,223.57) and (104.5,213.5) .. 
	(127.5,211.5) .. controls (150.5,209.5) and (175.02,200.09) .. 
	(189.5,200) node[above]{$\lambda_1$} .. controls (203.98,199.91) and (211,210) ..
	(220,230) ;

	%2nd band
	\draw[cyan,thick]    
	(100,190) .. controls (106.84,188.97) and (135.14,187.24) .. 
	(157.5,179) .. controls (179.86,170.76) and (176.5,160) .. 
	(190.5,160) node[above]{$\lambda_2$} .. controls (204.5,160) and (200.5,178) ..
	(220,190) ;

	%3rd band
	\draw[teal,thick]  (100,130) .. controls (114.5,130) and (111.5,170) .. 
	(125.5,170) .. controls (139.5,170) and (135.5,140) .. 
	(149.5,140) node[above]{$\lambda_3$} .. controls (163.5,140) and (166.5,162) ..
	(175.5,146)  .. controls (184.5,130) and (212.5,143) .. 
	(220,130) ;

	%4th band
	\draw[olive,thick]   
	(100,90) .. controls (115.5,84) and (107.5,94) .. 
	(114.5,94) .. controls (121.5,94) and (119.5,81) .. 
	(129.5,83)  .. controls (139.5,83) and (135.5,87) .. 
	(143.5,91) node[below]{$\lambda_4$} .. controls (151.5,98) and (146.55,83.87) .. 
	(155.5,85) .. controls (164.45,86.13) and (155.5,101) .. 
	(165.5,100)  .. controls (175.5,101) and (172.5,81) .. 
	(190,80)  .. controls (208.5,81) and (205.5,105) .. 
	(220,90) ;

	%axes
	\draw[<-,>=latex] (100,70) node[left]{$\lambda$} -- (100,240) -- (165,240) node[below]{$k \in \Omega^*$} -- (220,240) -- (220,70) ;

	%spectrum
	\draw[<-,>=latex] (250,70) node[right]{$\sigma(- \Delta/2 + V)$} -- (250,240)  ;

	\draw[very thick,blue] (250,230) -- (250,200);
	\draw[very thick,cyan] (249.6,190) -- (249.6,160);
	\draw[very thick,teal] (250.4,170) -- (250.4,130);
	\draw[very thick,olive] (250,100) -- (250,80);

	%band gaps
	\draw[<->] (260,100) -- (260,115) node[right]{Bandgap} -- (260,130);
	\draw[<->] (260,200) -- (260,195) node[right]{Bandgap} -- (260,190);




""",width="17cm",options="xscale=0.05,yscale=-0.05",preamble=raw"\usepackage{amsfonts}")

# ╔═╡ 2cf96d2f-8a37-44e8-8b54-5e2e35542fcc
md"""
---
"""

# ╔═╡ c24c96be-57d6-48dc-b79e-c9a13b9493b3
TableOfContents()

# ╔═╡ 7a952990-b524-43b6-a01b-cf6b6c6eaa14
let
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	
	Sidebar(Markdown.parse(read("sidebar.md", String)), 330)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DFTK = "acf6eb54-70d9-11e9-0013-234b7a5f5337"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"

[compat]
DFTK = "~0.7.18"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.7"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.74"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "60f5893cceb4916d3b04a92be4744d85d58ba4c8"

[[deps.ADTypes]]
git-tree-sha1 = "27cecae79e5cc9935255f90c53bb831cc3c870d7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.18.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AtomsBase]]
deps = ["LinearAlgebra", "PeriodicTable", "Preferences", "Printf", "Requires", "StaticArrays", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "922c2469c526996566dbabd273d15701ed2aacfe"
uuid = "a963bdd2-2df7-4f54-a1ee-49d51e6be12a"
version = "0.5.2"

    [deps.AtomsBase.extensions]
    AtomsBaseAtomsViewExt = "AtomsView"

    [deps.AtomsBase.weakdeps]
    AtomsView = "ee286e10-dd2d-4ff2-afcb-0a3cd50c8041"

[[deps.AtomsCalculators]]
deps = ["AtomsBase", "Compat", "StaticArrays", "Test", "Unitful"]
git-tree-sha1 = "99c477979252acdb189f3e4161970e42634e426a"
uuid = "a3e0e189-c65a-42c1-833c-339540406eb1"
version = "0.2.3"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bravais]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a39abce7826834ba1a517e44e9ed673e6621192a"
uuid = "ada6cbde-b013-4edf-aa94-f6abe8bd6e6b"
version = "0.2.5"

[[deps.Brillouin]]
deps = ["Bravais", "DirectQhull", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Reexport", "Requires", "StaticArrays"]
git-tree-sha1 = "c1cb59741b310573c43406c652031e3228d2f187"
uuid = "23470ee3-d0df-4052-8b1a-8cbd6363e7f0"
version = "0.5.31"

    [deps.Brillouin.extensions]
    BrillouinMakieExt = "Makie"
    BrillouinPlotlyJSExt = "PlotlyJS"
    BrillouinSpglibExt = "Spglib"

    [deps.Brillouin.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"
    Spglib = "f761d5c5-86db-4880-b97f-9680a7cccfb5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CUDA_Driver_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2023be0b10c56d259ea84a94dbfc021aa452f2c6"
uuid = "4ee394cb-3365-5eb0-8335-949819d2adfc"
version = "13.0.2+0"

[[deps.CUDA_Runtime_jll]]
deps = ["Artifacts", "CUDA_Driver_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "92cd84e2b760e471d647153ea5efc5789fc5e8b2"
uuid = "76a88914-d11a-5bdc-97e0-2f5a05c973a2"
version = "0.19.2+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ComponentArrays]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "ConstructionBase", "Functors", "LinearAlgebra", "StaticArrayInterface", "StaticArraysCore"]
git-tree-sha1 = "29dfb059630454c0282779d58e1a8539573b5945"
uuid = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
version = "0.15.30"

    [deps.ComponentArrays.extensions]
    ComponentArraysGPUArraysExt = "GPUArrays"
    ComponentArraysKernelAbstractionsExt = "KernelAbstractions"
    ComponentArraysOptimisersExt = "Optimisers"
    ComponentArraysReactantExt = "Reactant"
    ComponentArraysRecursiveArrayToolsExt = "RecursiveArrayTools"
    ComponentArraysReverseDiffExt = "ReverseDiff"
    ComponentArraysSciMLBaseExt = "SciMLBase"
    ComponentArraysTrackerExt = "Tracker"
    ComponentArraysZygoteExt = "Zygote"

    [deps.ComponentArrays.weakdeps]
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Optimisers = "3bd65402-5787-11e9-1adc-39752487f4e2"
    Reactant = "3c362404-f566-11ee-1572-e11a4b42c853"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SciMLBase = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CrystallographyCore]]
deps = ["LinearAlgebra", "StaticArrays", "StructEquality"]
git-tree-sha1 = "521a8ed44455592672d887632eef056c1ba4678c"
uuid = "80545937-1184-4bc9-b283-396e91386b5c"
version = "0.3.3"

[[deps.DFTK]]
deps = ["AbstractFFTs", "Artifacts", "AtomsBase", "AtomsCalculators", "Brillouin", "Dates", "DftFunctionals", "DifferentiationInterface", "DocStringExtensions", "FFTW", "ForwardDiff", "GPUArraysCore", "Interpolations", "IterTools", "KrylovKit", "Libxc", "LineSearches", "LinearAlgebra", "LinearMaps", "MPI", "Markdown", "Optim", "PeriodicTable", "PkgVersion", "Polynomials", "PrecompileTools", "Preferences", "Printf", "PseudoPotentialData", "PseudoPotentialIO", "Random", "Roots", "SparseArrays", "SpecialFunctions", "Spglib", "StaticArrays", "Statistics", "TimerOutputs", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "6533c3c3740166c56bccc9e6260ab794c147acf7"
uuid = "acf6eb54-70d9-11e9-0013-234b7a5f5337"
version = "0.7.18"

    [deps.DFTK.extensions]
    DFTKAMDGPUExt = "AMDGPU"
    DFTKCUDAExt = "CUDA"
    DFTKGenericLinearAlgebraExt = "GenericLinearAlgebra"
    DFTKGeometryOptimizationExt = "GeometryOptimization"
    DFTKIntervalArithmeticExt = "IntervalArithmetic"
    DFTKJLD2Ext = "JLD2"
    DFTKJSON3Ext = "JSON3"
    DFTKPlotsExt = "Plots"
    DFTKWannier90Ext = "wannier90_jll"
    DFTKWannierExt = "Wannier"
    DFTKWriteVTKExt = "WriteVTK"

    [deps.DFTK.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GenericLinearAlgebra = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
    GeometryOptimization = "673bf261-a53d-43b9-876f-d3c1fc8329c2"
    IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
    JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
    Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
    Wannier = "2b19380a-1f7e-4d7d-b1b8-8aa60b3321c9"
    WriteVTK = "64499a7a-5c06-52f2-abe2-ccb03c286192"
    wannier90_jll = "c5400fa0-8d08-52c2-913f-1e3f656c1ce9"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DftFunctionals]]
deps = ["ComponentArrays", "DiffResults", "ForwardDiff"]
git-tree-sha1 = "b904a2c78f0e768187a13f425719557b2d1db47e"
uuid = "6bd331d2-b28d-4fd3-880e-1a1c7f37947f"
version = "0.3.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "6d5153dc500d644d4d672723aa27a614ee84ab3b"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.11"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DirectQhull]]
deps = ["Qhull_jll"]
git-tree-sha1 = "49951cee00be263d4dc187f38aa96333d74f4d1c"
uuid = "c3f9d41a-afcb-471e-bc58-0b8d83bd86f4"
version = "0.2.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "7ea1aa5869e2626ccae84480e4f37185bc6f41d3"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.2.3"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "5bfcd42851cf2f1b303f51525a54dc5e98d408a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.15.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "9340ca07ca27093ff68418b7558ca37b05f8aeb1"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.29.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cd33c7538e68650bd0ddbb3f5bd50a4a0fa95b50"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.0"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "8e2d86e06ceb4580110d9e716be26658effc5bfd"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "da121cbdc95b065da07fbb93638367737969693f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.8+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_ICU_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "6ccbc4fdf65c8197738c2d68cc55b74b19c97ac2"
uuid = "655565e8-fb53-5cb3-b0cd-aec1ca0647ea"
version = "2.8.1+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "3d468106a05408f9f7b6f161d9e7715159af247b"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.12.2+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "20b6765a3016e1fca0c9c93c80d50061b94218b7"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "69.1.0+0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "5b6bb73f555bc753a6153deec3717b8904f5551c"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.3.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "PackageExtensionCompat", "Printf", "Random", "VectorInterface"]
git-tree-sha1 = "6dcba71deb016d646f1c1bcfcaacc764a198b8e6"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.10.2"
weakdeps = ["ChainRulesCore"]

    [deps.KrylovKit.extensions]
    KrylovKitChainRulesCoreExt = "ChainRulesCore"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd10d2cc78d34c0e2a3a36420ab607b611debfbb"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.7"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.Libxc]]
deps = ["Libxc_GPU_jll", "Libxc_jll"]
git-tree-sha1 = "31acc2c3c93b8d3c4bb43c816e082fcf971c43d3"
uuid = "66e17ffc-8502-11e9-23b5-c9248d0eb96d"
version = "0.3.20"

    [deps.Libxc.extensions]
    LibxcCudaExt = "CUDA"

    [deps.Libxc.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.Libxc_GPU_jll]]
deps = ["Artifacts", "CUDA_Runtime_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "25bde0bdc076ff720ce0b2383e044f0b62048a03"
uuid = "25af9330-9b41-55d4-a324-1a83c0a0a1ac"
version = "7.0.0+2"

[[deps.Libxc_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "94f7feffa99688628909494a890581fe3b6148dd"
uuid = "a56a6d9d-ad03-58af-ab61-878bf78270d6"
version = "7.0.0+2"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "4adee99b7262ad2a1a4bbbc59d993d24e55ea96f"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.4.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearMaps]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7f6be2e4cdaaf558623d93113d6ddade7b916209"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.11.4"
weakdeps = ["ChainRulesCore", "SparseArrays", "Statistics"]

    [deps.LinearMaps.extensions]
    LinearMapsChainRulesCoreExt = "ChainRulesCore"
    LinearMapsSparseArraysExt = "SparseArrays"
    LinearMapsStatisticsExt = "Statistics"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MPI]]
deps = ["Distributed", "DocStringExtensions", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "PkgVersion", "PrecompileTools", "Requires", "Serialization", "Sockets"]
git-tree-sha1 = "a61ecf714d71064b766d481ef43c094d4c6e3c52"
uuid = "da04e1cc-30fd-572f-bb4f-1f8673147195"
version = "0.20.23"

    [deps.MPI.extensions]
    AMDGPUExt = "AMDGPU"
    CUDAExt = "CUDA"

    [deps.MPI.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "9341048b9f723f2ae2a72a5269ac2f15f80534dc"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.3.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "e214f2a20bdd64c04cd3e4ff62d3c9be7e969a59"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.4+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "ec764453819f802fc1e144bfe750c454181bd66d"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.8+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "386b47442468acfb1add94bf2d85365dea10cbab"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ad31332567b189f508a3ea8957a2640b1147ab00"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "61942645c38dd2b5b78e2082c9b51ab315315d10"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.2"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9dd97171646850ee607593965ce1f55063d8d3f9"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.PeriodicTable]]
deps = ["Base64", "Unitful"]
git-tree-sha1 = "238aa6298007565529f911b734e18addd56985e1"
uuid = "7b2266bf-644c-5ea3-82d8-af4bbd25a884"
version = "1.2.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "f202a1ca4f6e165238d8175df63a7e26a51e04dc"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.7"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "dacc8be63916b078b592806acd13bb5e5137d7e9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5b2a443bc3ac0c5c1e40dc39fdd332c49b412833"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.74"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.0"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "02148a0cb2532f22c0589ceb75c110e168fb3d1f"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "21.9.0+0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PseudoPotentialData]]
deps = ["Artifacts", "Compat", "DocStringExtensions", "LazyArtifacts", "TOML"]
git-tree-sha1 = "754ae22669bfcb466fd4b0b3cacc1060fc0a7c57"
uuid = "5751a51d-ac76-4487-a056-413ecf6fbe19"
version = "0.2.4"
weakdeps = ["AtomsBase"]

    [deps.PseudoPotentialData.extensions]
    PseudoPotentialDataAtomsBaseExt = "AtomsBase"

[[deps.PseudoPotentialIO]]
deps = ["EzXML", "LinearAlgebra"]
git-tree-sha1 = "88cf9598d70015889c99920ff3dacca0eb26ae90"
uuid = "cb339c56-07fa-4cb2-923a-142469552264"
version = "0.1.1"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6f3ac0623e1173c006cc7377798ec3fb33fa504"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1004+0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "8a433b1ede5e9be9a7ba5b1cc6698daa8d718f1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Spglib]]
deps = ["CrystallographyCore", "StaticArraysCore", "StructEquality", "spglib_jll"]
git-tree-sha1 = "3260e525f1aa94d4b420b5ff845c672190c25dc5"
uuid = "f761d5c5-86db-4880-b97f-9680a7cccfb5"
version = "0.9.7"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "49440414711eddc7227724ae6e570c7d5559a086"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "064b532283c97daae49e544bb9cb413c26511f8c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.8"

[[deps.StructEquality]]
deps = ["Compat"]
git-tree-sha1 = "192a9f1de3cfef80ab1a4ba7b150bb0e11ceedcf"
uuid = "6ec83bb0-ed9f-11e9-3b4c-2b04cb4e219c"
version = "2.1.0"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "79529b493a44927dd5b13dde1c7ce957c2d049e4"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.0"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "tectonic_jll"]
git-tree-sha1 = "79e2d29b216ef24a0f4f905532b900dcf529aa06"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "6258d453843c466d84c17a58732dda5deeb8d3af"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.24.0"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

[[deps.UnitfulAtomic]]
deps = ["Unitful"]
git-tree-sha1 = "903be579194534af1c4b4778d1ace676ca042238"
uuid = "a7773ee8-282e-5fa2-be4e-bd808c38a91a"
version = "1.0.0"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9166406dedd38c111a6574e9814be83d267f8aec"
uuid = "409d34a3-91d5-4945-b6ec-7529ddf182d8"
version = "0.5.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.spglib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "bc328924cf4975fe49e6416f7e1622e8ceda55e8"
uuid = "ac4a9f1e-bdb2-5204-990c-47c8b2f70d4e"
version = "2.1.0+0"

[[deps.tectonic_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "Graphite2_jll", "HarfBuzz_ICU_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "54867b00af20c70b52a1f9c00043864d8b926a21"
uuid = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"
version = "0.13.1+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "fbf139bce07a534df0e699dbb5f5cc9346f95cc1"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.9.2+0"
"""

# ╔═╡ Cell order:
# ╟─e747d030-7487-4598-99d6-93afefd425e1
# ╟─a0407e4d-3edc-48ea-a335-b815b796e378
# ╠═c7351c0e-c5d5-43e8-b0db-24c7b1985ae1
# ╟─1441bdfc-715e-46da-a4d2-8b788c71853f
# ╟─7971ffe3-39f0-4cf8-9cee-fd55483a8ed9
# ╟─8d9fe76f-dfb3-4683-83c2-5d7cd2aff712
# ╟─df23c23e-b509-4f48-9ab9-dd2b59834d01
# ╟─9b17dde4-6652-4800-a568-1244cf519593
# ╟─fbbe239e-c06e-41b1-b1e3-857585f994bd
# ╟─ea689cad-2859-41b0-aef3-36651a8d166b
# ╟─6ba1ce9a-6ced-4a82-9b5a-6b07a25c534b
# ╟─6c5c1231-9e16-4a43-8dec-68df585d78b9
# ╟─3563cb4b-1e0d-47fd-8546-7a3885d397d1
# ╟─b90fe635-247e-4842-acbc-2a07204c89c8
# ╟─20029c12-e608-4235-96c6-19fab70f02f3
# ╟─9cf1f790-835a-438d-be09-2b4eea978351
# ╟─7c6dcdc6-692d-4a47-85ac-b9baee3a1be5
# ╟─696d9342-d8fd-442f-82ed-81032b6fcba2
# ╟─fec13421-c63f-4e96-9d0d-4d4f43ab172d
# ╟─2386c0dd-0b29-4168-9f25-2187d43b1461
# ╠═3ceea351-5110-4a6d-a525-fdf259c74675
# ╠═2a765fbb-759b-4480-913c-e8f01ef321a0
# ╠═5aeaa6b8-50b9-4c5f-bcce-592965b68353
# ╟─700ac798-5f2f-456e-8be8-cd40c76e6c96
# ╠═15fb9074-edbf-4458-9a41-bd3d8707a501
# ╟─32086662-d60d-47cd-acb6-0f1a4f18e351
# ╟─1c07a977-68b7-44e1-9adf-ad1af1b44ebc
# ╟─e2b26548-3f7d-4950-b22c-e7106912df49
# ╟─d7644757-dcee-4414-9dda-e09476779ab9
# ╠═0d419731-99e9-4380-a667-d2f935044f17
# ╟─0b85ae26-7e44-4506-a096-71d2a49def36
# ╟─da661541-2007-4681-baea-36164b43efc7
# ╟─0309d400-92ed-46ad-afce-a516a7f1380b
# ╟─a790a782-2252-4650-933a-4c002bd3e943
# ╟─3de086c4-877d-4d46-825a-376a7555bc7c
# ╟─c0894d52-b2e1-4015-a56b-1f447099dbda
# ╠═47d5c6a3-294f-4fd6-be95-72209103db60
# ╟─c1f717e0-fa07-4cbe-9105-a8745ffd053c
# ╠═53aed70a-bdd0-4c28-8755-31ae96daf413
# ╟─f5c11992-a773-4420-9abd-9681e7ca9eed
# ╟─4779dfeb-c180-4d98-b006-571ca7adadcb
# ╠═ed486e03-f003-4267-851c-f8ce9713dfb2
# ╟─9864b4d7-baef-4798-982d-5a75cb3a2793
# ╠═3e5db56a-9b1c-4e15-b5ad-5ddd2f8c8dd8
# ╟─762968dd-1d23-4ada-8f90-c6a98667efff
# ╟─d2f9823d-0b80-400b-8ddd-bc029a24b987
# ╠═0c996cdc-b54d-4167-b9e9-714ee4911025
# ╟─a3f148b9-911e-4fc7-8de3-2c4ffcf97756
# ╠═61fed22d-8bed-4694-a74f-68b025b45d50
# ╟─053648f0-d079-48f3-bc74-264738693d1f
# ╠═8308ea7e-0f61-4630-ac03-e7c71c94fbd7
# ╟─f0f4d26f-610c-4658-91ae-59a3b589b159
# ╟─6092a639-f4ba-4833-9f2a-b9a513e30bb1
# ╟─6bd7bc71-b827-4f33-8302-83520ba7cd7b
# ╟─083e1f7a-e09e-479a-9cd9-83203e8c23dc
# ╟─f8c35c68-0e60-45d7-ab93-5aace2bd1589
# ╟─db338525-379e-44e5-b3eb-8d98e33d171b
# ╟─9a6793e3-a46a-4736-9521-c5c6a1636aa3
# ╟─17d99340-383a-4b29-a1a2-244024c71155
# ╟─fcf275be-b4c2-4d06-b502-a0c52e74b533
# ╟─9d4a8178-52f4-4846-94d6-74b8d68abbb4
# ╟─dede1f95-42ac-4325-ad20-07f884b58e2e
# ╟─21523a64-b1a9-422d-8dc0-ebb324f35380
# ╟─62e5c394-dad2-4a62-a7ae-d4999f741505
# ╟─6da749bc-0181-49c2-8a7a-9c180930dce6
# ╟─08b20ea7-f3c3-4302-bfc3-a13bbd57930a
# ╟─2cf96d2f-8a37-44e8-8b54-5e2e35542fcc
# ╟─c24c96be-57d6-48dc-b79e-c9a13b9493b3
# ╟─7a952990-b524-43b6-a01b-cf6b6c6eaa14
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
