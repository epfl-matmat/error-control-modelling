### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 52b4d68c-d6f1-11ee-18b4-db386a5995ad
begin
	using HypertextLiteral
	using PlutoUI
	using PlutoTeachingTools
	
	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ b944f68e-85a6-4613-8c5e-0fe3d47ff0d4
md"""
# Introduction
"""

# ╔═╡ 53e49738-8673-4ecd-9e86-b1431200b132
md"""
## Opening remarks

```math
\text{\textit{All models are wrong, only some models are useful.}}
```

**[Slides of the presentation](/)**: These summarise the following paragraphs and the content of the course in a short presentation.
$\textcolor{red}{\text{todo: add slide link here}}$
"""

# ╔═╡ aa175689-f915-42a1-b472-19deb5811f46
md"""
Constructing suitable **approximations to model the world** around us
sits at the **heart of all sciences**. As the above famous quote implies,
modelling a scientific problem inevitable comes with the introduction
of errors. Even the most accurate models of physics, e.g. the standard model
of particle physics or relativistic quantum mechanics have their limitations.
Moreover these models are typically too involved to be applied to
all but the simplest systems. In particular for "real" systems of macroscopic
size additional approximations.
Typically this not only involves the selection of **more approximate**,
thus simpler,
**physical models**, but also employing **numerical** instead of an analytical techniques
for solving the underlying equations.
"""

# ╔═╡ 7229e2b6-cf77-4b32-b52d-1fc7d218b007
md"""
Which physical effects can be neglected and which numerical techniques are most
appropriate is often application-specific and **many alternatives exist**.
Understanding the **dominant sources of error** in a simulation is thus
a crucial component in understanding whether each respective **approximation choice** has been appropriate and --- most importantly --- which of them **needs to be refined** in order to improve the result.
"""

# ╔═╡ bd045f8e-7baa-4fbf-94bc-22e39cac9a0d
md"""
The goal of this lecture is to discuss mathematical techniques
which aid with the estimation of error in scientific simulations.
Within the lecture we will distinguish various forms of error and consider (some) of them in separate discussions throughout the semester. For example we will distinguish the **model error**, i.e. the error introduced the selected physical model compared to the ground truth,
as well as the **numerical error**,
i.e. all additional error from the numerical procedure employed for
solving the chosen physical model.
Of particular focus in the lecture will be understanding the various aspects
of the numerical error,
that is in particular the errors due to the chosen **discretisation basis**,
the **algorithms** to solve the numerical problem as well as the chosen
**floating-point arithmetic**.
"""

# ╔═╡ cd9950c8-1b3d-4a1b-8606-9fcf4f39599f
md"""
Looking back into the **history of scientific modelling**
a number of incidents with **hundreds of millions of dollars
of damage** could have been **prevented** if better numerical error control
measures would have been implemented.
One example is the 
[sinking of the Sleipner A oil rig](https://en.wikipedia.org/wiki/Sleipner_A)
platform
due to an insufficient finite element discretisation.
A second example is the
[self-destruction of an Ariane 5 rocket](https://en.wikipedia.org/wiki/Ariane_5)
in its test flight
due to an error in the performed floating-point arithmetic.
In this lecture we will discuss techniques how such errors could have been prevented.
Some guiding questions for our lecture are:
- What can go wrong in a simulation?
- How can we detect this and what can we do about it?
- If we know errors how can this be used to make simulations faster and more reliable?
"""

# ╔═╡ d6db19e2-78fb-4619-b2c7-c2784461c583
md"""
## Eigenvalue problems in scientific modelling

Eigenvalue problems occur naturally in many domains of science.
In physics and engineering problems they are frequently related to vibrations,
e.g. bridges or buildings can swing under wind load.
A famous example is the Tacoma narrows bridge, which collapsed in 1940 just half
a year after its opening due to strong winds bringing the bridge into strong vibrations.
Similarly the London millennium foot bridge had to be closed due to its strong wobbeling
by the passing pedestrians.

In my own field, atomistic modelling, the underlying workhorse theory is the Schrödinger equation.
In quantum mechanics the eigenpairs (states) of an operator describing the physical system under
study form the key quantity of interest in simulations.

Moreover, as we will discuss in the next examples in more detail,
there is a natural connections from solving partial
differential equations~(PDE) to eigenvalue problems involving the PDE operator.

Another natural connection arises when analysing the convergence of
iterative numerical schemes,
such as multidimensional fixed-point problems,
iterative schemes like Newton or the conjugate-gradient algorithm.
The convergence of CG, for example,
is related to the condition number of the underlying linear system,
which in turn are related to singular values of the system matrix
--- an Quantity, which can be obtained by solving an eigenvalue problem.
"""

# ╔═╡ 189b99a5-776f-4ad6-b05a-66b9801a97c4
md"""
## The connection between linear and eigenvalue problems : the heat equation

Let $\Omega \subset \mathbb R^n$ with a Lipshitz (smooth) boundary.
We seek $u : \Omega \times \mathbb R^+ → \mathbb R$ such that 
```math
\begin{align}
\frac{∂ u}{∂ t} - \laplacian u &= 0 && \text{in } \Omega \times \mathbb  R^+ \tag{Heat equation}
\\
u(\cdot,0) &= u_0 && \text{in } \Omega \tag{Initial conditions}
\\
u &= 0 && \text{on } ∂ \Omega \times \mathbb R^+ \tag{Dirichlet b.c.}
\end{align}
```
Where the initial conditions $u_0 : \Omega → \mathbb R$ are of appropriate regularity.

"""

# ╔═╡ e52eb9de-c5b3-424e-b2d5-e1b926f76fd4
md"""
To solve this problem, we perform separation of variables
```math
u(x,t) = \sum_{n = 1}^{\infty} c_n(t) v_n(x).
```
We will see later that this postulate is justified.

To start, let's just consider one term of the sum, i.e. $u(x,t) = c(t) v(x)$. 

- From the heat equation
  ```math
  	c'(t) v(x) - c(t) \laplacian v(x) = 0
  ```

- Assuming that $v(x)$ and $c(t)$ are non-zero, we can divide
  ```math
  \begin{align}
  \frac{c'(t)}{c(t)} = \frac{\laplacian v(x)}{v(x)} && ∀ \ t>0, x ∈ \Omega.
  \end{align}
  ```
- The left only depends on $t$, and the right only on $x$.
  As a result, we need
  ```math
  	\frac{c'(t)}{c(t)} = \frac{\laplacian v(x)}{v(x)}  = μ \in \mathbb R
  ```
  for some $μ$.
"""

# ╔═╡ d5581136-7b5e-4987-a67d-241286d424c8
md"""
- Consider the time-dependant part ($c$)
  ```math
  \begin{align}
    c'(t) &= μ c(t)
  \\
  	\Rightarrow c(t) &= α \exp(μ t).
  \end{align}
  ```
- Separately, consider the space-dependant part ($v$)
  ```math
  \begin{align}
  \laplacian v &= μ v && \text{in } \Omega
  \\
  v &= 0 && \text{on } ∂\Omega
  \end{align}
  ```
  *which is the Dirchlet-Laplace eigenvalue problem* !
"""

# ╔═╡ 4db1319a-0b7c-42d5-ba03-00f106542c80
md"""
- One can show that the operator $- \laplacian$ with homogeneous Dirichlet boundary conditions has a sequence of eigenvalues $0 < λ_1 ≤ λ_2 ≤ … → ∞$ and orthonormal eigenfunctions $\{ \phi_n \}_{n=1}^∞$, i.e.
  ```math
  - \laplacian \phi_n = λ_n \phi_n
  ```
  with
  ```math
  \langle \phi_n, \phi_m \rangle = ∫_\Omega \phi_n(x) \phi_m(x) dx = δ_{mn}.
  ```

- These form a complete basis of $H^1_0 (\Omega)$, the relevant Hilbert space for this problem.

This justifies the ansatz
```math
u(x,t) = \sum_{n=1}^∞ c_n(t) \phi_n(x) 
```
we made before, where $\phi_n$ are now eigenfunctions of the Laplace operator.
"""

# ╔═╡ 17dd79e4-c998-48cc-854e-34e20f3338f9
md"""
To fully translate the problem to the eigenfunction setting, we need to consider the initial condition $u_0.$

- On the one hand, in the eigenbasis,
  ```math
  u_0 = \sum_{i=1}^∞ \underbrace{\langle \phi_i, u_0 \rangle}_{u_{0,i}} \phi_i.
  ```

- On the other hand,
  ```math
  c_i(t) = α_i \exp(- λ_i t).
  ```
  which yields
  ```math
  u(x,t) = \sum_{i=1}^∞ α_i \exp(- λ_i t) \phi_i(x).
  ```

- Taking $t → 0$, we obtain
  ```math
  \sum_{i=1}^\infty u_{0,i} \phi_i(x) = u_0(x) = u(0,x) = \sum_{i=1}^\infty \alpha_i \phi_i(x).
  ```

- Since $\phi_i$ is a basis, we deduce $\alpha_i = u_{0,i}$ and finally obtain
  ```math
  u(x,t) = \sum_{i=1}^\infty u_{0,i} \exp(- \lambda_i t) \phi_i(x)
  ```

We observe the following :
- The solution decays to $u_\infty = 0$ as $t \to \infty$.
- The decay rate depends on $\lambda_i$.
"""

# ╔═╡ 9f07f419-fbe9-4ff5-8cbb-2bd056868b61
md"""
More importantly for this course, this shows that eigenvalue problems, which we will focus on, are key in understanding linear problems, even though they might seem unrelated at a first glance.
One of these linear problems of particular relevance is solving Schrödinger equation, which we treat in the next section.
"""

# ╔═╡ b102f351-f5b1-47c0-8737-f813b5b745f6
md"""
## The Schrödinger equation and quantum mechanics

At the microscopic level, the physics and chemistry of materials is governed by the interaction of electrons and nuclei. 
At this scale, the regime of quantum mechanics applies.

In quantum mechanics the state of a system is described the complex-valued square-integrable wave function $\Psi : (x,t) ↦ \Psi(x,t) \in \mathbb C$, where $x$ corresponds to the degrees of freedom (e.g. position or spin) of the particles (e.g. electrons) in the system.

For simplicity, take $x \in \mathbb R^3$, i.e. one particle in 3D.
Part of the meaning of $\Psi$ appears through its modulus squared. 
Indeed, for all $t$, $|\Psi(x,t)|^2$ corresponds to the probability distribution of finding the particle at position $x$ at time $t$.

Further information appears through its Fourier transform, as $|\hat \Psi (p,t)|^2$ corresponds to the probability distribution of finding the particle at momentum $p$ at time $t$, where the Fourier transform is given by
```math
\hat \Psi (p,t) = \int_{\mathbb R^3} \Psi (x,t) e^{-2 \pi i \ p \cdot x}dx.
```
For those already familiar with quantum mechanics, note that we set $\hbar = 1$.
"""

# ╔═╡ e1132fa4-c9ff-487e-b716-0b3157c2cdd0
md"""
In analogy to classical mechanics, one may obtain the dynamics of the particle by investigating its total energy (considering the mass of the particle to be 1).

- The total energy is given by
  ```math
  \begin{align}
    E(t) &= \text{kinetic energy} + \text{potential energy}
    \\
    &= \int_{\mathbb R^3} \frac{|p|^2}{2} | \hat \Psi (p,t) |^2 dp + \int_{\mathbb   R^3} V(x) | \Psi (x,t) |^2 dx
  \end{align}
  ```
  where $V : \mathbb R^2 \to \mathbb R$ is the potential inducing the particle dynamics.

- From the Fourier identity $p \hat \Psi (p) = -i \widehat{(\nabla \Psi)}(p)$, we can develop the first term into
  ```math
  \begin{align}
  	\int_{\mathbb R^3} |p|^2 | \hat \Psi (p,t) |^2 dp &= \int_{\mathbb R^3} |p \hat \Psi (p,t)|^2 dp
  \\
  	&= \int_{\mathbb R^3} \lvert -i \widehat{(\nabla \Psi)}(p,t)|^2 dp
  	\\
  	&= \int_{\mathbb R^3} \lvert  \widehat{(\nabla \Psi)}(p,t)|^2 dp
  	\\
  	&\hspace{-0.9em} \stackrel{\text{Parceval}}{=} \int_{\mathbb R^3} | \nabla \Psi(x,t) |^2 dx
  \end{align}
  ```
  where we used Parceval's theorem to go from the third to the fourth line.
"""

# ╔═╡ 33c563af-3436-4275-becd-5333d480b4f4
md"""
- With this we obtain,
  ```math
  \begin{align}
  E(t) &= \frac1{2} \int_{\mathbb R^3} | \nabla \Psi(x,t) |^2 dx + \int_{\mathbb R^3}   V(x) |\Psi(x,t)|^2 dx
  \\
  & \stackrel{IBP}{=} - \frac1{2} \int_{\mathbb R^3} \overline{\Psi(x,t)} \laplacian    \Psi(x,t)  dx + \int_{\mathbb R^3} \overline{\Psi(x,t)} V(x) \Psi(x,t) dx
  \\
  &= \int_{\mathbb R^3} \overline{\Psi(x,t)} \ \opH \ \Psi(x,t) dx
  \\
  &= \langle \Psi, \opH \Psi \rangle_{L^2(\mathbb R^3)}
  \end{align}
  ```
  Where we used integration by parts to go from the first to the second line. 
  
- In this, the *Hamiltonian*
  ```math
    \opH \coloneqq - \frac1{2} \laplacian + V(x)
  ```
  gives the total energy of the system when applied to the wavefunction $\Psi$, where
  -  $- \laplacian / 2$ is the kinetic energy operator
  -  $V(x)$ is the potential energy, applied point-wise in real space.
"""

# ╔═╡ 8e7daba1-277b-4632-b995-5be76616d0d3
md"""
The dynamics of the system are now given by the time-dependant Schrödinger equation (TDSE) 
```math
i \frac{\partial \Psi}{\partial t} = \opH \Psi \tag{TDSE}
```
In general, solving this equation is quite complicated.
However, since it is linear, we know that superpositions of multiple solutions are also a solution. 
This motivates a separation of variables similar to the heat equation.

- Inserting $\Psi(x,t) = c(t) \varphi(x)$ into the TDSE yields 
  ```math
  i \frac{c'(t)}{c(t)} = \frac{\opH \varphi(x)}{\varphi(x)} = E = \text{cst}
  ```

- For the time-dependence we obtain the ODE
  ```math
  \begin{align}
  c'(t) &= - i E c(t)
  \\
  \Rightarrow c(t) &= \exp(- i Et)
  \end{align}
  ```
- The space dependence yields the time-independent Schrödinger equation (TISE)
  ```math
  \opH \varphi(x) = E \varphi(x) \tag{TISE}
  ```
"""

# ╔═╡ 94c2368c-243e-4e29-ba9b-a6c9253fa74f
md"""
- As before, if we are able to find (orthonormal) eigenpairs $(E_n,\phi_n)$, then the TDSE can be solved as
  ```math
  \Psi(x,t) = \sum_{n =1}^\infty c_n e^{- i E_n t} \phi_n(x)
  ```
  with $c_n = \langle \phi_n , \Psi \rangle$.

- Of key importance in quantum mechanics in particular is the computation of the eigenpair corresponding to the smallest eigenvalue of $\opH$, called the *ground state*.
"""

# ╔═╡ d26799e2-bd9a-4ff5-89b2-7d4a8cc61b74
TableOfContents()

# ╔═╡ 8993ae47-1173-46c1-8f1e-55412215b61f
let
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	
	Sidebar(Markdown.parse(read("sidebar.md", String)), 265)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
HypertextLiteral = "~0.9.5"
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "72c07b07bb225d1401f37584678b084190a9d3b1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "c0216e792f518b39b22212127d4a84dc31e4e386"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.5"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

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

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "a6adc2dcfe4187c40dc7c2c9d2128e326360e90a"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.32"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "5b0d630f3020b82c0775a51d05895852f8506f50"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.4"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "eeaedcf337f33c039f9f3a209a8db992deefd7e9"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.4.8"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "85ddd93ea15dcd8493400600e09104a9e94bb18d"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.15"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─b944f68e-85a6-4613-8c5e-0fe3d47ff0d4
# ╟─52b4d68c-d6f1-11ee-18b4-db386a5995ad
# ╟─53e49738-8673-4ecd-9e86-b1431200b132
# ╟─aa175689-f915-42a1-b472-19deb5811f46
# ╟─7229e2b6-cf77-4b32-b52d-1fc7d218b007
# ╟─bd045f8e-7baa-4fbf-94bc-22e39cac9a0d
# ╟─cd9950c8-1b3d-4a1b-8606-9fcf4f39599f
# ╟─d6db19e2-78fb-4619-b2c7-c2784461c583
# ╟─189b99a5-776f-4ad6-b05a-66b9801a97c4
# ╟─e52eb9de-c5b3-424e-b2d5-e1b926f76fd4
# ╟─d5581136-7b5e-4987-a67d-241286d424c8
# ╟─4db1319a-0b7c-42d5-ba03-00f106542c80
# ╟─17dd79e4-c998-48cc-854e-34e20f3338f9
# ╟─9f07f419-fbe9-4ff5-8cbb-2bd056868b61
# ╟─b102f351-f5b1-47c0-8737-f813b5b745f6
# ╟─e1132fa4-c9ff-487e-b716-0b3157c2cdd0
# ╟─33c563af-3436-4275-becd-5333d480b4f4
# ╟─8e7daba1-277b-4632-b995-5be76616d0d3
# ╟─94c2368c-243e-4e29-ba9b-a6c9253fa74f
# ╟─d26799e2-bd9a-4ff5-89b2-7d4a8cc61b74
# ╟─8993ae47-1173-46c1-8f1e-55412215b61f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
