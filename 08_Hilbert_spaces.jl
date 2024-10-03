### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 46704182-ef55-11ee-31c7-db589297db34
begin
	import TikzPictures.TikzPicture
	using LaTeXStrings
	using PlutoTeachingTools
	using PlutoUI
	using HypertextLiteral

	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ 7d7ac845-c19c-4aff-84ed-8a486eb1bc56
md"""
# Hilbert Spaces
"""

# ╔═╡ 49f2bf4a-8f1a-4122-aad6-8b50db93ac1f
md"""
In our introductory discussion abut quantum mechanics we already discussed the importance of the probability density $|\Psi(x)|^{2}$.
Clearly for this to make sense the wavefunction $\Psi: \mathbb{R}^{n} \rightarrow \mathbb{C}$ needs to be normalized, i.e. the integral
```math
\int \overline{\Psi(x)} \Psi(x) d x
```
needs to be un-infinite.
Such an integral is *not* well-defined for ordinary functions $\mathbb{R}^{n} \rightarrow \mathbb{C}$ as not all functions are square-integrable over $\mathbb{R}^{n}$ (e.g. a non-zero constant is not).
This motivates the study of function spaces equipped with norms and the structure the norms impose.

Typical norms for function spaces are the $L^p(\Omega)$-norms (or just $L^{p}$-norms when $\Omega \subset \mathbb{R}^{n}$ is clear from context) :
```math
\| f \|_{p} \equiv\left(\int_{\Omega}|f|^{p}\right)^{1 / p}.
```
With this motivation in mind we will now study the structure of normed vector spaces and their companions such as Hilbert and Sobolev spaces.

"""

# ╔═╡ 526efb20-8bfd-4060-9171-979c43b92435
md"""

## Completeness

!!! note "Definition (Completeness)"
	A normed vector space $(V , \| \cdot \|)$ is called **complete** if every Cauchy sequence of vectors in $V$ converges to on element in $V$.

!!! tip "Remark (Cauchy sequence)"
	Let us recall that a sequence $\left(x_{n}\right)_{n \in \mathbb{N}}$ of elements $x_{n} \in V$ is Cauchy if
	```math
	\forall \varepsilon>0 \quad \exists M=M(\varepsilon) \in \mathbb{N}: \quad\left\|x_{n}-x_{m}\right\|<\varepsilon \quad \forall n, m \geq M(\varepsilon) .
	```
	What the definition says is that for such sequences there is a unique $x_{*} \in V$ with
	```math
	\forall \varepsilon>0 \quad \exists N=N(\varepsilon) \in \mathbb{N}: \quad\left\|x_{n}-x_{*}\right\|<\varepsilon \quad \forall x_{n} \geq N(\varepsilon) .
	```
	
"""

# ╔═╡ 338c4b80-c710-4e8a-8dcd-e11ff464cbc8
md"""
To make completeness as a concept more clear consider a counter-example. 

- Let's consider the rational numbers $\mathbb{Q}$
  (which can be seen as a vector space over the field $\mathbb{Q}$ itself)
  and ask whether it is complete with respect to $| \cdot |$.

- Consider the sequence
  ```math
  x_{n}=\sum_{k=0}^{n} \frac{1}{k !} \in \mathbb{Q}
  ```
  of rational numbers. 
  It is well known that
  ```math
  \lim _{n \rightarrow \infty} x_{n}=e \notin \mathbb{Q} .
  ```
  Thus, $\mathbb{Q}$ is not complete.

- One may, however, build the **completion** of $\mathbb Q$ by including all possible limit points of all sequences with elements in $\mathbb Q$. This is one way to define the real numbers.

A subtle point about completeness is that it depends on the norm which is used to determine whether a sequence is Cauchy or not. 
In other words, a vector space may be complete with respect to one norm but not another. 
Similarly, the completion of a space with respect to different norms may yield different spaces.

In practice the choice of norm is only important for infinite-dimensional vector spaces as for finite-dimensional real/complex vector spaces all norms are equivalent.

"""

# ╔═╡ 09efc623-3c15-4b91-8c3d-bda48baa8727
md"""
!!! note "Definition (Banach space)"
	A normed vector space which is complete is called a **Banach space**.

!!! warning "Example 1 (Rⁿ)" 
	 $\mathbb{R}^{n}$ with any $p$-norm $\|\cdot\|_{p}$.

!!! warning "Example 2 (Lᵖ)"
	For $1 \leq p<\infty$ the $L^{p}$ - spaces
	```math
		L^{p}(\Omega)=\left\{ f: \Omega \rightarrow \mathbb{C}  \ \middle | \ \int_{\Omega}  |f(x)|^{p} d x<\infty \right\} 
	```
	with $\Omega \subset \mathbb R^d$.

	
"""

# ╔═╡ c2b6dc90-b1ca-4641-9bfe-44d1558574f4
md"""
!!! note "Definition (Hilbert space)"
	A Banach space where the norm is induced by an inner product is a **Hilbert space**.

!!! warning "Example 3 (Square integrable functions)"
	An important Hilbert space is the space of square integrable functions
	```math
		L^{2}\left(\mathbb{R}^{d}\right)=\left\{\psi: \mathbb{R}^{d} \rightarrow \mathbb{C} \ \middle | \ \int_{\mathbb{R}^{d}} | \psi(x) |^{2} d x<\infty\right\}
	```
	with inner product
	```math
		 \langle\psi,  \varphi \rangle_{L^{2}}=\int_{\mathbb{R}^{d}} \overline{\psi(x)} \varphi(x) d x .
	```
"""

# ╔═╡ 1100206d-b44d-4625-b425-f506cc258337
md"""

## Compactness

We will now discuss compactness, a notion which we will meet in various places in the rest of the lecture. 
Here we need it to introduce two important variants of $L^p$ spaces.

!!! note "Definition (Compactness)"
	Let $V$ be a normed vector space. 
	A subset $K$ of $V$ is called compact if every sequence $(x_{n} )_n \subset K$ has a converging subsequence whose limit is an element of $K$.


A colloquial way of stating this is :
!!! tip ""
	On compact sets Bolzano-Weierstrass works.

The extraction of subspaces is central in analysis. 
For example the *extremal value theorem*, which we used to prove the existence of eigenvalues, is based on this argument and can be generalized to compact sets.

!!! tip "Remark (Topological compactness)"
	Equivalently to this definition in metric spaces is a coverage-based notion of compactness ($K \subseteq V$ is compact if every coverage of open sets admits a coverage using only a finite subcollection) due to the Heine-Borel theorem.

!!! tip "Remark (Compactness in finite dimensions)"
	In finite dimensions :
	```math
		\text{Compact } \iff \text{ Closed and bounded}
	```

In infinite dimensions things are trickier as we will see below, pointing as to the fact that our eigenvalue existence proof (Lemma 2.2) will not go through for general operators.
"""

# ╔═╡ 4b412cc4-514c-4f30-b70d-babf8e58979f
md"""
!!! note "Theorem 1"
	Let $V$ be a normed vector space for which every bounded sequence in $V$ contains a convergent subsequence.
	Then, $\dim (V)<\infty$.

> *Proof* by contradiction. 
> Assume $\dim (V)=\infty$. We construct a bounded sequence that has no convergent subsequence:
>
> - Take $x_{1} \in V$ such that $\|x_{1} \|=1$
> - Choose $x_{2} \in V$ such that $\|x_{2} \|=1$ and $\|x_{1}-x_{2}\| \geq 1 / 2$.
> - By the Riesz lemma we can choose $x_{k} \in V$ such that $\left\|x_{k}\right\|=1$ and $\left\|x_{k}-x_{i}\right\| \geq 1 / 2$ $\forall i=1, \dots, k-1 .$
> 
> This constructs a sequence $\left(x_{n}\right)_{n} \subset V$, such that $\left\|x_{n}\right\|=1$ (i.e. it is bounded), but $\left\|x_{n}-x_{m}\right\| \geqslant 1 \quad \forall m \neq n$. 
> As the Cauchy criterion is never satisfied we cannot extract a convergent subsequence.
> $\hspace{11cm} \square$

"""

# ╔═╡ c14b2c6a-b7d7-465d-a591-474c703545d5
md"""
## Weak derivatives & Sobolev spaces

As was the case for matrices, computing Rayleigh quotients
```math
R_{\opH}(\psi)=\frac{\langle\psi,  \opH \psi\rangle_{L^{2}}}{\langle\psi, \psi\rangle_{L^{2}}} \tag{1}
```
for operators $\opH$ as well as the minimisation of these objects will be central when approximating spectra $\sigma(\opH)$.
We must therefore ensure that this quantity is well-defined.

- All operators $\opH$ which we have studied so far and pretty much all Hamiltonians in quantum mechanics involve the Laplace operator (Kinetic energy). 
  However, if $\psi \in L^{2} (\mathbb{R}^{d} )$ there is no guarantee that derivatives of $\psi$ are still in $L^{2} (\mathbb{R}^{d} )$.

- Thus, taking wavefunctions to be $L^{2}$ is *not* sufficient to ensure the $L^{2}$ inner product in (1) even makes sense.
- Therefore, we seek Hilbert spaces with a more restrictive structure, which ensures derivatives to be $L^{2}$ - integrable as well.

Since we are only interested in integrating over derivatives, it turns out that we can allow ourselves a weaker notion of differentiability.
To define this, we first need to introduce two new function spaces :
-  $C_{0}^{\infty}(\Omega)$ is the space of infinitely differentiable functions $\Omega \rightarrow \mathbb{C}$ with compact support.
- The space of locally integrable functions $L_{loc}^p$ :  

!!! note "Definition (Locally integrable functions)"
	For $\Omega \subset \mathbb{R}^{d}$, we define
	```math
		L_{l o c}^{p}(\Omega)=\left \{ 
		f : \Omega \rightarrow \mathbb{C} \ \middle \vert  \  
		f\rvert_K \in L^{p}(K) \quad \forall K \in \Omega, K \text { compact} \right \}
	```

We have $L^{p}(\Omega) \subset L_{loc}^{p}(\Omega)$, since this is a less strict criterion : the growth of the function towards the boundary $\partial \Omega$ - or towards $\infty$ - is not limited.

!!! tip "Remark"
	Note $L^{p}(\Omega) \subset L_{loc}^{p}(\Omega) \subset L_{loc}^{1}(\Omega)$ for $p \geq 1$ and $\Omega \subset \mathbb{R}^{d}$ potentially unbounded, while $L^p(\Omega) \subset L^{1}(\Omega)$ only if $\Omega \subset \mathbb{R}^{d}$ is open and bounded.

With this in mind, we can introduce the weak derivative :
"""

# ╔═╡ 374523ab-decb-42e9-88bf-3420f7552d39
md"""
!!! note "Definition (Weak derivative)"
	The function $u \in L_{loc}^{2}(\Omega)$ with $\Omega \in \mathbb{R}^{d}$ open has a **weak derivative** along the coordinate $x_i$ in $L^{2}(\Omega)$ 
	if there exists a $g_{i} \in L^{2}(\Omega)$ with

	```math
	\begin{align}
		\forall \varphi \in C_{0}^{\infty}(\Omega) 
		&&
		\left\langle g_{i}, \varphi \right\rangle_{L^{2}(\Omega)}
		=
		-\left\langle u, \frac{\partial \varphi}{\partial x_{i}}\right\rangle_{L^{2}(\Omega)}.
	\end{align}
	```
	We then usually employ the same notation as for strong derivatives and e.g. denote $\frac{\partial u}{\partial x_{i}} = g_i$ **in the weak sense**.
	
	Defining for $\alpha= (\alpha_{1}, \ldots, \alpha_{n} ) \in \mathbb{N}_{0}^{n}$ the notation
	```math
	\begin{align}
		\|\alpha\|_{1}=\sum_{i=1}^{n} \alpha_{i} 
		&&
		D ^{\alpha}=\frac{\partial^{\|\alpha\|_{1}}}{\partial x_{1}^{\alpha_{1}} \ldots \partial x_{n}^{\alpha_{n}}} .
	\end{align}
	```
	we analogously define higher weak derivatives $D^{\alpha} u \in L^{2}(\Omega)$ if there exists $g_{\alpha} \in L^{2}(\Omega)$ with
	
	```math
	\begin{align}
		\forall \varphi \in C_{0}^{\infty}(\Omega)
		&&
		\left\langle g_{\alpha}, \varphi \right\rangle_{L^{2}(\Omega)}
		=
		(-1)^{\|\alpha\|_{1}}\left\langle u, D^{\alpha} \varphi \right\rangle_{L^{2}(\Omega)}
	\end{align}
	```
"""

# ╔═╡ 972451e6-4e52-4d60-ac58-744e87d752d5
md"""

The motivation for this definition is partial integration, where the boundary terms vanish since $\varphi$ has compact support.
Explicitly, for $\Omega = (a,b)$, we have
```math
	\int_a^b u(x) \frac{\partial \varphi(x)}{\partial x} dx = \underbrace{\bigg [ u(x) \varphi(x)  \bigg ]^b_a}_{0 \text{ ($\varphi$ has compact support)}} - \int_a^b \frac{\partial u(x)}{\partial x} \varphi(x) dx 
```


If the strong derivative exists, the weak one agrees with it. 
The weak derivative is well-defined.
"""

# ╔═╡ a84d9ce7-b0d9-4c88-a433-dcab67459a11
md"""
This notion allows us to construct new Hilbert spaces :

!!! note "Definition (Sobolev spaces)"
	Let $\Omega \subset \mathbb{R}^{d} $ open. 
	The Sobolev space $H^n(\Omega), n \in \mathbb{N}$
	```math
		H^{n}(\Omega)= \{\psi \in L^{2}(\Omega) \mid D^{\alpha} \psi \in L^{2}(\Omega) \quad \forall \alpha,\|\alpha\|_{1} \leq n  \}
	```
	with inner product
	```math
		\langle f,  g\rangle_{H^{n}}=\sum_{ \|  \alpha \|_{1} \leq n}\left\langle D^{\alpha} f,  D^{\alpha} g\right\rangle_{L^{2}}
	```
	and induced norm
	```math
		\|f\|_{H^{n}}=\sum_{\| \alpha \|_{1} \leq n}\left\|D^{\alpha} f\right\|_{L^{2}}
	```
	is a *Hilbert space.*


!!! tip "Remark"
	The Fourier transform is a useful tool to classify Sobolev spaces
	```math
		f \in H^{n} (\mathbb{R}^{d} ) 
		\iff 
		\int_{\mathbb{R}^{d}}\left(1+|p|^{2 n}\right) |\hat{f}(p) |^{2} d p<\infty
	```
"""

# ╔═╡ c0d4d23b-5b3e-441d-8161-59f8f78fbb52
md"""
## Denseness and separability

Before returning our focus to operators, the final concepts we need to discuss are denseness and separability.

!!! note "Definition (Dense subspaces)"
	A subspace $S \subset V$ is **dense** on $V$ if each vector $x_{*} \in V$ is either (1) also a member of $V$ or (2) one can find a Cauchy sequence $\left(x_{n}\right)_{n} \subset S$ which converges to $x_{*}$.


!!! tip "Remark" 
	 $\mathbb{Q}$ is dense in $\mathbb{R}$.

If $S \subset V$ is dense this means that, using only elements from $S$, we can construct an approximation to any $x_{*} \in V$ to any desired accuracy. 
Clearly, finding good dense subspaces
is key to do numerical computations involving
infinite dimensional spaces like $H^{n}(\Omega)$.
In particular the ones with a countable number of basis functions are useful as this provides a natural way to start with a crude approximation (use few basis functions) and then keep refining until the desired accuracy is reached.


Fortunately,
"""

# ╔═╡ 04821b47-064b-407e-b6a6-aa8439b0f0f6
md"""
!!! note "Definition (Separability)"
	A Hilbert space is separable if it has a dense countable subset.

A Hilbert space $\hilbert$ is separable if and only if it admits a countable orthonormal basis $\{\varphi_{\mu} \}_{\mu=1}^{\infty} \subset \hilbert.$ 
With this each element $\psi \in \hilbert$ can be identified as
```math
	\psi=\sum_{\mu=1}^{\infty} c_{\mu} \varphi_{\mu}
```
to a square-summable infinite sequence $\left(c_{\mu}\right)_{\mu} \in \mathbb{C}$. 
Every infinite-dimensional separable Hilbert space is thus isometrically isomorphic to the sequence space
```math
	\ell^{2}(\mathbb{C})=
	\left\{a: \mathbb N \rightarrow \mathbb{C} \ \middle |  \ \sum_{n=1}^{\infty} | a(n) |^{2}<\infty\right\}.

```

Separability is also closely related to being able to approximate the infinite-dimensional space numerically.
"""

# ╔═╡ a6aafce2-6039-499e-a281-59b0483d8c4b
md"""

!!! warning "Example 4 (Examples of seperable Hilbert spaces)"
	Examples of separable Hilbert spaces :
	-  $L^p(\Omega)$ for $1 \leq p<\infty$ if $\Omega \subset \mathbb{R}^{d}$ is open.
	-  $H^{n}(\Omega)$ as subspaces of $L^{2}(\Omega)$

!!! warning "Example 5 (Classic counterexample)"
	 $L^{\infty}([0,1])$ is *not* separable.
	Consider the family of characteristic functions $f_{t}= \indicator_{[0, t]}$ with $0<t \leq 1$. 
	Clearly,
	```math
		\forall \delta \text{ s.t. } 0<\delta<t \leq 1:\left\|f_{\delta}-f_{t}\right\|_{L^{\infty}}=\| \indicator_{( \delta, t]} \|_{L^{\infty}}=1 \tag{2}
	```
	Let $(g_k)_{k \in \mathbb N}$ be a sequence whose elements form a dense countable subset of $V=L^{\infty}([0,1])$. 
	For such a sequence one can show that for each element $x \in V$ of the vector space and for all $\alpha > 0$ there exists a $k \in \mathbb{N}$ such that $\left\|g_{k}-x\right\|<\alpha$. 
	Therefore, there is a $k \in \mathbb{N}$ for each $0<t \leq 1$ such that
	```math
	\left\|g_{k}-f_{t}\right\|<1 / 2
	```
	Because of (2) this can only be true for a single $t=t(k)$. 
	This results in a surjective map $\mathbb N \supset K \ni k \mapsto \delta(k) \in(0,1]$, which is a contradiction since $(0,1]$ is not countable.
"""

# ╔═╡ d17f5c96-645f-447a-ac3d-24cfeaed2e35
md"""
## Summary of concepts

- **Completeness** : Limits of converging sequences remain in normed vector space/Banach space
- **Compactness** : Bounded sequences of compact normed vector spaces admit a converging subsequence.
- **Compact subspace of Banach space** : Useful in iterative procedures and optimisation. Roughly, "boundedness $\Rightarrow$ convergence", which is key in proving the existence of eigenvalues in some operators.
- **Sobolev spaces** $H^{n}$ : Hilbert space more regular than $L^{2}$, where weak derivatives remain $L^{2}$-integrable.
- **Separability/dense subspaces** : Approximation by finite dimensional subspaces is possible as it admits a countable basis.

"""

# ╔═╡ b0b223a5-4b57-4a9f-a525-b0eb3796cd22


# ╔═╡ 32be4960-6c49-4e73-8172-75786ae42ad3
TableOfContents()

# ╔═╡ b9d75bc0-5e4a-4585-b4c0-327c485ad26f
let
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	
	Sidebar(Markdown.parse(read("sidebar.md", String)), 305)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.58"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "c52e33438f230a8285860e7d80a61412a18560fd"

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

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

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

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

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
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "5d3a5a206297af3868151bb4a2cf27ebce46f16d"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.33"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "0b898aba6cb0b01fb96245fa5375accb651a241a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.0.0"

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

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

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

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "02148a0cb2532f22c0589ceb75c110e168fb3d1f"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "21.9.0+0"

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
git-tree-sha1 = "677b65e17aeb6b4a0be1982e281ec03b0f55155c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.16"

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

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "tectonic_jll"]
git-tree-sha1 = "79e2d29b216ef24a0f4f905532b900dcf529aa06"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.0"

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

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.tectonic_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "Graphite2_jll", "HarfBuzz_ICU_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "54867b00af20c70b52a1f9c00043864d8b926a21"
uuid = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"
version = "0.13.1+0"
"""

# ╔═╡ Cell order:
# ╟─7d7ac845-c19c-4aff-84ed-8a486eb1bc56
# ╟─46704182-ef55-11ee-31c7-db589297db34
# ╟─49f2bf4a-8f1a-4122-aad6-8b50db93ac1f
# ╟─526efb20-8bfd-4060-9171-979c43b92435
# ╟─338c4b80-c710-4e8a-8dcd-e11ff464cbc8
# ╟─09efc623-3c15-4b91-8c3d-bda48baa8727
# ╟─c2b6dc90-b1ca-4641-9bfe-44d1558574f4
# ╟─1100206d-b44d-4625-b425-f506cc258337
# ╟─4b412cc4-514c-4f30-b70d-babf8e58979f
# ╟─c14b2c6a-b7d7-465d-a591-474c703545d5
# ╟─374523ab-decb-42e9-88bf-3420f7552d39
# ╟─972451e6-4e52-4d60-ac58-744e87d752d5
# ╟─a84d9ce7-b0d9-4c88-a433-dcab67459a11
# ╟─c0d4d23b-5b3e-441d-8161-59f8f78fbb52
# ╟─04821b47-064b-407e-b6a6-aa8439b0f0f6
# ╟─a6aafce2-6039-499e-a281-59b0483d8c4b
# ╟─d17f5c96-645f-447a-ac3d-24cfeaed2e35
# ╟─b0b223a5-4b57-4a9f-a525-b0eb3796cd22
# ╟─32be4960-6c49-4e73-8172-75786ae42ad3
# ╟─b9d75bc0-5e4a-4585-b4c0-327c485ad26f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
