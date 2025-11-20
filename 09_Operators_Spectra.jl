### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ e80f5b7f-dfa0-4785-86c9-fdc6c387f3ca
begin
	import TikzPictures.TikzPicture
	using HypertextLiteral
	using PlutoUI
	using PlutoTeachingTools
	using LaTeXStrings

	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ fd7a9505-fa0b-4f52-ba35-e476f093a7ed
begin
	using QuadGK
	using LinearAlgebra
	
	# basis functions and their Laplacian
	basis_function(x,n) = sin(x * n) * sqrt(2/π)
	laplacian(x,n) = - n^2 * sin(x * n) * sqrt(2/π)

	# potential
	function gaussian(x)
		A = -1000
		sigma = π/16
		return A * exp(-((x - π/2)/sigma)^2)
	end

	# matrix elements
	function Mij(i,j,V)
		# This is definitely not the most efficient approach from a computational standpoint, but it works
		quad = quadgk(x -> conj(basis_function(x,i)) * (- laplacian(x,j) / 2 + V(x) * basis_function(x,j)), 0, π, atol=1e-6)
		return quad[1] 
	end

	# reference values
	eig = let
		n = 256
		M = zeros(ComplexF64,n,n)
			for i in 1:n
				for j in i:n
					M[i,j] = Mij(i,j,gaussian)
				end
			end
		M = Hermitian(M)
		eig = eigen(M)
	end

	# varying basis size
	n_range = 2:2:80
	eigenvalues = ones(Union{Missing,Float64},maximum(n_range),length(n_range))
	eigenvalues .= missing
	for (i,n) in enumerate(n_range)
		M = zeros(ComplexF64,n,n)
		for i in 1:n
			for j in i:n
				M[i,j] = Mij(i,j,gaussian)
			end
		end
		M = Hermitian(M)
		eigenvalues[1:n,i] = eigen(M).values
	end
end

# ╔═╡ 1f6efd53-2429-46e7-a579-fedd920250ff
md"""
# Operators and their spectra
"""

# ╔═╡ 041e2cfd-0ec9-4766-98dd-7a90b2a4da84
md"""
Previously we discussed the essential properties of function spaces, which can be seen as infinite-dimensional generalization of Euclidean space where functions play the role of vectors.

We noted a number of desirable properties (compactness, completeness, separability), which in infinite dimensions are considerably harder to obtain or verify then in finite dimensions.

Along these lines, linear operators, i.e. objects mapping functions to functions, are the natural generalization of matrices, where again the infinite-dimensional setting leads to a number of subtle and surprising differences in the context of computing approximate spactra or eigenpairs, of which we want to sketch a few aspects on a high level.

First we formally define:
"""

# ╔═╡ 897dfb72-71d7-4ee3-a874-8f932fffd261
md"""
!!! note "Definition (Operator)"
	Let $V, W$ be $\mathbb{C}$ vector spaces. 
	A map $\opA : V \rightarrow W$ is called **linear operator** or just **operator** if
	```math
	\forall x, y \in V \quad \forall \alpha, \beta \in \mathbb{C} \quad \opA (\alpha x+\beta y)=\alpha \opA (x)+\beta \opA (y)
	```
	
	If $W=\mathbb C$ we can also call the operator a **linear functional**.
"""

# ╔═╡ d3ad8a0a-fa8a-4864-ada3-fe237e65c088
Foldable("Examples of operators", md"""
!!! warning ""
	- The zero operator $\opA : V \rightarrow W: x \mapsto 0_{W}$
	- The identity operator $\operatorname{id} : V \rightarrow V: x \mapsto x$
	- Matrices, e.g.
	  ```math
	  A=\left(\begin{array}{lll}
	  1 & 2 & 3 \\
	  4 & 5 & 6
	  \end{array}\right)
	  ```
	  defines the operator $\opA : \mathbb{R}^{3} \rightarrow \mathbb{R}^{2}$ :
	  ```math
	  \left(\begin{array}{l}
	  x_{1} \\
	  x_{2} \\
	  x_{3}
	  \end{array}\right)=x \mapsto A x=\left(\begin{array}{l}
	  x_{1}+2 x_{2}+3 x_{3} \\
	  4 x_1 +5 x_{2}+6 x_{3}
	  \end{array}\right)
	  ```
	
	- Differentiation: with $V=C^{1}([a, b])$ and $W=C^{0}([a, b])$, $\opA (f)=f^{\prime}$ is a linear operator.

""")

# ╔═╡ a8ec0c38-3f60-4713-a1f7-9a8e054a2533
md"""
Guided by our discussions in the matrix setting, a natural idea is to seek approximate eigenpairs of $\mathcal{A}$ following the Courant-Fisher theorem,
that is following the idea of **projection methods**,
which we introduced in the context of numerically computing eigenpairs.

Roughly speaking our goal is to find finite-dimensional spaces $S$,
such that applying Courant-Fisher to the projection of $\mathcal{A}$ into these spaces yields an approximation of the first few eigenpairs.
Or in other words we would like the approximations
```math 
		\mu_k (\opA) = \inf_{\substack{s \subset S \\ \dim (s) = k}} \max_{0 \neq \varphi \in s} R_{\mathcal{A}}(\varphi)
		= \inf_{\substack{s \subset S \\ \dim (s) = k}} \max_{0 \neq \varphi \in s} \frac{\langle φ | \mathcal{A} φ \rangle_\hilbert}{\|φ\|_\hilbert}
```
to converge to the smallest members of the spectrum $\sigma(\mathcal{A})$
(see [below](#Spectra-of-operators) a formal definition of spectrum for operators). 
"""

# ╔═╡ c61e470f-cdf7-4c44-bdf2-c16cad11856c
md"""
Under some conditions the $\mu_k(\opA)$ indeed converge to the smallest elements of $\sigma(\mathcal{A})$, but there are also counterxamples.
For example, for $\opA = - \Delta$ defined on all of $\mathbb{R}^3)$
(our [Dirichlet Laplacian example](https://teaching.matmat.org/error-control/02_Laplace_error_sources.html)), *all* $\mu_k \to 0$ as we make $S$ appropriately larger. In particular even in the limit of perfect discretisation we would thus only approximate the single element $0 \in \sigma(-\Delta)$ of the spectrum. We will provide some ideas why this is.
"""

# ╔═╡ b86c0bf2-cea9-4f6d-89bd-30634e0a2523
md"""
A subtle point in the definition of operators is that the **choice of the operator domain** $V$
(for a given "recipe" to evaluate the operator such as $\opA = - \Delta$)
has a crucial influence on the obtained spectrum $\sigma(\mathcal{A})$.
For example even in the case of physical operators such as
the Schrödinger operators $\opH = - Δ / 2 + v$ with smooth potential $v$,
taking the domain *too small* (e.g. $V = C^\infty_0$)
leads to $\sigma(\mathcal{A}) = \mathbb{C}$.
This is in contrast to the physical regime where
the spectra of such operators should be entirely real.
From the modelling perspective this means that if we target
the wrong domain $V$ with our approximation spaces $S$,
then we may obtain numerical results that are unphysical
and/or converge to an unphysical solution when taking
the limit towards a perfect discretisation.

"""

# ╔═╡ b64cdb0e-3f9a-48d7-a9c6-afebf057e6ef
md"""
While it is out of scope to discuss all details and subtleties
about the spectral approximation of operators in this course,
this chapter aims to provide some high-level overview of the
differences between matrices and operators alongside optional,
more in-depth discussions in a long appendix.

We also want to remark that some properties from the matrix
world do indeed carry forward to operators. A non-exhaustive list
is given in the foldable below:
"""

# ╔═╡ a18579a4-bf14-42ed-b16c-b01f8a8e88de
Foldable("Properties of matrices that hold for operators", md"""
- Zeros are mapped onto each other, since
  ```math
  0_{W}=0_{\mathbb{C}} \opA(x)=\opA\left(0_{\mathbb{C}} x\right)=\opA   \left(0_{V}\right).
  ```

- For a subspace $\tilde{V} \leq V$ we have that $\opA (\tilde{V}) \leq W$ is a subspace of $W$ since, for all $\opA(x), \opA(y) \in \opA(\tilde{V})$, we have
  ```math
  \alpha \opA(x)+\beta \opA(y)= \opA \underbrace{(\alpha x+\beta y)}_{\in   \tilde{V}} \in A(\tilde{V})
  ```
  which proves closure.

- By a similar argument one proves that linear operators conserve linear combinations
  ```math
  \opA \left(\sum_{i=1}^{n} \alpha_{i} x_{i}\right)=\sum_{i=1}^{n} \alpha_{i}  \opA \left(x_{i}\right)
  ```
  i.e. 
  ```math
  \opA  (\operatorname{span} \{x_{1},\dots,x_{n} \} )=\operatorname{span} ( \{\opA (x_{1} ), \dots, \opA(x_{n} ) \} ).
  ```

- Concatenations of two linear operators $\opA_{1}, \opA_{2}$ remain linear :
  ```math
  \begin{align}
   (\opA_{1} \circ \opA_{2} )(\alpha x+\beta y) & =\opA_{1} (\opA_{2}(\alpha x+\beta y) ) \\
  & = \opA_{1} (\alpha \opA_{2}(x)+\beta \opA_{2}(y) ) \\
  & =\alpha (\opA_1 \circ \opA_{2} )(x)+\beta (\opA_{1} \circ \opA_{2} )(y)
  \end{align}
  ```
- Kernel and image are subspaces :
  ```math
  \begin{align}
  \ker(A) &= \{x \in V: \opA (x)=0\} \leq V \\
  \im (A)  &= \{y \in W: \exists x \in V \text{ s.t. } \opA (x)=y\} \leq W
  \end{align}
  ```
""")

# ╔═╡ 678cfb59-33b7-46ea-bd66-09ebb20e09f1
md"""
## Domains of operators

Picking up on the subtle point of operator domains, let's first establish some standard choices.
"""

# ╔═╡ 5f979847-44af-4457-9220-6840b3b1bd11
md"""
In the following we will consider the specific case of linear operators on a Hilbert space $\hilbert$, i.e. operators of the form $\opA: D(\opA) \rightarrow \hilbert$. The **domain** $D(\opA) \leq \hilbert$ is a dense subspace of $\hilbert$. Unless otherwise noted the standard choice of the domain is
```math
	D(\mathcal{A})=\{f \in \hilbert \mid \opA f \in  \hilbert \},
```
i.e. the largest possible subspace of $\hilbert$, such that applying the operator does not take us out of $\hilbert$. Unless otherwise noted we will always employ this domain.

- Other choices of the domain are well possible and sometimes useful. Moreover and perhaps surprisingly it has a decisive influence on the properties of an operator. For example it is related to the operator being self-adjoint or not --- which is in turn related to the question whether the operator spectrum is physical or not, as we will discuss. 

- To illustrate the point of allowing a domain to be smaller than or different from the Hilbert space, we consider the Schrödinger operator of an isolated system, which has the form $\opH=-\Delta+ V$, where $V$ is a potential. 
  The natural Hilbert space for this setting is $L^{2}(\mathbb{R}^{d})$. 
  Thus to ensure $\opH f \in L^{2} (\mathbb{R}^{d})$ $\forall f \in D(\opH)$, we need to take $D(\opH)=H^{2} (\mathbb{R}^{d} )$, one of the Sobolev spaces we introduced last time. 
  Clearly $D(\opH)= \hilbert =L^{2}\left(\mathbb{R}^{d}\right)$ is *not* reasonable as here the Laplacian reduces the regularity.
"""

# ╔═╡ 15012412-d229-45e4-8b12-0dc89d3baaa2
md"""
Similar to matrices we can also define norms of operators making reference to the norm of $\hilbert$:

!!! note "Definition (Operator norm)"
	Let $\left(V,\|\cdot\|_{V}\right)$ and $\left(W,\|\cdot\|_{W}\right)$ denote two Banach spaces. The corresponding operator norm is
	```math
	\|\opA \|_{\boundedoperators (V, W)} \coloneqq \sup _{0 \neq x \in V} \frac{\| \opA (v)\|_{W}}{\| v \|_{V}}.
	```
	When dealing with an operator $\opA: D(\opA) \rightarrow \hilbert$
	and when the underlying Hilbert space is clear, we frequently denote this as
	```math
	\| \opA \|_\text{op} \equiv \sup _{0 \neq f \in D(\opA)} \frac{\| \opA f\|_{\hilbert}}{\|f\|_{\hilbert}},
	```
	i.e. the same setup as the matrix norm, but with the Hilbert space norm and domain inserted.
"""

# ╔═╡ 99968d73-d35d-4d8f-b851-1b468c865a31
md"""
## Bounded operators and continuity

To define the spectrum $\sigma(\opA)$ in infinite dimensions we need one further
property for operators:
"""

# ╔═╡ b44590f7-adcd-46d8-8889-71aa0f774116
md"""
!!! note "Definition (bounded operators)"
	An operator $\opA : V \to W$ is called bounded if $\| \opA \|_{\boundedoperators (V, W)} < \infty$.
	The set of all bounded operators is a Banach space.
"""

# ╔═╡ 9e4cdc94-8ffa-42de-b541-c8beddec8e50
Foldable("Choice of norm and boundedness (optional)",
md"""
Similar to the case of completeness, the choice of norm is essential for boundedness.
To illustrate this, we consider the boundedness of the derivative operator.

- Consider the operator $\partial: X \rightarrow Y$ mapping between the normed spaces $X=C^{1}\left([0,1], \| \cdot \|_{\infty}\right)$ and $Y=C^{0}\left([0,1],\|\cdot\|_{\infty}\right)$, i.e. where the supremum norm is used to measure distances. 
  The action of $\partial$ is $\partial f=f^{\prime}$ for $f \in X$. 
  In this case, $\partial$ is unbounded. 
  To see that consider the sequence $\left\{f_{n}(x)\right\}_{n=1}^{\infty} \subset X$ with $f_{n}(x)=\sin(2 \pi n x)$. 
  Clearly
  ```math
  \begin{align}
  	\sup _{0 \neq f \in C^{1}([0,1])} \frac{\|\partial f\|_{\infty}}{\|f\|_{\infty}} 
  	& \geq \sup_{n=1, \ldots, \infty} \frac{\left\|\partial f_{n}\right\|_{\infty}}{\left\|f_{n}\right\|_{\infty}} 
  	\\
  	& = \sup_{n=1, \dots, \infty} \frac{2 \pi n}{1}=\infty,
  \end{align}
  ```
  which implies the unboundedness of $\partial$.

- Now we consider $\partial: \tilde{X} \rightarrow Y$ with $\tilde{X}=C^{1}( [ 0,1], \| \cdot \|_{C^1})$ and $Y=C^{0}\left([0,1],\|\cdot\|_{\infty}\right)$, where $\tilde{X}$ is now equipped with the stronger norm
  ```math
  	\|f\|_{C^{1}}=\|f\|_{\infty}+\left\|f^{\prime}\right\|_{\infty} .
  ```
  Then $\partial$ is bounded, since
  ```math
  \begin{align}
  	\sup_{0 \neq f \in C^1([0,1])} \frac{\| \partial f \|_\infty}{\| f \|_{C^1}}
  	&= \sup_{0 \neq f \in C^1([0,1])} \frac{\| f' \| _\infty}{\| f \|_\infty + \| f' \|_\infty}
  	\\
  	& \leq \sup _{0 \neq f \in \tilde{x}} \frac{\left\|f^{\prime}\right\|_{\infty}}{\left\|f^{\prime}\right\|_{\infty}}=1 .
  \end{align}
  ```
""")

# ╔═╡ 321383d4-7b51-4b39-91cc-97c87c465e1a
md"""
For linear operators boundedness turns out to be equivalent to demanding them to be continuous, a key result in functional analysis. Details can be found in the foldable below.
"""

# ╔═╡ 4ea66162-8da9-4da9-8cae-14322da43c98
Foldable("Theorem 1: Boundedness ⇔ Continuous for linear operators (optional)",
md"""
	 
!!! note "Definition (Continuity)"
	Let $\left(V,\|\cdot\|_{V}\right)$ and $\left(W,\|\cdot\|_{W}\right)$ be two Banach spaces. 
	An operator $\opA: V \rightarrow W$ is called **continuous** if for every sequence $\{x_{n} \}_{n} \subset V$ with $x_{n} \rightarrow x \in V$ as $n \rightarrow \infty$ it holds
	```math
	\opA x_{n} \rightarrow \opA x \quad \text { as } \quad n \rightarrow \infty .
	```

!!! note "Theorem 1 (Boundedness-continuity equivalence)"
	Let $\left(V,\| \cdot \|_{V}\right)$ and $\left(W,\|\cdot\|_{W}\right)$ be two Banach spaces. 
	A linear operator $\opA : V \rightarrow W$ is bounded if and only if it is continuous.

> *Proof.*
> 
>  $\boxed \Rightarrow$ If $\opA$ is bounded then there exists a constant $c$ such that for all $x \in V$ we have $\| \opA x\left\|_{W} \leq c \right\| x \|_{V}$. Let $(x_{n})_n \subset V$ be a convergent sequence with limit $x \in V$. Then
> ```math
>  	\| \opA x_{n}- \opA x \|_{W} \leq c \|x_{n}-x \|_{V}
> ```
> by linearity of $\opA$. 
> Since $\left\|x_{n}-x\right\| \rightarrow 0$ as $n \rightarrow \infty$, $\opA x_{n} \to \opA x$ as $n \rightarrow \infty$, i.e. $\opA$ is continuous.
> 
> 
>  $\boxed \Leftarrow$ 
> If $\opA$ is continuous we can state an $\varepsilon$-$\delta$ criterion: For any $\varepsilon>0$, there is a $\delta>0$ such that 
> ```math
> 	\| \opA x_{1}-\opA x_{2} \|_{W} \leq \varepsilon \quad \forall x_{1}, x_{2} \in V \text{ with } \|x_{1}-x_{2} \|_{V} \leq \delta
> ```
> We set $\varepsilon=1$ and consider the case $x_{1}=0$ and $x_{2}=x \in V$ arbitrary. 
> Then there exists a $\delta>0$ such that
> ```math
> 	\| \opA x \|_{W} \leq 1 \quad \text { if } \quad \|x\|_{V} \leq \delta
> ```
> By linearity, we equivalently find, for all $\tilde{x} \in V$,
> ```math
> 	\|\tilde{x}\|_{V} \leq 1 \Rightarrow\| \opA \tilde{x}\|_{W} \leq 1 / \delta
> ```
> Therefore, we can define a constant $C>0$
> ```math
> 	C=\max_{\substack{ x \in V\\\| x \| \leq 1}} \frac{\| \opA x\|_{W}}{\|x\|_{V}} \leq 1 / \delta<\infty \tag{1}
> ```
> Denoting $x=\|x\|_V \hat{x}$ with $\| \hat x \|=1$, for all $x \in V$ we have
> ```math
>	\|\mathcal{A} x\|_{W}=\|x\|_{V}\| \opA \hat{x}\|_{W} \stackrel{(1)}{\leq}\|x\|_{V} C \|\hat{x}\|_{V}=C\|x\|_{V}
> ```
> which proves that $\opA$ is bounded.
> $\hspace{8cm} \square$ 
		 
""")

# ╔═╡ d332ddd1-942a-41c0-b8d1-cef005693f0d
md"""
!!! tip "Remarks"
	- In **finite dimensions** case **all linear operators** (matrices) **are   bounded** and thus continuous.

	- Typical **Schrödinger operators are not bounded**.

"""

# ╔═╡ 24f5dd6c-3fa5-4486-ad92-c1bbab44b4ce
md"""
## Spectra of operators

Extending our discussion of the finite-dimensional setting in the [Matrix perturbation theory chapter](https://teaching.matmat.org/error-control/07_Matrix_perturbation_theory.html) we can define:

"""

# ╔═╡ 59c06a2d-d980-458f-bf7a-44bb4f4a8a80
md"""
!!! note "Definition (Resolvent set)"
	Let $\hilbert$ be a separable Hilbert space, $\opA : D(\opA) \rightarrow \hilbert$. The **resolvent set** is
	```math
		\resolvent(\opA)=\{z \in \mathbb{C} \mid (\mathcal{A}-z) : D(\opA) \rightarrow \hilbert \text{ is invertible } {\color{noteblue} \underbrace{\color{black} \text{with bounded inverse}}_{\text{This is new in $\infty$-dimensions}}} \}
	```	
	where bounded inverse means that the operator
	```math
		(\opA -z)^{-1}: \hilbert \rightarrow D(\opA) \leq \hilbert
	```
	exists and there is a $C>0$ such that
	```math
	\left\|(\opA-z)^{-1} x\right\| \leq C \|x\| \quad \forall x \in \hilbert .
	```
	For $z \in \resolvent(x)$ the **resolvent** $R_{z}(\opA)=(\opA-z)^{-1}$ exists and is bounded.

We new that the new aspects in infinite dimensions is the additional requirement for $A - z$ to be bounded. A rationale for this requirement is given below:
"""

# ╔═╡ 625f255b-85b4-4ada-90bb-bb4c7d4ac597
Foldable("Rationale for the additional boundedness requirement", md"""
Having discussed the basic properties of operators, we now turn our attention towards their spectra.

- Generalizing from matrices, an eigenpair $(\lambda, \varphi) \in \mathbb{C} \times \hilbert$ of the operator $\opA$ satisfies
  ```math
  \tag{$\ast$}
  \begin{align}
  	\opA \varphi=\lambda \varphi & \Leftrightarrow(\opA-\lambda) \varphi=0 
  	\\
  	&\Leftrightarrow  \varphi \in \ker(\opA-\lambda)
  \end{align}
  ```
  where by $\opA-\lambda$ we understand the operator $\opA -\lambda \operatorname{id}_\hilbert$.

- Under a slight abuse of notation we can write
  $\ker(\opA - \lambda) = (\opA-λ)^{-1} \{0\}$,
  i.e. the solution set of applying the **resolvent**
  $R_{z}(\opA)=(\opA-z)^{-1}$ to zero.

- For eigenvalues $\opA-λ$ is clearly not invertible as $(\ast)$
  shows that $\opA-λ$ is not injective. However, let's assume for a second
  it was and that a unique non-zero solution $φ = (\opA-λ)^{-1} 0$ existed.
  Then we had
  ```math
  0 < \|φ\| = \left\|(\opA-\lambda)^{-1} 0\right\|.
  ```
  This implies that $(\opA-λ)^{-1}$ *cannot be bounded*, because if it was
  than there existed a constant $C>0$ such that
  ```math
  	\left\|(\opA-\lambda)^{-1} f \right\| \leq C\|f\| \qquad \forall f \in \hilbert
  ```
  and in particular
  ```math
  \left\|(\opA-\lambda)^{-1} 0\right\| = C \, \|0\| = 0,
  ```
  which is a contradiction to $φ$ being non-zero.

As in the finite-dimensional case we first construct the *resolvent set*, which includes all the points that cannot be eigenvalues, i.e. the ones where the resolvent exists ($\opA-z$ can be inverted) and also the above aspect of a non-bounded $(\opA-λ)^{-1}$ is excluded:
""")

# ╔═╡ 8db7d0bc-ff6b-4387-8009-d83e74faaad0
md"""
By construction the set $\resolvent(A)$ contains all $z \in \mathbb{C}$ for which $(\opA-z) x=y$ admits a unique solution $x \in D(\opA)$ for a given $y \in \hilbert$. 
Thus, for $(\opA-z) x=0$, only the trivial solution $x=0$ is possible. 

To obtain eigenvalues we thus have to study the complement, as before
"""

# ╔═╡ 86e5f05a-6562-48c4-80f4-10c7cee0698e
md"""
!!! note "Definition (Spectrum)"
	The **spectrum** is $\sigma (\opA) = \mathbb C \setminus \resolvent(\opA)$.
"""

# ╔═╡ 9823dc80-2adb-4e21-9588-fdd7dc1b3545
md"""
From the definition of $\resolvent(\opA)$ there can be three reasons for a value $\lambda \in \mathbb{C}$ to be in $\sigma(\opA)$, respectively not in $\rho(\opA)$, namely

- (1)  $(\opA-z)$ not injective $\qquad$ (*eigenvalues* or *point spectrum*, see below)

- (2)  $(\opA-z)$ not surjective.

- (3)  $(\opA-z)^{-1}$ not bounded.

Recall that injectivity and surjectivity together make ($\opA - z$) invertible.

For eigenpairs $(\lambda, \varphi) \in \mathbb{C} \times \hilbert$ one desirable property is that $(\lambda, C \varphi)$ for $0 \neq C \in \mathbb{C}$ is also an eigenpair. 
Thus, we want $(A-z) \varphi=0$ to admit multiple solutions, implying $(\opA-z)$ to be not injective.
"""

# ╔═╡ 1b98187a-2dde-46b4-8489-3f039d7a7362
md"""
!!! note "Definition (Eigenvalue, point spectrum, eigenspace)"
	The set 
	```math
		\sigma_{p}(\opA)=\{z \in \mathbb{C} \mid \opA-z \text{ is not injective}  \}
	```
	is the set of **eigenvalues** or **point spectrum** of $\opA$.
	To each eigenvalue we can associate an **eigenspace**
	```math
	\ker(\opA-\lambda)=\{\varphi \in D(A) \mid \opA \varphi=\lambda \varphi\} \neq\{0\} .
	```
"""

# ╔═╡ 3ae93aa5-4998-4d33-9126-9ce6c055d84a
md"""
!!! tip "Remark (Spectra in finite dimensions)"
	For $\hilbert=\mathbb C^{N}$, (3) is always given, and (1) and (2) are equivalent because
	```math
		\operatorname{dim}(\ker(\opA-\lambda))+\operatorname{dim}(\im (\opA-\lambda))=N .
	```
	Therefore $\sigma(\opA)=\sigma_p(\opA) \neq \varnothing$ (see Lemma 3 in [Matrix eigenproblems](https://teaching.matmat.org/error-control/03_Matrix_eigenproblems.html)).
"""

# ╔═╡ feb88cc7-66ac-463c-81b4-a223902c5897
md"""
In infinite dimensions one can easily construct examples where conditions (1), (2) or (3) fail separately, so the **spectrum can contain elements**, which are **not eigenvalues**. Similarly **$\sigma(\opA)$ can be empty**. See the Foldable below for explicit constructions of some of these cases.
"""

# ╔═╡ 8d7c4bd7-fc1e-4dc4-869c-ae33a27b015c
Foldable("Infinite dimensions: Non-surjective resolvent or empty spectrum (optional)", md"""
- In infinite dimensions injective **operators can fail to be surjective**.
  Therefore the **spectrum can contain elements**, which are **not eigenvalues**.
  To illustrate this, let $\hilbert=\ell^{2}(\mathbb{R})$ and define $\opA : \hilbert \to \hilbert$ by $\opA\left(\left(x_{n}\right)_{n}\right)=\left(0, x_{1}, x_{2}, \ldots\right)$ for $\left(x_{n}\right)_{n} \subset \mathbb{R}$ (i.e. we just insert a zero element at the front of the sequence). 
  $\opA$ is injective since $\opA \left(\left(x_{n}\right)_{n}\right)=0$ implies $\left(x_{n}\right)_{n}=0$, but is not surjective.
  So $0 \in \sigma(A)$, even though it is *not* an eigenvalue.


- A case of **surjective but not injective** can also be constructed.


- In infinite dimensions **$\sigma(\opA)$ can be empty**. 
  Consider $\op P_x f=-i \frac{\partial}{\partial x} f$ on $\hilbert=L^{2}(( 0,1) )$ with $D ( \op P_{x} )=\{f \in H^{2}((0,1)) \text{ with } f(0)=0 \}$. We define
  ```math
  \left(\op S_{z} f\right)(x)=i \int_{0}^{x} e^{i z(x-y)} f(y) d y
  ```
  as an operator on $\hilbert$. 
  We compute $1=(\op P_x -z) \op S_{z}= \op S_{z}( \op P _x -z)$ for all $z \in \mathbb{C}$. 
  So $\op P_x-z$ is invertible. 
  Still $\op S_{z}$ is bounded because $\left\|\op S_{z} f\right\|_{L^{2}} \leq \|f\|_{L^{2}}$, so have $z \in \resolvent(\opA)$.
""")

# ╔═╡ 871ab9cb-8fd3-407d-8a32-15d3c04c1cd5
md"""
A nice result is that **$\sigma(\opA)$ is closed** and **$\rho(\opA)$ is open**, which is a consequence of the resolvent $R_z(\opA)$ being analytic for a small open ball around every $z$, see the Lemma2 in the foldable for more details.
"""

# ╔═╡ 6caaef7c-085f-433e-9a08-16f1b8243759
Foldable("Lemma 2: Discussion Analyticity of resolvent (optional)", md"""
!!! note "Lemma 2"
	Let $\opA : D(\opA) \rightarrow \hilbert$ be a linear operator and $z \in \resolvent(\opA)$.
	Then, $B_r(z) \subset \resolvent(A)$ with radius
	```math
		r=\frac{1}{\|\left(\opA-z)^{-1} \|\right.}
	```
	and the resolvent can be expressed as the series 
	```math
		R_{t}(\opA)=(\opA-t)^{-1}=\sum_{n \geq 0}(t-z)^{n}(\opA-z)^{-(n+1)} \quad \forall t \in B_{r}(z)
	```
	Since, $\forall z \in \resolvent(A)$, this allows us to construct an open neighbourhood, which is also part of $\resolvent(A)$, we deduce that $\resolvent(\opA)$ is *open* and $\sigma(\opA)$ is *closed*.

> *Proof* (same as the argument as Proposition 3 of [perturbation theory](https://teaching.matmat.org/error-control/07_Matrix_perturbation_theory.html)).
> Let $z \in \rho(x)$ and $t \in B_{r}(z) \subset \rho(\opA)$. Then
> ```math
> \begin{align}
> 		(\opA-t)^{-1} &= [(\opA-z)-(t-z)]^{-1}  \\
>		&= {\color{gray} \underbrace{\color{prooftext} (\opA-z)^{-1}}_{\substack{\text{Inverse exists} \\ \text{because } t \in \resolvent(\opA)}}} 
>		\left [\operatorname{id}-(t-z)(\opA-z)^{-1}\right]^{-1}
> \end{align}
> \tag{2}
> ```
> Since $|t-z|<\frac{1}{\left\|(\opA-z)^{-1}\right\|}$, we have $\left\|(t-z)(\opA-z)^{-1}\right\|<1$, which admits a Neumann expansion in the last term of (2), giving the desired result.
> $\hspace{5cm} \square$
""")

# ╔═╡ b2af25d8-e00d-4478-8ddc-ca809478b4b1
md"""
## Self-adjoint operators

As discussed in principle multiple choices for the operator domain $D(\opA)$ can sometimes be reasonable and thus seem to leave considerable ambiguity.

- As the Schrödinger operator example shows, there is clearly some form of a maximal domain, since we clearly need $\varphi \in D(\opH)=H^{2} (\mathbb R^{d})$ to assure $\opH \varphi \in L^{2} (\mathbb{R}^{d} )$. 
  But what about choosing a smaller domain, e.g. $D(H)=C_{0}^{\infty} (\mathbb{R}^{d} )$ ?


- It turns out ([detailed discussion in the first appendix](#Appendix-:-In-depth-treatment-of-self-adjointness)) that the choice of a sufficiently large domain $D(\opA)$ is required for the operator to be *closed*, which is a slightly weaker form of being continuous. 
  Being closed is necessary for being self-adjoint and moreover if $\opA$ is not closed, then $\sigma(A)=\mathbb{C}$, which is not useful for a physical model.


- Roughly speaking we therefore want $D(\opA)$ to be exactly large enough to get self-adjointness.


Without going into all such detail we will now sketch the definition of a self-adjoint operator and refer to the [first appendix](#Appendix-:-In-depth-treatment-of-self-adjointness) for more details.
"""

# ╔═╡ e712629f-eaeb-4fd1-9aa0-fd3ec4765148
md"""
We first need

!!! note "Definition (Adjoint operators)"
    Given an operator $\opA: D(\opA) \rightarrow \hilbert$ 
    with $D(\opA) \subset \hilbert$, the adjoint operator
    $\opA ^*: D (\opA^{*} ) \rightarrow \hilbert$ is defined by $y \mapsto \opA^{*} y$,
    where $\opA^{*} y \in \hilbert$ is the unique element satisfying
    ```math
    \left\langle x, \opA^{*} y\right\rangle=\langle\opA x, y\rangle. \qquad \forall x \in D(\opA), \forall y \in D\left(\opA^{*}\right).
    ```
    Notably in this definition $D (\opA^{*} )$ is taken to be the
    *largest domain possible* for this equality to make sense.
    Uniqueness of $\opA^{*} y$ given $y$ is than guaranteed
    by the Riesz representation theorem.
"""

# ╔═╡ 974ce554-562c-4d32-aab6-29513f2fff5a
md"""
and

!!! note "Definition (Symmetric operators)"
	The operator $\opA$ on $\hilbert$ is called **symmetric** if, for all $x, y \in D(\opA) \times D(\opA)$,
	```math
		\langle\opA x, y\rangle= \langle x, \opA y\rangle.
	```

This is a first step towards self-adjoint operators and typically what is called "Hermitian operator" in physics texts.

Comparing this definition of symmetric operators and the adjoint definition,
given in the foldable, one notices 
above we notice the difference in the sets over which $y$ pans out:
**Symmetry and self-adjointness are not equivalent concepts**. 

"""

# ╔═╡ 25804d6a-9eb4-4c24-89cf-96042eac70d8
md"""
!!! note "Definition (Self-adjointness)"
	An operator is self-adjoint if it is symmetric with $D(\opA^*)=D(\opA)$, or simply an operator for which $\opA=\opA^{*}$.

Notably this implies that if a symmetric operator $\opA : \hilbert \to \hilbert$ is defined on the entire $\hilbert = D(\opA)$, then $\opA$ is self-adjoint.
"""

# ╔═╡ a15f5cc7-753d-4aa9-b8b5-6b701882d310
md"""
The following result provides the familiar result that **self-adjoint operators** have a real spectrum, i.e. a **physically meaningful spectrum**.
"""

# ╔═╡ 187f55b0-f71d-4f9d-b00b-e68966b7964d
md"""
!!! note "Theorem 3"
	Let $\opA$ be a symmetric operator with domain $D(\opA) \subset \hilbert$,
	that is 
	```math 
	\begin{align}
	    \langle \varphi, \opA \varphi \rangle & = \langle \opA \varphi, \varphi \rangle && \forall \varphi \in D(\opA) 
	\end{align}
	```
	 Then the following statements are equivalent
	
	1.  $\opA$ is self-adjoint, i.e. $D(\opA^*) = D(\opA)$
	
	2.  $\varnothing \neq \sigma(\opA) \subset \mathbb R$
	
	3.  $\exists \lambda \in \mathbb C$ such that $\opA - \lambda$ and
	    $\opA - \bar \lambda$ are both surjective from $D(\opA) \to \hilbert$.
"""

# ╔═╡ 392444bd-c154-4088-aa00-2b32018a90d6
md"""
### Summary of operator properties

- For infinite dimensions $D(\opA)=\hilbert$ in rarely feasible. 
  Choosing $D(\opA)$ sensibly (e.g. Sobolev spaces) is crucial to get self-adjointness and physically meaningful spectra.


- Many the special properties we discussed are automatically given in finite dimensions or can be easily checked (e.g. symmetry).


- We have the following implication relationships :
```math
\boxed{
\begin{align}
	\text{self-adjoint} &\Rightarrow \text{closed} \\
	\text{self-adjoint} &\Rightarrow \text{symmetric} \\
	\text{self-adjoint} &\Rightarrow \sigma(\opA) \subset \mathbb{R}
\end{align}}
```

```math
\boxed{
\begin{align}
	\text{(bounded} &\Leftrightarrow \text{continuous}) \Rightarrow \text{closed}
	\\
	&\Uparrow
	\\
	\dim(\hilbert) &< \infty \Rightarrow {\color{gray} \underbrace{\color{black} \sigma(\opA) = \sigma_p(\opA)}_{\text{spectrum only has eigenvalues}}}
	\\
	&\Downarrow
	\\
	D(\opA) &= \hilbert \Rightarrow (\text{symmetric} \Leftrightarrow \text{self-adjoint})
\end{align}}
```

"""

# ╔═╡ 8cdc2867-259f-4242-88ab-d17948cfc1bb
md"""
### Examples of self-adjoint operators

We provide some examples of self-adjoint operators on $\hilbert=L^{2}(\mathbb R^{d})$ and state the spectra of some of them (Proofs are given
[in the second appendix on Weyl sequences](#Example-spectra-of-operators)).
We will note that the previously introduced **Sobolev spaces arise naturally** as the correct **domain to ensure self-adjointness** of physical operators.
"""

# ╔═╡ a53f5693-6a82-481c-a1ba-cb10652e96ae
md"""
!!! warning "Identity"
	 $\mathop{\mathrm{id}}: \hilbert \to \hilbert$ is self-adjoint with
	$D(\mathop{\mathrm{id}}) = \hilbert$. As expected $\sigma(\mathop{\mathrm{id}}) = \{ 1\}$.
"""

# ╔═╡ b0798ff6-c3e8-4cd4-bff7-c7de1faf1c97
md"""
!!! warning "Laplace operator"
	 $\op A=-\laplacian$ with $D(-\Delta)=H^{2}(\mathbb{R}^{d})$ is self-adjoint, but with $D(-\Delta)=C_{0}^{\infty} (\mathbb R^{d})$ it is not.

	Using $D(-\Delta)=H^{2}(\mathbb{R}^{d})$ its spectrum is $\sigma(-\Delta) = [0, \infty)$.
"""

# ╔═╡ 5f17289c-0234-444e-9e7b-8d39c2851b96
md"""
!!! warning "Momentum operator"
	 $\op P_{j}=-i \frac{\partial}{\partial x_{j}}$ with 
	$D(\op P_{j})=H^{2}(\mathbb{R}^{d})$ is self-adjoint.
"""

# ╔═╡ 47c76553-a775-4e89-a6bf-abba188e792c
md"""
!!! warning "Schrödinger operators"
	The Schrödinger operators $\opH=-\Delta+V$ in $\hilbert=L^{2}\left(\mathbb{R}^{3}\right)$ are self-adjoint for $V \in L^{2}(\mathbb{R}^{3})+L^{\infty}(\mathbb{R}^{3})$ and $D(\opH)=H^{2}(\mathbb R^{3})$. The condition on the potential is to be understood as $V=V_{S}+V_{L}$ with $V_{S} \in L^{2}(\mathbb{R}^{3})$ (short-range component) and $V_{L} \in L^{\infty}(\mathbb{R}^{3})$ (long-range component). 
	This is satisfied for
	- Any bounded potential $V= \color{warnyellow} \underbrace{\color{black} 0}_{\in L^2(\mathbb R^3)} {\color{black}+} \underbrace{\color{black} V}_{\in L^\infty(\mathbb R^3)}$
	- The Coulomb potential
	```math
		V=\frac{1}{|x|} \quad \text { with } \quad V_{S}= 
		{\color{warnyellow} \underbrace{\color{black} \frac{e^{-|x|}}{|x|}}_{\in L^{2}(\mathbb R^{3})} }
		\text{ and } V_{L}  = 
		{\color{warnyellow} \underbrace{\color{black} \frac{1-e^{-|x|}}{|x|}}_{\in L^\infty(\mathbb R^3)}}
	```
"""

# ╔═╡ 6a10bd68-78fd-4323-a225-f545a0ee6b31
md"""
## Spectral characterization

Not all members of the spectrum have the same mathematical properties. In particular when it comes to approximating the spectrum, some parts of it are easier to obtain numerically than others.

"""

# ╔═╡ 60f38f19-df1f-4134-9b08-bb17d278ea4b
md"""
We already introduced the concept of an eigenvalue as well as the point spectrum $\sigma_p$ :
```math 
\begin{align}
    \sigma_p(\opA) &= \{ \lambda \in \mathbb C \mid \lambda \text{ is an eigenvalue of } \opA \}
	\\
	&= \{ \lambda \in \mathbb C \mid \mathrm{Ker}(\opA - \lambda) \neq \{0\} \}
\end{align}
```
A natural characterization of the spectrum is to spit $\sigma (\opA) = \sigma_p (\opA) \ \dot \cup \ \sigma_{\text{cont}} (\opA)$ where 

!!! note "Definition 3 (Continuous spectrum)"
	Let $\opA : D(\opA) \to \hilbert$ be a self-adjoint operator.
	We call the *continuous spectrum* $\sigma _{\text{cont}}$
	```math 
	\begin{align}
	        \sigma_{\text{cont}} (\opA) = \sigma \setminus \sigma_p (\opA)
	    
	\end{align}
	```

This classification turns out to be not very useful for considering spectral approximations as small perturbation of $\opA$ can easily mix $\sigma_p$ and $\sigma_{\text{cont}}$.
Whe therefore use an alternate :

"""

# ╔═╡ 3c0a48e2-fe1b-4150-8435-9f1a28f6340a
md"""
!!! note "Definition 4 (Isolated point)"
	Let $\opA : D(\opA) \to \hilbert$ be a self-adjoint operator.
	An entry
	$\lambda \in \sigma(\opA)$ is called an *isolated point* if there exists $\varepsilon > 0$ such that
	```math 
	\begin{align}
	[\lambda - \varepsilon, \lambda + \varepsilon] \cap \sigma(\opA) = \{ \lambda \}   
	\end{align}
	```

!!! tip "Isolated points are eigenvalues"
    One can show that every isolated point of $\sigma(\opA)$ is an eigenvalue.


They are thus interesting entries, motivating the following spectrum :
"""

# ╔═╡ afe2868e-56ed-4e96-bd1a-3f6546cdd2b8
md"""
!!! note "Definition 5 (Discrete and essential spectrum)"
	Let $\opA : D(\opA) \to \hilbert$ be a self-adjoint operator.
	```math
		\sigma_{\text{disc}} = \{ \lambda \in \sigma(\opA) \mid \lambda \text{ is an eigenvalue and has finite multiplicity} \}
	```
	is called the *discrete spectrum*.

	The complement
	```math
		\sigma_{\text{ess}} = \sigma(\opA) \setminus \sigma_{\text{disc}} (\opA)
	```
	is called the *essential spectrum*.


"""

# ╔═╡ 34e2d8bb-1a7f-48fa-b17d-c4a84f76cea0
md"""
### Typical spectra of Schrödinger operators

Let us illustrate these concepts by sketching the spectrum of Schrödinger operators $\opH = - \laplacian / 2 + V$ on $\hilbert = L^2(\mathbb R^3)$.
We take $V \in L^2(\mathbb R^3)  +  L^\infty_\varepsilon(\mathbb R^3)$, that is to say that the potential needs to be split into a square integrable ($L^2$) short-range part and a bounded long-range part that also is negligible at infinity ($L^\infty_\varepsilon$).
For example, the Coulomb potential satisfies this.
Since $L^\infty_\varepsilon \subset L^\infty$, $\opH$ is self-adjoint using $D(\opH) = H^2(\mathbb R^3)$.

Then one can show 
```math
	\sigma_{\text{ess}} (\opH) = [0, \infty).
```
Therefore, the spectrum looks like this :

---
"""


# ╔═╡ 6535886e-94e9-47e3-aeea-d3db829e7566
TikzPicture(L"""
% %real line
    % \draw[>=latex,->] (-5,-5) -- (5,0) node[right]{$\mathbb R$};
    % \draw (0,0) node{$|$} ;
    % \draw (0,0.1) node[above]{0};

    %sigma_p
    \foreach \time [evaluate=\time as \energy using {- 5 / \time^1 + 5/30} ] in {1,...,30}{
        \draw (\energy,-1.5) node{$\times$};
        }

    \foreach \energy  in {2,2.5,3.5,0}{
        \draw (\energy,-1.5) node[color=blue]{$\times$};
        }   

    \draw (5,-1.5) node[right]{$\sigma_p$};

    \draw[->] (0,-1) node[above]{$\infty$ multiplicity} -- (0,-1.35);

	\foreach \energy  in {2,2.5,3.5}{
        \draw[->,blue] (\energy,-1) -- (\energy,-1.35);
        }   

	\draw[blue] (3,-1) node[above,blue]{Non-isolated eigenvalues} ;

    %sigma_cont

    \draw (0,-2) node[blue]{$($} -- (5,-2) node[right]{$\sigma_{\textrm{cont}}$};

    \foreach \energy  in {2,2.5,3.5}{
        \draw (\energy,-2) node[color=white]{$\bullet$} node[blue]{$\circ$} ;
        }   

    %sigma_disc
    \foreach \time [evaluate=\time as \energy using {- 5 / \time^1 + 5/30} ] in {1,...,30}{
        \draw (\energy,-3) node{$\times$};
        }

    \draw (5,-3) node[right]{$\sigma_{\textrm{disc}}$};

    %sigma_ess
    \draw (0,-3.5) node[blue]{$[$} -- (5,-3.5) node[right]{$\sigma_{\textrm{ess}}$};
""",width="20cm",options="scale=1.1")

# ╔═╡ aefb411d-53f2-4571-a8b7-170e5b09240e
md"""
---

Notice how eigenvalues (members of $\sigma_p$) may be embedded inside the continuous spectrum.
These are very unstable to perturbations (so-called resonance states) and, as we will see now, hard to approximate.

Notice that, unlike $\sigma_{\text{cont}}$, the essential spectrum $\sigma_{ess}$ is always closed for a self-adjoint operator.


"""

# ╔═╡ f33be757-b599-4356-b713-41281124529d
md"""
## Courant-Fisher

To finish up our discussion we will now work towards a generalization of the Courant-Fisher characterization.
For this, we once again consider the quadratic form $q_\opA(u)$
```math
    q_\opA (u) = \langle u, \opA u \rangle
```
and the sequilinear form
```math
	    a_\opA (u,v) = \langle u, \opA v \rangle
```
which we had already encountered in the context of matrices.

A priori, these are only defined for $u,v \in D(\opA)$.
However, considering a completion
in the norm induced by the quadratic form, uniquely leads
to a **larger space $Q(\opA)$**, called the **form domain** of $\opA$
on which $q_\opA$ can be defined. Notice that we can find a subspace inclusion
```math 
\begin{align}
    D(\opA) \leqslant  Q(\opA) \leqslant \hilbert.
\end{align}
```
"""

# ╔═╡ 86d8f7b6-68f6-42fd-9564-9c8a82a35e70
md"""
To make this explicit, consider an example. We take $\opA = - Δ$ on $L^2(\mathbb{R}^d)$.
Using partial integration we rewrite (e.g. for sufficiently regular functions $u$ and $v$):
```math
\langle u, -Δ v \rangle = \int_{\mathbb{R}^d} u (-Δ v) = \int_{\mathbb{R}^d} ∇u ⋅ ∇v.
```
While the LHS is only valid for $u, v \in H^2(\mathbb{R}^d)$,
the RHS can be easily extended to all $u, v \in H^1(\mathbb{R}^d)$
and indeed a rigorous construction finds $Q(-\Delta) = H^1((\mathbb{R}^d)$
as its form domain.

The following foldable provides some more details on the construction of
the form domain.
"""

# ╔═╡ e13f6917-81b8-4b4f-bad9-9b30cbc1db51
Foldable("Details on the form domain (optional)", md"""
If $\opA$ is bounded from below then there exists $\alpha > -\infty$ such that
```math
	q_\opA (u) \geq \alpha \| u \| ^2
```
Note that is the case for all Schrödinger operators we consider.
For simplicity, we take $\alpha > 0$ (otherwise consider a shifted operator $\tilde \opA = \opA - \alpha + \varepsilon$ with $\varepsilon >0$).
Then 
```math
	\begin{align}
    \langle u,v \rangle_\opA \equiv \langle u, \opA v \rangle && u,v \in D(\opA)
	\end{align}
```
is an inner product with induced norm
```math 
\begin{align}
    \| u \|_A \equiv q_\opA(u) > 0 && u \in D(\opA) \tag{Energy norm}
\end{align}
```
 $D(\opA)$ is not complete with respect to $\| \cdot \|_\opA$.
Therefore, the idea is now to consider sequences $(u_n)_n \subset D(\opA)$, which are clearly Cauchy with respect to $\| \cdot \|_\opA$, and extend the definition, i.e. if $u_n \to u$ then we demand $q_\opA (u_n) \to q_\opA (u)$.
In this way, we can extend the definition *uniquely* to $a_\opA (u,v)$, $q_\opA(u)$ with $u,v \in Q(\opA)$, called the *form domain* of $\opA$, and which satisfies 
```math 
\begin{align}
    D(\opA) \leqslant  Q(\opA) \leqslant \hilbert
\end{align}
```
Where in each case the embedding is dense and continuous (Lewin [*Théorie spectrale et mécanique quantique*](https://doi.org/10.1007/978-3-030-93436-1) 2022, Theorem 3.10).
""")

# ╔═╡ 7bd639d0-0918-4872-b2ef-b23da592e2b8
md"""
!!! tip "Remark"
	Typical form domains use again Sobolev spaces. For example 
	```math
	\begin{align}
		Q(- \laplacian) &= H^1(\mathbb R^d) && \text{on } L^2(\mathbb R^d)
		\\
		Q(- \laplacian + V) &= H^1(\mathbb R^3) && \text{on } L^2(\mathbb R^3) \text{ with } V \text{ appropriate}
	\end{align}
	```

"""

# ╔═╡ 2ae96252-9aa1-4c55-aff5-e74a90472c2c
md"""For us the main purpose of form domains is that they are the "right domain" to use to formulate the Courant-Fisher theorem for operators. Some more ideas why form domains are useful are given in the Foldable below.
"""

# ╔═╡ bc90f5ca-2c95-4345-a315-f73f2fbfe48c
Foldable("Why do we even need a form domain ? (optional)", md"""

One importance of the form domain is that is allows for a weak formulation :

!!! note "Theorem 5"
	Let $\opA$ be self-adjoint. The following are equivalent
	
	1.  $\varphi \in D(\opA)$ and $\opA \varphi = \lambda \varphi$. (*Strong formulation*)
	
	2.  $\varphi \in Q(\opA)$ and $\exists \lambda \in \mathbb R$ s.t.
	    $a_A(f,\varphi) = \lambda \langle f, \varphi \rangle$ for all $f \in Q(\opA)$. (*Weak
	    formulation*)

Another reason is that it makes working with singular 
potentials $V$ more convenient. E.g. consider a potential $V$
with singularities even worse than standard Coulomb
(e.g. $V = \frac1{|r|^α}$ and $0 < α < 1$) the appropriate
domain to make the Schrödinger operator $(-Δ + V)$ self-adjoint can be
very subtle, e.g. we at least need $(-Δ + V)ψ \in L^2$ for all $ψ \in D(-Δ+V)$,
such that $D(-Δ+V)$ is not simply $H^2$.
However, one can show that in such cases we still have $Q(-Δ+V) = H^1$),
see chapter 3.3 in Lewin [*Théorie spectrale et mécanique quantique*](https://doi.org/10.1007/978-3-030-93436-1) 2022.

""")

# ╔═╡ 2d67bf62-0d3a-4b66-8028-53661602719d
md"""
Finally, we can state the equivalent the min-max principle that we already saw for matrices, extended to the case of operators. 

!!! note "Theorem 7 (Courant-Fisher)"
	Let $\opA$ be self-adjoint and bounded
	from below with form domain $Q(\opA) \leqslant \hilbert$. 
	Then
	```math 
	\begin{align}
	        \mu_k (\opA) \coloneqq \inf_{\substack{W \subset Q(\opA) \\ \dim (W) = k}} \max_{\substack{\varphi \in W \\ \| \varphi \|_\hilbert = 1}} q_\opA(\varphi) = \inf_{\substack{W \subset Q(\opA) \\ \dim (W) = k}} \max_{0 \neq \varphi \in W} R_{\mathcal{A}}(\varphi)
	\end{align}
	```
	 is equal to
	
	1.  The k$^{\textsf{th}}$ eigenvalue (counting multiplicities) of $\opA$ if $\opA$ has at
	    least $k$ eigenvalues below $\Sigma(\opA).$
	
	2.  Otherwise, $\mu_k(\opA) = \Sigma(\opA)$.

    In this we used the *bottom of the essential spectrum*
    ```math
    \Sigma(\opA) \coloneqq \left\{ \begin{array}{ll}
		\infty & \text{if } \sigma_{\text{ess}}(\opA) = \varnothing \\
		\min \sigma_{\text{ess}} (\opA) &  \text{else}
	\end{array}\right.
    ```
	and re-introduced the *Rayleigh quotient*, this time for operators 
	```math
		R_\opA (\varphi ) \coloneqq \frac{q_\opA (\varphi)}{\| \varphi \|^2_\hilbert }.
	```
"""

# ╔═╡ ab658994-7d4d-45f7-94a5-ebf1773fcc56
md"""
One way to check for the existence of $k$ eigenvalues below $\Sigma(\opA)$ is :

!!! tip ""
	If there is a subspace $W \leq Q(\opA)$ with $\text{dim} W = k$, such that 
	```math 
	\max_{0 \neq \varphi \in W} R_\opA (\varphi) < \Sigma(\opA)
	```
	then there are at least k eigenvalues.
"""

# ╔═╡ 6d7cf76a-cfc2-4477-ac8b-3a70204286ac
md"""
In other words, as long as $\mu_k(\opA) < \Sigma(\opA)$, the situation is as for Hermitian matrices.
As soon as $k$ eigenvalues have been found, we can no longer learn anything about $\sigma(\opA)$ by increasing the subspace size.

This leads to a severe restriction for the numerical approximation of eigenspectra, as we are in most cases restricted to below $\Sigma(\opA)$.
However, in this case we are able to estimate eigenspectra, and obtain error bounds in a manner similar to that of matrices.
	
"""

# ╔═╡ fa883a72-9e03-4a66-8625-89273f19073b
md"""
## Approximation of $\lambda_k < \Sigma(\opA)$ 

A straightforward technique to approximate the eigenspace of operators is to employ a simple *projection technique*, i.e. to employ a *finite-dimensional subspace* $S \subset Q(\opA)$ in the infimum.
The idea is the same as we discussed in the context of diagonalization algorithms.
*Assuming* $\mu_k(\opA) < \Sigma(\opA)$, i.e. that $\opA$ has $k$ eigenvalues $\lambda_1, \dots, \lambda_k$ (counting multiplicities) below $\Sigma(\opA)$ yields
```math 
\begin{align}
\tag{2}
    \lambda_k = \mu_k (\opA) = \inf_{\substack{W \subset {\color{blue}Q(\opA)} \\ \dim (W) = k}} \max_{0 \neq \varphi \in W} R_A(\varphi) \leq \inf_{\substack{W \subset {\color{blue} S} \\ \dim (W) = k}} \max_{0 \neq \varphi \in W} R_A(\varphi)
\end{align}
```


- Assuming $S$ to be $d$-dimensional, and taking a basis $\mathbb B = (\chi_1, \dots, \chi_d)$, we realize that the right-hand side of (2) to be nothing else than the $k$-th eigenvalue of the matrix
  ```math 
  \begin{align}
      (M_S^\opA)_{ij} = q_\opA(\chi_i, \chi_j) \cong \langle \chi_i , \opA \chi_j \rangle
  \end{align}
  ```
  Note that the second equality is only fully correct for $S \subset D(\opA)$. In other words 
  ```math 
  \begin{align}
      \lambda_k (\opA) = \mu_k(\opA) \leq \lambda_k (M_S^\opA)
  \end{align}
  ```
"""

# ╔═╡ 4c2952d8-2957-42ca-b05a-068941cf5250
md"""
- In practice one then keeps augmenting the subspace $S$, e.g. by increasing the approximation basis $\mathbb B$ following the principle 
  ```math
  	\lambda_k(\opA) = \mu_k(\opA) = \inf_{\substack{S \subset Q(\opA) \\ \dim Q(\opA) \geq k}} \lambda_k (M_S^\opA)
  ```

- Notice, however, how this technique is *unable to approximate above $\Sigma(\opA)$.*
  In other words, *some eigenvalues* of $\opA$, those embedded in $\sigma_{\text{ess}} (\opA)$ *cannot be found* using this projection technique.
"""

# ╔═╡ 523daf58-de5b-4392-bb7f-94c1be73136f
md"""
---
#### Illustration

Eigenenergies of the 1D gaussian potential well, i.e. the spectrum of the following Schrödinger operator
```math
\mathcal H = - \frac1{2} \frac{d^2}{dx^2} - A e^{-(x/\sigma)^2}
```
"""

# ╔═╡ 3d8f2b32-046e-446f-b28e-feecc53d02cd
md"""
y axis limits (slide to zoom in/out) $\qquad$ $(@bind ylim PlutoUI.Slider([2000,1000,500,200,100,50,25]; default=1000, show_value=true))
"""

# ╔═╡ 685f4ded-2a50-4ea6-af79-04192afa2de7
begin
	using Plots
	
	plot(ylims=(-ylim,ylim),ylabel="Calculated eigenvalues",xlabel="Dimension of discrete basis")
	
	hline!([0],color=:black,label=nothing)
	for i in 1:maximum(n_range)
		ydata = eigenvalues[i,:]
		color = :gray
		if minimum(skipmissing(ydata)) < 0
			color = i
			hline!([eig.values[i]],color=i,linestyle=:dash,label=nothing,linewidth=0.5)
		end
		plot!(n_range,ydata,marker=:x,label=nothing,color=color,linewidth=0.5)
	end
	plot!()
end

# ╔═╡ 58b80515-233c-4ad5-8650-7841c381e186
md"""
---
"""

# ╔═╡ 9a1daee7-9ae8-4f19-bb97-5a9ab392c9a4
md"""
## Compact operators and compact resolvents

Two types of "special operators" are very important and will be defined in passing.

!!! info "Definition (Compact operator)"
	A self-adjoint operator $\opA : D(\opA) \to \hilbert$ with domain $D(\opA) \subset \hilbert$ is called **compact** if the image of the unit ball
	is a compact set.

	Equivalently this means that $\opA v_n \to 0$ for every
	sequence $v_n \rightharpoonup 0$, where $\rightharpoonup$ denotes
	week convergence (see [the second appendix for a definition of weak convergence](#Appendix:-Spectral-characterisation-based-on-Weyl-sequences)).
"""

# ╔═╡ d535075f-f5fa-4c9a-97cb-fda30ce42173
md"""
Not many operators are compact, but their spectral properties are very nice:

!!! info "Theorem 8 (Spectra of compact operators)"
	Let $\hilbert$ be an infinite-dimensional Hilbert space and
	$\opA$ be a self-adjoint operator. Then $\opA$ is compact if and only if
	$\sigma_\text{ess}(\opA) = \{0\}$. A self-adjoint, compact operator
	can be diagonalised in an orthonormal basis of $\hilbert$ formed by
	taking the eigenvectors corresponding to the eigenvalues.
	The eigenvalues of $\opA$ accumulate at zero ($λ_n \to 0$).
	All non-zero eigenvalues are of finite multiplcity.

For **compact operators** we thus have $\Sigma(\opA) = 0$, meaning
that we can apply projection techniques to **approximate all eigenvalues
less than $0$**.
"""

# ╔═╡ c67c8802-cf11-440a-9afe-48f137777f34
md"""
More common and still rather nice are the related class of operators:

!!! info "Definition (Operators with compact resolvent)"
	We call a self-adjoint operator $\opA$ an **operator with compact resolvent**
	(o.c.r.) if $(\opA + i)^{-1}$ is a compact operator.

	Equivalently one can show that $\opA$ is an o.c.r.
	if for all $z \in \resolvent(\opA)$ the operator $(\opA - z)^{-1}$
	is compact.
"""

# ╔═╡ 8f57a7e2-3449-4034-9309-6233c414ac0f
md"""
!!! info "Theorem 9 (Spectra of o.c.r.)"
	Let $\hilbert$ be an infinite-dimensional Hilbert space and
	$\opA$ be a self-adjoint operator. Then $\opA$ is with compact resolvent
	if and only if $\sigma_\text{ess}(\opA) = \{\}$.
	The spectrum of an operator with compact resolvent $\opA$
	consists only of eigenvalues, each with finite multiplicity.
	The sequence $\lambda_n$ of all eigenvalues in the
	spectrum of $\opA$ satisfies
	```math
	\lim_{n\to\infty} |λ_n| = +\infty.
	```
	Further, $\opA$ can be diagonalised in an orthonormal basis of $\hilbert$.

For **operators with compact resolvents** we thus find that there is **no essential spectrum**, so we can apply projection techniques to **approximate all eigenvalues**,
just like in the case of matrices.
"""

# ╔═╡ 0007d42e-5a2b-45fe-9bcc-603e54704d1e
md"""
## Kato-Temple bound

It turns out with minor modifications the Bauer-Fike and Kato-Temple bounds
still hold for operators and eigenvalues below the essential spectrum.
We formulate the Kato-Temple case explicitly:

"""

# ╔═╡ aeb26e12-b13d-4aeb-ab6a-75e7885fdb03
md"""
!!! note "Theorem 11 (Kato-Temple bound)"
	Let $\opA$ be a self-adjoint operator,
	$\tilde \varphi \in D(\opA)$ with $\| \tilde \varphi \| = 1$, $\tilde \lambda = q_\opA(\tilde \varphi)$,
	and $r = \opA \tilde \varphi - \tilde \lambda \tilde \varphi$. Let $\lambda$ be
	the eigenvalue closest to $\tilde \lambda$ and 
	```math 
	\begin{align}
	        \delta = \min_{s \in \sigma(\opA) \setminus \{ \lambda \} } | s - \tilde \lambda |
	    
	\end{align}
	```
	 such that
	$\sigma(\opA) \cap [ \tilde \lambda - \delta, \tilde \lambda + \delta ] = \{ \lambda \}$.
	Then 
	```math 
	\begin{align}
	        | \tilde \lambda - \lambda | \leq \frac{\| r \|^2}{\delta}
	    
	\end{align}
	```
"""

# ╔═╡ e07daf16-dfcd-49eb-b821-f7b8a7a5760e
md"""The outline of a proof is given in the [third appendix](#Appendix:-Proof-outline-for-Kato-Temple-bound)."""

# ╔═╡ cda93cc1-e9ba-43c8-9cdc-18b7de1315dc
md"""
# Appendix
"""

# ╔═╡ 23d46169-107b-41a7-8ef1-410f63c046bc
md"""
## Appendix : In-depth treatment of self-adjointness 

To arrive at the definition of self-adjoint operators we will first look at the graph of operators, which might seem trivial aud unrelated, but is closely related to the question of choosing a good domain $D(\opA)$ for an operator, and the relation between $D(\opA)$, the spectrum $\sigma (\opA)$, and self-adjointness.

"""

# ╔═╡ 34fafa7e-3a64-433f-b8a8-7bfaf7228dba
md"""
### Graph and closure

!!! note "Definition (Graph)"
	Let $A: D(\opA) \rightarrow \hilbert$. The **graph** of the operator $\opA$ is
	```math
	\graph(\opA)=\{(\varphi, \opA \varphi) \mid \varphi \in D(\opA )\} \subset D(\opA) \times \hilbert
	```
"""

# ╔═╡ ff0ad3c9-f430-4f99-8df5-bec5c6c06e3f
md"""
Not all subsets of $\hilbert \times \hilbert$ are graphs of operators as the following lemma shows.

!!! note "Lemma 12"
	 A set $\graph \subset \hilbert \times \hilbert$ is the graph of an operator $\opA$ if and only if

	- (1)  $\graph$ is a vector subspace of $\hilbert \times \hilbert$

	- (2)  $(0, y) \in \graph$ implies $y=0$ 

	- (3) The projection $D =\{x \in \hilbert \mid \exists y \in \hilbert$ with $(x, y) \in \graph\}$ is dense in $\hilbert$.
"""

# ╔═╡ 5f3954d3-bbe5-474a-834c-42e8603ff4bf
md"""
> *Proof*.
>
>  $\boxed{\Rightarrow}$ : (1) - (3) follow directly since $\opA$ is a linear operator and $D(\opA)$ is dense in $\hilbert$.
>
>  $\boxed{\Leftarrow}$ : (2) ensures that each $x$ is assigned a unique image $y = \opA x$, since due to (1) $\graph$ is a vector space.
> Therefore, the application of $\opA$ is well-defined and linear.
> (3) finally ensures the density of $D(\opA)$ in $\hilbert$.
"""

# ╔═╡ a41650f2-2e4d-4c76-abe7-65e8867065e5
md"""
!!! note "Definition (Extension and restriction)"
	Let $\opA$ and $\op B$ be two operators on $\hilbert$. 
	$\op B$ is called on **extension** of $\opA$ if $D(\opA) \subset D(\op B)$ and if $\forall x \in D(\op A): \op A x= \op B x$. 
	In that case we also call $\opA$ a **restriction** of $\op B$.

!!! warning "Example"
	One can define the derivative operator $f \mapsto f^{\prime}$ on $\hilbert =L^{2}(\mathbb{R})$. The domain $D(\op B)=C_{0}^{\infty}(\mathbb{R})$ (the space of compact supported smooth functions) is a valid choice since $C_{0}^{\infty} \subset L^{2}(\mathbb{R})$ leading to an operator $\opA$. Alternatively we could define $\op B$ with domain $D(\op B)= H^{1}(\mathbb{R})$, which is the largest possible domain for this operator. Since $C_{0}^{\infty} \subset H^1(\mathbb{R})$, $\op B$ is an extension of $\opA$ and $\opA$ is a restriction of $B$.
"""

# ╔═╡ 3c24ac4f-6903-4cfd-93f3-be7af58f1448
md"""
!!! note "Definition (Closed operators)"
	An operator $\opA$ with domain $D(\op A)$ is called **closed** if its graph $\graph(\opA)$ is closed in $\hilbert \times \hilbert$.
	That is, if for each converging sequence $(x_{n} ) \subset D(\opA)$ with $x_{n} \rightarrow x$ and $\opA x_n \to y$, it holds that $x \in D(\opA)$ and $\opA x= y$.


!!! tip "Remark"
	Being closed can be seen as a slightly weaker form of being continuous. 
	The latter requires that every converging sequence $x_{n} \to x$ *always* leads to a converging sequence $\opA x_{n} \to \opA x$. 
	For a closed operator, not all sequences $\left(\opA x_{n}\right)_{n}$ need to converge, and only for those that do we have the requirement $\opA x_{n} \rightarrow \opA x$.

	A consequence is 
	```math
	\begin{align}
		\textsf{bounded or continuous} \Rightarrow \textsf{closed}
	\end{align}
	```

"""

# ╔═╡ 3b36fc40-9c77-49a1-bc9a-abc5752a0daf
md"""
!!! note "Proposition 13"
	If $\opA$ with domain $D(\opA)$ is not closed then $\sigma(\opA)=\mathbb{C}$.

> *Proof* by contradiction.
> We assume $\sigma(\opA) \neq \mathbb{C}$ and show that $\opA$ is necessarily closed. 
> In this case, there exists a complex $z \notin \sigma(\opA)$, and we consider a sequence $(x_{n}) \subset D(\opA)$ such that $x_{n} \rightarrow x$ and $\op A x_{n} \rightarrow y$. 
> Since $\opA-z$ is invertible with bounded (and thus continuous) inverse we have
>```math
>\begin{aligned}
>	(\opA-z)^{-1} \opA x_{n} & =(\opA-z)^{-1}(\opA-z+z) x_{n} \\
>	& =x_{n}+z(\opA-z)^{-1} x_{n}
>\end{aligned}
>```
>Going to the limit we obtain
>```math
> 	(\opA-z)^{-1} y=x+z(\opA-z)^{-1} x
>```
> making use of $(\opA-z)^{-1}$ being continuous. 
> Therefore $x=(\opA-z)^{-1}(y-z x)$, which is in $D(\opA)$ since $(\opA-z)^{-1}$ maps $\hilbert \rightarrow D(\opA)$. 
> Self-multiplying by $(\opA-z)$ we obtain
>```math
> 	(\opA x-z) x=y-z x  \iff  \opA x=y 
> ```
> This contradicts our original statement. $\hspace{8cm} \square$
"""

# ╔═╡ 353bee69-8a2e-4da2-83b6-6a2f30621fd6
md"""
Clearly, this result is in stark contrast with physical intuition, where we seek operators with real spectra, ideally on isolated countable points.

We now illustrate the point that non-closed operators can arise if $D(\opA)$ is chosen too small.


- We consider the momentum operator on $\hilbert=L^{2}(\mathbb{R})$, i.e. $\op P : f \mapsto -i f^{\prime}$, and set $D(\op P)=C_{0}^{\infty}(\mathbb{R}).$ 
  We claim $\op P$ is not closed. 
  To see this consider the Gaussian $f=e^{-x^2} \in H^{1}(\mathbb{R}) \setminus C_{0}^{\infty}(\mathbb{R})$. 
  By density of $C_{0}^{\infty}(\mathbb{R})$ in $H^{1}(\mathbb{R})$ we can construct a sequence $(f_{n})_n \subset C_{0}^{\infty}(\mathbb{R})$ such that $f_{n} \rightarrow f$ and $f_{n}^{\prime} \rightarrow f^{\prime}$. 
  While $(f_{n},-i f_{n}') \in \graph(\op P)$, the limiting pair $\left(f, -i f^{\prime}\right)$ is not. 
  Thus $\op P$ is not closed.


- A natural idea is to attempt the "closure" of on operator by adding elements to $\graph(\opA)$ until it is a closed set.


- Since the set should still be the graph of an operator we need to ensure that the properties of Lemma 12 hold. (1) and (3) are easily kept in such a closure process. 
  However, (2) is less easily conserved.


- In particular, there are operators which cannot be closed. 
  Due to proposition 13 they are in general not of physical interest, and we can disregard such peculiarities.
"""

# ╔═╡ b0072ed0-6d00-457e-94e7-f6a400125d23
md"""
!!! note "Definition (Closable operators)"
	Let $\opA$ be on operator with domain $D(\opA)$. 
	One calls $\opA$ **closable** if it admits at least one closed extension. 
	In that case the closed extension with the smallest domain $D(\bar \opA)$ is denoted $\bar \opA$ and called the **closure** of $\opA$.
"""

# ╔═╡ c2a9e93c-e3ed-4040-a703-5bddedf50cb6
md"""
!!! note "Theorem 14 (Closure of ∂ and Δ in Rᵈ)"
	In $\hilbert =L^{2}  (\mathbb{R}^{d})$, let $\op P_{j}^\text{min}$ be the momentum operator defined by $\op P_{j}^\text{min} f=-i \partial_{x_j} f$ with domain $D(\op P_{j}^\text{min})=C_{0}^{\infty} (\mathbb R^{d})$ for $j=1, \ldots d$.
	$\op P_{j}^{\text{min}}$ is closable with closure $\op P_{j} \equiv \overline{\op P_{j}^\text{min}}$ with domain $D(\op P_{j})=H^{1}(\mathbb{R}^{d})$.

	Let further $\op{L}^\text{min} f=-\laplacian f$ in the domain $D(\op L^\text{min})=C_{0}^{\infty} (\mathbb R^d)$, $\op L^\text{min}$ is a again closable with closure $\op L \equiv \overline{\op L^\text{min}}$ of domain $D(\op {L})=H^{2}(\mathbb{R}^{d} )$.
"""

# ╔═╡ d2b6c0c1-aaf8-439c-a02c-707079f9ece5
md"""

> *Proof.* 
> We start by proving that $\op P_{j}$ is closed. 
> Let thus $\left(f_{n}, \op P_{j} f_{n}\right)=\left(f_{n},-i \partial_{x_{j}} f_{n}\right)$ be a sequence of points on the graph of $\op P_{j}$, which converges  to $(f, g)$ in $\hilbert \times \hilbert$. 
> If we show $f \in D (\op P_j )$ and $g= \op P_{j} f$, then $\op P_{j}$ is indeed closed. 
> For this, we use integration by parts against a test function $\varphi \in C_{0}^{\infty} (\mathbb{R}^{d} )$ and convergence in $L^{2} (\mathbb{R}^{d} )$ :
>```math
>\begin{align}
>i \int_{\mathbb{R}^{d}} \varphi(x) g(x) d x 
>& =\lim _{x \rightarrow \infty} \int_{\mathbb{R}^{d}} \varphi(x) \partial_{x_{j}} f_{n}(x) d x 
> \\
>& =-\lim _{n \rightarrow \infty} \int_{\mathbb{R}^{d}} \partial_{x_{j}} \varphi(x) f_{n}(x) d x 
> \\
>& =- \int_{\mathbb{R}^{d}} \partial_{x_{j}} \varphi(x) f(x) d x .
>\end{align}
>```
>This shows that $i g =\partial_{x_j} f$ in the sense of distributions and in particular (since $g \in L^{2}\left(\mathbb{R}^{d}\right)$) that $\partial_{x_j} f \in L^{2} (\mathbb{R}^{d} )$.
>Therefore $f \in H^{1} (\mathbb{R}^{d} )=D\left(\op P_{j}\right)$ and $g=\op P_{j} f_{n}$ which shows that $\op P_{j}$ is closed. 
>Thus, due to $\op P_{j}^\text{min} \subset \op P_{j}$, $\op P_{j}^\text{min}$ is closable.
>An immediate consequence is that $\overline{\graph(\op P_j^\text{min})} \subset \graph(\op P_j)$.
>
> To show $\overline{\op P_j^\text{min}} = \op P_j$ it remains to show $\graph (\op P_j) \subset \overline{\graph(\op P_j^\text{min})}$.
> For this we rely on the denseness of $C_{0}^{\infty} (\mathbb{R}^{2} )$ in $D (\op P_{j})=H^{1}(\mathbb{R}^{d})$ and in $L^{2}\left(\mathbb{R}^{d}\right)$ : this ensures that for each $(f, -\partial_{x_{j}} f) \in L^{2}(\mathbb{R}^{d}) \times L^{2}(\mathbb{R}^{d})$ one can find a sequence $(f_{n}) \subset C_{0}^{\infty}(\mathbb{R}^{d})$ such that $f_{n} \rightarrow f$ and $\partial_{x_{j}} f_{n} \rightarrow \partial_{x_j} f$ strongly in $L^{2}(\mathbb{R}^{d})$. 
> Thus $\graph (\op P_{j}) \subset \overline{\graph (\op P_{j}^\text{min})}$ and $\overline{\op P_{j}^\text{min}}=\op P_{j}$.
>
> The arguments for $\overline{\op L^\text{min}} = \op L$ are similar, relying again on the denseness of $C_0^\infty(\mathbb R^d)$ in $H^2(\mathbb R^d)$.
> $\hspace{13cm} \square$
"""

# ╔═╡ 8ed5e1b6-a7be-4b59-9921-238a3212cad7
md"""
!!! tip "Remark"
	We notice the natural appearance of Sobolev spaces as the domains of standard operators of quantum mechanics in $L^{2} (\mathbb{R}^{d} )$.
"""

# ╔═╡ d0e94968-d053-454c-9a50-6d1097c264a8
md"""
### Adjoint operators and graph

Looking at the graph of an operator $\opA$ also allows for a different point of view on adjoints.

- For the adjoint $\opA^*$ we would like to have 
  ```math
  	\langle y, \op A x\rangle=\left\langle \op  A^{*} y, x\right\rangle \qquad \forall x \in D(\opA), y \in D(\opA^*) \tag{3}
  ```


- The idea is to define $\opA^*$ with the largest $D\left(\op A^{*}\right)$ possible to satisfy the above equality.


- Introducing the inner product
  ```math
  	\langle (x_{1}, x_{2} ),(y_{1}, y_{2} ) \rangle_{\hilbert \times \hilbert}= \langle x_{1}, y_{1} \rangle + \langle x_{2}, y_{2} \rangle
  ```
  on the space $\hilbert \times \hilbert$, Equation (3) can also be written as
  ```math
  	\langle (y, \opA^{*} y ), (\opA x, -x ) \rangle_{\hilbert \times \hilbert}=0,
  ```
  which suggests defining the graph of $\opA^*$ as the orthogonal complement of the flipped graph of $\graph (\opA)$
  ```math
  	\graph (\opA^{*} )= \{ (\opA x, -x ) \mid x \in D(\opA) \}^{\perp}
  ```
  where by $\perp$ we consider orthogonality with respect to $\langle \cdot, \cdot \rangle_{\hilbert \times \hilbert}$.
"""

# ╔═╡ df3f76f1-116d-4863-8d6a-c32076c01d6b
md"""
- Again we use Lemma 12 to test whether this set is indeed the graph of an operator $\opA^*$ with domain $D(\opA^*)$.

  (i) is satisfied since the orthogonal complement of a set is always a vector space, thus $\graph (\opA^*)$ is necessarily a vector space.

  For (ii) we employ the denseness of $D(\opA)$ in $\hilbert$. 
  Note that $(0, y) \in \graph(\opA^{*} )$ if and only if $\langle x, y\rangle=0$ for all $x \in D(\opA)$ by definition. 
  Therefore consider
  ```math
  	D(\opA)^{\perp}= \{y \in \hilbert \mid (0, y) \in \graph (\opA^{*} ) \} .
  ```
  Since $D(\opA)$ is dense, this implies $D(\opA)^\perp = \{ 0 \}$ and in turn $y=0$ as required.

  The tricky part is (iii) : the denseness of $D (\opA ^{*} )$. 
  Consider $D (\opA^{*} ) ^\perp$. 
  We have $y \in D (\opA^{*} ) ^\perp$ if and only if
  ```math
  \begin{align}
  	{\color{gray}
  	\underbrace{\color{black} (y, 0) \in \graph (\opA^{*} )^{\perp}}
  		_{\substack{
  			\text{Because Eq. (3) is exactly \emph{not} satisfied} \\
  			\text{for those $x$, i.e. we are $\hilbert \times \hilbert$} \\
  			\text{orthogonal to $\graph(\opA^*)$}
  		}}
  	} 
  	&= ( \{ (\opA x, -x ) \mid x \in D(\opA) \}^{\perp} )^{\perp} 
  	\\
  	&=\overline{ \{ (\opA x, -x \mid x \in D(\opA) \} }.
  \end{align}
  ```
  or equivalently if $(0,y) \in \overline{\graph(\opA)}$. In other words we have  found that
  ```math
  	D (\opA^{*} )^{\perp}=\{x \in \hilbert \mid (0, x) \in \overline{\graph(\opA)}\} .
  ```
  That is to say $D (\opA^{*} )$ is dense if and only if $\overline{\graph(\opA)}$ satisfies property (ii) of Lemma 12., i.e. if $\opA$ is closable.
"""

# ╔═╡ 1d249ea8-9b1a-4a47-bb08-9fdb2525ad3c
md"""
To summarise the discussion:

- The adjoint $\opA^{*}$ with dense domain $D (\opA^{*} )$ is only well-defined if $\opA$ is closable.
-  $\opA^{*}$ is always closed.

The next result is no surprise after this discussion:

!!! note "Lemma 15"
	Let $\opA: D(\opA)  \to \hilbert$ be a closable operator. 
	Then it holds that 
	```math
		\left ( \opA^* \right )^* = \overline{\opA}
	```
	and
	```math
		\left ( \left ( \opA ^* \right )^* \right )^* = \opA^* = \overline{\opA^*}.
	```
"""

# ╔═╡ a8ad5b99-db46-423f-8541-755c4a3099c8
md"""
### Symmetry, graphs, and related extensions

- Recall an operator $\opA$ defined on $D(\opA)$ is called **symmetric** if for all  $x, y \in D(\opA)$ we have $\langle x, \opA y\rangle=\langle \opA x, y\rangle$.

- Equivalently one could state that $\opA^*$ needs to be an extension of $\opA$, i.e. $\opA \subset \opA^*$ or yet equivalently $\graph(\opA) \subset \graph (\opA^{*} )$.

- Since $\opA^{*}$ is always closed, $\opA \subset \opA^{*}$ implies that a symmetric $\opA$ is always closable with $\overline \opA= (\opA ^{*} )^{*}$. For a symmetric $\opA$ the operators $\opA ^{*}$ and $(\opA ^{*} )^{*}$ are thus two important closed extensions it always possesses.

- If $\op A \subset \op B$ with $\op B$ a symmetric extension of $\opA$, then
  ```math
  	A \subset \overline{\opA} \subset \op B \subset \op B^{*} \subset \op A^{*}
  ```
  since $\op A \subset \op B$ implies $\op B^{*} \subset \opA^{*}$ by definition of the adjoint.

- To summarize
  - All symmetric extensions of $\opA$ are between $\opA$ and $\opA^*$.
  - All closed and symmetric extensions of $\opA$ are between $\overline \opA$ and $\opA^*$.
"""

# ╔═╡ 77cde55f-7b10-4587-aebd-289eeff3d06d
md"""
The spectrum of a symmetric operator is already rather restricted.
Only 4 cases can arise :

!!! note "Theorem 16 (Spectrum of symmetric operators)"
	Let $\opA$ with $D(\opA)$ be a symmetric, closed operator. 
    Its spectrum is one of 

	- the entire complex plane $\sigma(\opA)=\mathbb{C}$
	- the upper half plane $\sigma(\opA)=\mathbb{C}_{+}= \{ z \in \mathbb{C} \mid \im (z) \geq 0\}$
	- the lower half plane $\sigma(\op{A})=\mathbb{C}_{-}=\{z \in \mathbb{C} \mid \im (z) \leq 0\}$
	- non-empty and included in $\mathbb{R}$ : $\varnothing \neq \sigma(\opA) \subset \mathbb{R}$

	In any case the spectrum never has eigenvalues in $\mathbb{C} \setminus \mathbb{R}$, meaning $\ker (\opA-z)=\{0\}$ if $\im (z) \neq 0$.
	
	Furthermore, if $\opA$ with $D(\opA)$ is symmetric and $z \in \mathbb{C} \setminus \mathbb{R}$ is such that $\opA - z$ is surjective then $\opA$ is closed and $z \in \resolvent(\opA)$, meaning that $\sigma(\opA)$ is restricted to the half-plane not containing $z$.

"""

# ╔═╡ 0d1ee873-cb69-4e9b-ad7c-ccdf36ce4dad
md"""
- Just having the information that $\opA - z$ is surjective for a single $z$ with $\im (z) \neq 0$ is sufficient to exclude a half plane as the possible spectrum and implies $\opA = \bar \opA$.


- Further if $\opA$ symmetric and bounded over all $D(\opA)=\hilbert$ (e.g. finite dimensions), its spectrum is necessarily included in a ball of radius $\| \opA \|$. 
  To see that consider a $z \in \mathbb{C}$ with $|z|>\| \opA \|$ :
  ```math
      \left \| 1+ \frac{1}{z} (\opA - z ) \right \| = \frac{\| \opA \|}{|z|} <1
  ```
  So $-\big(\frac1z (\opA - z )\big)^{-1}= \left[1- \big(1+ \frac1z (\opA - z )\big) \right]^{-1}$ can be
  computed by a convergent Neumann series making $-\frac1z (\opA-z)$ and thus $(\opA-z)$ invertible with bounded inverse, and thus $z \in \resolvent(\opA)$. 
  Therefore, $\sigma(\opA)$ is contained in a ball of radius $\| \opA \|$ leaving only the fourth case, i.e. $\sigma(\opA) \subset(-\| \opA \|,\| \opA \|) \subset \mathbb{R}$.
"""

# ╔═╡ 6dd9f4ed-6c02-4966-8397-5faad13f17d1
md"""
- The proof of Theorem 16 is based on a simple equality : for all $x \in D(\opA)$ and $a, b \in \mathbb{R}$ we have 
  ```math 
  \begin{align}
     \|( \opA -a-i b) x \|^{2} & = \langle (\opA -a-i b) x,  (\opA -a-i b) x\rangle 
     \\
     &  =\| (\opA -a) x \|^{2}+b^{2}\|x\|^{2}-2 b \im ( \langle (\opA-a) x,  x \rangle ) 
     \\
     & =\| (\opA-a) x \|^{2}+b^{2} \| x \|^{2}
  \end{align}
  ```
  where we used that for symmetric operators
  ```math
   \langle( \opA-a) x, x\rangle=\langle x,(\opA-a) x\rangle=\overline{\langle(\opA-a) x, x \rangle}
  ```
  is real.
  
  The resulting relation 
  ```math 
        \|( \opA -a-i b) x\|^{2}
        =\|(\opA -a) x \|^{2}+b^{2}\|x\|^{2} 
        \geq b^{2}\|x\|^{2}
  ```
  yields that, for $b \neq 0$, $(\opA -a-i b) x=0$ implies $x=0$, i.e. that there *cannot be non-real eigenvalues* for a symmetric operator.

  Furthermore, if the inverse $(\opA-a-i b)^{-1}$ exits it is automaticity bounded by
  ```math
  \|(A-a-i b)^{-1} \| \leq \frac{1}{|b|}
  \tag{4}
  ```

  Therefore, for a symmetric operator, the only reason for $a+i b \in \sigma(\opA)$ is that $\opA -a-i b$ is not surjective.

  For the full proof of Theorem 16, see Lewin ([*Thérorie spectrale et mécanique quantique*](https://doi.org/10.1007/978-3-030-93436-1) 2022), Theorem 2.23.
"""

# ╔═╡ 8f7acddd-9192-4b8b-a269-756194e46087
md"""
### Self-adjointess and graph

- We recall that an operator $\opA$ with $D(\opA) \subset \hilbert$ is *self-adjoint* if $\opA=\opA^{*}$, i.e. if $\opA$ is symmetric (implying $\opA \subset \op{A}^{*}$) and if $D (\opA^{*} )=D(\opA)$


- Recall that for symmetric operators we have the graph inclusion
  ```math
      \graph(\opA)
      = 
      \{(x, \opA x) \in D(\opA) \times \hilbert \} 
      \subset 
      \{(\opA y,-y) \in \hilbert \times D (\opA) \}^{\perp}
      =
      G (\opA^{*} ).
  ```
  With self-adjointness, the inclusion becomes an equality : $G(\opA)=G (\opA^{*})$.


- As a result a symmetric operator is self-adjoint if and only if
  ```math 
  \begin{align}
      \langle x, \opA z\rangle
      &=
      \langle y, z \rangle && &&  \forall z \in D(\opA)  
      \\
      \Rightarrow  x &\in D(\opA) && \text{ and } && \opA x=y.
  \end{align}
  ```

- Note that a self-adjoint operator is never the extension or restriction of a self-adjoint operator. 
  This goes back to our earlier realisation that $\op A \subset \op B$ implies $\op B^{*} \subset \op A^{*}$ thus $\op A \subset \op B \subset \op B^{*} \subset \op A^{*}$. 
  Therefore, $\opA=\opA^{*}$ implies $\opA=\op B$.
  Crucially, one **cannot grow or shrink $D(\opA)$ preserving self-adjointness**.
"""

# ╔═╡ ea2e50cf-fe18-4f3d-aaab-f5bdc2fed9ae
md"""
## Appendix: Spectral characterisation based on Weyl sequences
"""

# ╔═╡ 7d4a2aef-1545-4d1e-94fb-30bbef09f439
md"""
### Definition of Weyl sequences
Another important characterisation method of the members of $\sigma(\opA)$ is based on the convergence of bounded sequences, where we already noted differences between the finite and infinite dimensional case.
To fully appreciate the details we need a few more definitions.
"""

# ╔═╡ a2cd3c3a-86c7-426b-af4c-dc853f769961
md"""

!!! note "Definition (Weak Convergence)"
	Let $\hilbert$ be a Hilbert
	space. A sequence $(\phi_n) \in \hilbert$ is said to *converge weakly*
	against a $\varphi \in \hilbert$ if 
	```math 
	\begin{align}
	    \lim_{n \to \infty} \langle \phi_n , f \rangle &= \langle \varphi , f \rangle  && \forall f \in \mathcal \hilbert.
	\end{align}
	```
	 In this case we also write
	$\phi_n \rightharpoonup \varphi$.


!!! tip "Remark (Strong convergence implies weak convergence)"
	Let
	$\phi_n \to \varphi$ strongly in $\hilbert$.
	Then, it holds
	$\| \phi_n - \varphi \| \to 0$ as $n \to \infty$. 
	Let further $f \in \hilbert$. We obtain 
	```math 
	\lim_{n \to \infty} | \langle \phi_n - \varphi, f \rangle | \leq \lim_{n \to \infty} \| \phi_n - \varphi \| \cdot \| f \| = 0
	```
	Hence, strong convergence implies weak convergence.


!!! tip "Remark"
	If $\phi_n \rightharpoonup \varphi$ weakly and
	$\| \phi_n \| \to \| \varphi \|$ strongly, then $\phi_n \to \varphi$ strongly
	because 
	```math 
	\begin{align}
	        \| \phi_n - \varphi \| ^2 &= \| \phi_n \| ^2 - 2 \langle \varphi, \phi_n \rangle + \| \varphi \|^2
	    
	\end{align}
	```
	 Which gives, as $n \to \infty$ 
	```math 
	\begin{align}
	        2 \| \varphi \| ^2 - 2 \langle \varphi, \varphi \rangle & = 0
	    
	\end{align}
	```




"""

# ╔═╡ 45bb6731-1b91-49e2-9863-8614ccf3d238
md"""
In infinite dimensions closed and bounded sets are no longer necessarily compact.
Thus bounded sequences may not have strongly converging subsequences.
However,

!!! note "Theorem 18" 
	Let $\hilbert$ be a Hilbert space and
	$(\phi_n) \subset \hilbert$ a bounded sequence. Then, there always exists a
	*weakly* convergent subsequence.


!!! warning "Example 1"
	Consider the sequence $(e_i) = (0,0,\dots,0,1,0,0,\dots)$, the
	sequence of unit vectors in $\ell^2 (\mathbb R)$.   
	The
	sequence is bounded, as each element is inside the infinite dimensional unit ball ($\| e_i \| = 1$), and no strongly convergent subsequence exists. 
	However,
	$e_i \rightharpoonup 0$.

With this in mind we return to the characterization of the spectra of self-adjoint operators.

!!! note "Definition (Weyl Sequence)"
	A sequence $(\phi_n) \subset D(\opA)$
	with $\| \phi_n \| = 1$, which satisfies
	$\| (\opA - \lambda) \phi_n \| \to 0$ for a $\lambda \in \mathbb R$ is called a
	*Weyl sequence*.



"""

# ╔═╡ acb94f49-ba8a-4e2d-bf5e-8f85b101c836
md"### Weyl sequences and spectra"

# ╔═╡ 1cbfdb25-84a6-4b83-89a8-574338ef0127
md"""
!!! note "Theorem 19"
	Let $\opA$ be self-adjoint with $D(\opA) \subset \hilbert$ and
	$\lambda \in \mathbb R$. The following are equivalent

	1.  $\lambda \in \sigma(\opA)$
	
	2.  $\inf_{\varphi \in D(\opA), \| \varphi \| = 1} \| (\opA - \lambda) \varphi \| = 0$
	
	3.  There exists a Weyl sequence for $\lambda$.
"""

# ╔═╡ d1a9808e-f197-457e-bb63-aa824ad6733b
md"""
!!! tip "Remark"
	The third statement explains nicely a key difference between finite and infinite dimensions.

	- **In finite dimensions**, the unit sphere is compact, since it is closed and bounded. 
	  Therefore, for each Weyl sequence $\phi_n$ with $\| \phi_n \| = 1$, we can extract a convergent subsequence
	  $\phi_{n_k} \to \varphi \in \hilbert$.

	  Since all operators are bounded / continuous in finite dimensions, $\opA \phi_{n_k} \to \opA \varphi$. 
	  Thus, $\| (\opA-\lambda) \phi_{n_k} \| \to 0$ implies $\opA \varphi = \lambda \varphi$.
	
	  Furthermore, from the triangle inequality 
	  ```math 
	  \begin{align}
	  		\bigl \vert \| \phi_{n_k} \| - \| \phi_{n_k} - \varphi \| \bigr | \leq \| \varphi \| \leq \| \phi_{n_k} \| + \| \phi_{n_k} - \varphi \|
	  	
	  \end{align}
	  ```
	  In addition, since $\phi_{n_k} \to \varphi$ strongly and $\| \phi_n \| = 1$ these three statements imply $\| \varphi \| = 1$. 

	  A Weyl sequence for $\lambda$ thus always yields an eigenpair and $\sigma(\opA)$ consists only of eigenvalues.
	  Moreover *Weyl sequences converge (strongly) to eigenvalues.*


	- **In infinite dimensions**, $\| \phi_n \| = 1$ implies that we have a bounded sequence. 
	  Because of Theorem 2, there exists a bounded subsequence $\phi_{n_k}$ with $\phi_{n_k} \rightharpoonup \varphi$ weakly.
	  Further, the strong convergence $\| (\opA - \lambda) \phi_n \| \to 0$ implies weak convergence. Thus, for any $f \in D(\opA)$, it holds
	
	  ```math 
	  \begin{align}
	          0 &= \lim_{n\to \infty} \langle f, (\opA - \lambda) \phi_n \rangle \\
	          &= \lim_{n \to \infty} \langle (\opA - \lambda) f , \phi_n \rangle \\
	          &= \langle (\opA - \lambda) f, \varphi \rangle
	    
	  \end{align}
	  ```
	  where we used the symmetry of $\opA - \lambda$ and the weak convergence of $\phi_n$. 
	  Therefore, 
	  ```math 
	  \begin{align}
	        \tag{1}
	        \langle \opA f , \varphi \rangle &= \lambda \langle f, \varphi \rangle &&   \forall f \in D(\opA)
	  \end{align}
	  ```
	
	  Next, we argue based on the graph of $\opA$
	  ```math 
	  \begin{align}
	        G(\opA) &= \{ ( \varphi, \opA \varphi ) \in \hilbert \times \hilbert \mid   \varphi \in D(\opA) \}
	  \end{align}
	  ```
	  and the graph of the adjoint, which can be written
	
	  ```math 
	  \begin{align}
	          G(\opA^*) &= \{ ( g, \opA^* g ) \in \hilbert \times \hilbert \mid g \in D(\opA^*) = D(\opA) \} \\
	        &= \{ (g,h) \in \hilbert \times \hilbert \mid \langle \opA f, g \rangle = \langle f, h \rangle \ \forall f \in D(\opA) \}.
	  \end{align}
	  ```
	  Employing (1), we deduce by comparing the expressions that $( \varphi, \lambda \varphi ) \in G(\opA^*)$. 
	  The only element for which this is possible is $(g , \opA ^* g ) = ( \varphi, \lambda \varphi )$. 

	  Hence,
	  ```math 
	  \begin{align}
	  	\lambda \varphi = \opA^* \varphi = \opA \varphi
	  \end{align}
	  ```
	   using the self-adjointness of $\opA$.

	  Consider the case where $\lambda \in \sigma(\opA)$, but $\lambda$ is *not* an eigenvalue. 
	  Then $\mathop{\mathrm{Ker}}(\opA - \lambda) = \{ 0 \}$ and it must hold $\phi_n \rightharpoonup \varphi = 0$. 
	  *Weyl sequences converge weakly to zero.*

"""

# ╔═╡ c5db4017-6cca-4063-b4d7-75f11c198a92
md"""
One consequence of Theorem 19 is :

!!! note "Theorem 20"
	Let $\opA$ be a self-adjoint operator on
	$D(\opA) \subset \hilbert$. Then 
	```math 
	\begin{align}
	        \inf \sigma (\opA) &= \inf_{0 \neq \varphi \in D(\opA)} \frac{\langle \varphi, \opA \varphi \rangle}{\langle \varphi, \varphi \rangle}
	        \\
	        \sup \sigma(\opA) &= \sup_{0 \neq \varphi \in D(\opA)} \frac{\langle \varphi, \opA \varphi \rangle}{\langle \varphi, \varphi \rangle}
	    
	\end{align}
	```

In particular, $\sigma(\opA) \subset [a, \infty )$ with
$a > - \infty$ if and only if
$\langle \varphi, \opA \varphi \rangle \geq a \| \varphi \| ^2$ for all $\varphi \in D(\opA)$.
Such operators are called *bounded from below*.
Similarly, upper semi-bounded operators satisfy $\langle \varphi, \op B  \varphi \rangle \leq b \| \varphi \| ^2$ and satisfy $\sigma (\op B) \subset (- \infty , b]$

"""

# ╔═╡ 7a7d4e4d-b432-4615-9d0a-a7c308d8053d
md"""
Weyl sequences turn out to yield a useful characterization of the spectra of self-adjoint operators, and overview of which is shown here :

Spectrum | Weyl sequence for spectral characterization
---|:---
$\lambda \in \sigma(\opA)$ | $\exists \text{ a Weyl sequence } (\phi_n) \subset D(\opA) \text{ s.t. } \| \phi_n \| = 1 \text{ and } \| (\opA - \lambda) \phi_n \| \to 0$
$\lambda \in \sigma_{\text{ess}} (\opA)$ | $\exists \text{ a Weyl sequence s.t. } \phi_n \rightharpoonup 0 \text{ (weakly)}$
$\lambda \in \sigma_{\text{disc}} (\opA)$ | $\text{\emph{All} Weyl sequences have subsequences } \phi_{n_k} \to \varphi \text{ (strongly)}$
$\lambda \in \sigma_{\text{cont}} (\opA)$ | $\text{\emph{All} Weyl sequences verify } \phi_n \rightharpoonup 0 \text{ (weakly)}$
$\lambda \in \sigma_{\text{p}}(\opA)$ | $\exists \text{ a Weyl sequence with a weak limit different from 0, i.e. } \mathop{\mathrm{Ker}}(\opA - \lambda) \neq 0.$
"""

# ╔═╡ 0d6276aa-d1d2-4047-8c9c-db1a3508e5bd
md"""
Based on the Weyl sequences we can also obtain the alternative definition for the bottom of the essential spectrum:
"""

# ╔═╡ e771a411-d86f-4fae-9415-dba270369e49
md"""
!!! note "Theorem 21"
	Let $\opA$ be a self-adjoint operator with form domain $Q(\opA) \subset \hilbert$, then
	the *bottom of the essential
	spectrum* $\Sigma(\opA)$ is uniquely defined by 
	```math 
	\begin{align}
	        \Sigma(\opA) \coloneqq \min \sigma_{\text{ess}} (\opA) = \min_{\substack{(v_n) \in Q(\opA)^\mathbb N \\ \| v_n \| = 1 \\ v_n \rightharpoonup 0}} \liminf_{n \to \infty} q_\opA(v_n)
	    
	\end{align}
	```
	with the convention $\Sigma(\opA) = \infty$ if and only if $\sigma_{\text{ess}}(\opA) = \varnothing$.
"""

# ╔═╡ f29b8d33-be83-4579-bff9-5d385c2f6a5f
md"""
### Example spectra of operators

Let us use Theorem 19 and Theorem 20 to deduce the spectra for a few self-adjoint operators on $\hilbert = L^2(\mathbb R ^d)$ 

!!! warning "Example 2 (Identity operator)"
	 $\mathop{\mathrm{id}}: \hilbert \to \hilbert$ (i.e. $D(\mathop{\mathrm{id}}) = \hilbert$) is clearly bounded from above and below by 1. 
	Thus, $\sigma(\mathop{\mathrm{id}}) = \{ 1\}$
"""

# ╔═╡ 9b7f9af8-21cb-47f1-bb08-a5c02cd0a21e
md"""
!!! warning "Example 3 (Multiplication by a continuous function)"
	Let
	$V : \mathbb R^d \to \mathbb R, \op V : D(\op V) \to \hilbert$ where $D(\op V) = \{ f \in L^2 (\mathbb R^2) \mid Vf \in L^2 (\mathbb R^2) \}$.
	Again, $\op V$ is bounded from below by $\inf_{x \in \mathbb R^d} V(x)$ and from above by
	$\sup_{x \in \mathbb R^d} V(x)$. Therefore, 
	```math 
	\begin{align}
	        \sigma(V) = \overline{\mathop{\mathrm{range}}(V)} = \left [ \inf_{x\in \mathbb R^3} V, \sup_{x \in \mathbb R^3} V \right ] 
	    
	\end{align}
	```
"""

# ╔═╡ 477bb397-e65c-4dc5-ac92-194f21dfe3eb
md"""
!!! warning "Example 4 (Laplace operator - Δ)"
	We want to show explicitly that $\sigma ( - \laplacian) = [ 0 , \infty )$ for the Laplace operator $- \laplacian$ with
	$D(- \laplacian) = H^2 (\mathbb R^d)$. 
	We already stated this operator to be self-adjoint.
	Using integration by parts, we find
	```math 
	\begin{align}
	        \forall \varphi \in H^2(\mathbb R^d) && \langle \varphi, - \laplacian \varphi \rangle = \langle \nabla \varphi, \nabla \varphi \rangle \geq 0.
	\end{align}
	```
	Thus $\sigma(- \laplacian) \subset [0, \infty)$.

	To show the reverse we take $k_0 \in \mathbb R^d, f \in H^2(\mathbb R^d)$ and define a Weyl sequence 
	```math
		f_n(x) = n^{-d/2} f(x/n) e^{i x \cdot k_0}
	```
	with Fourier transform 
	```math
		\hat f_n(k) = n^{d/2} \hat f (n (k - k_0)).
	```
	Using $| \cdot |$ to denote norms in $\mathbb R^d$ :
	```math
	\begin{align}
	\| (- \laplacian - |k_0| )^2 f_n \| ^2 &= \int_{\mathbb R^d} \bigl | |k|^2 - |k_0|^2 \bigr | ^2 |\hat f_n (k) |^2 dk 
	\\
	&= \int_{\mathbb R^d} \left | \left | k_0 + \frac{p}{n} \right |^2 - |k_0|^2 \right  |^2 |\hat f (p) |^2 dp
	\\
	&=\frac1{n} \int_{\mathbb R^d} \left | 2 p \cdot k_0 + \frac{|p|^2}{n} \right |^2 |\hat f (p) |^2 dp
	\end{align}
	```
	which converges to 0 as $n \to \infty$.
	Therefore (Theorem 3) $|k_0|^2 \in \sigma(- \laplacian).$
	As $k_0$ spans $\mathbb R^d$, $|k_0|^2$ spans $[0, \infty)$, so $[0, \infty) \subset \sigma(- \laplacian)$. 

	Therefore, we have $[0,\infty) \subset \sigma(-\laplacian) \subset [0, \infty)$.
	Hence, $[0, \infty) = \sigma(-\laplacian)$.
"""

# ╔═╡ 66e35560-eaed-4c50-bb66-4974f711c0fe
md"""
## Appendix: Proof outline for Kato-Temple bound

Here we want to provide an outline of the proof of the [Kato-Temple bound for operators, Theorem 11](#Kato-Temple-bound).

Note first, that (as for Matrices) Theorem 11 follows automatically if we can proof the following:
"""

# ╔═╡ da97ef7b-e068-4e56-a375-d1008bd242d5
md"""
!!! note "Theorem 22 (Temple's Inequality)"
	Let
	$\opA, \tilde \varphi, \tilde \lambda, r$ as in
	Theorem 8. Suppose $\alpha, \beta \in \mathbb R$ with
	$\alpha < \tilde \lambda < \beta$ and
	$(\alpha, \beta ) \cap \sigma(\opA) = \{\lambda \}$. Then 
	```math 
	\begin{align}
	        \tilde \lambda - \frac{\| r\|^2}{\beta - \tilde \lambda} \leq \lambda \leq \tilde \lambda + \frac{\| r \|^2}{\alpha - \tilde \lambda}
	    
	\end{align}
	```
"""

# ╔═╡ fdb9a4cd-abcf-4fad-a47c-2f8d5cdd85c5
md"""

To prove this, we need a bit of spectral calculus.
We consider Cauchy's Formula in infinite dimensions:
"""

# ╔═╡ 98c8ec12-86dd-4835-a5a8-cf7cec4adb51
md"""
!!! note "Theorem 23 (Cauchy's Formula)"
	Let $\opA$ with
	$D(\opA) \leq \hilbert$ be a self-adjoint operator. Let $a,b \in \resolvent(\opA) \cap \mathbb R$
	with $a<b$. Then 
	```math 
	\begin{align}
	        \mathbf 1 _{(a,b)} (\opA) = \mathbf 1 _{[a,b]} (\opA) = - \frac1{2 \pi i} \oint_C (\opA-z)^{-1} dz
	    
	\end{align}
	```
	 For all positively oriented (counter-clockwise) contours $C$
	enclosing $[a,b]$ and crossing the real axis at $a$ and $b$.

---
"""

# ╔═╡ 39fb2f57-5338-44ff-8e7d-84d5edfcde18
TikzPicture(L"""
        %real line
        \draw[>=latex,->] (-5,0) -- (5,0) node[above right,blue]{$\sigma(A)$} ;
        \draw (5,0) node[right]{$\mathbb{R}$};

    
        %sigma_p

    
        \foreach \energy  in {-4,-1,2,2.5,3.5}{
            \draw[blue] (\energy,0) node{$\times$};
            }   
    
        \draw[ultra thick,color=blue] (-3,0) -- (-2,0);
        \draw[ultra thick,color=blue] (0,0) -- (1.5,0);

		\draw[rounded corners,purple] (-0.5,0.5) rectangle (2.25,-0.5) {};


        \draw[purple]  -- (0.875,0.5) node{$<$} -- (2.25,0.5) -- (2.25,-0.5) 		node[below]{$b$} --  (0.875,-0.5) node{$>$} -- (-0.5,-0.5) node[below]{$a$};

""",width="20cm",options="scale=1",preamble=raw"\usepackage{amsfonts}")

# ╔═╡ 91281544-95d5-44d4-87d3-4eb4295da72c
md"""

Illustration of Cauchy's formula.

---
"""

# ╔═╡ 0c9917f1-236a-491e-82ed-e30668e44978
md"""
!!! tip "Remarks"
	*  $a,b \notin \sigma(\opA)$
	* Similar to our discussion in the finite-dimensional case $\mathbf 1_{[a,b]} (\opA)$ is a projector.


Of special importance are the *spectral projectors*
```math 
\begin{align}
    P^\opA (\lambda) = \mathbf 1 _{(-\infty , \lambda ]} (\opA)
\end{align}.
```

To gather some understanding, we first study the finite-dimensional case :

"""

# ╔═╡ 169a8014-4e59-491b-99b8-9fd7417ec1f1
md"""
##### Finite dimensions

We study the finite dimensional case with a matrix
$M \in \mathbb R^{d \times d}$ with distinct eigenvalues
$\lambda_1 < \dots < \lambda_m ( m \leq d)$.
Note that we assume single eigenvalues for simplicity.
In this setting
```math 
\begin{align}
    P^M(\lambda) = \mathbf 1 _{(-\infty,\lambda]} (M) = \bigoplus_{\lambda_i \leq \lambda} \mathbf 1 _{\{λ_i\}} (M)
\end{align}
```
where
```math
\im\Big(\mathbf 1 _{\{λ_i\}} (M)\Big) = \mathop{\mathrm{Ker}}(M - \lambda_i)
```
is the eigenspace of eigenvalue $\lambda_i$
and $\mathbf 1 _{\{λ_i\}} (M)$ is the projector into the eigenspace
of $λ_i$.

We notice that $P^M (\lambda)$ is piecewise constant, and $P^M (\lambda) \in \mathbb R^{d \times d}$. 
Therefore, we can compute
its derivatives with respect to $\lambda$ (in the distributional sense) :
```math 
\begin{align}
    \frac{dP^M(\lambda)}{d \lambda} = \sum_{i=1}^d \delta (\lambda_i - \lambda_j) \mathbf 1 _{\{ \lambda_j \}} (M)
\end{align}
```
where $\mathbf 1 _{\{ \lambda_j \}} (M)$ is the eigenspace of eigenvalue $\lambda_j$.
By integration, we obtain in particular
```math
\int_\mathbb R \lambda dP^M(\lambda) = \int_\mathbb R \lambda \frac{dP^M(\lambda)}{d \lambda} d\lambda = \sum_{i=1}^d  \lambda_j \mathbf 1 _{\{ \lambda_j \}} (M) = M.
```

Now, we want to write something similar for $\opA$ in infinite dimensions.
"""

# ╔═╡ f31ef6df-5e6c-4753-8184-09e9cc245016
md"""
##### Infinite dimensions

We fix $v \in D(\opA)$. 
Now the function
```math 
\begin{align}
    P_v(\lambda) = \langle v, \mathbf 1 _{(-\infty, \lambda]} (\opA) v \rangle = \langle v, P^\opA(\lambda) v \rangle
\end{align}
```
turns out to be bounded and increasing, making its
distributional derivative a measure, called the *spectral measure*. With
this, we can write 
```math 
\begin{align}
    \langle v, Av \rangle = \int_\mathbb R \lambda \ d \langle v, P^\opA (\lambda) v \rangle.
\end{align}
```
which is the infinite dimensional form of spectral resolution.
More generally, for any (measurable) function $p$
and $v,w \in D(\opA)$ 
```math 
\begin{align}
    \langle v , p(\opA) w \rangle = \int_\mathbb R p(\lambda) \ d \langle v, P^\opA (\lambda) w \rangle.
\end{align}
```
which can in turn be used to define a self-adjoint operator $p(\opA)$ (functional calculus).
Without going into details with this the following result is plausible :
"""

# ╔═╡ d70acb78-64ae-4d34-88f6-2190b5272996
md"""
!!! note "Lemma 24" 
	If $f$ is a polynomial and $\opA$ is self-adjoint, then 
	```math
		\sigma(f(\opA)) = \{ f(\lambda) \mid \lambda \in \sigma(\opA) \}.
	```
"""

# ╔═╡ 73f548f8-f6e2-4c4d-9f9b-349fceb3993e
md"Based on this Lemme we are finally in a position to give the"

# ╔═╡ 089cbbb0-a73e-48a9-b8ef-b6cd56c1c71c
md"""
> **Proof of Theorem 22.** 
> We note that our assumptions imply $(\lambda , \beta) \cap \sigma(\opA) = \varnothing$ such that, using Lemma 10 and the polynomial $f(x) = (x - \beta) (x - \lambda)$, we obtain that $(\opA - \beta) (\opA - \lambda) \geq 0$, i.e. that the operator $(\opA - \beta) (\opA - \lambda)$ only has non-negative spectrum.
> Therefore 
> ```math
> \begin{align}
> 	0 &\leq \langle \tilde \varphi, (\opA - \beta) (\opA - \lambda) \tilde \varphi \rangle
> 	\\
> 	&= \langle \tilde \varphi, (A - \tilde \lambda  + \tilde \lambda - \beta) (\opA - \tilde \lambda + \tilde \lambda - \lambda) \tilde \varphi \rangle
> 	\\
> 	&= \| r \|^2 + \langle \tilde \varphi, (\opA - \tilde λ) \tilde \varphi \rangle (\tilde λ - λ)
> 	\\
> 	& \qquad + (\tilde λ  + \beta ) \langle \tilde \varphi , (\opA - \tilde λ) \tilde \varphi \rangle
> 	\\
> 	& \qquad + (\tilde λ - \beta ) (\tilde λ - λ)
> 	\\
> 	& = \| r \|^2 + (\tilde λ - \beta ) (\tilde λ - λ )
> \end{align}
> ```
> Dividing by $(\tilde λ - \beta)$ and rearranging yields the first inequality.
> Similarly, $(\opA - \alpha) (\opA -\lambda) \geq 0$ yields the second.
> $\hspace{8cm} \square$
"""

# ╔═╡ 3a4bcd93-36dd-4751-902c-f15b63e25cf9
md"----"

# ╔═╡ c1408ac8-a9c0-4834-bd7e-bb23dab01488
TableOfContents()

# ╔═╡ 848e2e37-3605-4feb-9424-cdb76da54957
let
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	
	Sidebar(Markdown.parse(read("sidebar.md", String)), 690)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.7"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.73"
QuadGK = "~2.11.2"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "90c0329afd19b3381ba51a19a10b0ef777201995"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

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

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

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

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

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

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

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
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

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

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

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

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9dd97171646850ee607593965ce1f55063d8d3f9"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

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
git-tree-sha1 = "3faff84e6f97a7f18e0dd24373daa229fd358db5"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.73"

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
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
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

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

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
git-tree-sha1 = "a136f98cefaf3e2924a66bd75173d1c891ab7453"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.7"

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

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "85c7811eddec9e7f22615371c3cc81a504c508ee"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+2"

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

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

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
# ╟─1f6efd53-2429-46e7-a579-fedd920250ff
# ╟─e80f5b7f-dfa0-4785-86c9-fdc6c387f3ca
# ╟─041e2cfd-0ec9-4766-98dd-7a90b2a4da84
# ╟─897dfb72-71d7-4ee3-a874-8f932fffd261
# ╟─d3ad8a0a-fa8a-4864-ada3-fe237e65c088
# ╟─a8ec0c38-3f60-4713-a1f7-9a8e054a2533
# ╟─c61e470f-cdf7-4c44-bdf2-c16cad11856c
# ╟─b86c0bf2-cea9-4f6d-89bd-30634e0a2523
# ╟─b64cdb0e-3f9a-48d7-a9c6-afebf057e6ef
# ╟─a18579a4-bf14-42ed-b16c-b01f8a8e88de
# ╟─678cfb59-33b7-46ea-bd66-09ebb20e09f1
# ╟─5f979847-44af-4457-9220-6840b3b1bd11
# ╟─15012412-d229-45e4-8b12-0dc89d3baaa2
# ╟─99968d73-d35d-4d8f-b851-1b468c865a31
# ╟─b44590f7-adcd-46d8-8889-71aa0f774116
# ╟─9e4cdc94-8ffa-42de-b541-c8beddec8e50
# ╟─321383d4-7b51-4b39-91cc-97c87c465e1a
# ╟─4ea66162-8da9-4da9-8cae-14322da43c98
# ╟─d332ddd1-942a-41c0-b8d1-cef005693f0d
# ╟─24f5dd6c-3fa5-4486-ad92-c1bbab44b4ce
# ╟─59c06a2d-d980-458f-bf7a-44bb4f4a8a80
# ╟─625f255b-85b4-4ada-90bb-bb4c7d4ac597
# ╟─8db7d0bc-ff6b-4387-8009-d83e74faaad0
# ╟─86e5f05a-6562-48c4-80f4-10c7cee0698e
# ╟─9823dc80-2adb-4e21-9588-fdd7dc1b3545
# ╟─1b98187a-2dde-46b4-8489-3f039d7a7362
# ╟─3ae93aa5-4998-4d33-9126-9ce6c055d84a
# ╟─feb88cc7-66ac-463c-81b4-a223902c5897
# ╟─8d7c4bd7-fc1e-4dc4-869c-ae33a27b015c
# ╟─871ab9cb-8fd3-407d-8a32-15d3c04c1cd5
# ╟─6caaef7c-085f-433e-9a08-16f1b8243759
# ╟─b2af25d8-e00d-4478-8ddc-ca809478b4b1
# ╟─e712629f-eaeb-4fd1-9aa0-fd3ec4765148
# ╟─974ce554-562c-4d32-aab6-29513f2fff5a
# ╟─25804d6a-9eb4-4c24-89cf-96042eac70d8
# ╟─a15f5cc7-753d-4aa9-b8b5-6b701882d310
# ╟─187f55b0-f71d-4f9d-b00b-e68966b7964d
# ╟─392444bd-c154-4088-aa00-2b32018a90d6
# ╟─8cdc2867-259f-4242-88ab-d17948cfc1bb
# ╟─a53f5693-6a82-481c-a1ba-cb10652e96ae
# ╟─b0798ff6-c3e8-4cd4-bff7-c7de1faf1c97
# ╟─5f17289c-0234-444e-9e7b-8d39c2851b96
# ╟─47c76553-a775-4e89-a6bf-abba188e792c
# ╟─6a10bd68-78fd-4323-a225-f545a0ee6b31
# ╟─60f38f19-df1f-4134-9b08-bb17d278ea4b
# ╟─3c0a48e2-fe1b-4150-8435-9f1a28f6340a
# ╟─afe2868e-56ed-4e96-bd1a-3f6546cdd2b8
# ╟─34e2d8bb-1a7f-48fa-b17d-c4a84f76cea0
# ╟─6535886e-94e9-47e3-aeea-d3db829e7566
# ╟─aefb411d-53f2-4571-a8b7-170e5b09240e
# ╟─f33be757-b599-4356-b713-41281124529d
# ╟─86d8f7b6-68f6-42fd-9564-9c8a82a35e70
# ╟─e13f6917-81b8-4b4f-bad9-9b30cbc1db51
# ╟─7bd639d0-0918-4872-b2ef-b23da592e2b8
# ╟─2ae96252-9aa1-4c55-aff5-e74a90472c2c
# ╟─bc90f5ca-2c95-4345-a315-f73f2fbfe48c
# ╟─2d67bf62-0d3a-4b66-8028-53661602719d
# ╟─ab658994-7d4d-45f7-94a5-ebf1773fcc56
# ╟─6d7cf76a-cfc2-4477-ac8b-3a70204286ac
# ╟─fa883a72-9e03-4a66-8625-89273f19073b
# ╟─4c2952d8-2957-42ca-b05a-068941cf5250
# ╟─523daf58-de5b-4392-bb7f-94c1be73136f
# ╟─fd7a9505-fa0b-4f52-ba35-e476f093a7ed
# ╟─3d8f2b32-046e-446f-b28e-feecc53d02cd
# ╟─685f4ded-2a50-4ea6-af79-04192afa2de7
# ╟─58b80515-233c-4ad5-8650-7841c381e186
# ╟─9a1daee7-9ae8-4f19-bb97-5a9ab392c9a4
# ╟─d535075f-f5fa-4c9a-97cb-fda30ce42173
# ╟─c67c8802-cf11-440a-9afe-48f137777f34
# ╟─8f57a7e2-3449-4034-9309-6233c414ac0f
# ╟─0007d42e-5a2b-45fe-9bcc-603e54704d1e
# ╟─aeb26e12-b13d-4aeb-ab6a-75e7885fdb03
# ╟─e07daf16-dfcd-49eb-b821-f7b8a7a5760e
# ╟─cda93cc1-e9ba-43c8-9cdc-18b7de1315dc
# ╟─23d46169-107b-41a7-8ef1-410f63c046bc
# ╟─34fafa7e-3a64-433f-b8a8-7bfaf7228dba
# ╟─ff0ad3c9-f430-4f99-8df5-bec5c6c06e3f
# ╟─5f3954d3-bbe5-474a-834c-42e8603ff4bf
# ╟─a41650f2-2e4d-4c76-abe7-65e8867065e5
# ╟─3c24ac4f-6903-4cfd-93f3-be7af58f1448
# ╟─3b36fc40-9c77-49a1-bc9a-abc5752a0daf
# ╟─353bee69-8a2e-4da2-83b6-6a2f30621fd6
# ╟─b0072ed0-6d00-457e-94e7-f6a400125d23
# ╟─c2a9e93c-e3ed-4040-a703-5bddedf50cb6
# ╟─d2b6c0c1-aaf8-439c-a02c-707079f9ece5
# ╟─8ed5e1b6-a7be-4b59-9921-238a3212cad7
# ╟─d0e94968-d053-454c-9a50-6d1097c264a8
# ╟─df3f76f1-116d-4863-8d6a-c32076c01d6b
# ╟─1d249ea8-9b1a-4a47-bb08-9fdb2525ad3c
# ╟─a8ad5b99-db46-423f-8541-755c4a3099c8
# ╟─77cde55f-7b10-4587-aebd-289eeff3d06d
# ╟─0d1ee873-cb69-4e9b-ad7c-ccdf36ce4dad
# ╟─6dd9f4ed-6c02-4966-8397-5faad13f17d1
# ╟─8f7acddd-9192-4b8b-a269-756194e46087
# ╟─ea2e50cf-fe18-4f3d-aaab-f5bdc2fed9ae
# ╟─7d4a2aef-1545-4d1e-94fb-30bbef09f439
# ╟─a2cd3c3a-86c7-426b-af4c-dc853f769961
# ╟─45bb6731-1b91-49e2-9863-8614ccf3d238
# ╟─acb94f49-ba8a-4e2d-bf5e-8f85b101c836
# ╟─1cbfdb25-84a6-4b83-89a8-574338ef0127
# ╟─d1a9808e-f197-457e-bb63-aa824ad6733b
# ╟─c5db4017-6cca-4063-b4d7-75f11c198a92
# ╟─7a7d4e4d-b432-4615-9d0a-a7c308d8053d
# ╟─0d6276aa-d1d2-4047-8c9c-db1a3508e5bd
# ╟─e771a411-d86f-4fae-9415-dba270369e49
# ╟─f29b8d33-be83-4579-bff9-5d385c2f6a5f
# ╟─9b7f9af8-21cb-47f1-bb08-a5c02cd0a21e
# ╟─477bb397-e65c-4dc5-ac92-194f21dfe3eb
# ╟─66e35560-eaed-4c50-bb66-4974f711c0fe
# ╟─da97ef7b-e068-4e56-a375-d1008bd242d5
# ╟─fdb9a4cd-abcf-4fad-a47c-2f8d5cdd85c5
# ╟─98c8ec12-86dd-4835-a5a8-cf7cec4adb51
# ╟─39fb2f57-5338-44ff-8e7d-84d5edfcde18
# ╟─91281544-95d5-44d4-87d3-4eb4295da72c
# ╟─0c9917f1-236a-491e-82ed-e30668e44978
# ╟─169a8014-4e59-491b-99b8-9fd7417ec1f1
# ╟─f31ef6df-5e6c-4753-8184-09e9cc245016
# ╟─d70acb78-64ae-4d34-88f6-2190b5272996
# ╟─73f548f8-f6e2-4c4d-9f9b-349fceb3993e
# ╟─089cbbb0-a73e-48a9-b8ef-b6cd56c1c71c
# ╟─3a4bcd93-36dd-4751-902c-f15b63e25cf9
# ╟─c1408ac8-a9c0-4834-bd7e-bb23dab01488
# ╟─848e2e37-3605-4feb-9424-cdb76da54957
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
