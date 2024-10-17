### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ e80f5b7f-dfa0-4785-86c9-fdc6c387f3ca
begin
	using HypertextLiteral
	using PlutoUI
	using PlutoTeachingTools

	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ 1f6efd53-2429-46e7-a579-fedd920250ff
md"""
# Operators
"""

# ╔═╡ 041e2cfd-0ec9-4766-98dd-7a90b2a4da84
md"""
Previously we discussed the essential properties of function spaces, which can be seen as infinite-dimensional generalization of Euclidean space where functions play the role of vectors.

We noted a number of important properties (compactness, completeness, separability), which in infinite dimensions are considerably harder to obtain or verify then in finite dimensions.

Along these lines, linear operators, i.e. objects mapping functions to functions, are the natural generalization of matrices.
Conversely, we can consider matrices as operators between finite dimensional spaces.
Once again, these objects highlight some subtle differences between infinite and finite dimensional spaces.
"""

# ╔═╡ c3fdd2eb-f9c0-49cc-8c56-274971a1fc3b
md"""

## Linear Operators

!!! note "Definition (Operator)"
	Let $V, W$ be $\mathbb{C}$ vector spaces. 
	A map $\opA : V \rightarrow W$ is called **linear operator** or just **operator** if
	```math
	\forall x, y \in V \quad \forall \alpha, \beta \in \mathbb{C} \quad \opA (\alpha x+\beta y)=\alpha \opA (x)+\beta \opA (y)
	```
	
	If $W=\mathbb C$ we can also call the operator a **linear functional**.
"""

# ╔═╡ d3ad8a0a-fa8a-4864-ada3-fe237e65c088
md"""
!!! warning "Examples of operators"
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


"""

# ╔═╡ c5c84bdb-957d-4909-be69-e8bf31051fb6
md"""
!!! tip "Remarks (Properties of matrices that keep holding)"
	Zeros are mapped onto each other, since
	```math
	0_{W}=0_{\mathbb{C}} \opA(x)=\opA\left(0_{\mathbb{C}} x\right)=\opA \left(0_{V}\right) .
	```

!!! tip ""
	For a subspace $\tilde{V} \leq V$ we have that $\opA (\tilde{V}) \leq W$ is a subspace of $W$ since, for all $\opA(x), \opA(y) \in \opA(\tilde{V})$, we have
	```math
	\alpha \opA(x)+\beta \opA(y)= \opA \underbrace{(\alpha x+\beta y)}_{\in \tilde{V}} \in A(\tilde{V})
	```
	which proves closure.

!!! tip ""
	By a similar argument one proves that linear operators conserve linear combinations
	```math
		\opA \left(\sum_{i=1}^{n} \alpha_{i} x_{i}\right)=\sum_{i=1}^{n} \alpha_{i}  \opA \left(x_{i}\right)
	```
	i.e. 
	```math
		\opA  (\operatorname{span} \{x_{1},\dots,x_{n} \} )=\operatorname{span} ( \{\opA (x_{1} ), \dots, \opA(x_{n} ) \} ).
	```

!!! tip ""
	Concatenations of two linear operators $\opA_{1}, \opA_{2}$ remain linear :

	```math
	\begin{align}
	 (\opA_{1} \circ \opA_{2} )(\alpha x+\beta y) & =\opA_{1} (\opA_{2}(\alpha x+\beta y) ) \\
	& = \opA_{1} (\alpha \opA_{2}(x)+\beta \opA_{2}(y) ) \\
	& =\alpha (\opA_1 \circ \opA_{2} )(x)+\beta (\opA_{1} \circ \opA_{2} )(y)
	\end{align}
	```

!!! tip ""
	Kernel and image are subspaces :
	```math
	\begin{align}
	\ker(A) &= \{x \in V: \opA (x)=0\} \leq V \\
	\im (A)  &= \{y \in W: \exists x \in V \text{ s.t. } \opA (x)=y\} \leq W
	\end{align}
	```
"""

# ╔═╡ 8ce37b41-90c1-4a85-b6eb-814e0e51ca6f
md"""
!!! note "Definition (Operator norm and bounded operators)"
	Let $\left(V,\|\cdot\|_{V}\right)$ and $\left(W,\|\cdot\|_{W}\right)$ denote two Banach spaces. 
	A linear operator $\opA: V \rightarrow W$ is called **bounded** if the norm
	```math
	\|\opA \|_{\boundedoperators (V, W)} \coloneqq \sup _{0 \neq x \in V} \frac{\| \opA (v)\|_{W}}{\| v \|_{V}}<\infty
	```
	is bounded. 
	The set $\boundedoperators (V,W)$ of all bounded operators with this norm is a Banach space.

"""

# ╔═╡ 9e4cdc94-8ffa-42de-b541-c8beddec8e50
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

"""

# ╔═╡ 20da7e84-32e7-49d3-97a1-2e0793878ef5
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
"""

# ╔═╡ 0816a190-9344-4a39-941b-5b9b375dc7b4
md"""
> *Proof.*
> 
>  $\boxed \Rightarrow$ If $\opA$ is bounded then for all $x \in V$ there exists a constant $c$ such that $\| \opA x\left\|_{W}=c \right\| x \|_{V}$. Let $(x_{n})_n \subset V$ be a convergent sequence with limit $x \in V$. Then
> ```math
>  	\| \opA x_{n}- \opA x \|_{W}=c \|x_{n}-x \|_{V}
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
> $\hspace{10cm} \square$
"""

# ╔═╡ d332ddd1-942a-41c0-b8d1-cef005693f0d
md"""
!!! tip "Remark (Boundedness in finite dimensions)"
	In the finite dimensional case all linear operators (matrices) are bounded and thus continuous.

!!! tip "Remark (Boundedness of Schrödinger operators)"
	In passing we note that typical Schrödinger operators are not bounded.

"""

# ╔═╡ a2e86a6f-88e7-4cbd-90c6-6d0e20af7869
md"""
In the following we will consider the specific case of linear operators on a Hilbert space $\hilbert$, i.e. operators of the form $\opA: D(\opA) \rightarrow \hilbert$. The **domain** $D(\opA) \leq \hilbert$ is a dense subspace of $V$. Unless otherwise noted the standard choice of the domain is
```math
	D(x)=\{f \in \hilbert \mid \opA f \in  \hilbert \},
```
i.e. the largest possible subspace of $\hilbert$, such that applying the operator does not take us ot of $\hilbert$. Unless otherwise noted we will always employ this domain.

- Other choices of the domain are well possible and sometimes useful. Moreover and perhaps surprisingly it has a decisive influence on the properties of an operator. For example it is related to the operator being self-adjoint or not, as we will discuss. 


- To illustrate the point of allowing a domain to be smaller than or different from the Hilbert space, we consider the Schrödinger operator of an isolated system, which has the form $\opH=-\Delta+ V$, where $V$ is a potential. 
  The natural Hilbert space for this setting is $L^{2}(\mathbb{R}^{d})$. 
  Thus to ensure $\opH f \in L^{2} (\mathbb{R}^{d})$ $\forall f \in D(\opH)$, we need to take $D(\opH)=H^{2} (\mathbb{R}^{d} )$, one of the Sobolev spaces we introduced last time. 
  Clearly both $D(\opH)= \hilbert =L^{2}\left(\mathbb{R}^{d}\right)$ or $D(H)=\hilbert=H^{2}\left(\mathbb{R}^{d}\right)$ are *not* reasonable as in both cases the Laplacian reduces the regularity.


- For convenience and when the underlying Hilbert space is clear, we frequently denote the operator norm as
  ```math
  	\| \opA \|_{op} \equiv \sup _{0 \neq f \in D(\opA)} \frac{\| \opA f\|_{\hilbert}}{\|f\|_{\hilbert}},
  ```
  i.e. the same setup as the matrix norm, but with the Hilbert space norm and domain inserted.
"""

# ╔═╡ 85df6f49-6935-45f0-9526-0b2b4075f3af
md"""
## Spectra of operators

Having discussed the basic properties of operators, we now turn our attention towards their spectra.

- Generalizing from matrices, an eigenpair $(\lambda, \varphi) \in \mathbb{C} \times \hilbert$ of the operator $\opA$ satisfies
  ```math
  \begin{align}
  	\opA \varphi=\lambda \varphi & \Leftrightarrow(\opA-\lambda) \varphi=0 
  	\\
  	&\Leftrightarrow  \varphi \in \color{gray} \underbrace{\color{black} (\opA-\lambda)^{-1}\{0\}}_{\ker(\opA-\lambda)}
  \end{align}
  ```
  Where by $\opA-\lambda$ we understand the operator $\opA -\lambda \operatorname{id}_\hilbert$. 

  Again, the **resolvent** $R_{z}(\opA)=(\opA-z)^{-1}$ arises naturally, and its study is closely linked to the spectrum of $\opA$.


- First we note that for eigenvalues $(\opA-\lambda)^{-1}$ *cannot be a bounded operator*. 
  To show this, assume it was bounded, i.e. we had a $C>0$ such that
  ```math
  	\left\|(\opA-\lambda)^{-1} x \right\| \leq C\|x\| \qquad \forall x \in \hilbert
  ```
  in particular for one non-zero $\varphi \in \ker(\opA-\lambda) \leq \hilbert$
  ```math
  	0<\|\varphi\|=\left\|(\opA-\lambda)^{-1} 0\right\| \leqslant C \| 0 \|=0,
  ```
  which is a contradiction.

  This motivates the next definition
"""

# ╔═╡ 59c06a2d-d980-458f-bf7a-44bb4f4a8a80
md"""
!!! note "Definition (Resolvent set)"
	Let $\hilbert$ be a separable Hilbert space, $\opA : D(\opA) \rightarrow \opA$. The **resolvent set** is
	```math
		\resolvent(\opA)=\{z \in \mathbb{C} \mid (A-z) : D(\opA) \rightarrow \hilbert \text{ is invertible } {\color{noteblue} \underbrace{\color{black} \text{with bounded inverse}}_{\text{This is new in $\infty$-dimensions}}} \}
	```	
	where bounded inverse means that the operator
	```math
		(\opA -z)^{-1}: \hilbert \rightarrow D(\opA) \leq \hilbert
	```
	exists and there is a $C>0$ such that
	```math
	\left\|(\opA-z)^{-1} x\right\| \leq C \|x\| \quad \forall x \in \hilbert .
	```
	For $z \in \resolvent(x)$ the resolvent $R_{z}(\opA)=(\opA-z)^{-1}$ exists and is bounded.

"""

# ╔═╡ 10aeb8f3-ecdb-4e7c-baad-caa3e1f8994b
md"""
By construction the set $\resolvent(A)$ contains all $z \in \mathbb{C}$ for which $(\opA-z) x=y$ admits a unique solution $x \in D(\opA)$ for a given $y \in \hilbert$. 
Thus, for $(\opA-z) x=0$, only the trivial solution $x=0$ is possible. 
To obtain eigenvalues we thus have to study the complement. As before

!!! note "Definition (Spectrum)"
	The **spectrum** is $\sigma (\opA) = \mathbb C \setminus \resolvent(\opA)$.
"""

# ╔═╡ 9823dc80-2adb-4e21-9588-fdd7dc1b3545
md"""
From the definition of $\resolvent(\opA)$ there can be three reasons for a value $\lambda \in \mathbb{C}$ to not be in $\rho(\opA)$, namely

1.  $(\opA-z)$ not injective.

2.  $(\opA-z)$ not surjective.

3.  $(\opA-z)^{-1}$ not bounded.

Recall that injectivity and surjectivity together make ($\opA - z$) invertible.

For eigenpairs $(\lambda, \varphi) \in \mathbb{C} \times \hilbert$ one desirable property is that $(\lambda, C \varphi)$ for $0 \neq C \in \mathbb{C}$ is also an eigenpair. 
Thus, we want $(A-z) \varphi=0$ to admit multiple solutions, implying $(\opA-z)$ to be not injective.
"""

# ╔═╡ 5827ce09-7199-4357-a5f4-8ec7ae396fca
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

!!! tip "Remark (Spectra in finite dimensions)"
	For $\hilbert=\mathbb C^{N}$, (3) is always given, and 1. and 2. are equivalent because
	```math
		\operatorname{dim}(\ker(\opA-\lambda))+\operatorname{dim}(\im (\opA-\lambda))=N .
	```
	Therefore $\sigma(\opA)=\sigma_p(\opA) \neq \varnothing$ (see Lemma [Matrix Eigenproblems].3).
"""

# ╔═╡ 6b9eb655-048a-42f6-b50c-d928071636bb
md"""
!!! tip "Remark (Behaviour in infinite dimensions)" 
	- In infinite dimensions injective operators can fail to be surjective. 
	  To illustrate this, let $\hilbert=\ell^{2}(\mathbb{R})$ and define $\opA : \hilbert \to \hilbert$ by $\opA\left(\left(x_{n}\right)_{n}\right)=\left(0, x_{1}, x_{2}, \ldots\right)$ for $\left(x_{n}\right)_{n} \subset \mathbb{R}$ (i.e. we just insert a zero element at the front of the sequence). 
	  $\opA$ is injective since $\opA \left(\left(x_{n}\right)_{n}\right)=0$ implies $\left(x_{n}\right)_{n}=0$, but is not surjective.
	  So $0 \in \sigma(A)$, even though it is *not* an eigenvalue.


	- The reverse (surjective but not injective) can also be constructed.


	- In infinite dimensions $\sigma(\opA)$ can be empty. 
	  Consider $\op P_x f=-i \frac{\partial}{\partial x} f$ on $\hilbert=L^{2}(( 0,1) )$ with $D ( \op P_{x} )=\{f \in H^{2}((0,1)) \text{ with } f(0)=0 \}$. We define
	  ```math
	  	\left(\op S_{z} f\right)(x)=i \int_{0}^{x} e^{i z(x-y)} f(y) d y
	  ```
	  as an operator on $\hilbert$. 
	  We compute $1=(\op P_x -z) \op S_{z}= \op S_{z}( \op P _x -z)$ for all $z \in \mathbb{C}$. 
	  So $\op P_x-z$ is invertible. 
	  Still $\op S_{z}$ is bounded because $\left\|\op S_{z} f\right\|_{L^{2}} \leq \|f\|_{L^{2}}$, so have $z \in \resolvent(\opA)$.

To conclude we state a generalisation of an earlier statement about the series expansion of the resolvent.
"""

# ╔═╡ 07bf8bf2-d2b1-4ed8-b95c-f6c7d4cc6e45
md"""
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
"""

# ╔═╡ 12967e6f-5cda-4e16-95ab-2b6ae841d7bc
md"""
> *Proof* (same as the argument before Theorem [pertubation theory].1).
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

"""

# ╔═╡ b2af25d8-e00d-4478-8ddc-ca809478b4b1
md"""
## Self-adjoint operators

As discussed in principle multiple choices for the operator domain $D(\opA)$ can sometimes be reasonable and thus seems to leave considerable ambiguity.

- As the Schrödinger operator example shows, there is clearly some form of a maximal domain, since we clearly need $\varphi \in D(\opH)=H^{2} (\mathbb R^{d})$ to assure $\opH \varphi \in L^{2} (\mathbb{R}^{d} )$. 
  But what about choosing a smaller domain, e.g. $D(H)=C_{0}^{\infty} (\mathbb{R}^{d} )$ ?


- It turns out (detailed discussion in the appendix) that the choice of a sufficiently large domain $D(\opA)$ is required for the operator to be *closed*, which is a slightly weaker form of being continuous. 
  Being closed is necessary for being self-adjoint and moreover if $\opA$ is not closed, then $\sigma(A)=\mathbb{C}$, which is not useful for a physical model.


- Roughly speaking we therefore want $D(\opA)$ to be exactly large enough to get self-adjointness.


Without going into all such detail we will now sketch the definition of a self-adjoint operator and refer to the appendix for more details.
"""

# ╔═╡ 9dd77c43-3d86-4ab0-9ca3-6389a32ea739
md"""

We first state the definition of the (Hilbert-)**adjoint operator**. 
Given $\opA: D(\opA) \rightarrow \hilbert$ with $D(\opA) \subset \hilbert$, the adjoint operator $\opA ^*: D (\opA^{*} ) \rightarrow \hilbert$ is defined by $y \mapsto \opA^{*} y$, where $\opA^{*} y \in \hilbert$ is the unique element satisfying
```math
\left\langle x, \opA^{*} y\right\rangle=\langle\opA x, y\rangle \qquad \forall x \in D(\opA^*), \forall y \in D\left(\opA^{*}\right)
```
Uniqueness is guaranteed by the Riesz representation theorem.

!!! note "Definition (Symmetry)"
	The operator $\opA$ on $\hilbert$ is called **symmetric** if, for all $x, y \in D(\opA) \times D(\opA)$,
	```math
		\langle\opA x, y\rangle= \langle x, \opA y\rangle.
	```

This is a first step towards self-adjoint operators and typically what is called "Hermitian operator" in physics texts.

Comparing this definition of symmetric operators and the adjoint definition above we notice the difference in the sets over which $y$ pans out.
Crucially, *symmetry and self-adjointness are not equivalent concepts*. 
However,

"""

# ╔═╡ 25804d6a-9eb4-4c24-89cf-96042eac70d8
md"""
!!! note "Definition (Self-adjointness)"
	 An operator is self-adjoint if it is symmetric with $D(\opA^*)=D(\opA)$, or simply an operator for which $\opA=\opA^{*}$.

The choice of domain is crucial for self-adjointness. 
One can show that choosing both a smaller or a larger domain $D(\opA)$ destroys the self-adjointness of $\opA$.
"""

# ╔═╡ 8e9a3665-775b-4097-bb00-b6a337176cab
md"""
We now illustrate this concept by providing some examples of self-adjoint operators on $\hilbert=L^{2}(\mathbb R^{d})$. 
We will note that the previously introduced Sobolev spaces arise naturally on the correct domain to ensure self-adjointness of physical operators.

!!! warning "Laplace operator"
	 $\op A=-\laplacian$ with $D(-\Delta)=H^{2}(\mathbb{R}^{d})$ is self-adjoint, but with $D(-\Delta)=C_{0}^{\infty} (\mathbb R^{d})$ it is not.

!!! warning "Momentum operator"
	 $\op P_{j}=-i \frac{\partial}{\partial x_{j}}$ with 
	$D(\op P_{j})=H^{2}(\mathbb{R}^{d})$ is self-adjoint.

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

# ╔═╡ 5f500170-8d73-4d2b-87cb-c706c4966ba5
md"""
We arrive at the following statement, which illustrates once again the difference between finite and infinite dimensions :

!!! note "Proposition 3"
	If a symmetric operator $\opA : \hilbert \to \hilbert$ is defined on the entire $\hilbert = D(\opA)$, then $\opA$ is self-adjoint and bounded. 
"""

# ╔═╡ 392444bd-c154-4088-aa00-2b32018a90d6
md"""
## Summary

- For infinite dimensions $D(\opA)=\hilbert$ in rarely feasible. 
  Choosing $D(\opA)$ sensibly (e.g. Sobolev spaces) is crucial to get self-adjointness and, as we will see, real and physically meaningful spectra.


- Many the special properties we discussed are automatically given in finite dimensions or can be easily checked (e.g. symmetry).


- We have the following implication relationships :
```math
\boxed{
\begin{align}
	\text{self-adjoint} &\Rightarrow \text{closed}
	\\
	& \Rightarrow \text{symmetric}
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

# ╔═╡ 23d46169-107b-41a7-8ef1-410f63c046bc
md"""
## Appendix : In-depth treatment of self-adjointness 

To arrive at the definition of self-adjoint operators we will first look at the graph of operators, which might seem trivial aud unrelated, but is closely related to the question of choosing a good domain $D(\opA)$ for an operator, and the relation between $D(\opA)$, the spectrum $\sigma (\opA)$, and self-adjointness.

"""

# ╔═╡ d6cd60e6-a1d0-402b-98e2-7cb75588ea0e
md"""
### Graph and closure

!!! note "Definition (Graph)"
	Let $A: D(\opA) \rightarrow \hilbert$. The **graph** of the operator $\opA$ is
	```math
	\graph(\opA)=\{(\varphi, \opA \varphi) \mid \varphi \in D(\opA )\} \subset D(\opA) \times \hilbert
	```

Not all subsets of $\hilbert \times \hilbert$ are graphs of operators as the following lemma shows.

!!! note "Lemma 4"
	 A set $\graph \subset \hilbert \times \hilbert$ is the graph of an operator if and only if

	1.  $\graph$ is a vector subspace of $\hilbert \times \hilbert$

	2.  $(0, y) \in \graph$ implies $y=0$ 

	3. The projection $D =\{x \in \hilbert \mid \exists y \in \hilbert$ with $(x, y) \in \graph\}$ is dense in $\hilbert$.
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
		\textsf{bounded / continuous} \Rightarrow \textsf{closed}
	\end{align}
	```

"""

# ╔═╡ 3b36fc40-9c77-49a1-bc9a-abc5752a0daf
md"""
!!! note "Proposition 5"
	If $\opA$ with domain $D(\opA)$ is not closed then $\sigma(\opA)=\mathbb{C}$.

> *Proof* by contradiction.
> We assume $\sigma(\opA) \neq \mathbb{C}$ and show that $\opA$ is necessarily closed. 
> In this case, there exists a $z \notin \sigma(\opA)$, and we consider a sequence $(x_{n}) \subset D(\opA)$ such that $x_{n} \rightarrow x$ and $\op A x_{n} \rightarrow y$. 
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
  By density of $C_{0}^{\infty}(\mathbb{R})$ in $H^{2}(\mathbb{R})$ we can construct a sequence $(f_{n})_n \subset C_{0}^{\infty}(\mathbb{R})$ such that $f_{n} \rightarrow f$ and $f_{n}^{\prime} \rightarrow f^{\prime}$. 
  While $(f_{n},-i f_{n}') \in \graph(\op P)$, the limiting pair $\left(f, -i f^{\prime}\right)$ is not. 
  Thus $\op P$ is not closed.


- A natural idea is to attempt the "closure" of on operator by adding elements to $\graph(\opA)$ until it is a closed set.


- Since the set should still be the graph of an operator we need to ensure that the properties of Lemma 4 hold. (1) and (3) are easily kept in such a closure process. 
  However, (2) is less easily conserved.


- In particular, there are operators which cannot be closed. 
  Due to proposition 5 they are in general not of physical interest, and we can disregard such peculiarities.
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
!!! note "Theorem 6 (Closure of ∂ and Δ in Rᵈ)"
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
- Again we use Lemma 4 to test whether this set is indeed the graph of an operator $\opA^*$ with domain $D(\opA^*)$.

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
  That is to say $D (\opA^{*} )$ is dense if and only if $\overline{\graph(\opA)}$ satisfies property (ii) of Lemma 4., i.e. if $\opA$ is closable.
"""

# ╔═╡ 1d249ea8-9b1a-4a47-bb08-9fdb2525ad3c
md"""
To summarise the discussion:

- The adjoint $\opA^{*}$ with dense domain $D (\opA^{*} )$ is only well-defined if $\opA$ is closable.
-  $\opA^{*}$ is always closed.

The next result is no surprise after this discussion:

!!! note "Lemma 7"
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
  - All closed and symmetric extensions of $\opA$ are between $\overline \opA$ and $\opA$.
"""

# ╔═╡ 77cde55f-7b10-4587-aebd-289eeff3d06d
md"""
The spectrum of a symmetric operator is already rather restricted.
Only 4 cases can arise :

!!! note "Theorem 8 (Spectrum of symmetric operators)"
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
- The proof of Theorem 8 is based on a simple equality : for all $x \in D(\opA)$ and $a, b \in \mathbb{R}$ we have 
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

  For the full proof of Theorem 8, see Lewin ([*Thérorie spectrale et mécanique quantique*](https://doi.org/10.1007/978-3-030-93436-1) 2022), Theorem 2.23.
"""

# ╔═╡ 8f7acddd-9192-4b8b-a269-756194e46087
md"""
### Self-Adjointess and Graph

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

"""

# ╔═╡ ca779a00-60d2-434f-bca6-a387d9566783
md"""
- Note that a self-adjoint operator is never the extension or restriction of a self-adjoint operator. 
  This goes back to our earlier realisation that $\op A \subset \op B$ implies $\op B^{*} \subset \op A^{*}$ thus $\op A \subset \op B \subset \op B^{*} \subset \op A^{*}$. 
  Therefore, $\opA=\opA^{*}$ implies $\opA=\op B$.
  Crucially, one can therefore not *grow or shrink $D(\opA)$ preserving self-adjointness*.
"""

# ╔═╡ a31a2162-c299-4f2f-890e-cc450ebccf8f


# ╔═╡ c1408ac8-a9c0-4834-bd7e-bb23dab01488
TableOfContents()

# ╔═╡ 848e2e37-3605-4feb-9424-cdb76da54957
let
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	
	Sidebar(Markdown.parse(read("sidebar.md", String)), 395)
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

julia_version = "1.10.4"
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
git-tree-sha1 = "5d3a5a206297af3868151bb4a2cf27ebce46f16d"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.33"

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
version = "5.8.0+1"

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
# ╟─1f6efd53-2429-46e7-a579-fedd920250ff
# ╟─e80f5b7f-dfa0-4785-86c9-fdc6c387f3ca
# ╟─041e2cfd-0ec9-4766-98dd-7a90b2a4da84
# ╟─c3fdd2eb-f9c0-49cc-8c56-274971a1fc3b
# ╟─d3ad8a0a-fa8a-4864-ada3-fe237e65c088
# ╟─c5c84bdb-957d-4909-be69-e8bf31051fb6
# ╟─8ce37b41-90c1-4a85-b6eb-814e0e51ca6f
# ╟─9e4cdc94-8ffa-42de-b541-c8beddec8e50
# ╟─20da7e84-32e7-49d3-97a1-2e0793878ef5
# ╟─0816a190-9344-4a39-941b-5b9b375dc7b4
# ╟─d332ddd1-942a-41c0-b8d1-cef005693f0d
# ╟─a2e86a6f-88e7-4cbd-90c6-6d0e20af7869
# ╟─85df6f49-6935-45f0-9526-0b2b4075f3af
# ╟─59c06a2d-d980-458f-bf7a-44bb4f4a8a80
# ╟─10aeb8f3-ecdb-4e7c-baad-caa3e1f8994b
# ╟─9823dc80-2adb-4e21-9588-fdd7dc1b3545
# ╟─5827ce09-7199-4357-a5f4-8ec7ae396fca
# ╟─6b9eb655-048a-42f6-b50c-d928071636bb
# ╟─07bf8bf2-d2b1-4ed8-b95c-f6c7d4cc6e45
# ╟─12967e6f-5cda-4e16-95ab-2b6ae841d7bc
# ╟─b2af25d8-e00d-4478-8ddc-ca809478b4b1
# ╟─9dd77c43-3d86-4ab0-9ca3-6389a32ea739
# ╟─25804d6a-9eb4-4c24-89cf-96042eac70d8
# ╟─8e9a3665-775b-4097-bb00-b6a337176cab
# ╟─5f500170-8d73-4d2b-87cb-c706c4966ba5
# ╟─392444bd-c154-4088-aa00-2b32018a90d6
# ╟─23d46169-107b-41a7-8ef1-410f63c046bc
# ╟─d6cd60e6-a1d0-402b-98e2-7cb75588ea0e
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
# ╟─ca779a00-60d2-434f-bca6-a387d9566783
# ╟─a31a2162-c299-4f2f-890e-cc450ebccf8f
# ╟─c1408ac8-a9c0-4834-bd7e-bb23dab01488
# ╟─848e2e37-3605-4feb-9424-cdb76da54957
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
