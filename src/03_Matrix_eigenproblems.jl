### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 317ff75e-d718-11ee-1bbb-1374836729ec
begin
	using PlutoUI
	using HypertextLiteral
	using PlutoTeachingTools
	
	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ cec2f53f-c4a3-4bb8-bc67-99243a3e9d4f
md"""
# Matrix Eigenproblems
"""

# ╔═╡ 9bf047e1-624e-4656-9c4d-d1021617ee9d
md"""
## Preliminaries
(We follow Saad 2011)

!!! note "Definitions (Eigenvalue/vector/space, spectrum, spectral radius)"
	Let $A \in \mathbb C^{n \times n}$ be a square matrix.

	-  $λ \in \mathbb C$ is an *eigenvalue* of $A$ if and only if $A - \lambda I$ is *not* invertible, i.e. if $\ker(A - \lambda I) \neq \{0\}$. 
	
	- A non-zero element $v \in \ker (A- \lambda I)$, i.e. $A v = \lambda v$, is called an *eigenvector*.

	- The space $\eigenspace_A(\lambda) = \ker (A - \lambda I)$ is called *eigenspace*.

	- The set of all eigenvalues $\sigma(A)$ is called the *spectrum*.

	- We can further define the *spectral radius* $\spectralradius(A)$ as $\max_{\lambda \in \sigma(A)} | \lambda |$.

"""

# ╔═╡ 0d33d753-d27b-43a5-aef2-6a020fb7be5b
md"""
!!! note "Definitions (Trace, kernel, image)"
	Let $A \in \mathbb C^{n \times n}$. We further define the following quantities
	```math
	\begin{align}
	&\text{Trace} & \text{tr} (A) &\coloneqq \sum_i A_{ii}
	\\
	&\text{Kernel} & \ker (A) &\coloneqq \{ v \in \mathbb C^n \mid Av = 0\} 
	\\
	&\text{Image} & \im (A) &\coloneqq \{ w \in \mathbb C^n \mid \exists v \in \mathbb C^n \text{ s.t. } Av = w \}
	\end{align}
	```

!!! tip "Remark"
	 $\ker(A)$ and $\im (A)$ are subspaces of $\mathbb C^n$, and
	```math
	\dim ( \ker(A)) + \dim (\im(A)) = n.
	```
"""

# ╔═╡ cac1fb8f-f495-4ddb-ace0-58c17e39ebc3
md"""
!!! note "Theorem 1"
	Let $A \in \mathbb C^{n \times n}$. $A$ is injective if and only if $\ker (A) = \{ 0 \}$.

> *Proof.*
>
>  $\boxed \Rightarrow$ : Let $v \in \ker (A)$. Then
> ```math 
> \begin{align}
> Av = 0 &\Rightarrow A(-v) = - Av = 0 
> \\ &\hspace{-0.8em} \stackrel{\text{injective}}{\Rightarrow} v = -v
> \\ &\Rightarrow v = 0
> \end{align}
> ```
>
>  $\boxed \Leftarrow$ : Let $A v_1 = A v_2$. Then
> ```math
> \begin{align}
> A(v_1 - v_2) = A v_1 - A v_2 = 0 &\Rightarrow v_1 - v_2 \in \ker (A)
> \\&\hspace{-1.4em}  \stackrel{\ker (A) = \{ 0\}}{\Rightarrow} v_1 - v_2 = 0
> \\ &\Rightarrow v_1 = v_2
> \\ & && & \square
> \end{align}
> ```

"""

# ╔═╡ b1507c1f-1f1d-41f9-8b89-5d59a9a90e75
md"""
## Matrix norms

Clearly, $A \in \mathbb C^{n \times n}$ is a vector space.
What are good norms for it ?


"""

# ╔═╡ 88e4ba37-7420-4170-b9bc-3229cb257335
md"""

- Recall the conditions a norm should satisfy, for all $A \in \mathbb C^{n\times n}$ :
  ```math
  \begin{align}
  \| A \| &\geq 0 && \| A \| = 0 \iff A = 0 
  \\
  \| \alpha A \| &= | \alpha| \| A \| && \forall \alpha \in \mathbb C 
  \\
  \| A + B \| &\leq \| A \| + \| B \| && \forall B \in \mathbb C^{n \times n}     \tag{$\triangle$}
  \end{align}
  ```
  The third condition is known as the triangle inequality and sometimes denoted as $\triangle$.


- The $p$-norms of $\mathbb C^n$ ($\| v \|_p = (\sum_{i=1}^n |v_i|^p)^{1/p} \quad \forall v \in \mathbb C^n$) inherit a matrix norm for $p \geq 1$ or $p = \infty$ as 
  ```math
  \| A \|_p = \max_{0 \neq v \in \mathbb C^n} \frac{\| A v \|_p}{\|v\|_p},
  ```
  the so-called *induced matrix norm*.

- They are always consistent with the associated vector norm, i.e. 
  ```math
  \| A v \|_p = \| v \|_p \frac{\| A v \|_p }{\|v\|_p} \leq \| v\|_p \max_{0 \neq v   \in \mathbb C^n} \frac{\| A v \|_p}{\| v \|_p} = \| A \|_p \| v \|_p.
  ```
  Further, we obtain
  ```math
  	\| A B \|_p = \max_{x \neq v \in \mathbb C^n} \frac{\| A B v \|_p}{\| v \|_p} \leq \| A \|_p \max_{0 \neq v \in \mathbb C^n } \frac{\| B v \|_p}{\| v\|_p} = \| A \|_p \| B \|_p
  ```
  Thus,
  ```math
  \begin{align}
  \| A ^2 \|_p &\leq \| A \|_p^2.
  \end{align}
  ```
- More generally,
  ```math
  \begin{align}
  \| A ^k \|_p &\leq \| A \|_p^k.
  \end{align}
  ```
"""

# ╔═╡ e10a7029-411a-4c06-abd7-dc5f835fe7db
md"""
Here are some explicit definitions for induced matrix norms 
```math
\begin{align}
\| A \|_1 &= \max_{j = 1, \dots, n} \sum_{i=1}^n | A_{ij} | \tag{Column-sum norm}
\\
 \| A \|_\infty &= \max_{i = 1, \dots, n} \sum_{j=1}^n |A_{ij}| \tag{Row-sum norm}
\\
\| A\|_2 &= \sqrt{\spectralradius (A A^H)} = \sqrt{\spectralradius(A^H A)} 
\end{align}
```
Here, $A^H$ denotes the Hermitian conjugate, i.e. $A^H = \overline{ A ^T} = (\bar A ) ^T$.
For those more familiar with linear algebra over the reals, Hermitian matrices ($A^H = A$) generalize many properties of real symmetric matrices, which are themselves all Hermitian.
!!! tip "Remark"
	 $\| A \|_2$ is also equal to the largest singular value of $A$. 
	In addition, since $A^HA$ (and $A A ^H$) is a Hermitian positive semi-definite matrix, the spectral radius of $A ^H A$ (and $A A ^H$) is simply the largest eigenvalue of $A^HA$, and we can write 
	```math
	\| A \|_2 = \sqrt{\lambda_{max} ( A ^H A )} = \sqrt{\lambda_{max} ( A A ^H)}
	```
"""

# ╔═╡ 9e95ff3b-ca5e-4a88-a466-bf59b0578bd1
md"""
Another important matrix norm is the *Frobenius norm* 
```math
\| A \|_F = \sqrt{\mathrm{tr} ( A A ^H)} = \sqrt{\sum_{i,j = 1}^n | A _{ij}^2|}
```
Its relevance will become clear once we consider the equivalence of matrix norms, which we treat in the next section.
"""

# ╔═╡ 6aef964d-650c-460c-807d-9afcf6d7475e
md"""
## Equivalence of norms

!!! note "Definition (Norm equivalence)"
	The norms $\| \cdot \|_p$ and $\| \cdot \|_q$ over a space $V$ are equivalent if there exists positive constants $\alpha, \beta > 0$ such that for all $A \in V$ 
	```math
		\alpha \| A \|_p \leq \| A \|_q \leq \beta \| A \|_p
	```

!!! note "Proposition 2 (Equivalence of Frobenius norm and 2-norm)"
	We consider the vector space of matrices $A \in \mathbb C^{n \times n}$.
	The Frobenius norm $\| \cdot \|_F$ and the spectral norm $\| \cdot \|_2$ are equivalent.
	Moreover,
	```math
		\| A \|_2 \stackrel{(a)}{\leq} \| A \|_F \stackrel{(b)}{\leq} \sqrt{n} \| A \|_2
	```

This relationship is relevant as $\| A \|_2$ proves useful in theory but is hard to compute, unlike $\| A \|_F$ for which computation is quite straightforward.

> *Proof.*
> (a) is left as an exercise. For (b), we have
> ```math
> \begin{align}
> \| A \|_F^2 = \text{tr} ( A ^H A ) = \sum_{i=1}^n \lambda_i( A ^H A) \leq \sum_{i=1}^n \lambda_{max} ( A ^H A ) = n \| A \|_2^2
> && \square
> \end{align}
> ```
"""

# ╔═╡ bc5beb41-903f-439f-ac2f-40b523945563
md"""
## Similarity

!!! note "Definition (Similarity and diagonalisability)"
	 $A \in \mathbb C^{n \times n}$ is *similar* to $B \in \mathbb C^{n \times n}$ if there exists an invertible $S \in \mathbb C^{n \times n}$ such that $S^{-1} A S = B$.

	If $S$ is *unitary*, i.e. $S^{-1} = S^H$, then $A$ and $B$ are *unitarily similar*.

	 $A \in \mathbb C^{n \times n}$ is *diagonalisable* if it is similar to a diagonal matrix.

!!! tip "Remark (Eigenvalues of similar matrices)"
	Similar matrices have equal eigenvalues of equal multiplicities.
	If $(\lambda, v)$ is an eigenpair of $A$ and $B = S^{-1} A S$, then $(\lambda, S^{-1} v)$ is an eigenpair of $B$.
	Indeed, since
	```math
		B S^{-1} v = S^{-1} A S S^{-1} v = S^{-1} A v = \lambda S^{-1} v.
	```
	Thus, the eigenvalues of $A$ and $B$ are equal and their multiplicities the same.

!!! tip "Remark"
	If $A$ and $B$ are unitarily similar and $A$ is Hermitian, so is $B$, since
	```math
		B = U^{-1} A U = U^H A U.
	```

!!! tip "Remark"
	 $A \in \mathbb C^{n \times n}$ is diagonalizable if and only is there exists $n$ linearly independent eigenvectors $v_i, i = 1, \dots n$.
	To see this, we collect all eigenvectors as columns
	```math
		S = \begin{bmatrix} v_1 & v_2 & \dots & v_n \\ \downarrow & \downarrow & & \downarrow \end{bmatrix}
	```
	implying $AS = S \Lambda$ with $\Lambda = \text{diag} (\lambda_1, \dots, \lambda_n)$.
	Since $v_i$ are linearly independent, $S$ is invertible, thus $S^{-1} A S = \Lambda$.

!!! tip "Remark (Invariance of trace)"
	The trace is invariant under a similarity transform.
	Indeed, noticing $\text{tr} (ABC) = \text{tr} (BCA)$, we obtain
	```math
		\text{tr} ( S ^{-1} A S ) = \text{tr} (A S S^{-1}) = \text{tr} (A)
	```
"""

# ╔═╡ 4b0bef6c-39da-458a-bd72-5db98b972435
md"""
## Quadratic forms and the spectral theorem

"""

# ╔═╡ 8e87aa27-131a-4fc2-b837-fce855fdc8d4
md"""
!!! note "Definition (Quadratic form)"
	The polynomial $q_A : \mathbb C^n \to \mathbb C$
	```math
	q_A(x) = \sum_{i,j=1}^n \overline{x_i} A_{ij} x_j = \langle x , A x \rangle
	```
	associated to a matrix $A$ is called *quadratic form* of $A$. 
"""

# ╔═╡ 82281d67-5824-4974-a4ae-504f5f91c338
md"""
Based on the quadratic form we can show

!!! note "Lemma 3"
	Every matrix $A \in \mathbb C^{n \times n}$ has at least one eigenvalue $\lambda \in \mathbb C$.

> *Proof.*
> - Consider
> ```math
>  \min_{x \in \mathbb C^n,\, g(x) = 0} q_A(x) \quad \text{with}\quad g(x) = \langle x,x \rangle -1.
> ```
> 
> -  $q_A$ is a polynomial, and thus continuous. 
> - Let $S = \{ x \in \mathbb C^n \mid g(x) = 0 \} \neq \varnothing$, which is closed and bounded.
> - Therefore, $q_A$ must have a minimum over $S$ at $v \in S$ by the extreme value theorem (proved in the appendix).
> - Further, $q_A$ and $g$ are complex (Wirtinger) differentiable. 
> - Thus, by the Lagrange multiplier theorem there exists $\lambda \in \mathbb C$ such that
> ```math
> \begin{align}
> 		&\nabla q_A(v) = \lambda \nabla g(v)
> 		\\
> 		&(\ast) \begin{cases}
> 			v^H A  = \lambda v^H & \qquad \text{(Wirtinger derivative wrt }x) \\ A v = \lambda v & \qquad \text{(Wirtinger derivative wrt } \bar x)
> 		\end{cases}
> \end{align}
> ```
> - Hence, $\lambda$ is an eigenvalue of A. $\hspace{9.0cm} \square$
"""

# ╔═╡ abae949c-3d70-4d16-ac3b-18397f12fa51
md"""
Our special interest are Hermitian matrices, particularly due to the importance of self-adjoint operators in quantum mechanics and other physically relevant problems.
As a reminder, Hermitian matrices satisfy $A^H = A$, i.e. $\langle u, Au \rangle = \langle Au, u \rangle$.
For these, we have two key results :

!!! note "Corollary 4"
	Eigenvalues of Hermitian matrices are real.

> *Proof.* 
> If $A = A^H$ we can rewrite the first condition $(\ast)$ in the above proof as
> ```math
> v^H A^H = λ v^H \quad \Longleftrightarrow \quad A v = \bar \lambda v
> ```
> which in combination with the second condition of $(\ast)$ implies $\lambda = \bar \lambda$.
"""

# ╔═╡ b7ad00f2-4112-458c-bf79-6c78a722fde2
md"""
!!! note "Lemma 5"
	Let $A \in \mathbb C^{n \times n}$ be Hermitian and $W$ an invariant subspace of $A$, i.e. $\forall x \in W, A x \in W.$
	Then, $W^\perp = \{ x \mid \langle x,y \rangle = 0 \quad \forall y \in W\}$ is also an invariant subspace, i.e. if $x \in W^\perp$, then $A x \in W^\perp$.

> *Proof.* 
> Consider $x \in W^\perp$. For all $y \in W$, we have
> $\langle y , A x \rangle = \langle A y , x \rangle = 0$ since $A y \in W$.
> Thus, $A x \in W^\perp$.

"""

# ╔═╡ 7b4df739-334e-4c30-b42e-38395848a927
md"""
!!! note "Theorem 6 (Spectral theorem)"
	If $A \in \mathbb C^{n \times n}$ is a Hermitian matrix, then there exists a decomposition $A = U \Lambda U^H$ into a real diagonal matrix $\Lambda = \text{diag} (\lambda_1, \dots, \lambda_n)$ and a unitary matrix $U = (u_1, \dots, u_n)$ with $A u_i = \lambda_i u_i \quad \forall i = 1, \dots, n.$ 

> *Proof* by induction over $n$.
> 
>  $n = 1$ : Trivial, just take $\Lambda_{11} = A_{11}$ (the only element of $A$), $u_{11} = 1$.
>
>  $n \geq 1$ :
>  - By Lemma 3 and Theorem 4, there exists $\lambda_1 \in \mathbb R$ and a $v_1 \in \mathbb C^n$ such that $A v_1 = \lambda_1 v_1$.
> - Let $W = \text{span} \{v_1\}$ and $\{ v_2 ,\dots, v_n \}$ be an orthonormal basis for $W^\perp$, and collect $V = (v_1, \dots, v_n)$.
> - Then,
> ```math
> V^H A V = \begin{pmatrix} \lambda_1 & C \\ 0 & B \end{pmatrix}
> ```
> where $B \in \mathbb C^{(n-1) \times (n-1)}, C \in \mathbb C^{1 \times (n-1)}$.
> - Since $A$ is Hermitian and $V$ unitary, $V^H A V$ is Hermitian, which implies $B^H = B$ and $C = 0$.
> - By induction hypothesis $B = Q D Q^T$ with $D$ real and diagonal and $Q$ unitary.
>   The decomposition is thus obtained as 
> ```math
> \begin{align}
> 		U &= V \begin{pmatrix} 1 & 0 \\ 0 & Q \end{pmatrix}
> 		&\qquad
> 		\Lambda &= \begin{pmatrix} \lambda & 0 \\ 0 & D \end{pmatrix} && && \square
> \end{align}
> ```
"""

# ╔═╡ 7f03c100-5a14-4a66-b28b-3de2907775d1
md"""
!!! tip "Remark (Normal matrices)"
	The spectral theorem generalizes to normal matrices (with $A A^H = A^H A$), where the key point is that in this case the spectral theorem can be applied to the Hermitian matrix $A A^H$.

!!! tip "Remark (Eigenvectors of real symmetric matrices)"
	The eigenvectors of *real* symmetric matrices may be chosen to be real as well.

"""

# ╔═╡ f03007bf-14b7-4d2f-98cd-cfd17e3a97c8
md"""
## Rayleigh quotient and the min-max principle.

For the existence of eigenvalues in matrix $A$, the quadratic form $q _A(x) = \langle x , Ax \rangle$ played a crucial role. 
Essentially its extrema for $\| x \| = 1$ gave us the existence of eigenvalues (proof of Lemma 3).

Another important and related object is the *Rayleigh quotient* 
```math
R_A(x) = \frac{q_A(x)}{\| x \|^2} = \frac{\langle x , Ax \rangle}{\langle x,x \rangle}.
```
In this section, we will characterize its extrema.

!!! tip "Remark"
	 $R_A(x)$ is bounded in $\mathbb C^n$ since $| R_A(x) | \leq \| A \|_2$.

"""

# ╔═╡ a1caae5f-5ecd-4177-bd39-2c606d068629
md"""
We focus on Hermitian matrices $A \in \mathbb C^{n \times n}.$

- For these, $q_A(x)$ is always real, since $\langle x, Ax \rangle = \overline{\langle A x,x \rangle} = \overline{\langle x , Ax \rangle}$.
  Therefore, $R_A(x)$ is real as well.
- The converse is also true and is the subject of an exercise.

- We also notice that the Rayleigh quotient of an eigenvector $v_k$ gives the corresponding eigenvalue $\lambda_k$. 
  Indeed,
  ```math
  \frac{\langle v_k, A v_k \rangle}{\langle v_k, v_k \rangle} = \frac{\langle v_k,   \lambda_k v_k \rangle}{\langle v_k, v_k \rangle} = \lambda_k \frac{\langle v_k, v_k \rangle}{\langle v_k, v_k \rangle} = \lambda_k.
  ```

- As the eigenvalues of $A$ are real, we can order them ascendingly as
  ```math
  \lambda_1 \leq \lambda_2 \leq \dots \leq \lambda_n.
  ```
  This is the convention we will usually employ.
  Note that in this convention we enumerate all $n$ eigenvalues, even if they have the same value, i.e. $\lambda_i = \lambda_{i+1}$.
  The number of repeated eigenvalues gives their multiplicity.

- This convention admits the following important characterization :
"""

# ╔═╡ 109fb7ba-a0fd-4926-9b42-151cd98bf731
md"""
!!! note "Theorem 7 (Courant-Fisher min-max principle)"
	The $k$-th eigenvalue (in ascending order, counting multiplicities) of a Hermitian matrix $A \in \mathbb C^{n \times n}$ is characterized by 
	```math
	\begin{align}
		\lambda_k &= \max_{\substack{S \subset \mathbb C^n \\ \dim (S) = n - k + 1}} \min_{0 \neq v \in S} R_A(v)
		\\
		\lambda_k &= \min_{\substack{S \subset \mathbb C^n \\ \dim (S) = k}} \max_{0 \neq v \in S} R_A(v)
	\end{align}
	```
	where $S$ is an $n -k +1$-dimensional or a $k$-dimensional subspace of $\mathbb C^n$, respectively.
"""

# ╔═╡ cf4669a6-5347-4c3b-b284-28254259032e
md"""
> *Proof.*
> We start by proving the first statement.
> Since $A$ is Hermitian it is unitarily similar to a diagonal matrix of its eigenvalues 
> ```math 
>	\Lambda = \textrm{diag} (\lambda_1, \dots, \lambda_n) = U^{-1} A U
> ```
> with $U$ containing the eigenvectors as columns.
> 
> Let $S_k$ denote the span of the first $k$ eigenvectors in $U$, i.e. $S_k = \textrm{span} (v_1, \dots, v_k)$.
> Since the eigenvectors are linearly independent, $\dim S_k = k$.
>
> We consider an *arbitrary* $(n-k+1)$-dimensional subspace $S$.
> Since $\dim (S) = n-k+1$ and $\dim S_k = k$ thus the sum exceeding $n$,
> their intersection cannot reduce to $\{ 0 \}$, i.e. $S \cap S_k \neq \{ 0 \}$.
> Therefore, there exists $x \in S \cap S_k$ with $x \neq 0$.
> We write
> ```math
> x = \sum_{i=1}^k c_i v_i
> ```
> such that
> ```math
> 	R_A(x) = \frac{\langle x, Ax \rangle}{\langle x,x \rangle} = \frac{\sum_{i=1}^k |c_i|^2 \lambda_i}{\sum_{i=1}^k |c_i|^2} \leq \lambda_k \frac{\sum_{i=1}^k |c_i|^2}{\sum_{i=1}^k |c_i|^2}  = \lambda_k
> ```
> since by construction $\lambda_i \leq \lambda_k$ for $i = 1, \dots, k$.
> Thus, denoting the minimum of $R_A$ over subspace $S$ as $\mu(S)$, we obtain 
> ```math
> \tag{a}
> \mu(S) \coloneqq \min_{0 \neq v \in S} R_A(v) \leq R_A(x) \leq \lambda_k.
> ```
>
> On the other hand, now consider the *particular* $(n-k+1)$-dimensional subspace $\tilde S = \textrm{span} (v_k, \dots, v_n)$.
> Similarly to our previous construction, for any $\tilde x = \sum_{i=k}^n \tilde c_i v_i \in \tilde S$ it holds
> ```math
> 		R_A(\tilde x) = \frac{\sum_{i=k}^n \tilde |c_i|^2 \lambda_i}{\sum_{i=k}^n \tilde |c_i|^2} \geq \lambda_k
> ```
> and thus
> ```math
>   \tag{b}
>   \mu(\tilde S) = \min_{0 \neq \tilde x \in \tilde S} R_A(\tilde x) \geq \lambda_k.
> ```
>
> Therefore according to (a) *any* $(n-k+1)$-dimensional subspace $S$ has $\mu(S) \leq \lambda_k$, while according to (b) the special subspace $\tilde{S}$ of this size even has $\mu(\tilde{S}) \geq \lambda_k$, such that
> ```math
> \max_{\substack{S \subset \mathbb C^n \\ \dim (S) = n - k + 1}} \mu(S) = \lambda_k,
> ```
> which proves the first result.
> The second can be proven using an analogous strategy. $\qquad \square$
"""

# ╔═╡ a79a3ad9-f7ec-43f6-ae77-70d18d3bbc40
md"""
!!! tip "Remark (λ₁ and λₙ)"
	As a direct consequence we obtain 
	```math
	\begin{align}
		\lambda_1 &= \max_{\substack{S \subset \mathbb C^n \\ \dim (S) = n}} \min_{0 \neq v \in S} R_A(v) = \min_{0 \neq v \in \mathbb C^n} R_A(v)
		\\
		\lambda_n &= \min_{\substack{S \subset \mathbb C^n \\ \dim (S) = n}} \max_{0 \neq v \in S} R_A(v) = \max_{0 \neq v \in S} R_A(v)
	\end{align}
	```
"""

# ╔═╡ 2dbc7c07-3c25-4968-be48-4bd8a2b0fb8b
md"""
A further characterization is obtained as a corollary

!!! note "Corollary 8"
	The $k$-th eigenvalue (in ascending order, counting multiplicities) of a Hermitian matrix $A \in \mathbb C^{n \times n}$ is given by 
	```math
		\lambda_k = \frac{\langle v_k, A v_k \rangle}{\langle v_k, v_k \rangle} = \min_{\substack{0 \neq x \in \mathbb C^n \\ x \perp S_{k-1}}} R_A(x)
	```
	where $S_{k-1} = \textrm{span} (v_1, \dots, v_{k-1})$.
	The $k$-th eigenvector can be further obtained as the maximizer of $R_A$, orthogonal to all following eigenvectors (with larger eigenvalues).
"""

# ╔═╡ 36589e2a-2a7d-46f0-ba0b-1a6285940941
md"""
Overall, Theorem 7 and the Rayleigh quotient provide the theoretical framework to compute appropriate eigenpairs by minimization in subspaces.
We will explain this algorithmically in the next lecture.

The study of the Rayleigh quotient also provides direct insight into matrix properties, e.g. for a Hermitian matrix to be positive definite it is necessary and sufficient that $R_A(x) > 0$ for all $x \in \mathbb C^n$ with $x \neq 0$.
Similarly, to be positive semidefinite we need $R_A(x) \geq 0$.
"""

# ╔═╡ f9ed8b9f-1be9-465f-98e6-9d1bdd0627b1
md"""
## Appendix : Proof of the extreme value theorem
"""

# ╔═╡ 38c3ca34-7837-4c14-8ea4-a8a5a3daa415
md"""
We start by proving the boundedness theorem.

!!! note "Theorem 9 (Boundedness theorem)"
	If $f : [a,b] \to \mathbb R$ is continuous on $[a,b]$, then it is bounded on $[a,b]$.

> *Proof.* 
> - Suppose $f$ is not bounded by above on $[a,b]$.
> - Then, for every $n \in \mathbb N$ there exists $x_n \in [a,b]$ with $f(x_n) > n$.
>   This defines a sequence $(x_n)$.
> - As $[a,b]$ is bounded, there must be a convening subsequence $(x_{n_k})$ of $(x_n)$ with limit $x$ by the Bolzano-Weirestrass theorem.
> - As $[a,b]$ is closed, we must have $x \in [a,b]$.
> - Further, $f$ is continuous, thus $f(x_{n_k}) \to f(x)$ as $x_{n_k} \to x$.
> - But, $f(x_{n_k}) > n \geq k$, therefore $f(x_{n_k}) \to \infty$ as $x_{n_k} \to x$.
> - This contradicts our initial statement. 
> 
> We can consider an analogous argument for the lower bound. $\hspace{4cm}  \square$
"""

# ╔═╡ 7aaa09c1-db9f-4791-807e-50cc19ffaf8a
md"""
!!! note "Theorem 10 (Extreme value theorem)"
	A continuous real-valued function $f : [a,b] \to \mathbb R$ has a minimal and a maximal value.

> *Proof.*
> - By the boundedness theorem, $f$ is bounded from above.
> - By the completeness of real numbers, a supremum $s$ (least upper bound) of $f$ over $[a,b]$ exists.
> - We want to find $d \in [a,b]$ such that $f(d) = s$.
> 
> - As $s$ is the supremum, $s - 1/n$ with $n \in \mathbb N$ is not an upper bound for $f$, thus there exists $d_n \in [a,b]$ with $s - 1/n < f(d_n)$.
>   This defines a sequence $(d_n)$.
>
> - By construction 
>   ```math
>   		s - \frac1{n} < f(d_n) \leq s.
>   ```
>   Hence, $f(d_n) \to s$ as $n \to \infty$.
>
> - As $(d_n)$ is bounded, we obtain from the Bolzano-Weirestrass theorem and the fact that $[a,b]$ is closed, that there exists a converging subsequence $(d_{n_k}) \to d \in [a,b]$.
> - In addition, $f$ is continuous, so we have $f(d_{n_k}) \to f(d)$ as $k \to \infty$.
>
> - As $f(d_{n_k})$ is a subsequence of $f(d_n)$ and $f(d_n) \to s$ as $n \to \infty$, we must have $s = f(d).$ 
>
> Again, we can construct an analogous argument for the minimum $\hspace{3.5cm} \square$
"""

# ╔═╡ 3e81997b-65d5-4e2b-b2c2-95aec0ea10bb
TableOfContents()

# ╔═╡ 2fbc04ff-66e5-41ef-b366-2e54882ee74a
let
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	
	Sidebar(Markdown.parse(read("sidebar.md", String)), 330)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
HypertextLiteral = "~0.9.5"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.71"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "9cca7d1e0cdc9047e55a067f6d6e20c53054481f"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

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
version = "1.11.0"

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

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

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

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "dacc8be63916b078b592806acd13bb5e5137d7e9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

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

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

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
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─cec2f53f-c4a3-4bb8-bc67-99243a3e9d4f
# ╟─317ff75e-d718-11ee-1bbb-1374836729ec
# ╟─9bf047e1-624e-4656-9c4d-d1021617ee9d
# ╟─0d33d753-d27b-43a5-aef2-6a020fb7be5b
# ╟─cac1fb8f-f495-4ddb-ace0-58c17e39ebc3
# ╟─b1507c1f-1f1d-41f9-8b89-5d59a9a90e75
# ╟─88e4ba37-7420-4170-b9bc-3229cb257335
# ╟─e10a7029-411a-4c06-abd7-dc5f835fe7db
# ╟─9e95ff3b-ca5e-4a88-a466-bf59b0578bd1
# ╟─6aef964d-650c-460c-807d-9afcf6d7475e
# ╟─bc5beb41-903f-439f-ac2f-40b523945563
# ╟─4b0bef6c-39da-458a-bd72-5db98b972435
# ╟─8e87aa27-131a-4fc2-b837-fce855fdc8d4
# ╟─82281d67-5824-4974-a4ae-504f5f91c338
# ╟─abae949c-3d70-4d16-ac3b-18397f12fa51
# ╟─b7ad00f2-4112-458c-bf79-6c78a722fde2
# ╟─7b4df739-334e-4c30-b42e-38395848a927
# ╟─7f03c100-5a14-4a66-b28b-3de2907775d1
# ╟─f03007bf-14b7-4d2f-98cd-cfd17e3a97c8
# ╟─a1caae5f-5ecd-4177-bd39-2c606d068629
# ╟─109fb7ba-a0fd-4926-9b42-151cd98bf731
# ╟─cf4669a6-5347-4c3b-b284-28254259032e
# ╟─a79a3ad9-f7ec-43f6-ae77-70d18d3bbc40
# ╟─2dbc7c07-3c25-4968-be48-4bd8a2b0fb8b
# ╟─36589e2a-2a7d-46f0-ba0b-1a6285940941
# ╟─f9ed8b9f-1be9-465f-98e6-9d1bdd0627b1
# ╟─38c3ca34-7837-4c14-8ea4-a8a5a3daa415
# ╟─7aaa09c1-db9f-4791-807e-50cc19ffaf8a
# ╟─3e81997b-65d5-4e2b-b2c2-95aec0ea10bb
# ╟─2fbc04ff-66e5-41ef-b366-2e54882ee74a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
