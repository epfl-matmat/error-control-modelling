### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 8dab1cdf-f700-4dd9-a839-58fe8cfd6c4a
begin
	using PlutoUI
	import TikzPictures.TikzPicture
	using LaTeXStrings
	using HypertextLiteral

	toc = 
# sidebar --- DO NOT TOUCH THIS LINE
Markdown.parse( "**Error control in scientific modeling** 
" * read("sidebar.md",String)) 
# sidebar --- DO NOT TOUCH THIS LINE

# latex macros --- DO NOT TOUCH THIS LINE
include("latex_macros.jl") 
# latex macros --- DO NOT TOUCH THIS LINE
end

# ╔═╡ 1df1157f-86ec-405d-9d5d-d626da09c152
md"""
# Bounds on Eigenvalues

Suppose we have computed an approximate eigenpair $(\tilde \lambda, \tilde u)$ of a Hermitian matrix $A \in \mathbb C^{n \times n}$.
- How do we know our computation is correct, in particular if finite-precision arithmetic is employed ? 
- How do we know how far away we are from the exact solution ?

In one of the exercises we already saw that 
```math
\begin{align}
\lambda_i \leq \| A \|_F  && \forall i = 1, \dots, n
\end{align}
```
but we also noted it to be a rather crude bound.
This does, however, point to the fact that the matrix entries have something to say about the matrix eigenvalues.




"""

# ╔═╡ 7eba1e23-9f38-4478-84f4-e81ef653ba2c
md"""
## Gerschgorin circles

A first tighter bound is given by

"""

# ╔═╡ 1309b18d-e2f0-46fa-8e0c-3a2730d69d20
md"""
!!! note "Theorem 1 (Gerschgorin circles)"
	Any eigenvalue $\lambda$ of a matrix $A$ is located in one of the closed disks of the complex plane centred at $A_{i i}$ and having the radius 
	```math
		\sum_{\substack{j=1 \\ j \neq i}}^{j=n} \left|A_{i j}\right| \text {. }
	```
	In other words,
	```math
		\forall \lambda \in \sigma(A) \quad \exists i \text{ s.t. } \left|\lambda-A_{i i}\right| \leq \sum_{\substack{j=1 \\ j \neq i}}^{j=n}\left|A_{i j}\right| \tag{1}
	```



"""

# ╔═╡ 118355a0-2104-4216-9d11-53079fa2024d
md"""
> *Proof* by contradiction. Assume (1) does not hold. Then there is an eigenvalue $\lambda$ such that, for all $i=1, \dots, n$,
> ```math
>	\left|\lambda-A_{i i}\right| > \sum_{\substack{j=1 \\ j \neq i}}^{j=n} \left|A_{i j}\right| \text {. }  \tag{2}
> ```
> We then write $A-\lambda I=D-\lambda I+H$, with $D=\operatorname{diag}\left(A_{11}, \ldots, A_{i i}, \ldots, A_{n n}\right)$ and $H=A-D$ (i.e. zero on the diagonal). Due to (2), $D - \lambda I$ is invertible, thus
> ```math
> 	A-\lambda I=(D-\lambda I)(I+\underbrace{(D-\lambda I)^{-1} H}_{M})
> ```
> The elements of $M$ are
> ```math
> 	M_{i j}= \begin{cases}0 & \text{if } i=j \\ \frac{A_{i j}}{A_{i i}-\lambda} & \text {otherwise }\end{cases}
> ```
> Due to (2), we thus have $\sum_{j=1}^{n}\left|M_{i j}\right|<1$ $\forall i$, therefore $\|M\|_{\infty}< 1$, and $\|M\|_{2}<\|M\|_{\infty}<1$.
> Since $\|M\|_{2}$ bounds the modulus of the eigenvalues of $M$, $I+M$ is non-singular, which implies $A-\lambda I$ is non-singular, meaning that $\lambda$ is not an eigenvalue of $A$.
> This contradicts our initial statement.
>  $\hspace{8cm} \square$

"""

# ╔═╡ 1bae37c2-946f-467e-9b1d-c3ca5a129b98
md"""
!!! tip "Remark"
	Since the same result holds for $A^T$, we can formulate a version in terms of column sums instead of row sums.
	```math
	\forall \lambda \in \sigma(A) \quad \exists  j  \text { s.t. } \left|\lambda-A_{j j}\right| \leq \sum_{\substack{i=1 \\ i \neq j}}^{i = n} \left|A_{i j}\right|.
	```
"""

# ╔═╡ 9c189c48-497e-43c1-a50e-fc1ebe7f0714
md"""
The disks defined in this theorem are called **Gerschgorin disks**. There are $n$ disks and their union contains the spectrum of $A$.
Gerschgorin disks are particularly useful when the matrix is almost diagonal, i.e. when a diagonalization algorithm is close to convergence.
"""

# ╔═╡ 371399ec-a571-4aed-b910-7307d90edad1
md"""
## The residual & Bauer-Fike bound

Suppose now we have obtained an approximate eigenpair $(\tilde{\lambda}, \tilde{v})$ to $A$. 
We can define the **residual**
```math
r=A \tilde{v}-\tilde{\lambda} \tilde{v}
```
which can be computed and employed as a stopping criterion in iterations.

Our goal now is to relate the residual $r$ to the error on the eigenvalue, $|\lambda-\tilde{\lambda}|$.
"""

# ╔═╡ 5466b137-915b-4118-9dbf-f97750303784
md"""
!!! note "Theorem 2 (Bauer-Fike)"
	Let $A \in \mathbb{C}^{n \times n}$ Hermitian and let $(\tilde{\lambda}, \tilde{v})$ be an approximate eigenpair with $\|\tilde{v}\|_{2}=1$ and residual $r=A \tilde{v}-\tilde{\lambda} \tilde{v}$.
	Then, there exists an eigenvalue $\lambda$ of $A$ such that
	```math
		\left|\lambda- \tilde \lambda \right| \leq\|r\|_{2} .
	```

> *Proof.* 
> If $\tilde{\lambda} \in \sigma(A)$ the result is trivial. 
> Suppose $\tilde{\lambda}$ is not on eigenvalue of $A$. 
> Then $A-\tilde{\lambda} I$ is invertible, thus we write
> ```math
>	\begin{align}
>	\tilde{v} & =(A-\tilde{\lambda} I)^{-1} r \\
>	& = U (D-\tilde{\lambda} I)^{-1} U^{-1} r  \tag{3}
>	\end{align}
> ```
> where in the last step we used that $A$ is Hermitian, thus it can be diagonalised as $A=U D U^{-1}$ with $U$ unitary.
> Taking the 2-norm on both sids of (3) yields
> ```math
>	\begin{aligned}
>	1 =\|\tilde{v}\|_2 & =\left\|U (D-\tilde{\lambda} I)^{-1} U^{-1} r\right\|_2 \\
>	& \leq  \underbrace{\|U\|_{2}}_{=1}\| D -\tilde{\lambda} I ^{-1}\|_{2} 
> {\underbrace{\left\|U^{-1}\right\|_{2}}_{=1}} \| r \|_{2}  \\
>	& =\max _{i=1,\dots,n} \left|\lambda_{i}-\tilde{\lambda}\right|^{-1} \|r\|_{2}
>	\end{aligned}
> ```
> Since $\min _{i} x_{i}=\max x_i^{-1}$, we obtain
> ```math
> 	\min _{i=1, \ldots, n}\left|\lambda_{i}- \tilde \lambda \right| \leq\|r\|_{2}
>```
>as desired. 
> $\hspace{13cm} \square$
"""

# ╔═╡ fd5bcfb5-86c5-4097-a8b0-3f12ef8ed752
md"""
This is a simple way to get general error bound by establishing a **residual-error relationship**, i.e. a relation between the residual as a computable check for convergence and the error of our quantity of interest against the exact result. 
We will note that there is no need to know the exact result !

!!! tip "Remark"
	In general in *a posteriori* error analysis we want to establish relationship
	```math
	\|e\|_{p} \leq  C \|r\|_{q}
	```
	where $e$ is the error against the exact answer, $r$ the residual, and $C$ a known and computable constant. 
	Which norms $p$ and $q$ are the best choice *depends on context*. 
	For example, for measuring the error in the eigenvector we might choose the $\infty$-norm if we are interested in entry-wise error, or the 2-norm if we are interested in the error natural to the vector space $\mathbb{C}^{n}$.
	Note that there is no reason for $q$ or $C$ to be identical in both cases.

This suggests the following important point: error-residual relationships are not unique.
In fact for our case a better bound is the Kato-Temple bound, which we will derive next.

## Kato-Temple bound
"""

# ╔═╡ 54fea52c-6adb-4c83-b380-416bfdaeacec
md"""
!!! note "Lemma 3"
	Let $\tilde{v}$ be an approximate eigenvector of a Hermitian matrix $A$ and $\tilde \lambda$ the associated approximated eigenvalue, calculated from $\tilde{\lambda}=\langle\tilde{v}, A \tilde{v}\rangle=R_{A}(\tilde{v})$. Let $(\alpha, \beta)$ be an interval that contains no eigenvalue of $A$, but let $\tilde{\lambda} \in(\alpha, \beta)$. Then
	```math
	(\beta-\tilde{\lambda})(\tilde{\lambda}-\alpha) \leq\|r\|_{2}^{2} .
	```
"""

# ╔═╡ f172a676-599b-4be7-9b5d-ceba6b1794dc
md"""
> *Proof.*
> First notice $r \perp \tilde{v}$ :
>```math
>	\begin{align*}
>	\langle\tilde{v}, r\rangle & =\langle\tilde{v}, A \tilde{v}-\tilde{\lambda} \tilde{v}\rangle \\
>	& =\langle\tilde{v}, A \tilde{v}-\langle\tilde{v}, A \tilde{v}\rangle \tilde{v}\rangle \\
>  & =\langle\tilde{v}, A \tilde{v}\rangle-\langle\tilde{v}, A \tilde{v}\rangle=0  \tag{4}
>	\end{align*}
>```
> Thus,
>```math
> \begin{align}
> \left\langle(A-\alpha I) \tilde{v}, (A-\beta I) \tilde{v}\right\rangle
> &=\langle(A-\tilde{\lambda} I) \tilde{v} +(\tilde{\lambda}-\alpha) \tilde v, \tilde{v} 
> (A-\tilde{\lambda} I) \tilde{v}+(\tilde{\lambda}-\beta) \tilde{v} \rangle 
> \\
>& =\langle r+(\tilde{\lambda}-\alpha) \tilde{v}, r+(\tilde{\lambda}-\beta) \tilde{v}\rangle 
> \\
> & \stackrel{(4)}{=}\|r\|_{2}^{2}+(\tilde{\lambda}-\alpha)(\tilde{\lambda}-\beta). \tag{5}
>\end{align}
> ```
> Now expand $\tilde{v}$ in an eigenbasis of $A$, i.e. $\tilde{v}=\sum_{i=1}^{n} c_{i} v_{i}$ where $v_i$ is an eigenvector of $A$ with associated eigenvalue $\lambda_i$, which yields for the left-hand side of (5)
>```math
>	\begin{align}
>	\left\langle(A-\alpha I) \tilde v,(A-\beta I) \tilde{v}\right\rangle & =\sum_{i=1}^{n}\left|c_{i}\right|^{2}(\lambda_i -\alpha)(\lambda_i -\beta) \\
>	& \geq 0
>	\end{align}
>```
> since $(\alpha, \beta)$ contains no eigenpairs. 
> Considering the right-hand side of (5), we obtain
> ```math
> \|r\|_{2}^{2}+(\tilde{\lambda}-\alpha)(\tilde{\lambda}-\beta) \geq 0
> ```
> which is the desired result. $\hspace{11cm} \square$

"""

# ╔═╡ 3ec0e31c-3d09-43d1-9d9a-506320f7d964
md"""
!!! note "Theorem 4 (Kato-Temple)"
	Let $\tilde{v}$ be an approximate eigenvector to $A$, $\|\tilde{v}\|=1, \tilde{\lambda}=\langle\tilde{v}, A \tilde{v}\rangle$, and $r=A \tilde{v}-\tilde{\lambda} \tilde{v}.$ Assume we know an interval $(a, b)$ with $\tilde{\lambda} \in(a, b)$ and where $\lambda$ is the only eigenvalue of $A$ in $(a, b)$. 
	Then,
	```math
		\frac{\|r\|_{2}^{2}}{a-\tilde{\lambda}} \leq \tilde{\lambda}-\lambda \leq \frac{\|r\|_{2}^{2}}{b-\tilde{\lambda}}
	```
"""

# ╔═╡ a3e89cc9-55b6-423f-874a-b923026847ea
md"""
> *Proof.*
> Let $\lambda$ be the closest eigenvalue to $\tilde \lambda$. 
> If $\lambda<\tilde{\lambda}$, take $\alpha=\lambda$ and $\beta=b$ in Lemma 3.3 to yield
> ```math
> 	\begin{align}
> 	0 &\leq(b-\tilde{\lambda})(\tilde{\lambda}-\lambda) \leq\|r\|^{2}_2 \\
>  \Rightarrow \quad 0 &\leq \tilde{\lambda}-\lambda \leq \frac{\|r\|^{2}_2}{b-	\tilde{\lambda}}
> \end{align}
> ```
> Otherwise, set $\alpha=a$ and $\beta=\lambda$ to obtain
> ```math
> 	0 \leq \lambda-\tilde{\lambda} \leq \frac{\| r \|^{2}_2}{\tilde{\lambda}-a}
> ```
> Combining both results completes the proof. 
> $\hspace{7cm} \square$

"""

# ╔═╡ f9cdf1eb-e803-48be-ab67-cd81fc33be3e
md"""
!!! note "Corollary 5 (Symmetric version of Kato-Temple)"
	Let $\tilde{v}$ be an approximate eigenvector to $A$, $\|\tilde{v}\|=1, \tilde{\lambda}=\langle\tilde{v}, A \tilde{v}\rangle$ and $r=A \tilde{v}-\tilde{\lambda} \tilde{v}$. Let $\lambda$ be the eigenvalue closest to $\tilde \lambda$. 
	We define the distance of $\lambda$ to the rest of the spectrum as
	```math
		\delta=\min _{i}\left\{ |\lambda_{i}-\tilde{\lambda} |, \lambda_{i} \neq \lambda\right\} .
	```
	 $\delta$ is also sometimes called the **gap**.
	Then,
	```math
		|\tilde{\lambda}-\lambda| \leq \frac{\|r\|_{2}^{2}}{\delta}.
	```

> *Proof.*
> Theorem 3.4 with $a=\tilde{\lambda}-\delta$ and $b=\tilde{\lambda}+\delta$.
"""

# ╔═╡ 4eb5e944-000f-4755-b039-63591d1dab8a
md"""
Kato-Temple type bounds are considerably shaper than Bauer-Fike bounds.
However, $\delta$ is in general not directly computable. 
Instead, we usually seek a computable approximation for it.

- Let us assume a diagonalization routine yields approximate eigenvectors $\tilde{v}_{1}, \ldots, \tilde{v}_{n}$ with eigenvalue approximations $\tilde{\lambda}_{1}, \ldots, \tilde{\lambda}_{n}$ and residuals $\tilde{r}_{1}, \ldots, \tilde{r}_{n}$.
- A natural idea to approximate $\delta$ when targeting the error estimate of $\lambda_{i}$ is to take
  ```math
  \delta_{\text {est}}=\min \left(\left|\tilde{\lambda}_{i-1}-  \tilde{\lambda}_{i}\right|,\left|\tilde{\lambda}_{i}-\tilde{\lambda}_{i+1}\right|\right)
  ```
  This can be a good approximation, however it is not guaranteed that the error $|\lambda_{i}-\tilde{\lambda}_{i} |$ is smaller $\left\|r_{i}\right\|^{2} / \delta_{\text{est}}$ since it may happen that $\delta_{\text {est}} > \delta$.

- Pictorially consider the following situation where we want to bound $|\lambda-\tilde{\lambda}_{i} |$ by Kato-Temple, and where $\delta=\min \left(\left|\tilde{\lambda}_{i}-\lambda_{i+1}\right|,\left|\tilde{\lambda}_{i}-\lambda_{i-1}\right|\right)$.
"""

# ╔═╡ f6e33422-d41e-4444-b08b-b281f5eb0db3
TikzPicture(L"""
        \draw[>=latex, ->] (0,0) -- (10,0) node[right]{$\sigma(A)$};

		\draw[red, dotted] (3,0) -- (3,-1) (5.1,0) -- (5.1,-2) (6.3,-2) -- (6.3,0) ;
		\draw[purple, dotted] (3.4,0) -- (3.4,1) (5.1,0) -- (5.1,1) ;


		\draw (0.5,0) node{$\times$} node[above]{$\lambda_{i-3}$};
		\draw[blue] (0.5 + 0.2,0) node{$\times$} node[below]{$\tilde \lambda_{i-3}$};

		\draw (1.5,0) node{$\times$} node[above]{$\lambda_{i-2}$};
		\draw[blue] (1.5 + 0.2,0) node{$\times$} node[below]{$\tilde \lambda_{i-2}$};

		\draw (3,0) node{$\times$} node[above]{$\lambda_{i-1}$};
		\draw[blue] (3 + 0.4,0) node{$\times$} node[below]{$\tilde \lambda_{i-1}$};

		\draw (4.5,0) node{$\times$} node[above]{$\lambda_{i}$};
		\draw[blue] (4.5 + 0.6,0) node{$\times$} node[below]{$\tilde \lambda_{i}$};

		\draw (6.3,0) node{$\times$} node[above]{$\lambda_{i+1}$};
		\draw[blue] (6.3 + 1.1,0) node{$\times$} node[below]{$\tilde \lambda_{i+1}$};

		\draw (8.5,0) node{$\times$} node[above]{$\lambda_{i+2}$};
		\draw[blue] (8.5 + 0.4,0) node{$\times$} node[below]{$\tilde \lambda_{i+2}$};

		\draw[red] (3,-1) node{|}  -- node[below]{$\tilde \lambda_i - \lambda_{i-1}$} (5.1,-1) node{|} ;

		\draw[red] (5.1,-2) node{|}  -- node[below]{$\tilde \lambda_i - \lambda_{i+1} = \delta$} (6.3,-2) node{|} ;

		\draw[purple] (5.1,1) node{|} -- node[above]{$\tilde \lambda_i - \tilde \lambda_{i-1} = \delta_{\rm{est}} > \delta $} (3.4,1) node{|};


""",width="23cm",options="scale=1",preamble=raw"\usepackage{amsfonts}")

# ╔═╡ c7c14c64-133f-44e2-83ea-abde368de440
md"""
- Instead, we can use Theorem 2 (Bauer-Fike) to save the situation and get a guaranteed error bound !
  For example : 
  ```math
  	\begin{align}
  	\left|\tilde{\lambda}_{i}-\lambda_{i+1}\right| & =\left|\tilde{\lambda}_{i}-\tilde{\lambda}_{i+1}+\tilde{\lambda}_{i+1}-\lambda_{i+1}\right| \\
  	& \geq\left|\tilde{\lambda}_{i}-\tilde{\lambda}_{i+1}\right|-\left|\tilde{\lambda}_{i+1}-\lambda_{i+1}\right| \\
  	& \hspace{-0.5em} \stackrel{\text{Thm 2}}{=}\left|\tilde{\lambda}_{i}-\tilde{\lambda}_{i+1}\right|-\left\|r_{i+1}\right\|_{2}
  	\end{align}
  ```
  and similarly for $\left|\tilde{\lambda}_{i}-\lambda_{i-1}\right|$. 

- With this we obtain a computable lower bound to $\delta$, and thus a guaranteed upper bound to $\left\|r_{i}\right\|^{2} / \delta$.
"""

# ╔═╡ 1b5e1cf6-e041-4311-9d9a-1db6cb628a3e
md"""
!!! note "Theorem 6"
	Let $\tilde{v}$ be an approximate eigenvector to $A$, $\|\tilde{v}\|=1, \tilde{\lambda}=\langle\tilde{v}, A \tilde{v}\rangle$, and $r=A \tilde{v}-\tilde{\lambda} \tilde{v}$. Let $\lambda$ be the eigenvalue closest to $\tilde{\lambda}$ and $\delta=\min _{i} \{|\lambda_{i}-\tilde{\lambda}|, \lambda_{i} \neq \lambda \}$ be the distance from the spectrum (gap). 
	Let $v$ be an eigenvector of $A$ associated with $\lambda$. 
	Recall that $\theta(x, y)$, the angle between two vectors, is defined as 
	```math
		\cos \theta(x, y)=\frac{|\langle x, y\rangle|}{\|x\|\|y\|}.
	```
	Then,
	```math
		\sin \theta(\tilde{v}, v) \leq \frac{\|r\|_{2}}{\delta}
	```
"""

# ╔═╡ 3c3fabfc-f7e1-4e93-acfa-48c7448f968a
md"""
> *Proof.*
> We set $\theta \equiv \theta(v, \tilde{v})$ and write $\tilde{v}=v \cos \theta+w \sin \theta$, where $w \perp v$ is chosen appropriately. 
> We have
> ```math
>	\begin{align}
>	(A-\tilde{\lambda} I) \tilde{v} & =\cos \theta(A-\tilde{\lambda} I) v +\sin \theta(A-\tilde{\lambda} I) w \\
>	& =\cos \theta({ \lambda} -\tilde{\lambda}) v+\sin \theta(A-\tilde{\lambda} I) w
>	\end{align}
> ```
> Note that
>```math
> 	\begin{align}
>	\langle v, (A-\tilde{\lambda} I) w \rangle & =\langle(A-\tilde{\lambda} I) v, w\rangle \\
>	& =(\lambda-\tilde{\lambda})\langle v, w\rangle=0,
>	\end{align}
> ```
> i.e. the two vectors in the previous sum are orthogonal. 
> Therefore
> ```math
>	\begin{align}
>	\|r \|_{2}^{2} & = \| (A-\tilde{\lambda} I ) \tilde{v} \|_{2}^{2} \\
>	& =\cos ^{2} \theta \ |\lambda-\tilde{\lambda}|^{2}+\sin ^{2} \theta \ \|(A-\tilde{\lambda} I) w\|_{2}^{2}
>	\end{align}
> ```
> Hence,
> ```math
> 	\sin ^{2} \theta\|(A-\tilde{\lambda} I) w\|_{2}^{2} \leq\|r\|_{2}^{2}
> ```
> Since $w \perp v$, 
>```math
>	\begin{aligned}
>	\|(A-\tilde{\lambda} I) w\|_{2} & =\left\|\sum_{\substack{i=1 \\
>	\lambda_{i} \neq \lambda}}^{n}\left(\lambda_{i}-\tilde{\lambda}\right) v_{i} v_{i}^{H} w\right\|_{2} \\
>	& \geq \min _{\substack{i=1, \dots, n \\
>	\lambda_i \neq \lambda}}\left|\lambda_{i}-\tilde{\lambda}\right|=\delta
>	\end{aligned}
>```
> This concludes the proof. $\hspace{11cm} \square$
"""

# ╔═╡ 760385db-cc30-40c1-8f66-5692028a96e3
TableOfContents()

# ╔═╡ e5c087db-1df5-4654-a5f3-d4c27a5a9782
begin
	Sidebar(elts...; location="upper right") = @htl("""
	<aside class="sidebar" style='top: 255px;right: 17px;'>$elts</aside>
	
	<style>
	aside.sidebar {
		position: fixed;
		max-width: min(30%, 300px, calc(100vw - 750px));
		padding: 0.4rem;
		border-radius: 10px;
		max-height: calc(100vh - 320px);
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
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
PlutoUI = "~0.7.58"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "ed636307911dd852a4499e803d4a734588679e50"

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

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

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

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

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
version = "0.3.23+2"

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
git-tree-sha1 = "52ff2af32e591541550bd753c0da8b9bc92bb9d9"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.7+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

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
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

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
# ╟─8dab1cdf-f700-4dd9-a839-58fe8cfd6c4a
# ╟─1df1157f-86ec-405d-9d5d-d626da09c152
# ╟─7eba1e23-9f38-4478-84f4-e81ef653ba2c
# ╟─1309b18d-e2f0-46fa-8e0c-3a2730d69d20
# ╟─118355a0-2104-4216-9d11-53079fa2024d
# ╟─1bae37c2-946f-467e-9b1d-c3ca5a129b98
# ╟─9c189c48-497e-43c1-a50e-fc1ebe7f0714
# ╟─371399ec-a571-4aed-b910-7307d90edad1
# ╟─5466b137-915b-4118-9dbf-f97750303784
# ╟─fd5bcfb5-86c5-4097-a8b0-3f12ef8ed752
# ╟─54fea52c-6adb-4c83-b380-416bfdaeacec
# ╟─f172a676-599b-4be7-9b5d-ceba6b1794dc
# ╟─3ec0e31c-3d09-43d1-9d9a-506320f7d964
# ╟─a3e89cc9-55b6-423f-874a-b923026847ea
# ╟─f9cdf1eb-e803-48be-ab67-cd81fc33be3e
# ╟─4eb5e944-000f-4755-b039-63591d1dab8a
# ╟─f6e33422-d41e-4444-b08b-b281f5eb0db3
# ╟─c7c14c64-133f-44e2-83ea-abde368de440
# ╟─1b5e1cf6-e041-4311-9d9a-1db6cb628a3e
# ╟─3c3fabfc-f7e1-4e93-acfa-48c7448f968a
# ╟─760385db-cc30-40c1-8f66-5692028a96e3
# ╟─e5c087db-1df5-4654-a5f3-d4c27a5a9782
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
