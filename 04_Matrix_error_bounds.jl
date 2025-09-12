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

# ╔═╡ 8dab1cdf-f700-4dd9-a839-58fe8cfd6c4a
begin
	using PlutoUI
	import TikzPictures.TikzPicture
	using LaTeXStrings
	using LinearAlgebra
	using PlutoTeachingTools
	using Plots
	using HypertextLiteral

	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ c00bddac-a453-4831-81be-d5814a7e88b9
md"""
# Bounds on Eigenvalues
"""

# ╔═╡ 440fa8bb-c502-4c51-bd80-b0a86d3b8507
md"""
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
		\forall \lambda \in \sigma(A) \quad \exists\, i\, \text{ s.t. }\, \left|\lambda-A_{i i}\right| \leq \sum_{\substack{j=1 \\ j \neq i}}^{j=n}\left|A_{i j}\right| \tag{1}
	```



"""

# ╔═╡ 118355a0-2104-4216-9d11-53079fa2024d
md"""
> *Proof* by contradiction.
> - Assume (1) does not hold. Then there is an eigenvalue $\lambda$ such that, for all $i=1, \dots, n$,
> ```math
>	\left|\lambda-A_{i i}\right| > \sum_{\substack{j=1 \\ j \neq i}}^{j=n} \left|A_{i j}\right| \text {. }  \tag{2}
> ```
> - We then write
>   ```math
>   A-\lambda I=D-\lambda I+H,
>   ```
>   with $D=\operatorname{diag}\left(A_{11}, \ldots, A_{i i}, \ldots, A_{n n}\right)$ and $H=A-D$ (i.e. zero on the diagonal).
> - Due to (2), $D - \lambda I$ is invertible, thus
> ```math
> 	A-\lambda I=(D-\lambda I)\Big(I+\underbrace{(D-\lambda I)^{-1} H}_{M}\Big)
> ```
> - The elements of $M$ are
> ```math
> 	M_{i j}= \begin{cases}0 & \text{if } i=j \\ \frac{A_{i j}}{A_{i i}-\lambda} & \text {otherwise }\end{cases}
> ```
> - Due to (2), we thus have $\sum_{j=1}^{n}\left|M_{i j}\right|<1$ $\forall i$, therefore $\|M\|_{\infty}< 1$, and $\|M\|_{2}<\|M\|_{\infty}<1$.
> - Since $\|M\|_{2}$ bounds the modulus of the eigenvalues of $M$, we have that $I+M$ is non-singular, which implies $A-\lambda I$ is non-singular. Therefore $\lambda$ cannot be an eigenvalue of $A$.
> - This contradicts our initial statement.
>   $\hspace{7cm} \square$

"""

# ╔═╡ 1bae37c2-946f-467e-9b1d-c3ca5a129b98
md"""
!!! tip "Remark"
	Since the same result holds for $A^T$, we can formulate a version in terms of column sums instead of row sums.
	```math
	\forall \lambda \in \sigma(A) \quad \exists  \,j\,  \text { s.t. }\, \left|\lambda-A_{j j}\right| \leq \sum_{\substack{i=1 \\ i \neq j}}^{i = n} \left|A_{i j}\right|.
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
> - If $\tilde{\lambda} \in \sigma(A)$ the result is trivial. 
> - Suppose $\tilde{\lambda}$ is not on eigenvalue of $A$. 
>   Then $A-\tilde{\lambda} I$ is invertible, thus we write
>   ```math
>   \begin{align}
>   \tilde{v} & =(A-\tilde{\lambda} I)^{-1} r \\
>   & = U \, (D-\tilde{\lambda} I)^{-1} \, U^{-1} \, r  \tag{3}
>   \end{align}
>   ```
>   where in the last step we used that $A$ is Hermitian, thus it can be diagonalised as $A=U D U^{-1}$ with $U$ unitary.
> - Taking the 2-norm on both sids of (3) yields
>   ```math
>   \begin{aligned}
>   1 =\|\tilde{v}\|_2 & =\left\|U\,(D-\tilde{\lambda} I)^{-1}\, U^{-1}\, r\right\|_2 \\
>   & \leq  \underbrace{\|U\|_{2}}_{=1} \ \| (D -\tilde{\lambda} I)^{-1}\|_{2} \ 
>   {\underbrace{\left\|U^{-1}\right\|_{2}}_{=1}} \, \| r \|_{2}  \\
>   & =\max _{i=1,\dots,n} \left|\lambda_{i}-\tilde{\lambda}\right|^{-1} \|r\|_{2}
>   \end{aligned}
>   ```
> - Since $\text{argmin}_{i}\, x_{i}=\text{argmax}_i \, x_i^{-1}$, we obtain
>   ```math
>   \min _{i=1, \ldots, n}\left|\lambda_{i}- \tilde \lambda \right| \leq\|r\|_{2}
>   ```
>   as desired. 
>   $\hspace{12cm} \square$
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
	where $e$ is the **error against the exact answer**, $r$ the **residual**,
    and $C$ a known and computable constant. 
	**Which norms** $p$ and $q$ are the best choice *depends on context*. 

	For example, for measuring the error in the eigenvector we might choose the $\infty$-norm if we are interested in **entry-wise error**, or the 2-norm if we are interested in the error **natural to the vector space** $\mathbb{C}^{n}$.
	Note that there is no reason for $q$ or $C$ to be identical in both cases.

This suggests the following important point: error-residual relationships are not unique.
In fact for our case a better bound is the Kato-Temple bound, which we will derive next.

## Kato-Temple bound
"""

# ╔═╡ 54fea52c-6adb-4c83-b380-416bfdaeacec
md"""
!!! note "Lemma 3"
	Let $\tilde{v}$ be an approximate eigenvector of a Hermitian matrix $A$ with (for simplicity) $\|\tilde{v}\| = 1$. Let $\tilde \lambda$ be the associated approximated eigenvalue, calculated from $\tilde{\lambda}=\langle\tilde{v}, A \tilde{v}\rangle=R_{A}(\tilde{v})$. Further let $(\alpha, \beta)$ be an interval that contains no eigenvalue of $A$, but let $\tilde{\lambda} \in(\alpha, \beta)$. Then
	```math
	(\beta-\tilde{\lambda})(\tilde{\lambda}-\alpha) \leq\|r\|_{2}^{2} .
	```
"""

# ╔═╡ f172a676-599b-4be7-9b5d-ceba6b1794dc
md"""
> *Proof.*
> - First notice $r \perp \tilde{v}$ :
>   ```math
>   \begin{align*}
>   \langle\tilde{v}, r\rangle & =\left\langle\tilde{v}, A \tilde{v}-\tilde{\lambda}\, \tilde{v}\right\rangle \\
>   & =\Big\langle\tilde{v}, A \tilde{v}-\left\langle\tilde{v}, A \tilde{v}\right\rangle \tilde{v}\Big\rangle \\
>   & =\langle\tilde{v}, A \tilde{v}\rangle-\langle\tilde{v}, A \tilde{v}\rangle=0  \tag{4}
>   \end{align*}
>   ```
> - Using this result as well as the trick of adding and subtracting $\tilde{\lambda} \tilde{v}$:
>   ```math
>   \begin{align}
>   \Big\langle(A-\alpha I) \, \tilde{v}, (A-\beta I) \, \tilde{v}\Big\rangle
>   &=\Big\langle(A-\tilde{\lambda} I)\, \tilde{v} +(\tilde{\lambda}-\alpha)\, \tilde v,\\
>   &\hspace{50pt}
>   (A-\tilde{\lambda} I) \, \tilde{v}+(\tilde{\lambda}-\beta)\, \tilde{v} \Big\rangle 
>   \\
>   & =\Big\langle r+(\tilde{\lambda}-\alpha) \tilde{v},\ r+(\tilde{\lambda}-\beta) \tilde{v}\Big\rangle 
>   \\
>   & \stackrel{(4)}{=}\|r\|_{2}^{2}+(\tilde{\lambda}-\alpha)(\tilde{\lambda}-\beta). \tag{5}
>   \end{align}
>   ```
> - Now expand $\tilde{v}$ in an eigenbasis of $A$, i.e.
>   ```math
>   \tilde{v}=\sum_{i=1}^{n} c_{i} v_{i}
>   ```
>   where $v_i$ is an eigenvector of $A$ with associated eigenvalue $\lambda_i$.
> - This yields for the left-hand side of (5)
>   ```math
>   \begin{align}
>   \big\langle(A-\alpha I) \, \tilde v,\ (A-\beta I) \, \tilde{v}\big\rangle & =\sum_{i=1}^{n}\left|c_{i}\right|^{2}\,(\lambda_i -\alpha)(\lambda_i -\beta) \\
>   & \geq 0
>   \end{align}
>   ```
>   since the interval $(\alpha, \beta)$ contains no eigenpairs.
> - Considering the right-hand side of (5), we obtain
>   ```math
>   \|r\|_{2}^{2}+(\tilde{\lambda}-\alpha)(\tilde{\lambda}-\beta) \geq 0
>   ```
>   which is the desired result. $\hspace{10cm} \square$

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
	Let $\tilde{v}$ be an approximate eigenvector to $A$, $\|\tilde{v}\|=1, \tilde{\lambda}=\langle\tilde{v}, A \tilde{v}\rangle$ and $r=A \tilde{v}-\tilde{\lambda} \tilde{v}$. Let $\lambda_i$ be the eigenvalue closest to $\tilde \lambda$. 
	We define the distance of $\tilde{\lambda}$ to the rest of the spectrum as
	```math
		\delta=\min _{j,\, j\neq i} |\lambda_{j}-\tilde{\lambda} | .
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
Note that Corollary 5 does not give a useful bound if $\delta = 0$, which can happen for degenerate eigenvalues.

Given a sufficiently good eigenvector $\tilde{v}$ **Kato-Temple-type bounds**
are considerably **sharper** than Bauer-Fike bounds.
Their main complication is that the gap $\delta$ requires access to the **exact** eigenvalue and is thus usually not directly computable.
However, one can usually determine an approximate gap as we will detail below.

- Let us assume a diagonalization routine yields approximate eigenvectors $\tilde{v}_{1}, \ldots, \tilde{v}_{n}$ with eigenvalue approximations $\tilde{\lambda}_{1}, \ldots, \tilde{\lambda}_{n}$ and residuals $\tilde{r}_{1}, \ldots, \tilde{r}_{n}$. We further assume that we have not lost any eigenpair, i.e. $\tilde{λ}_i$ approximates $λ_i$ for all $i = 1, \ldots, n$.
- **A naive idea** to approximate the gap $\delta$ for eigenvalue $\lambda_{i}$ is to take
  ```math
  \delta_{\text {est}}=\min \left(\left|\tilde{\lambda}_{i-1}-  \tilde{\lambda}_{i}\right|,\left|\tilde{\lambda}_{i}-\tilde{\lambda}_{i+1}\right|\right)
  ```
  While this *can be* a good approximation, a disadvantage of this approach is
  that this expression **is not a guaranteed lower bound** to $\delta$.
  In particular it may happen that $\delta_{\text {est}} > \delta$, which implies
  that $\left\|r_{i}\right\|^{2} / \delta_{\text{est}}$
  may be smaller than the actual error ! The Kato-Temple estimate is **no longer a guaranteed bound**.

- To see this pictorially, let us assume for simplicity that all exact eigenvalues are simple. We want to bound $|\lambda-\tilde{\lambda}_{i} |$ by Kato-Temple, and where the exact gap $\delta=\min \left(\left|\tilde{\lambda}_{i}-\lambda_{i+1}\right|,\left|\tilde{\lambda}_{i}-\lambda_{i-1}\right|\right)$. Then it can happen that:
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
- To avoid this problem, we can use Theorem 2 (Bauer-Fike) to obtain a **guaranteed lower bound on the gap**. For example:
  ```math
  	\begin{align}
  	\left|\tilde{\lambda}_{i}-\lambda_{i+1}\right| & =\left|\tilde{\lambda}_{i}-\tilde{\lambda}_{i+1}+\tilde{\lambda}_{i+1}-\lambda_{i+1}\right| \\
  	& \stackrel{\text{rev.}~\Delta}{\geq}\left|\tilde{\lambda}_{i}-\tilde{\lambda}_{i+1}\right|-\left|\tilde{\lambda}_{i+1}-\lambda_{i+1}\right| \\
  	& \hspace{-0.5em} \stackrel{\text{Thm 2}}{\geq}\left|\tilde{\lambda}_{i}-\tilde{\lambda}_{i+1}\right|-\left\|r_{i+1}\right\|_{2}
  	\end{align}
  ```
  and similarly for $\left|\tilde{\lambda}_{i}-\lambda_{i-1}\right|$. 

- This results in computable lower bounds to $\left|\tilde{\lambda}_{i}-\lambda_{i+1}\right|$ and $\left|\tilde{\lambda}_{i}-\lambda_{i-1}\right|$, thus $\delta$.
  Using this gap estimate thus yields a guaranteed upper bound to $\left\|r_{i}\right\|^{2} / \delta$.

- Note that when $\lambda_i$ and $\lambda_{i+1}$ are too close, this estimate can be negative, thus not a valid lower bound for $\delta$.
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

# ╔═╡ 9107b5a1-5768-436d-b705-6b9215121aa0
TODO("Summary on other bounds")

# ╔═╡ bc4e51aa-674e-4a46-a20d-d9bfd82ff108
md"""
## Error bounds showcase
"""

# ╔═╡ 45f6befb-4518-4467-81f3-ecd6da88ccbb
md"""
- M12 = $(@bind M12 PlutoUI.Slider(-1.0:0.001:1.0; default=0.001, show_value=true))
- M13 = $(@bind M13 PlutoUI.Slider(-1.0:0.001:1.0; default=0.1, show_value=true))
- M14 = $(@bind M14 PlutoUI.Slider(-1.0:0.001:1.0; default=0.1, show_value=true))
- M23 = $(@bind M23 PlutoUI.Slider(-1.0:0.001:1.0; default=-0.05, show_value=true))
"""

# ╔═╡ 50024c47-f094-4fed-8ec1-1861a982a49a
M = [
	1.0  M12  M13  M14  0.0;
	M12  2.0  M23  0.0  -0.10;
	M13  M23  3.0  0.1  0.05;
	M14  0.0  0.1  4.0  0.0;
	0.0  -0.1  0.05  0.0  5.0
];

# ╔═╡ ef9e68d5-b5d5-456e-a13d-8c0149c6b2d6
md"""
Consider the near-diagonal matrix


M = $(latexify_md(M))


with some Sliders to tune the off-diagonal elements:
"""

# ╔═╡ a7e14252-02ee-49c1-9f46-5699cf590923
md"""
- Plot Geschgorin disks: $(@bind show_geschgorin PlutoUI.CheckBox(default=true))
- Plot Bauer-Fike estimate: $(@bind show_bauer_fike PlutoUI.CheckBox(default=true))
- Plot Kato-Temple estimate: $(@bind show_kato_temple PlutoUI.CheckBox(default=true))
"""

# ╔═╡ 2d7cbb67-5c7a-455c-8c19-e9b9e8f4b832
begin
	geschgorin_centres = diag(M)
	geschgorin_radii = [sum(abs, row) - abs(row[i]) for (i, row) in enumerate(eachrow(M))]
	computed_eigenvalues = diag(M)
	exact_eigenvalues = eigvals(M)
	computed_eigenvectors = collect(eachcol(Matrix{Float64}(I, 5, 5)))
	residuals = map(computed_eigenvalues, computed_eigenvectors) do λ, v
		M * v - λ * v
	end
	error_Bauer_Fike = norm.(residuals)
	δ = map(1:size(M, 2)) do i
		δ_left = δ_right = Inf
	
		λtilde = computed_eigenvalues
		if i > 1
			δ_left  = abs(λtilde[i] - λtilde[i-1]) - norm(residuals[i-1])
		end
		if i < size(M, 2)
			δ_right = abs(λtilde[i] - λtilde[i+1]) - norm(residuals[i+1])
		end
		max(0.0, min(δ_left, δ_right))
	end
	error_Kato_Temple = norm.(residuals).^2 ./ δ
end;

# ╔═╡ 08508856-57bd-4d17-b8c6-3b0150e5b9af
begin
	function make_circle(xmidpoint, radius; n=100)
		θs = (0:n) .* 2π ./ n
		map(radius .* cos.(θs), radius .* sin.(θs)) do x, y
			x + xmidpoint, y
		end
	end
	
	p = plot(; aspect_ratio=1.0, xlims=[0, 6], ylims=[-1.5, 1.5], legend=:bottomright)
	scatter!(p, exact_eigenvalues, zeros(size(M, 1));
		    label="eigenvalues", mark=:x, markersize=8)
	
	scatter!(p, computed_eigenvalues, zeros(size(M, 1));
			 label="computed", mark=:+, markersize=8)
	
	for i in 1:size(M, 1)
		if show_geschgorin
			label = i==1 ? "Gerschgorin" : ""
			circle = make_circle(geschgorin_centres[i], geschgorin_radii[i])
			plot!(p, circle; color=3, label, lw=2)
		end

		if show_bauer_fike
			label = i==1 ? "Bauer-Fike" : ""
			circle = make_circle(computed_eigenvalues[i], error_Bauer_Fike[i])
			plot!(p, circle; color=4, label, lw=2)
		end

		if show_kato_temple
			label = i==1 ? "Kato-Temple" : ""
			circle = make_circle(computed_eigenvalues[i], error_Kato_Temple[i])
			plot!(p, circle ; color=5, label, lw=2)
		end
	end
	p
end

# ╔═╡ 760385db-cc30-40c1-8f66-5692028a96e3
TableOfContents()

# ╔═╡ e5c087db-1df5-4654-a5f3-d4c27a5a9782
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
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Plots = "~1.40.7"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.58"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "bea0b2664a0216395b8d096b166a222104912631"

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
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

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
git-tree-sha1 = "6c72198e6a101cccdd4c9731d3985e904ba26037"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.1"

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
git-tree-sha1 = "7bb1361afdb33c7f2b085aa49ea8fe1b0fb14e58"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.1+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "83dc665d0312b41367b7263e8a4d172eac1897f4"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.4"

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
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

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
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

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
git-tree-sha1 = "706dfd3c0dd56ca090e86884db6eda70fa7dd4af"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d3c8af829abaeba27181db4acb485b18d15d89c6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.1+0"

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
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

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
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

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
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

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
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

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
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

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
git-tree-sha1 = "2c962245732371acd51700dbb268af311bddd719"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.6"

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
git-tree-sha1 = "59071150afa35787c1656ba234cf03fdf8e2603f"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.8+0"

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
git-tree-sha1 = "9caba99d38404b285db8801d5c45ef4f4f425a6d"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.1+0"

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
git-tree-sha1 = "4bba74fa59ab0755167ad24f98800fe5d727175b"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.12.1+0"

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
# ╟─c00bddac-a453-4831-81be-d5814a7e88b9
# ╟─8dab1cdf-f700-4dd9-a839-58fe8cfd6c4a
# ╟─440fa8bb-c502-4c51-bd80-b0a86d3b8507
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
# ╠═9107b5a1-5768-436d-b705-6b9215121aa0
# ╟─bc4e51aa-674e-4a46-a20d-d9bfd82ff108
# ╟─ef9e68d5-b5d5-456e-a13d-8c0149c6b2d6
# ╟─50024c47-f094-4fed-8ec1-1861a982a49a
# ╟─45f6befb-4518-4467-81f3-ecd6da88ccbb
# ╟─08508856-57bd-4d17-b8c6-3b0150e5b9af
# ╟─a7e14252-02ee-49c1-9f46-5699cf590923
# ╟─2d7cbb67-5c7a-455c-8c19-e9b9e8f4b832
# ╟─760385db-cc30-40c1-8f66-5692028a96e3
# ╟─e5c087db-1df5-4654-a5f3-d4c27a5a9782
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
