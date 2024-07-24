### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ a0407e4d-3edc-48ea-a335-b815b796e378
begin
	import TikzPictures.TikzPicture
	using LaTeXStrings
	using PlutoUI
	using PlutoTeachingTools
	using HypertextLiteral

	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ e747d030-7487-4598-99d6-93afefd425e1
md"""
# Periodic Problems
"""

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
    L^2_{per} (\Omega) = \{ f \in L^2_{loc} (\mathbb R^3) \mid f \text{ is } \mathbb L \text{ periodic} \}
  \end{align}
  ```
  with inner product 
  ```math
  \begin{align}
    \langle f, g \rangle _{L_{per}^2 (\Omega)} = \int_{\Omega} \overline{f(x)} g(x) dx.
  \end{align}
  ```


- A convenient basis for $L^2_{per} (\Omega)$ is the plane wave or Fourier basis 
  ```math
  \begin{align}
      \mathbb B = \left \{ e_G(x) \equiv \frac1{\sqrt{\Omega}} e^{i G x} \middle \vert G \in \mathbb L ^* \right \} 
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
- Finally, we introduce some concepts from Fourier analysis.
  For a function $f \in L^2_{per} (\Omega)$ the following expressions are frequently encountered :
  ```math
  \begin{align}
      f = \sum_{G \in L^*} \hat f_G e_G \tag{Fourier Series Expansion}
  \end{align}
  ```
  with 
  ```math
  \begin{align}
      \hat f_G \equiv \langle f | e_G \rangle = \int_{\Omega} f e^{- i G x} dx    \tag{Fourier Coefficients}
  \end{align}
  ```
  and 
```math
\begin{align}
      \| f \|_{L^2_{per} (\Omega)} = \int_\Omega | f |^2 dx = \sum_{G \in \mathbb L^*} |\hat f_G|^2 < \infty. \tag{Parseval's Identity}
\end{align}
```

"""

# ╔═╡ 6c5c1231-9e16-4a43-8dec-68df585d78b9
md"""
- Employing a Fourier basis, *periodic Sobolev spaces* can be easily characterized as 
  ```math
  \begin{align}
      H^S_{per} (\Omega) = \left \{ f \in L^2_{per} (\Omega)  \middle | \sum_{G \in   \mathbb L^*} (1 + |G|^2)^S 
    |\hat f_G|^2 < \infty \right  \} 
  \end{align}
  ```
  where $S>0$.
"""

# ╔═╡ 32086662-d60d-47cd-acb6-0f1a4f18e351
md"""
## Periodic operators and periodic schrödinger problems

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
	Let $V \in L^{3/2}_{per} (\Omega)$ for $\Omega \in \mathbb R^3$. The
	operator 
	```math
		\opH = - \frac1{2} \laplacian + V
	```
	is self-adjoint on
	$\hilbert = L^2(\mathbb R^3)$ with $D(\opH) = H^2(\mathbb R^3)$ and $Q(\opH) = H^1(\mathbb R^3)$.
	$\opH$ is also bounded from below.


A perhaps surprising result is that, while $\opH$ is $\mathbb L$-periodic,
not all eigenfunctions of $\opH$ are $\mathbb L$-periodic.


"""

# ╔═╡ 3a8ff44a-c150-483a-b2f1-4fa085dd9075
md"""
## Bloch-Floquet theory

Since $\widehat{\laplacian
 f} (p) = - |p|^2 \hat f (p)$, the Fourier
transform 
```math
\begin{align}
    \hat f(p) = \frac1{(2 \pi)^{\frac{3}{2}}} \int_{\mathbb R} f(x) e^{-i p \cdot x} dx
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
    &= \frac1{(2 \pi)^{\frac{3}{2}}} \sum_{G \in \mathbb L^*} \int_{\Omega^*} \hat f(G+k) e^{i (k+G)\cdot x} dk
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
      (\bloch f) (x) = f_k(x) &= \frac1{\sqrt{|\Omega|}} \sum_{G \in \mathbb L^*} \hat f (G+k) e^{i G\cdot x}
     \\
    &= \sum_{G \in \mathbb L^*} \hat f (k+G) e_G(x).
  \end{align}
  ```
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
  f_k(x+R) &= f_k(x) \tag{Periodicity}
  \\
  f_{k+G} (x) &= e^{- i G \cdot x} f_k(x). \tag{Born-von Karman condition}
  \end{align}
  ```
  
- As a result we can always restrict $f_k(x)$ to $k \in \Omega^*$ and $x \in \Omega$.

"""

# ╔═╡ 9a6793e3-a46a-4736-9521-c5c6a1636aa3
md"""
- An explicit construction of $f_k$ from $f$ is available via **Poisson's formula** : 
  ```math
  \begin{align}
    f_k(x) &= \frac1{\sqrt{|\Omega|}} \sum_{G \in \mathbb L^*} \hat f (G + k) e^{i G \cdot x} 
    \\
    &= \frac1{\sqrt{|\Omega^*|}} \sum_{R \in \mathbb L} f(x+R) e^{- i k \cdot (x +R)}
  \end{align}
  ```
  This can be seen from 
  ```math
  \begin{align}
    \hat f(G+k) &= \frac1{(2 \pi)^{\frac{3}{2}}} \int_{\mathbb R^3} f(x) e^{-i (G + k) \cdot x} dx
    \\
    &= \frac1{(2 \pi)^{\frac{3}{2}}} \sum_{R \in \mathbb L} \int_\Omega f(R+y) e^{-i (G+k) \cdot (R+y)} dy
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
	$L^2_{qp} (\Omega^*, L^2_{per}(\Omega))$ (defined below)
	with the property that each $f \in C_0^\infty (\mathbb R^d)$ is mapped to 
	```math
	\Omega^* \ni K \mapsto f_k = \sum_{R \in \mathbb L } u(\bullet + R) e^{- k \cdot (\bullet + R)} \in L^2_{per}(\Omega).
	```
	The extension to $L^2(\mathbb R^d)$ follows by density of $C_0^\infty$ in $L^2$.

	In this $L^2_{qp} (\Omega^*, L^2_{per}(\Omega))$ is the space of $\mathbb L^*$-*quasiperiodic* functions which are $L^2_{per}(\Omega)$-valued :
	```math
	\begin{align}
	        L^2_{qp} (\Omega^*, L^2_{per}(\Omega)) = 
	         \left \{ \mathbb R^d \ni k \mapsto u_k \in L^2_{per}(\Omega) \middle \vert \int_{\Omega^*} \| u_k \|^2_{L^2_{per}(\Omega)} \ dk < \infty \right .
	        \\
	        \text{where } u_{k+G} = u_k e^{-i G \cdot x}  \text{ for all } G \in\mathbb L^*
	\\
		\text{and almost every } k \in \mathbb R^d \bigg \}.
	    
	\end{align}
	```
	On this vector space we have the inner product 
	```math
	\langle f,g \rangle_{L^2_{qp} (\Omega^*, L^2_{per}(\Omega))} = \frac1{|\Omega^*|} \int_{\Omega^*} \langle f_k, g_k \rangle_{L^2_{per}(\Omega)} dk.
	```
"""

# ╔═╡ fcf275be-b4c2-4d06-b502-a0c52e74b533
md"""
!!! tip "Remark"
	By way of Parceval's identity we can also understand $L^2_{qp} (\Omega^*, L^2_{per}(\Omega)) \equiv L^2(\Omega^* \times \Omega)$, i.e. $\bloch$ maps single variable functions of $y \in \mathbb R^d$ to functions of two variables : $k \in \Omega^*$ and $x \in \Omega$.

	We can illustrate it as
	```math
	\begin{align}
		L^2(\mathbb R^d) && \stackrel{\bloch}{\longrightarrow} && L^2_{qp} (\Omega^*, L^2_{per}(\Omega)) = L^2(\Omega^* \times \Omega)
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
\bloch H^1(\mathbb R^3) = L^2_{qp} (\Omega^*, H^1_{per}(\Omega)).
```
Similarly
```math
\bloch (- \laplacian _x) f = (- \laplacian_x f)_k = \lvert -i \nabla_k + k |^2 f_k = \lvert -i \nabla_k + k |^2 \bloch f
```
suggesting
```math
\bloch H^2(\mathbb R^3) = L^2_{qp} (\Omega^*, H^2_{per} (\Omega))
```
i.e. that regularity in $f$ induces equal regularity in $f_k$ for all $k \in \Omega^*$.

Notably, this also suggests the Bloch transform of the kinetic energy operator $- \laplacian_x$ as
```math
\bloch (- \laplacian_x) \bloch^{-1} = \lvert - i ∇_x + k|^2
```
interpreted in the sense of distributions and defined on the Hilbert space $L^2_{qp} (\Omega^*, L^2_{per} (\Omega))$ with domain $D( \bloch (- \laplacian_x) \bloch^{-1} ) = L^2_{qp} (\Omega^*, H^2_{per} (\Omega))$.

This result is formalized in the following theorem.
"""

# ╔═╡ 62e5c394-dad2-4a62-a7ae-d4999f741505
md"""
Denoting $\mathscr L (V) = \{ f : V \to V | f \text{ linear} \}$ the space of linear maps between a function space and itself, we have :

!!! note "Theorem 3 (Bloch fibers)"
	Any bounded $\mathbb L$-periodic operator $\opA$ on $L^2(\mathbb R^3)$ is decomposed by the Bloch-Floquet transform in the following sense : there exists a function $k \mapsto \opA_k$ in $L^\infty_{qp} (\mathbb R^3, \mathscr L(L^2_{per} (\Omega))$ such that for any $u \in L^2(\mathbb R^3)$, any $G \in \mathbb L^*$, and almost any $k \in \Omega^*$ is holds that 
	```math
	\begin{align}
	(\opA u)_k = \opA_k u_k && \opA_{k+G} = e^{-i G \cdot x} \opA_k e^{i G \cdot x}
	\end{align}
	```
	where the operators $(\opA_k)_{k \in \Omega^*}$ on $\hilbert = L^2_{per} (\Omega)$ are called the *Bloch fibers* of $\opA$.

	The Bloch decomposition can be extended to unbounded operators on $L^2(\mathbb R^3)$, in which case $\opA_k$ are unbounded as well.

!!! warning "Example 1 (Laplace operator)"
	For $\opA = - \laplacian$, the Bloch fibers aro $\opA_k = \lvert - i \nabla_x + k|^2$.
"""

# ╔═╡ 6da749bc-0181-49c2-8a7a-9c180930dce6
md"""
## Fibers of periodic Schrödinger operators

We consider $\opH = - \frac1{2} \laplacian + V$ on $L^2(\mathbb R^3)$ with domain $D(\opH) = H^2(\mathbb R^3)$ and $V \in L^2_{per} (\mathbb R^3)$ such that $\opH$ is $\mathbb L$-periodic.
Recall this operator is self-adjoint and as *no* eigenvalues.

- The fibers of $\opH$ are $\opH_k = \lvert - i \nabla_x + k|^2 + V$, where we note
  -  $\opH_k$ is a self-adjoint and bounded from below operator on $L^2_{per}(\Omega)$ with domain $H^2_{per}(\Omega)$ and form domain $H^1_{per}(\Omega)$
  -  $\opH_k$ has a compact resolvent $\resolvent(\opH_k)$, such that each $\opH_k$ has a purely discrete spectrum with eigenvalues accumulating at $+ \infty$ and eigenfunctions forming an orthonormal basis for $L^2_{per}(\Omega)$.


- To understand the spectrum of $\opH$, it suffices completely to study the spectrum of all $\opH_k$. 
  Indeed
  ```math
  	\sigma(\opH) = \bigcup_{k \in \overline{\Omega^*}} \sigma(\opH_k).
  ```


- We thus seek an $L^2_{per} (\Omega)$ orthonormal basis $\{ (λ_{k,n}, \phi_{k,n} )\}_{n \in \mathbb N} \subset (\mathbb R \times L^2_{per}(\Omega))^{\mathbb N}$ of eigenstates of $\opH_k$ :
  ```math
  \begin{align}
  	\opH_k \phi_{k,n} &= λ_{k,n} \phi_{k,n}
  	\\
  	\langle \phi_{k,n} , \phi_{k,m} \rangle_{L^2_{per}(\Omega)} &= \delta_{mn} && \forall n,m \in \mathbb N.
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

# ╔═╡ 41c0f837-055c-4077-b7f8-966e84ee709f
md"""
---

## Numerical treatment of periodic treatments

From our treatment it is clear that the lower end of the spectrum of $\sigma(\opH)$ for periodic potentials can be approximated by considering the lowest few bands via computing the lowest few eigenpairs of the fibers $\opH_k, k \in \overline{\Omega^*}$.

- For each $k \in \overline{\Omega^*}$, $\opH_k$ has only eigenvalues.
  Projection techniques are thus well-suited to approximate its eigenpairs from above.
  Moreover, Temple's inequality can be used for an a posteriori error bound.


- In practice, we select first a subset $\mathbb K \subset \overline{\Omega^*}$, the so-called *$k$-grid* or *$k$-point mesh*.


- A popular basis for discretizing $\opH_k$ and solving for its lowest-energy eigenstates using a Ritz-Galerkin projection method are the *k-adapted plane waves*
  ```math
  	\mathbb B^{E_{crit}}_{k} = \left  \{ e_G \middle \vert G \in \mathbb L^* , \frac1{2} |G+k|^2 < E_{crit} \right \}
  ```
  Note that thus a **different basis set** for each $k \in \mathbb K$ is employed.
"""

# ╔═╡ 9b7c42b4-9592-4ae3-8a86-41b6013885f3


# ╔═╡ e59d4d62-8ceb-43ad-ab86-159eac59500f
TableOfContents()

# ╔═╡ 7a952990-b524-43b6-a01b-cf6b6c6eaa14
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
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.59"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "64f925facca6edc922f6dbf4fb9dd734096e4cd4"

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
# ╟─e747d030-7487-4598-99d6-93afefd425e1
# ╟─a0407e4d-3edc-48ea-a335-b815b796e378
# ╟─1441bdfc-715e-46da-a4d2-8b788c71853f
# ╟─7971ffe3-39f0-4cf8-9cee-fd55483a8ed9
# ╟─8d9fe76f-dfb3-4683-83c2-5d7cd2aff712
# ╟─df23c23e-b509-4f48-9ab9-dd2b59834d01
# ╟─9b17dde4-6652-4800-a568-1244cf519593
# ╟─fbbe239e-c06e-41b1-b1e3-857585f994bd
# ╟─ea689cad-2859-41b0-aef3-36651a8d166b
# ╟─6ba1ce9a-6ced-4a82-9b5a-6b07a25c534b
# ╟─6c5c1231-9e16-4a43-8dec-68df585d78b9
# ╟─32086662-d60d-47cd-acb6-0f1a4f18e351
# ╟─3a8ff44a-c150-483a-b2f1-4fa085dd9075
# ╟─083e1f7a-e09e-479a-9cd9-83203e8c23dc
# ╟─f8c35c68-0e60-45d7-ab93-5aace2bd1589
# ╟─db338525-379e-44e5-b3eb-8d98e33d171b
# ╟─9a6793e3-a46a-4736-9521-c5c6a1636aa3
# ╟─17d99340-383a-4b29-a1a2-244024c71155
# ╟─fcf275be-b4c2-4d06-b502-a0c52e74b533
# ╟─9d4a8178-52f4-4846-94d6-74b8d68abbb4
# ╟─21523a64-b1a9-422d-8dc0-ebb324f35380
# ╟─62e5c394-dad2-4a62-a7ae-d4999f741505
# ╟─6da749bc-0181-49c2-8a7a-9c180930dce6
# ╟─08b20ea7-f3c3-4302-bfc3-a13bbd57930a
# ╟─41c0f837-055c-4077-b7f8-966e84ee709f
# ╟─9b7c42b4-9592-4ae3-8a86-41b6013885f3
# ╟─e59d4d62-8ceb-43ad-ab86-159eac59500f
# ╟─7a952990-b524-43b6-a01b-cf6b6c6eaa14
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
