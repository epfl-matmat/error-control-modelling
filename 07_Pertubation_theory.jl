### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 32b2d18c-edc0-11ee-13d9-83a8a646b6bd
begin
	import TikzPictures.TikzPicture
	using LaTeXStrings
	using PlutoUI
	using HypertextLiteral
	using PlutoTeachingTools

	RobustLocalResource("https://teaching.matmat.org/error-control/latex_macros.md", "latex_macros.md")
	Markdown.parse(read("latex_macros.md", String))
end

# ╔═╡ 57032426-2ff0-442d-aa5b-37a7c621d83f
md"""
# Matrix Pertubation Theory
"""

# ╔═╡ 25758608-0d58-474e-b6b1-d3463938d1dc
md"""
We saw in the previous lectures how to find the eigenpairs $(\tilde{\lambda}_{i}, \tilde{x}_{i} )$ of a Hermitian matrix $\tilde{A} \in \mathbb{C}^{n \times n}$.
However, in our discussion about floating point numbers we already noted that sometimes $A$ itself cannot be represented exactly, such that effectively we find the eigenpairs of $\tilde{A}+\Delta A$, where $\Delta A \in \mathbb{C}^{n \times n}$ is a small error.
The same happens if we deliberately employ reduced precision or approximate the computation of $\tilde{A}$ itself.
A natural question is thus : 
- How do the eigenpairs change if we go from $A$ to $\tilde{A}+\Delta A$ ?
One way to track this change in $\tilde{\lambda}_{i}$ and $\tilde{x}_{i}$ approximately without recomputing the eigenpairs of $\tilde{A}+\Delta A$ from scratch is **perturbation theory**.

Our setup is as follows: 
- We have available, for the Hermitian matrix $\tilde{A} \in \mathbb{C}^{n \times n}$, an approximate set of eigenpairs $(\tilde{\lambda}_{i}, \tilde{x}_{i})$, and we want to find a correction to the lowest few eigenpairs, such that they are approximately eigenpairs of $A(t)=\tilde{A}+t \Delta A$ where $t \in \mathbb{C}$ and $\Delta A \in \mathbb{C}^{n \times n}$ Hermitian. 
- We care in particular about the regime $t \simeq 0$.
- Clearly the eigenpairs of $A(t)$ are also $t$-dependant with $\lambda_{i}(t=0)=\tilde \lambda_{i}, x_{i}(t=0)=\tilde x_{i}$. 
- It would be natural to just expand $\lambda_{i}$ and $x_{i}$ in a Maclaurin series of $t$. 
  However, it is not clear *a priori* whether this series converges, what its convergence radius is and thus if this would even make sense.
- We will thus first study linear pertubations $\tilde A + t\Delta  A$, for which we first need the notion of a resolvent.

"""

# ╔═╡ 490b2be8-78e5-46a3-b313-413f90f8265d
md"""
## Resolvent

!!! note "Definition (Resolvent set and resolvent)"
	Let $A \in \mathbb{C}^{n \times n}$ and $z \in \mathbb{C}$. The set
	```math
		\resolvent(A)=\{z \in \mathbb{C} \mid A-zI \text{ is invertible} \}
	```
	is called the **resolvent set** of $A$.
	For $z \in \rho(A)$ the transformation
	```math
		R_{z}(A)=(A-z I)^{-1}
	``` 
	is called the **resolvent** of $A$.

!!! tip "Remark"
	In order to distinguish the resolvent set from the spectral radius, both used in this chapter and denoted by "rho", we will use $\resolvent$ (`\rho`) for the resolvent set and $\spectralradius$ (`\varrho`) for the spectral radius.
"""

# ╔═╡ ad786ef1-0b62-4a9a-84ca-c776d95e7cbc
md"""
!!! tip "Remark"
	- The mapping $z \mapsto R_{z}(A)$ has singularities exactly at the eigenvalues of $A$. In fact
	  ```math
	  \sigma(A)=\mathbb{C} \setminus \resolvent(A)
	  ```
	  is one way to define the **spectrum** of $A$.
	
	- Anywhere else $z \mapsto R_{z}(t)$ is analytic. 
	  Indeed, we have for any $\varepsilon>0$ and any $z$ from an $\varepsilon$-ball around $z_0$ away from eigenvalues, i.e. $z \in B_{\varepsilon} (z_{0}) \subset \resolvent(A)$, that 
	  ```math
	  \begin{aligned}
	  R_{z}(A)=(A-z)^{-1} & =\left[\left(A-z_{0} I\right)-\left(z-z_{0}\right) I\right]^{-1} \\
	  & =R_{z_{0}}(A)\left[I-\left(z-z_{0}\right) R_{z_{0}}(A)\right]^{-1}
	  \end{aligned}
	  ```
	  Now $\left[I-\left(z-z_{0}\right) R_{z_{0}}(A)\right]^{-1}$ can be expanded into a Neumann series
	  ```math
	  \left[I-\left(z-z_{0}\right) R_{z_{0}}(A)\right]^{-1}=\sum_{k=0}^{\infty}\left[\left(z-z_{0}\right) R_{z_{0}}(A)\right]^{k}
	  ```
	  whenever the spectral radius of $\left(z-z_{0}\right) R_{z_{0}}(A)$ is less the unity. 
	  The Taylor expansion of $R_{z}(A)$ thus exists as long as $z \in B_{\varepsilon} (z_{0} )$ with $\varepsilon<1 / \spectralradius \left(R_{z_{0}}(A)\right)$ and has the form
	
	  ```math
	  R_{z}(A)=\sum_{k=0}^{\infty}\left(z-z_{0}\right)^{k} R_{z_{0}}(A)^{k+1}
	  ```
	- A consequence (use Cramer's rule) is that the singularities at $z=\lambda_{1}, \dots, \lambda_{n}$ are not essential. 
	  Thus, the resolvent is **mesomorphic** on $z$ (only has a countable number of singularities, which are removable).

"""

# ╔═╡ d5a28a85-68f8-45de-912b-6944edbf1de7
md"""

We also mention the resolvent identities, which will prove useful later :
```math
\begin{align}
	R_{z_{1}}(A)-R_{z_{2}}(A) &= \left(z_{1}-z_{2}\right) R_{z_{1}}(A) R_{z_{2}}(A)
	\tag{1st identity}
	\\
	R_{z}(A)-R_{z}(B)&=R_{z}(A)(B-A) R_{z}(B) 
	\tag{2nd identity}
\end{align}
```
The first identity is also known as the Hilbert identity.

In what follows we emphasize the use of the resolvent to obtain spectral information


"""

# ╔═╡ 9845a2db-1878-435d-bb4f-e8e77cbde85d
md"""
## Spectral projectors

Let $\lambda_{i}$ be $a_{n}$ eigenvalue of $A$
and let $\contour_{\lambda_{i}}$ be a clockwise contour enclosing no other eigenvalues of $A$, as illustrated below.
"""

# ╔═╡ c3fd5d2f-9907-4500-8f40-511752fc6238
TikzPicture(L"""
%Shape: Polygon Curved [id:ds06447437278331569] 
\draw[blue] (119,73.37) .. controls (143,57.37) and (154.03,80.73) .. (168.5,80.87) node{>} .. controls (182.97,81.01) and (196.5,58.87) .. (217.5,58.37) node{>} .. controls (238.5,57.87) and (242.5,70.37) .. (242,93.37) node{$
\wedge$} .. controls (241.5,116.37) and (216,128.87) .. (197,129.37) node{<} node[above]{$C_{\lambda_i}$} .. controls (178,129.87) and (159,125.37) .. (140,133.87) .. controls (121,142.37) and (91.83,115) .. (123,110.87) .. controls (154.17,106.73) and (152.62,102.92) .. (150.5,99.87) .. controls (148.38,96.81) and (95,89.37) .. (119,73.37) -- cycle ;

\draw[->] (50,100) -- node{$\times$} node[above]{$\lambda_{i-1}$} (200,100) node{$\times$} node[above]{$\lambda_{i}$} -- node{$\times$} node[above]{$\lambda_{i+1}$} (300,100) ;
""", width="23cm", options="scale=0.05", preamble=raw"\usepackage{amsfonts}")

# ╔═╡ 7f02b176-3057-4ed6-8ebe-31d645f16e48
md"""
We consider the integral
```math
P_{i}=-\frac{1}{2 \pi i} \oint_{\contour_{\lambda_{i}}} R_{z}(A) d z
```
and establish these are orthogonal projectors for $i \neq j$ :
"""

# ╔═╡ e2f62f58-6614-4b80-b225-9e3136027f8d
md"""
!!! note "Theorem 1"
	The linear transformations $P_{i}, i=1, \dots, n$ defined above associated with distinct eigenvalues $i=1, \dots, n$ are such that

	- (a) $P_{i}^{2}=P_{i}$
	
	- (b) $P_{i} P_{j}=P_{j} P_{i}=0$ if $i \neq j$

	- (c) $\sum_{i=1}^{n} P_{i}=I$
"""

# ╔═╡ 2210e9c0-53e3-472e-80e7-4151d14eec9d
md"""
> *Proof* of (a).
> Let $\contour_{\lambda_{i}}$ and $\contour_{\lambda_{i}}^{\prime}$ be two curves enclosing $\lambda_{i}$ with $\contour_{\lambda_{i}}^{\prime}$ enclosing $\contour_{\lambda_{i}}$ (illustrated below). Then
>```math
>	\begin{aligned}
>	(2 i \pi)^{2} P_{i}^{2} & =\oint_{\contour_{\lambda_{i}}} \oint_{\contour_{\lambda_{i}}^{\prime}} R_{z}(A) R_{z^{\prime}}(A) d z d z^{\prime} \\
>	& \hspace{-1.3em} \stackrel{\text{1st identity}}{=} \oint_{\contour_{\lambda_{i}}} \oint_{\contour_{\lambda_{i}^{\prime}}} \frac{1}{z^{\prime}-z}\left(R_{z^{\prime}}(A)-R_{z}(A)\right) d z dz'
>	\end{aligned}
>```
>Now we use Cauchy's integral formula to obtain
>```math
>	\oint_{\contour_{\lambda_{i}}} \frac{d z}{z^{\prime}-z}=0 
>```
>since $z'$ is outside of $\contour_{\lambda_i}$ and thus there are no singularities inside $\contour_{\lambda_i}$, and
>```math
>	\oint_{\contour_{i}^{\prime}} \frac{d_{z^{\prime}}}{z^{\prime}-z}=2 \pi i
>```
>since $z$ is inside $\contour_{\lambda_i}'$
>Therefore,
> ```math
>	\begin{align}
>	& \oint_{\contour_{\lambda_{i}}} \oint_{\contour_{\lambda_{i}}^{\prime}} \frac{R\left(z^{\prime}\right)}{z^{\prime}-z} d z^{\prime} d z=\oint_{\contour_{\lambda_{i}}^{\prime}} R\left(z^{\prime}\right)\left(\oint_{\contour_{\lambda_{i}}} \frac{d z}{z^{\prime}-z}\right) d z^{\prime}=0 
>	\\
>	& \oint_{\contour_{\lambda_{i}}} \oint_{\contour_{\lambda_i}^{\prime}} \frac{R(z)}{z^{\prime}-z} d z^{\prime} d z=\oint_{\contour_{\lambda_i}} R(z)\left(\oint_{\contour_{\lambda_{i}}^{\prime}} \frac{d z^{\prime}}{z^{\prime}-z}\right) d z = 2 \pi i \oint_{\contour_{\lambda_i}} R(z) dz
>	\end{align}
> ```
> Such that $P_{i}^{2}=P_{i}$. $\hspace{12cm} \square$
"""

# ╔═╡ 38b1496f-cc7f-47c7-b40f-119d2dd1e663
TikzPicture(L"""
%Rounded Rect [id:dp03377449697776913] 
\draw  [draw opacity=0][fill={rgb, 255:red, 242; green, 242; blue, 242 }  ,fill opacity=1 ] (20,45.5) .. controls (20,41.36) and (23.36,38) .. (27.5,38) -- (232.5,38) .. controls (236.64,38) and (240,41.36) .. (240,45.5) -- (240,147.37) .. controls (240,151.51) and (236.64,154.87) .. (232.5,154.87) -- (27.5,154.87) .. controls (23.36,154.87) and (20,151.51) .. (20,147.37) -- cycle ;

%Shape: Rectangle [id:dp6341955401888753] (buffer)
\draw  [draw opacity=0] (9.5,45.87) -- (20,45.87) -- (20,147.37) -- (9.5,147.37) -- cycle ;


%Shape: Polygon Curved [id:ds8789544884464016] 
\draw[blue]   (124,50.37) .. controls (158.5,50.37) and (187.5,70.87) .. (187,99.87)  node{$\wedge$} .. controls (186.5,128.87) and (166.5,137.37) .. (148,139.87) node{$\times$} node[above right]{$z'$} .. controls (129.5,142.37) and (83,145.5) .. (64.5,109.87) node[left]{$C_{\lambda_i}'$}  .. controls (46,74.23) and (89.5,50.37) .. (124,50.37) -- cycle ;

%Shape: Polygon Curved [id:ds35740336827854424] 
\draw[red]   (126.5,71.37) .. controls (149,71.87) and (168.5,83.37) .. (168,100.87) node{$\wedge$} .. controls (167.5,118.37) and (154,120.37) .. (142.5,123.37) node{$\times$} node[above right]{$z$} .. controls (131,126.37) and (102.5,129.37) .. (90.5,108.37)  node[left]{$C_{\lambda_i}$} .. controls (78.5,87.37) and (104,70.87) .. (126.5,71.37) -- cycle ;

\draw[->] (30,90) -- node{$\times$} node[above]{$\lambda_i$} (230,90);
""",width="20cm",options="scale=0.04",preamble=raw"\usepackage{amsfonts}")

# ╔═╡ aefedae9-b792-41db-bb52-c8540619bcf1
md"""
> *Proof* of (b).
> The proof follows a very similar argument since $\contour_{\lambda_{i}}$ and $\contour_{\lambda_{j}}$ can be taken to be non-intersecting.
"""

# ╔═╡ 36f4bb7f-e1d5-4587-b96a-574820601963
md"""
> *Proof* of (c).
>We consider
>```math
>	P=-\frac{1}{2 \pi i} \sum_{i=1}^{n} \oint_{\contour_{\lambda_{i}}} R(z) d z .
>```
>  $R(z)$ has no poles apart from $\lambda_{i}$, so we can replace the sum over integrals enclosing single eigenvalues by a single integral enclosing all eigenvalues. 
>We choose this curve to be a circle $C$ with radius $r$ about the origin. 
>Thus
>```math
>	P=-\frac{1}{2 \pi i} \oint_{C} R(z) d z
>```
>By a change of variables $t=1 / z$ :
>```math 
>	P=-\frac{1}{2 \pi i} \oint_{\contour_{-}^{\prime}}(A- \frac1{t} I)^{-1}\left(-\frac{d t}{-t^{2}}\right)=-\frac{1}{2 \pi i} \oint_{\contour_{+}^{\prime}}(t A-I)^{-1} \frac{d t}{t}
>```
>where $\contour_{+/-}^{\prime}$ are the circles about the origin with radius $1/r$ that run clockwise/counter-clockwise.
>Note that $r$ must be larger than the spectral radius of $A$, $\spectralradius(A)$, thus $\spectralradius(t A)<1$.
>Therefore, $(I-t A)^{-1}$ can be expanded into a Neumann series
>```math
>	(I-t A)^{-1}=\sum_{k=0}^{\infty}(t A)^{k}
>```
>which converges, from which
>```math
>	\begin{align}
>	P & =\frac{1}{2 \pi i} \oint_{\contour_{+}^{\prime}}\left[\sum_{k=0}^{\infty} t^{k-1} A^{k}\right] d t &&  \\
>	& \hspace{-1.5em} \stackrel{\text{Residue Thm}}{=} \frac{1}{2 \pi i} \oint_{\contour_{+}^{\prime}} \frac{I}{t} d t 
>	\\
>	& \hspace{-1.5em} \stackrel{\text{Residue Thm}}{=\def\laplacian{{\Delta}}} \frac{2 \pi i}{2 \pi i} I=I
>	&& \square
>	\end{align} 
>```
"""

# ╔═╡ 129c7fb1-479e-4e82-a17c-0b153d2d890d
md"""
!!! note "Theorem 2"
	Let $A \in \mathbb{C}^{n \times n}$ Hermitian. 
	We have
	```math
	\im \left(P_{i}\right)=\ker\left(A-\lambda_{i} I\right)=\eigenspace_A\left(\lambda_{i}\right)
	```
	i.e. the image of the projector $P_{i}$ associated to $\lambda_{i}$ is an invariant subspace of $A$, which is equal to the eigenspace associated to $\lambda_{i}$.
"""

# ╔═╡ fec3438a-8da1-4d10-9b71-7c7b87d08fdf
md"""
> *Proof.*
>Since $R_{z}(A)$ and $A$ commute we get by integration $P_{i} A=A P_{i}$, such that $\im \left(P_{i}\right)$ is invariant under $A$.
>
>Next we show $\im \left(P_{i}\right) \supseteq \ker\left(A-\lambda_{i} I\right)$.
>Consider an $x \in \ker\left(A-\lambda_{i} I\right)$ and expand
>```math
>	\begin{aligned}
>	R_{z}(A) x & =(A-z I)^{-1} x \\
>	& =\left[\left(A-\lambda_{i} I\right)-\left(z-\lambda_{i}\right) I\right]^{-1} x \\
>	& =-\frac{1}{z-\lambda_{i}}\left[I-\left(z-\lambda_{i}\right)^{-1}\left(A-\lambda_{i} I\right)\right]^{-1} x \\
>	& =-\frac{x}{z-\lambda_{i}} .
>	\end{aligned}
>```
>Therefore
>```math
>P_{i} \ x=\left[-\frac{1}{2 \pi i} \oint_{\contour_{i}} R_{z}(A) d z\right] x=\frac{x}{2 \pi i} \oint_{\contour_{i}} \frac{1}{z-\lambda_{i}} d z=x
>```
>and $x \in \im \left(P_{i}\right)$.
>
>Finally, we show $\im \left(P_{i}\right) \subseteq \ker\left(A-\lambda_{i} I\right)$. 
>With
> ```math
>\left(z-\lambda_{i}\right) R_{z}(A)=-I+\left(A-\lambda_{i} I\right) R_{z}(A)
>```
>we deduce
>```math
>	\begin{aligned}
>	-\frac{1}{2 \pi i} \oint_{\contour_{\lambda_{i}}}\left(z-\lambda_{i}\right) R_{z}(A) d z & =-\frac{1}{2 \pi i}\left(A-\lambda_{i} I\right) \oint_{\contour_{\lambda_i}} R_{z}(A) d z \\
>	& =\left(A-\lambda_{i} I\right) P_{i}
>	\end{aligned}
>```
>and more generally
>```math
>	-\frac{1}{2 \pi i} \oint_{\contour_{\lambda_{i}}}\left(z-\lambda_{i}\right)^{k} R_{z}(A) d z=\left(A-\lambda_{i} I\right)^{k} P_{i} . \tag{1}
>```
>Since the resolvent has no essential singularities, there exists a $k>0$ integer for which the left-hand side of (1) is zero. 
>Since $\forall x \in \im \left(P_{i}\right)$, $x=P_{i} \ x$ this shows $x \in \ker\left(A-\lambda_{i} I\right)^{k}.$ 
>For Hermitian, and thus diagonalisable matrices,
>```math
>	\ker\left(A-\lambda_{i} I\right)^{k}=\ker\left(A-\lambda_{i} I\right) \quad \forall k \geq 1
>```
> thus $\im \left(P_{i}\right) \subseteq  \ker \left(A-\lambda_{i} I\right)$. $\hspace{9cm} \square$
"""

# ╔═╡ f0bbb6e5-5bbf-4c5b-895c-6d4531b75138
md"""
!!! tip "Remark"
	Making use of the eigenbasis $\{x_{i}\}_{i}$ of $A$ we can also denote
	```math
		\indicator_{\{\lambda_{i}\}}(A) =
		\sum_{i=1}^{n} \indicator_{\{\lambda_{i}\}} (\lambda_{i}) x_{i} x_{i}^{H}
		=x_{i} x_i^H
		= \ker (A-\lambda_{i} I )=P_{i}
	```
This proves that the resolvent and associated spectral projectors contain all the information about the eigendecomposition of $A$. 

"""

# ╔═╡ cd7ed2fa-e9fa-45f9-ac6a-ba9a517e5593
md"""
## Analyticity

We proceed by considering the setting of our linear pertubation $A=A+t \Delta A$, which we study by considering the associated resolvent
```math
R_{z}(A(t))=(\tilde{A}+t \Delta A-z I)^{-1}
```
First we note

"""

# ╔═╡ 91ac9688-aa16-4d5c-8629-407141c9eee0
md"""
!!! note "Proposition 3"
	The resolvent $R_{z}(A(t))$ is analytic with respect to $t$ in the open disk of the complex plane with
	```math
		|t|< \frac1{\spectralradius( R_{z}(\tilde{A}) \Delta A )}
	```

> *Proof.*
> Note that
> ```math
>	\begin{aligned}
>	R_{z}(A(t)) & =(\tilde{A}+t \Delta A-z I)^{-1} \\
>	& =R_{z}(\tilde{A})\left[I+t R_{z}(\tilde{A}) \Delta A\right]^{-1}
>	\end{aligned}
> ```
> which is analytic in $t$ as long as $\spectralradius (t R_{z}(\tilde{A}) \Delta A )<1$ 
> $\hspace{6cm} \square$

"""

# ╔═╡ 4851946b-d584-495e-adc1-008e91555ec6
md"""
!!! note "Theorem 4"
	Let $\Gamma$ be a curve around one or a few eigenvalues of $\tilde{A}$ and let
	```math
	\spectralradius_a = \inf_{z \in \Gamma}\left[\spectralradius (R_{z}(\tilde{A}) \Delta A )\right]^{-1} \text {. }
	```
	
	Then $\spectralradius_{a}>0$ and the spectral projector
	
	```math
	P(t)=-\frac{1}{2 \pi i} \oint_{\Gamma} R_{z}(A(t)) d z
	```
	
	is analytic in the disk $|t|<\spectralradius_a$.
"""

# ╔═╡ 163c66ac-6c48-476d-a4d2-9f44ef54dbe9
md"""
> *Proof.*
> This result boils down to proving that under the given conditions the resolvent $R_{z}(A(t))$ is analytic for each $z \in \Gamma$, which is clearly the case if
>
> ```math
>	|t|< \inf _{z \in \Gamma} \frac1{\spectralradius (R_{z}(\tilde A) \Delta A )} =\spectralradius_a  .
> ```
>We now show that $\spectralradius_{a}>0$. Consider
>```math
>	\spectralradius (\Delta A R_{z}(\tilde{A}) ) 
>	\leq \|R_{z}(\tilde{A}) \Delta A \| 
>	\leq\|\Delta A\| \cdot \|R_{z}(\tilde{A}) \|.
>```
> Since $z \mapsto\left\|R_{z}(A)\right\|$ is continuous for $z \in \Gamma$, it reaches a maximum $z_{0}$ (extremal value theorem) for which
>```math
>	\begin{aligned}
>	\spectralradius (\Delta A R_{z}(\tilde{A}) ) 
>	& \leq\|\Delta A\| \|R_{z}(\tilde{A}) \| 
>	\\
>	& \leq\|\Delta A\| \|R_{z_{0}}(\tilde{A}) \| \equiv K
>	\end{aligned}
> ```
> Hence, $1 / \spectralradius_{a} \leq K \iff  \spectralradius_{a} \geq 1 / K>0$.
> $\hspace{8cm} \square$
"""

# ╔═╡ 2d9827a3-2752-4bac-9c33-c46fa93b52f9
md"""
An immediate corollary is that the rank of $P(t)$ stays constant as long as $|t|<\spectralradius_{a}$ :

!!! note "Corollary 5"
	The number $m$ of eigenvalues, counted with multiplicity, located in $\Gamma$ is constant given that $|t|<\spectralradius_a$.

Note, that this is usually too strong of a condition. 
The real condition is that $P(t)$ is analytic.

"""

# ╔═╡ 2c1b0372-5232-45bf-95d1-1a6641508ef0
md"""

!!! tip "Remark"
	- While individual eigenvalues may not be analytic is $t$, the average
	  ```math
	  	\hat{\lambda}(t)=\frac{1}{m} \sum_{i=1}^{m} \lambda_{i}(t)
	  ```
	  of the $m$ eigenvalues of $A(t)$, which are in $\Gamma$, usually is analytic.

	- Another way to obtain $\hat{\lambda}(t)$ is to consider the trace of the restriction of $A(t)$ to the invariant subspace $\im P(t)$, namely
	  ```math
	  \begin{align}
	  \operatorname{tr}\left(\left.A(t)\right|_{\im P(t)}\right) &   =\operatorname{tr}\left[ (P(t) A(t) P(t)) \middle \vert_{\im P(t)}\right] \\
	  & =\operatorname{tr}(P(t) A(t) P(t)) \\
	  & =\operatorname{tr}(A(t) P(t)) \\
	  & =\sum_{i=1}^{m} \lambda_{i}(t) \\
	  & =m \hat{\lambda}(t)
	  \end{align}
	  ```
"""

# ╔═╡ 653598ee-c9d7-4731-9c8f-687d5962c0d6
md"""
!!! note "Theorem 6"
	The linear transformation $A(t) P(t)$ and its weighted trace $\hat{\lambda}(t)$ are analytic in the disk $|t|<\spectralradius_{a}$.

> *Proof.* 
> $A(t) P(t)$ is analytic due to previous discussion. 
> $t \mapsto \operatorname{tr}(M(t))$ is analytic if $M(t)$ is analytic, which is the case.
> $\hspace{11cm} \square$

"""

# ╔═╡ bd9f3ea5-4fd8-4155-bc80-5d658b304abe
md"""
Because of Theorem 6, a *simple eigenvalue* is analytic around $t=0$ and stays simple. 
The vector $x_{i}(t)=P_{i}(t) \tilde{x}_{i}$ is an *eigenvector* of $A(t)$ associated to this simple eigenvalue if $\tilde{x}_{i}=x_{i}(t=0)$ is an eigenvector of $\tilde{A}=A(t=0)$. 
Clearly $x_{i}(t)$ is *analytic too* !

For a multiple eigenvalue the situation is more complex. 
In general if an eigenvalue is of multiplicity $m$ then it will split into $m$ distinct small branches $\lambda_{i}(t)$. 
Individually these are *not* analytic in $t$, but their average is.

This makes it important to identify groups of eigenvalues in practice.
"""

# ╔═╡ c5067a9b-e202-4b03-9f54-aa367b2f4ebd
md"""
A final nice result is

!!! note "Proposition 7"
	Any eigenvalue $\lambda_{i}(t)$ of $A(t)$ inside the curve $\Gamma$ satisfies
	```math
	\left|\lambda_{i}(t)-\tilde{\lambda}_{i}\right|= O(|t|) .
	```

> *Proof.*
> Due to Theorem 6, 
> $\im (\tilde{P}_{i} ) = \ker (\tilde{A}-\tilde{\lambda}_{i} I )$ 
> which implies $(\tilde A-\tilde{\lambda}_{i} I ) \tilde{P}_{i}=0$. For an eigenvector $x(t)$ with norm unity associated to $\lambda_{i}(t)$ we thus have
>```math
>	\begin{align}
>	\left(A(t)-\tilde{\lambda}_{i} I\right) P(t) u(t) & =\left(A(t)-\tilde{\lambda}_{i} I\right) u(t) \\
>	& =\left(\lambda(t)-\tilde{\lambda}_{i}\right) u(t)
>	\end{align}
> ```
> Taking norms on both sides :
>```math
>\begin{align}
> |\lambda(t)-\lambda_{i} | & = \| (A(t)-\tilde{\lambda}_{i} I ) P(t) u(t) \| \\
>& \leq \| (A(t)-\tilde{\lambda}_{i} I ) P(t) \| \\
>& = \| (A(t)-\tilde{\lambda}_{i} I ) P(t)-
> {\color{gray} \underbrace{\color{prooftext} (\tilde{A}- \tilde{\lambda}_{i} I ) \tilde{P}_{i} }_{=0} }  \|  \tag{2} \\
>& =O(|t|) \tag{3}
>\end{align}
>```
> where we went from (2) to (3) by considering that the (2) is analytic in $t$ and $\tilde A=A(0)$, $\tilde P_{i}= P(0)$, etc.
> $\hspace{10cm} \square$

"""

# ╔═╡ 22fe9a75-06bf-43b3-9cb0-92a9022ae653
md"""
## Pertubation theory setting

We return to the question of computing the lowest order variations of $\lambda(t)$ and $x(t)$ for $A(t)=\tilde{A}+t \Delta A$, having $(\tilde{\lambda}, \tilde{x})$ available. 
We assume $\lambda(t)$ to be simple and isolated with $\Gamma$ a contour around $\tilde{\lambda}$. 
Further we assume $|t|<\spectralradius_{a}=\inf _{z \in \Gamma} 1 / \spectralradius  (R_{z} (\tilde A) \Delta A)$ such that $\lambda(t)$ and $x(t)$ are analytic. 

- Thus, we may expand them in a Taylor series:
  ```math
  	\begin{align}
  	& \lambda(t)=\tilde{\lambda}+t \lambda^{\prime}(0)+\frac{1}{2} t^{2} \lambda^{\prime \prime}(0)+\cdots \\
  	& x(t)=\tilde{x}+t x^{\prime}(0)+\frac{1}{2} t^{2} x^{\prime \prime}(0)+\cdots
  	\end{align}
  ```

- Since $A(t) x(t)=\lambda(t) x(t)$ we have, to first order in $t$,
  ```math
  	\begin{align}
  	 \tilde{A} \tilde{x}+t\left(\Delta A \tilde{x}+\tilde{A} x^{\prime}(0)\right)
  	& =\tilde{\lambda} \tilde{x}+t\left(\lambda'(0) \tilde{x}+\tilde{\lambda} x^{\prime}(0)\right)+ O \left(| t |^{2}\right) 
  	\\
  	\Rightarrow \quad \Delta A \tilde{x}+\tilde{A} x^{\prime}(0)
  	& =\lambda^{\prime}(0) \tilde{x}+\tilde{\lambda} x^{\prime}(0)+O((t)) \tag{4}
  	\end{align}
  ```

- Now we denote by $P_{\tilde \lambda}$ the projector into the eigenspace of $\tilde{\lambda}$ (generated by $R_{z}(\tilde{A})$ via a contour integral) and denote $Q_{\tilde \lambda}=I-P_{\tilde \lambda}$.
"""

# ╔═╡ eb10b66a-c103-4964-83e4-e980fdcaa36b
md"""
- Projecting (4) onto $P_{\tilde \lambda}$, we obtain after a rearrangement
  ```math
  \begin{align}
  	P_{\tilde{\lambda}} \Delta A \tilde{x}-\lambda^{\prime}(0) {\color{gray} \underbrace{\color{black} P_{\tilde{\lambda}} \tilde{x}}_{= \tilde x}}
  	= {\color{gray} \underbrace{\color{black} (\tilde{A}-\tilde{\lambda} I) P_{\tilde{\lambda}}}
  	_{\substack{= 0 \text{ since} \\ \ker (\tilde A - \tilde \lambda I) = \im P_{\tilde \lambda} }}}
  	x^{\prime}(0) +O(|t|)
  \end{align}
  ```
  Multiplying by $\tilde x^H$, we obtain
  ```math
  \lambda^{\prime}(0)=\tilde{x}^{H} \Delta A \tilde{x}=\langle\tilde{x}, \Delta A   \tilde{x}\rangle
  ```

- Projecting (4) onto $Q_{\tilde \lambda}$, we obtain
  ```math
  Q_{\tilde{\lambda}}(\tilde{A}-\hat{\lambda} I) x^{\prime}(0)=-Q_{\tilde{\lambda}} \Delta A \tilde{x}+O(|t|)
  ```
  Using the invariance of $Q_{\tilde \lambda}$ with respect to $\tilde A - \tilde \lambda I$, we get
  ```math
  \begin{align}
   Q_{\tilde{\lambda}}(\tilde{A}-\tilde{\lambda} I) Q_{\tilde{\lambda}} x^{\prime}(0) &=-Q_{\tilde{\lambda}} \Delta A \tilde{x} +O(|t|) 
  \\
  \Rightarrow Q_{\tilde{\lambda}} x^{\prime}(0) &=-Q_{\tilde{\lambda}}(\tilde{A}-  \tilde{\lambda} I)^{-1} Q_{\tilde \lambda}  \Delta A \tilde{x} \\  
  &=-R_{\tilde{\lambda}}(\tilde{A}) Q_{\tilde{\lambda}} \Delta A \tilde x
  \end{align}
  ```
  which defines the first-order change *orthogonal* to the eigenspace of $\tilde   \lambda$.

"""

# ╔═╡ 645a99c9-dc40-484d-bfe1-ea5a9745c009
md"""
- To fix the remaining degree of freedom of $x^{\prime}(0)$, namely $\langle x^{\prime}(0), \tilde{x} \rangle$ we consider the normalisation of $x(t)$. 
  We require to first order
  ```math
  \begin{aligned}
  1 & =(x(t))^{H} x(t) \\
  & =\tilde{x}^H \tilde{x}+t \tilde{x}^{H} x^{\prime}(t)+t\left(x^{\prime}(t)\right)^{H} \tilde{x}+O\left(|t|^{2}\right) \\
  \Rightarrow 0 & =\operatorname{Re} \tilde{x}^{H} x^{\prime}(t)+O(|t|)
  \end{aligned}
  ```

- Additionally, we conventionally *choose* $\im x^H x^{\prime}(t)=0$. 
  This is *not* necessary, but reasonable to ensure $x(t)$ stays a real vector in case $\tilde{x}$ is real. 
  As a result $\tilde x$ and $x^{\prime}(t)$ are orthogonal.

The final first-order changes are
```math
\boxed{
\begin{array}{l}
\lambda'(0)=\langle\tilde{x}, \Delta A \tilde{x})  \\
x^{\prime}(0)=-(\tilde{A}-\tilde{\lambda} I)^{-1} Q_{\tilde{\lambda}} \Delta A \tilde{x}
\end{array}}
\tag{PT}
```

"""

# ╔═╡ bc34fa94-f155-41fb-92d1-3e51db702f26
md"""
## Interpretation and applications of this result

If a modelling error / error is the computation causes us to employ $\tilde A$ instead of $A$ for the computation, then this result gives is the leading forward error with $\Delta A=A-\tilde{A}$.

- We obtain the absolute forward error
  ```math
  	\begin{align}
  	\|\lambda'(0) \| & \leq \|\Delta A\|\|x\|=\|\Delta A\| 
  	\\
  	\|x'(0)       \| & \leq \|(\tilde{A}-\tilde{\lambda} I)^{-1} Q_{\tilde{\lambda}} \Delta A \| 
  	\\
  	& \leq \|(\tilde{A}-\tilde{\lambda} I)^{-1} Q_{\tilde{\lambda}} \|\| \Delta A\| 
  	\\
  	&= \frac{1}{\delta}\|\Delta A\| 
  	\end{align}
  ```
  where $\delta = \min _{\lambda \in \sigma(A)}|\lambda-\tilde{\lambda}|$ is the gap.
- Since the backward error is just $\| \Delta A\|$, this implies that the condition number is 1 for the eigenvalue computation and is $1 / \delta$ for the eigenvector computation.

"""

# ╔═╡ f3fe5786-f3d9-4c14-9f51-911732a9dfb0
md"""
#### Newton's method

With respect to the exact eigenproblem $A(1)=\tilde{A}+\Delta A$, the term $\Delta A \tilde{x}$ is exactly the residual of $(\tilde \lambda, \tilde{x} )$. 
We now consider one step of Newton's method for solving the system of equations:

```math
\begin{align}
(A-\lambda I) x &= 0 \\
\frac{1}{2} x^{H} x-1 &= 0
\end{align}
```

with $A=A(1), \lambda=\lambda(1), x=x(1)$, where we desire the unknown eigenpair $(\lambda, x)$ while having the approximation $(\tilde \lambda, \tilde{x} )$. 
One Newton step goes to
```math
	\begin{pmatrix}
		\tilde{x}_{\text {new}} \\
		\tilde{\lambda}_{\text {new}}
	\end{pmatrix} 
	=
	\begin{pmatrix}
		\tilde{x} \\
		\tilde{\lambda}
	\end{pmatrix}
	-
	\begin{pmatrix}
	A-\lambda I & -\tilde{x} \\
	\tilde{x}^H & 0
	\end{pmatrix}^{-1}
	\begin{pmatrix}
	r \\
	0
	\end{pmatrix}
```
where $r=A \tilde{x} - \tilde{\lambda} \tilde{x}=\Delta A \tilde{x}$. 
If we replace $A-\lambda I$ by the approximation $\tilde A-\tilde{\lambda} I$ we obtain a Quasi-Newton scheme, which under close inspection turns out to be exactly equivalent to (PT).

Note that the normalisation of $x$ in the correction with respect to $\tilde{x}$ is again enforced by the constraint $\frac{1}{2} x^{H} x-1=0$. 

Typically, $\tilde{x}_{\text {new}}-\tilde{x} \perp \tilde{x}$ is manually enforced by inverting $Q_{\tilde \lambda} (\tilde{A}-\tilde{\lambda} I) Q_{\tilde \lambda}$ instead of $(\tilde A - \tilde \lambda I)$.
"""

# ╔═╡ 1a28c717-c969-4ce9-be55-80b566f6436d


# ╔═╡ ed155363-2f4a-4aa8-9b41-79a5fa28c058
TableOfContents()

# ╔═╡ e5b16293-ca61-4fc6-b89e-82100792efc7
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
PlutoUI = "~0.7.58"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "925932f0baa3a074c1147a9fc8889e10396eea69"

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
# ╟─57032426-2ff0-442d-aa5b-37a7c621d83f
# ╟─32b2d18c-edc0-11ee-13d9-83a8a646b6bd
# ╟─25758608-0d58-474e-b6b1-d3463938d1dc
# ╟─490b2be8-78e5-46a3-b313-413f90f8265d
# ╟─ad786ef1-0b62-4a9a-84ca-c776d95e7cbc
# ╟─d5a28a85-68f8-45de-912b-6944edbf1de7
# ╟─9845a2db-1878-435d-bb4f-e8e77cbde85d
# ╟─c3fd5d2f-9907-4500-8f40-511752fc6238
# ╟─7f02b176-3057-4ed6-8ebe-31d645f16e48
# ╟─e2f62f58-6614-4b80-b225-9e3136027f8d
# ╟─2210e9c0-53e3-472e-80e7-4151d14eec9d
# ╟─38b1496f-cc7f-47c7-b40f-119d2dd1e663
# ╟─aefedae9-b792-41db-bb52-c8540619bcf1
# ╟─36f4bb7f-e1d5-4587-b96a-574820601963
# ╟─129c7fb1-479e-4e82-a17c-0b153d2d890d
# ╟─fec3438a-8da1-4d10-9b71-7c7b87d08fdf
# ╟─f0bbb6e5-5bbf-4c5b-895c-6d4531b75138
# ╟─cd7ed2fa-e9fa-45f9-ac6a-ba9a517e5593
# ╟─91ac9688-aa16-4d5c-8629-407141c9eee0
# ╟─4851946b-d584-495e-adc1-008e91555ec6
# ╟─163c66ac-6c48-476d-a4d2-9f44ef54dbe9
# ╟─2d9827a3-2752-4bac-9c33-c46fa93b52f9
# ╟─2c1b0372-5232-45bf-95d1-1a6641508ef0
# ╟─653598ee-c9d7-4731-9c8f-687d5962c0d6
# ╟─bd9f3ea5-4fd8-4155-bc80-5d658b304abe
# ╟─c5067a9b-e202-4b03-9f54-aa367b2f4ebd
# ╟─22fe9a75-06bf-43b3-9cb0-92a9022ae653
# ╟─eb10b66a-c103-4964-83e4-e980fdcaa36b
# ╟─645a99c9-dc40-484d-bfe1-ea5a9745c009
# ╟─bc34fa94-f155-41fb-92d1-3e51db702f26
# ╟─f3fe5786-f3d9-4c14-9f51-911732a9dfb0
# ╟─1a28c717-c969-4ce9-be55-80b566f6436d
# ╟─ed155363-2f4a-4aa8-9b41-79a5fa28c058
# ╟─e5b16293-ca61-4fc6-b89e-82100792efc7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
