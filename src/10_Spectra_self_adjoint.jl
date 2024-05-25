### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4004f2e8-5738-47e2-9d3c-8b68b896b9a5
begin
	import TikzPictures.TikzPicture
	using LaTeXStrings
	using PlutoUI
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

# ╔═╡ 565cdda1-d844-4bed-8f48-fe57c1428136
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

	#matrix elements
	function Mij(i,j,V)
		# I know that this is definitely not the most efficient approach from a computational standpoint, but it works
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

# ╔═╡ 1478c459-0933-4c3f-8e84-86d051bea4b3
md"# Spectra of Self-Adjoint Operators"


# ╔═╡ 21814beb-9d27-44c7-a41f-1191d856d1b7
md"""

## Self-adjoint operators

The following strong result is the justification for the elaborate construction of self-adjoint operators, where great care in constructing their domain is often needed.

!!! note "Theorem 1"
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

We thus obtain that the familiar result that self-adjoint operators have a real spectrum.

"""

# ╔═╡ d2b95a1b-5e4c-4e85-8a6e-d893c9c3cbde
md"""
## Weyl sequences for spectral characterization

Another important characterisation method is based on the convergence of bounded sequences, where we already noted differences between the finite and infinite dimensional case.
To fully appreciate the details we need a few more definitions.
"""

# ╔═╡ 206bf4ec-9bdd-4141-8ee4-570f077c06d0
md"""



!!! note "Definition 1 (Weak Convergence)"
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

# ╔═╡ 390974fb-36d3-44e1-8543-24a16e400512
md"""
In infinite dimensions closed and bounded sets are no longer necessarily compact.
Thus bounded sequences may not have strongly converging subsequences.
However,

!!! note "Theorem 2" 
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

!!! note "Definition 2 (Weyl Sequence)"
	A sequence $(\phi_n) \subset D(\opA)$
	with $\| \phi_n \| = 1$, which satisfies
	$\| (\opA - \lambda) \phi_n \| \to 0$ for a $\lambda \in \mathbb R$ is called a
	*Weyl sequence*.



"""

# ╔═╡ 95acb6bd-42d7-4ede-b425-a91ed9640d9a
md"""
!!! note "Theorem 3"
	Let $\opA$ be self-adjoint with $D(\opA) \subset \hilbert$ and
	$\lambda \in \mathbb R$. The following are equivalent

	1.  $\lambda \in \sigma(\opA)$
	
	2.  $\inf_{\varphi \in D(\opA), \| \varphi \| = 1} \| (\opA - \lambda) \varphi \| = 0$
	
	3.  There exists a Weyl sequence for $\lambda$.


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

# ╔═╡ 20958c09-7372-4c0b-bd3a-3011758cba2e
md"""
One consequence of Theorem 3 is :

!!! note "Theorem 4"
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
Similarly, upper semi-bounded operators satisfy $\langle \varphi, \op B  \varphi \rangle \leq b$ and satisfy $\sigma (\op B) \subset (- \infty , b]$

!!! tip "Remark"
	The striking similarity with Courant-Fisher for matrices already points to a possible generalization of the min-max statement statements for operators. 
	This we will pursue next.



"""

# ╔═╡ 2f377ee6-e23e-432f-aeb6-9097fa260027
md"""
Beforehand, let us pause for a second and use this result to deduce the spectra for a few self-adjoint operators on $\hilbert = L^2(\mathbb R ^d)$ 

!!! warning "Example 2 (Identity operator)"
	 $\mathop{\mathrm{id}}: \hilbert \to \hilbert$ (i.e. $D(\mathop{\mathrm{id}}) = \hilbert$) is clearly bounded from above and below by 1. 
	Thus, $\sigma(\mathop{\mathrm{id}}) = \{ 1\}$
	

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
	
	Now, we take $k_0 \in \mathbb R^d, f \in H^2(\mathbb R^d)$ and define a Weyl sequence 
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
	Hence, $[0, \infty) = \sigma(-\laplacian)$
"""

# ╔═╡ d2bf63fb-6035-49b2-9182-c606bb3643af
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

This classification turns out to be not very useful for considering spectral approximations as small pertubation of $\opA$ can easily mix $\sigma_p$ and $\sigma_{\text{cont}}$.
Whe therefore use an alternate :



"""

# ╔═╡ a803317a-07c5-4f95-bc45-4c24f3971385
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

One can show that every isolated point of $\sigma(\opA)$ is an eigenvalue.
They are thus interesting entries, motivating the following spectrum :


!!! note "Definition 5 (Discrete and essential spectrum)"
	Let $\opA : D(\opA) \to \hilbert$ be a self-adjoint operator.
	```math
		\sigma_{\text{disc}} = \{ \lambda \in \sigma(\opA) \mid \lambda \text{ is an isolated point and has finite multiplicity} \}
	```
	is called the *discrete spectrum*.

	The complement
	```math
		\sigma_{\text{ess}} = \sigma(\opA) \setminus \sigma_{\text{disc}} (\opA)
	```
	is called the *essential spectrum*.


"""

# ╔═╡ f2b8a3ba-71e9-48a9-b77d-91aa844f3da3
md"""
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


# ╔═╡ dda55c4f-1a71-477b-b460-0c67629c2904
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

# ╔═╡ 8e5dd330-181a-4279-b54b-47be263e13a9
md"""
---

Notice how eigenvalues (members of $\sigma_p$) may be embedded inside the continuous spectrum.
These are very unstable to perturbations (so-called resonance states) and, as we will see now, hard to approximate.

Notice that, unlike $\sigma_{\text{cont}}$, $\sigma_{ess}$ is always closed for a self-adjoint operator.


"""

# ╔═╡ e60a51b3-88ed-4e03-9937-4c86773e33e8
md"""
Overall, Weyl sequences turn out to yield a useful characterization of the spectra of self-adjoint operators, and overview of which is shown here :

Spectrum | Weyl sequence for spectral characterization
---|:---
$\lambda \in \sigma(\opA)$ | $\exists \text{ a Weyl sequence } (\phi_n) \subset D(\opA) \text{ s.t. } \| \phi_n \| = 1 \text{ and } \| (\opA - \lambda) \phi_n \| \to 0$
$\lambda \in \sigma_{\text{ess}} (\opA)$ | $\exists \text{ a Weyl sequence s.t. } \phi_n \rightharpoonup 0 \text{ (weakly)}$
$\lambda \in \sigma_{\text{disc}} (\opA)$ | $\text{\emph{All} Weyl sequences have subsequences } \phi_{n_k} \to \varphi \text{ (strongly)}$
$\lambda \in \sigma_{\text{cont}} (\opA)$ | $\text{\emph{All} Weyl sequences verify } \phi_n \rightharpoonup 0 \text{ (weakly)}$
$\lambda \in \sigma_{\text{p}}(\opA)$ | $\exists \text{ a Weyl sequence with a weak limit different from 0, i.e. } \mathop{\mathrm{Ker}}(\opA - \lambda) \neq 0.$
"""

# ╔═╡ 1b178edd-98f9-4533-b3df-38efcd026333
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
However, in the following paragraph we introduce the *form domain* $Q(\opA)$, an extension of $D(\opA)$

If $\opA$ is bounded from below then there exists $\alpha > -\infty$ such that
```math
	q_A (u) \geq \alpha \| u \| ^2
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
    \| u \|_A \equiv q_A(u) > 0 && u \in D(\opA) \tag{Energy norm}
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
Where in each case the embedding is dense and continuous (Lewin [*Thérorie spectrale et mécanique quantique*](https://doi.org/10.1007/978-3-030-93436-1) 2022, Theorem 3.10).

"""

# ╔═╡ 819cccb8-fe5b-4bf6-9998-33ce4ab30819
md"""

The importance of the form domain is that is allows for a weak formulation :

!!! note "Theorem 5"
	Let $\opA$ be self-adjoint. The following are equivalent
	
	1.  $\varphi \in D(\opA)$ and $\opA \varphi = \lambda \varphi$. (*Strong formulation*)
	
	2.  $\varphi \in Q(\opA)$ and $\exists \lambda \in \mathbb R$ s.t.
	    $a_A(f,\varphi) = \lambda \langle f, \varphi \rangle$ for all $f \in Q(\opA)$. (*Weak
	    formulation*)



"""

# ╔═╡ 68f786f1-45a4-475b-a0bb-37636819ddf2
md"""
Based on this formulation we obtain, again based on Weyl sequences :

!!! note "Theorem 6"
	Let $\opA$ be a self-adjoint operator with form domain $Q(\opA) \subset \hilbert$, then
	the *bottom of the essential
	spectrum* $\Sigma(\opA)$ is uniquely defined by 
	```math 
	\begin{align}
	        \Sigma(\opA) \coloneqq \min \sigma_{\text{ess}} (\opA) = \min_{\substack{(v_n) \in Q(\opA)^\mathbb N \\ \| v_n \| = 1 \\ v_n \rightharpoonup 0}} \liminf_{n \to \infty} q_A(v_n)
	    
	\end{align}
	```
	with the convention $\Sigma(\opA) = \infty$ if and only if $\sigma_{\text{ess}}(\opA) = \varnothing$.

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

# ╔═╡ f9f9e20a-0beb-46c6-863b-789d04d0b2af
md"""
Finally, we can state the equivalent the min-max principle that we already saw for matrices, extended to the case of operators. 

!!! note "Theorem 7 (Courant-Fisher)"
	Let $\opA$ be self-adjoint and bounded
	from below with form domain $Q(\opA) \leqslant \hilbert$. 
	Then
	```math 
	\begin{align}
	        \mu_k (\opA) \coloneqq \inf_{\substack{W \subset Q(\opA) \\ \dim (W) = k}} \max_{\substack{\varphi \in W \\ \| \varphi \|_\hilbert = 1}} q_A(\varphi) = \inf_{\substack{W \subset Q(\opA) \\ \dim (W) = k}} \max_{0 \neq \varphi \in W} R_A(\varphi)
	    
	\end{align}
	```
	 is equal to
	
	1.  The k$^{\textsf{th}}$ eigenvalue (counting multiplicities) of $\opA$ if $\opA$ has at
	    least $k$ eigenvalues below $\Sigma(\opA).$
	
	2.  Otherwise, $\mu_k(\opA) = \Sigma(\opA)$.

	We re-introduced the *Rayleigh quotient*, this time for operators 
	```math
		R_\opA (\varphi ) \coloneqq \frac{q_\opA (\varphi)}{\| \varphi \|^2_\hilbert }
	```

One way to check for the existence of $k$ eigenvalues below $\Sigma(\opA)$ is :

!!! tip ""
	If there is a subspace $W \leq Q(\opA)$ with $\text{dim} W = k$, such that 
	```math 
	\max_{0 \neq \varphi \in W} R_\opA (\varphi) < \Sigma(\opA)
	```
	then there are at least k eigenvalues.

In other words, as long as $\mu_k(\opA) < \Sigma(\opA)$, the situation is as for Hermitian matrices.
As soon as $k$ eigenvalues have been found, we can no longer learn anything about $\sigma(\opA)$ by increasing the subspace size.

This leads to a severe restriction for the numerical approximation of eigenspectra, as we are in most cases restricted to below $\Sigma(\opA)$.
However, in this case we are able to estimate eigenspectra, and obtain error bounds in a manner similar to that of matrices.
	
"""

# ╔═╡ 9ccec1c6-55e8-4db6-aa70-18dabc87eca8
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
      (M_S^\opA)_{ij} = q_A(\chi_i, \chi_j) \cong \langle \chi_i , \opA \chi_j \rangle
  \end{align}
  ```
  Note that the second equality is only fully correct for $S \subset D(\opA)$. In other words 
  ```math 
  \begin{align}
      \lambda_k (\opA) = \mu_k(\opA) \leq \lambda_k (M_S^\opA)
  \end{align}
  ```
"""

# ╔═╡ f743838b-4204-432e-a81e-77af58737cba
md"""
- In practice one then keeps augmenting the subspace $S$, e.g. by increasing the approximation basis $\mathbb B$ following the principle 
  ```math
  	\lambda_k(\opA) = \mu_k(\opA) = \inf_{\substack{S \subset Q(\opA) \\ \dim Q(\opA) \geq k}} \lambda_k (M_S^\opA)
  ```

- Notice, however, how this technique is *unable to approximate above $\Sigma(\opA)$.*
  In other words, *some eigenvalues* of $\opA$, those embedded in $\sigma_{\text{ess}} (\opA)$ *cannot be found* using this projection technique.
"""

# ╔═╡ e80b1728-592a-4a97-b4ab-20b3ef8604c4
md"""
---
#### Illustration

Eigenenergies of the 1D gaussian potential well, i.e. the spectrum of the following Schrödinger operator
```math
\mathcal H = - \frac1{2} \frac{d^2}{dx^2} - A e^{-(x/\sigma)^2}
```
"""

# ╔═╡ 843f48c0-b38f-4c94-8539-cb0ca1a85252
md"""
y axis limits (slide to zoom in/out) $\qquad$ $(@bind ylim PlutoUI.Slider([2000,1000,500,200,100,50,25]; default=1000, show_value=true))
"""

# ╔═╡ f41ef893-8910-4e6b-9f08-cccba464bac1
#plotting 
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

# ╔═╡ 59832711-e53c-405c-b83d-3dc4fc05fa53
md"""

---
## Kato-Temple bound

To close our discussion on the spectra of operators, we want to rationalize why the Kato-Temple bound can still be employed for the eigenpairs below the essential spectrum.
We restate :

"""

# ╔═╡ f950d40e-6ec7-4ee5-81c3-95314b52d029
md"""
!!! note "Theorem 8 (Kato-Temple bound)"
	Let $\opA$ be a self-adjoint operator,
	$\tilde \varphi \in D(\opA)$ with $\| \tilde \varphi \| = 1$, $\tilde \lambda = q_A(\tilde \varphi)$,
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

# ╔═╡ 3bee37f2-ec13-4799-b61e-af6879e3b7cc
md"""

Note that this results follows automatically if we can prove the following :

!!! note "Theorem 9 (Temple's Inequality)"
	Let
	$\opA, \tilde \varphi, \tilde \lambda, r$ as in
	Theorem 1. Suppose $\alpha, \beta \in \mathbb R$ with
	$\alpha < \tilde \lambda < \beta$ and
	$(\alpha, \beta ) \cap \sigma(\opA) = \{\lambda \}$. Then 
	```math 
	\begin{align}
	        \tilde \lambda - \frac{\| r\|^2}{\beta - \tilde \lambda} \leq \lambda \leq \tilde \lambda + \frac{\| r \|^2}{\alpha - \tilde \lambda}
	    
	\end{align}
	```
"""

# ╔═╡ 5731ef65-9c89-4859-b891-f5c2606c5772
md"""

To prove this, we need a bit of spectral calculus.
Like in the finite-dimensional case, we have :

!!! note "Theorem 10 (Cauchy's Formula)"
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

# ╔═╡ 6739d23c-1d13-4cad-93e9-b370b6da3270
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

# ╔═╡ 8c52943d-a2d6-4fa7-bae5-2c7c6fe54140
md"""

Illustration of Cauchy's formula.

---

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

# ╔═╡ 06fc9dd1-601c-4d4c-9f86-ddb9be7fc13f
md"""
##### Finite dimensions

We study the finite dimensional case with a matrix
$M \in \mathbb R^{d \times d}$ with distinct eigenvalues
$\lambda_1 < \dots < \lambda_m ( m \leq d)$.
Note that we assume single eigenvalues for simplicity.
In this setting
```math 
\begin{align}
    P^M(\lambda) = \mathbf 1 _{(-\infty,\lambda]} (M) = \bigoplus_{\lambda_i \leq \lambda} \mathop{\mathrm{Ker}}(M - \lambda_i)
\end{align}
```
where $\text{Ker}(M - \lambda_i)$ is th eigenspace of eigenvalue $\lambda_i$.
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

# ╔═╡ 66c953ac-d572-469a-94f1-67cc473c5704
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
More generally, for any (measurable) function $f$
and $v,w \in D(\opA)$ 
```math 
\begin{align}
    \langle v , f(\opA) w \rangle = \int_\mathbb R f(\lambda) \ d \langle v, P^\opA (\lambda) v \rangle.
\end{align}
```
which can in turn be used to define a self-adjoint operator $f(\opA)$ (functional calculus).
Without going into details with this the following result is plausible :
"""

# ╔═╡ 2de385c1-8ad1-4cbf-a997-473daf2bcc33
md"""
!!! note "Lemma 11" 
	If $f$ is a polynomial and $\opA$ is self-adjoint, then 
	```math
		\sigma(f(\opA)) = \{ f(\lambda) \mid \lambda \in \sigma(\opA) \}.
	```

With this result, we can finally prove
"""

# ╔═╡ 85c4f6c2-23dc-45b8-b733-f267b03cdde5
md"""
> **Proof of Theorem 9.** 
> We note that our assumptions imply $(\lambda , \beta) \cap \sigma(\opA) = \varnothing$ such that, using Lemma 11 and the polynomial $f(x) = (x - \beta) (x - \lambda)$, we obtain that $(\opA - \beta) (\opA - \lambda) \geq 0$, i.e. that the operator $(\opA - \beta) (\opA - \lambda)$ only has non-negative spectrum.
> Therefore 
> ```math
> \begin{align}
> 	0 &\geq \langle \tilde \varphi, (\opA - \beta) (\opA - \lambda) \tilde \varphi \rangle
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

# ╔═╡ 4d4d6fa7-99c5-4cb7-8d3e-830634a37daf


# ╔═╡ e1d81761-c3c4-4bef-9e14-96bdcf5c8eba
TableOfContents()

# ╔═╡ 5c5e2ead-9487-42a4-b049-ce7b90cc97d2
begin
	Sidebar(elts...; location="upper right") = @htl("""
	<aside class="sidebar" style='top: 305px;right: 17px;'>$elts</aside>
	
	<style>
	aside.sidebar {
		position: fixed;
		max-width: min(30%, 300px, calc(100vw - 750px));
		padding: 0.4rem;
		border-radius: 10px;
		max-height: calc(100vh - 370px);
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
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"

[compat]
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
Plots = "~1.40.4"
PlutoUI = "~0.7.55"
QuadGK = "~2.9.4"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "224cc730c0a759dc6dfd0d03398245e55847eb79"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

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
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "8e59b47b9dc525b70550ca082ce85bcd7f5477cd"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.5"

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

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

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
git-tree-sha1 = "60b1194df0a3298f460063de985eae7b01bc011a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.1+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

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

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

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

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "af81a32750ebc831ee28bdaaba6e1067decef51e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

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

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "442e1e7ac27dd5ff8825c3fa62fbd1e86397974b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.4"

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

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "68723afdb616445c6caaef6255067a8339f91325"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.55"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "02148a0cb2532f22c0589ceb75c110e168fb3d1f"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "21.9.0+0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

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
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "tectonic_jll"]
git-tree-sha1 = "79e2d29b216ef24a0f4f905532b900dcf529aa06"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "71509f04d045ec714c4748c785a59045c3736349"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.7"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

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

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

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

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

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
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─4004f2e8-5738-47e2-9d3c-8b68b896b9a5
# ╟─1478c459-0933-4c3f-8e84-86d051bea4b3
# ╟─21814beb-9d27-44c7-a41f-1191d856d1b7
# ╟─d2b95a1b-5e4c-4e85-8a6e-d893c9c3cbde
# ╟─206bf4ec-9bdd-4141-8ee4-570f077c06d0
# ╟─390974fb-36d3-44e1-8543-24a16e400512
# ╟─95acb6bd-42d7-4ede-b425-a91ed9640d9a
# ╟─20958c09-7372-4c0b-bd3a-3011758cba2e
# ╟─2f377ee6-e23e-432f-aeb6-9097fa260027
# ╟─d2bf63fb-6035-49b2-9182-c606bb3643af
# ╟─a803317a-07c5-4f95-bc45-4c24f3971385
# ╟─f2b8a3ba-71e9-48a9-b77d-91aa844f3da3
# ╟─dda55c4f-1a71-477b-b460-0c67629c2904
# ╟─8e5dd330-181a-4279-b54b-47be263e13a9
# ╟─e60a51b3-88ed-4e03-9937-4c86773e33e8
# ╟─1b178edd-98f9-4533-b3df-38efcd026333
# ╟─819cccb8-fe5b-4bf6-9998-33ce4ab30819
# ╟─68f786f1-45a4-475b-a0bb-37636819ddf2
# ╟─f9f9e20a-0beb-46c6-863b-789d04d0b2af
# ╟─9ccec1c6-55e8-4db6-aa70-18dabc87eca8
# ╟─f743838b-4204-432e-a81e-77af58737cba
# ╟─e80b1728-592a-4a97-b4ab-20b3ef8604c4
# ╟─565cdda1-d844-4bed-8f48-fe57c1428136
# ╟─843f48c0-b38f-4c94-8539-cb0ca1a85252
# ╟─f41ef893-8910-4e6b-9f08-cccba464bac1
# ╟─59832711-e53c-405c-b83d-3dc4fc05fa53
# ╟─f950d40e-6ec7-4ee5-81c3-95314b52d029
# ╟─3bee37f2-ec13-4799-b61e-af6879e3b7cc
# ╟─5731ef65-9c89-4859-b891-f5c2606c5772
# ╟─6739d23c-1d13-4cad-93e9-b370b6da3270
# ╟─8c52943d-a2d6-4fa7-bae5-2c7c6fe54140
# ╟─06fc9dd1-601c-4d4c-9f86-ddb9be7fc13f
# ╟─66c953ac-d572-469a-94f1-67cc473c5704
# ╟─2de385c1-8ad1-4cbf-a997-473daf2bcc33
# ╟─85c4f6c2-23dc-45b8-b733-f267b03cdde5
# ╟─4d4d6fa7-99c5-4cb7-8d3e-830634a37daf
# ╟─e1d81761-c3c4-4bef-9e14-96bdcf5c8eba
# ╟─5c5e2ead-9487-42a4-b049-ce7b90cc97d2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
