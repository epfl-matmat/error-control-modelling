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

# ╔═╡ 80ffd650-d95b-471a-9d2f-e2bc781c322e
begin
	using LinearAlgebra
	using LinearMaps
	using PlutoUI
	using PlutoTeachingTools
	using Printf
	using Plots
	using HypertextLiteral
end

# ╔═╡ 2a26de48-c3b6-405b-8c07-bae036dd805d
using DFTK

# ╔═╡ c498a8aa-df5c-4536-b54b-db5912fb411c
md"""
# 6 Overview of diagonalisation algorithms

This lecture will provide a selective overview of standard diagonalisation algorithms, contrasting their respective scope of applicability and the main ingredients. Our main angle will be to understand questions related to numerical stability in the key ingredients of the algorithms. Instead of taking a proof-guided approach, we will mostly follow a computational approach using the increased-precision and interval techniques we discussed in previous lectures.

A more comprehensive treatment of the topic can be found in the book [Numerical Methods for Large Eigenvalue Problems](https://epubs.siam.org/doi/book/10.1137/1.9781611970739) by Youssef Saad as well as the [Lecture notes on Large Scale Eigenvalue Problems](https://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf) by Peter Arbenz.

In our discussion we will always take $A \in \mathbb{C}^{N\times N}$ to be a Hermitian matrix and we will seek approximations to the eigenpairs $(\lambda_i, x_i) \in \mathbb{R} \times \mathbb{C}^N$, i.e.
```math
A x_i = \lambda_i x_i
```
where we impose an ordering $\lambda_1 \leq \lambda_2 \leq \cdots \leq \lambda_N$.
"""

# ╔═╡ 0f8dc304-6ade-4fbe-8b64-ca368defa2ba
TODO(
	"""
We should be careful to not make mistakes in statements with the multiple orderings for eigenvalues at hand (largest eigenvalue as in Courant-Fisher, vs largest-absolute eigenvalue as in Power method and related methods)
    - Apart from fixing typos around this, we might want to add a dedicated "warning" / reminder at the beginning of this lecture.

	"""
)

# ╔═╡ dddf152e-b1a3-470d-93b7-d213544f0215
TableOfContents()

# ╔═╡ 5729f2a4-7557-4e5c-b230-31faf6286ff3
md"""
## Single vector iteration

### Power method

We already discussed the **power method** in previous lectures and exercises as a way to obtain the eigenvector $x_J$ corresponding to the largest-abs eigenvalue. Via the Rayleigh quotient the corresponding eigenvalue $\lambda_J$ then becomes accessible as well.

**Remark**. Keep in mind that the largest-abs eigenvalue $\lambda_J$ is in general not the same as the largest eigenvalue $\lambda_N$, in case a negative eigenvalue is dominating.

For completeness, our algorithmic recipe to obtain the eigenvector to tolerance $\varepsilon$ was:

> 0. Start from a vector $x^{(0)}$, which has non-zero overlap with $x_n$.
> 1. Compute the matrix-vector product: $\tilde{x}^{(i)} \leftarrow A x^{(i)}$
> 2. Normalise: $x^{(i+1)} \leftarrow \tilde{x}^{(i)} / \| \tilde{x}^{(i)} \|$
> 3. If $\| x^{(i+1)} - x^{(i)} \| < \varepsilon$, exit the procedure, else  $i \leftarrow i+1$ and return to 1.

Below we provide a slightly extended iteration, which additionally employs the Rayleigh quotient to estimate the eigenvalue and the residual
```math
r^{(i)} = A x^{(i)} - \lambda^{(i)} x^{(i)}
```
in each iteration,
which is used to control convergence instead of the difference of eigenvectors.
"""

# ╔═╡ f8827f52-cba1-433f-9c14-9a0b621b2c35
function power_method(A; x=randn(eltype(A), size(A, 2)), tol=1e-6, maxiter=100, verbose=true)
	T = real(eltype(A))

	eigenvalues    = T[]
	residual_norms = T[]
	λ = NaN
	for i in 1:maxiter
		normalize!(x)
		Ax = A * x
		λ = x'Ax
		norm_r = norm(Ax - λ * x)

		verbose && @printf "%3i %8.4g %8.4g\n" i λ norm_r
		push!(eigenvalues, λ)
		push!(residual_norms, norm_r)
		norm_r < tol && break

		x = Ax
	end
	converged = residual_norms[end] < tol
	!converged && verbose && @warn "Power method not converged."
	
	(; λ, x, converged, eigenvalues, residual_norms)
end

# ╔═╡ eb66e66c-278f-465e-a56d-cac7432723f7
@bind logδ PlutoUI.Slider(-1:0.1:2; show_value=true, default=1)

# ╔═╡ 23f9fd84-0b58-4ade-a72c-b3e303a79480
Apower = Diagonal([1, 30 - 10^logδ, 30 + 10^logδ])

# ╔═╡ bfff1b7f-bf26-4373-a499-014b871067ad
let
	(; residual_norms) = power_method(Apower, verbose=false, tol=1e-6, maxiter=200)
	plot(residual_norms, yaxis=:log, ylims=[1e-6, 1e2], label="residual_norm", xlims=[0, 200], lw=2)

	λs = sort(eigvals(Apower))
	plot!(1:200, abs(λs[end-1] / λs[end]).^(1:200); label="Expected rate", ls=:dash, c=1, lw=2)
end

# ╔═╡ 72264fdd-d151-41b3-abd7-3fd6b0b4fb84
md"""
From playing with the slider we notice that when the two largest-abs eigenvalues get too close, the convergence deteriorates. This can be understood using the following theorem:

**Theorem 5.1** Assume that $|\lambda_J| \neq |\lambda_{J-1}|$, i.e. that the largest eigenvalue is simple. Then the initial vector $x^{(0)}$ has either no overlap with the eigenvector $x_J$ or the sequence $x^{(t)}$ generated by the power method converges to $x_J$.

This is a well-known result for which proofs can be found e.g. in the book by Saad or [on wikipedia](https://en.wikipedia.org/wiki/Power_iteration). The key outcome of the proof is that the convergence rate depends on the ratio
```math
\frac{|λ_{J-1}|}{|λ_{J}|}
```
thus deteriorates as $|λ_{J-1}|$ and $|λ_{J}|$ are too close, e.g. if $λ_J$ is a multiple eigenvalue.
"""

# ╔═╡ fa98f7d8-bc72-4ac8-b46a-2456e37e42de
md"""
### Spectral transformations

A natural question to ask now how this method can be extended, e.g. to compute the smallest eigenvalue or eigenvalues in the middle of the spectrum of a matrix. The answer to this provide the following well-known **spectral transformations**:

**Theorem 5.2** (Spectral transformations)**.** Assume $(λ_i, x_i)$ is an eigenpair of $A$, then
1. If $A$ is invertible, then $(\frac{1}{λ_i}, x_i)$ is an eigenpair of $A^{-1}$.
2. For every $σ \in \mathbb{R}$ $(λ_i - σ, x_i)$ is an eigenpair of $A - σ I$.
3. If $A - σ I$ is invertible, then $(\frac{1}{λ_i - σ}, x_i)$ is an eigenpair of $(A - σI)^{-1}$.

**Proof:** 1. and 2. follow from a simple direct calculation. 3. is a corollary of 1 & 2.
"""

# ╔═╡ 291cb850-5089-480e-b803-1e2a20852b95
md"""
Part 1 of this theorem is the justification for using the inverse power method to approximate the smallest eigenpair of $A$ as done already in one of the previous exercises. Part 3 of this theorem, called **shift and invert** provides us with a recipe to compute an arbitrary isolated eigenvalue.

To illustrate this we consider the matrix
"""

# ╔═╡ 98c91264-19cc-41e5-951a-a201601fbfc5
Ashift = [0.4 -0.6 0.2; -0.3 0.7 -0.4; -0.1 -0.4 0.5]

# ╔═╡ aea3267b-1164-418c-98a8-70cc37d1c959
md"""which has eigenvalues"""

# ╔═╡ e33b4164-b980-446f-9398-9ec68cf224f9
eigvals(Ashift)

# ╔═╡ 7661b0fa-839a-47a1-a0ea-a89e8d4abc3a
md"""Suppose we don't know the exact eigenvalues, but we know that the eigenpair $(λ_i, x_i)$ we are interested in has an eigenvalue around $0.4$. Then we set $σ = 0.4$ and employ the matrix $(A - σ I)^{-1}$ in the `power_method`, which as per our previous discussion will converge to the eigenvalue of  $(A - σ I)^{-1}$. The eigenvalues of the latter matrix are $\left\{\frac{1}{λ_i - σ} \, \middle| \, i =1,\ldots,3 \right\}$. Clearly, the largest-abs one of these eigenvalues is attained for the $λ_i$ closest to $σ$.
From Theorem 5.2 we additionally know that the eigenvector of $(A - σ I)^{-1}$ associated to $\frac{1}{λ_i - σ}$ is exactly the eigenvector of $A$ associated to $λ_i$. Therefore, we have found a way to compute the *exact* eigenpair $(λ_i, x_i)$ based on the approximate eigenvalue $σ = 0.4$.

Computationally we do this as such:
"""

# ╔═╡ f923e4ba-1833-4b2e-968d-8527fc70720d
let
	λ₂ = eigvals(Ashift)[2]  # Exact targeted eigenvalue
	
	σ = 0.4  # Our approximation to the eigenvalue of interest
	shifted  = Ashift - σ * I        # Shift the matrix
	factorised = factorize(shifted)  # Factorise to obtain fast \ operation

	# Use a LinearMap to represent the inverse lazily
	# (i.e. `inverted * v = factorised \ v`)
	shift_and_invert = InverseMap(factorised)

	# Run power method on (A - σI)⁻¹
	(; λ, x) = power_method(shift_and_invert)

	# Compute eigenvalue using the Rayleigh quotient of the original matrix
	λ_rayleigh = dot(x, Ashift, x)

	println()
	println("Power method converged to eigenvalue: $λ")
	println("Exact value for 1/(λ₂ - σ) is:        $(1/(λ₂ - σ))")
	println("Computed eigenvalue of Ashift:        $λ_rayleigh")
	println("Exact eigenvalue using eigvals:       $λ₂")
end

# ╔═╡ 6a049358-0e6c-4ea9-8950-52f5bdda287e
md"""
## Rayleigh quotient iterations

We already noted earlier that the Rayleigh quotient is a natural approximation for the eigenvalue corresponding most closely to an approximate eigenvector $v$. A natural idea is thus to replace the shift $σ$ in the shift-and-invert strategy by the Rayleigh quotient $R_A(x^{(i)})$ of the current iterate.

This is theoretically backed up by the following result.

**Proposition 5.3.** Let $A \in \mathbb{C}^{N\times N}$ be a Hermitian matrix and $v\in \mathbb{C}^N$. Then
```math
R_A(v) = \frac{v^H A v}{v^H v} = \text{argmin}_{α \in \mathbb{C}} \ \| A v - α v \|.
```

**Proof:** We denote $ρ = R_A(v)$ and consider the function
```math
f(z) = \|A v - (ρ + z) v\|^2
```
where $z \in \mathbb{C}$. We observe that
```math
\begin{aligned}
f(z) &= v^H A^2 v - (2ρ + z + \overline{z})\, v^H A v + |ρ + z|^2 v^H v \\
&= v^H A^2 v - 2\big(ρ + \text{Re}(z)\big) \, v^H A v + \big( ρ^2 + 2ρ \text{Re}(z) + |z|^2 \big) v^H v\\
&= v^H A^2 v + \big(-2ρ^2 - 2ρ\text{Re}(z)\big)v^H v + \big( ρ^2 + 2ρ \text{Re}(z) + |z|^2 \big) v^H v\\
&= v^H A^2 v + (- ρ^2 + |z|^2) v^H v
\end{aligned}
```
Since $f(z) \geq 0$ by construction and $v^H v \geq 0$,
this function is minimal for $z = 0$,
which concludes the proof. 

With this setup we obtain Rayleigh-quotient iteration (RQI):
"""

# ╔═╡ b0e6abc4-7581-4a7b-9160-792d2824a1cf
function rayleigh_quotient_iteration(A; x=randn(eltype(A), size(A, 2)),
                                     tol=1e-6, maxiter=100, verbose=true)
	T = real(eltype(A))

	eigenvalues    = T[]
	residual_norms = T[]
	λ = NaN
	for i in 1:maxiter
		normalize!(x)
		Ax = A * x
		λ = x'Ax
		norm_r = norm(Ax - λ * x)

		verbose && @printf "%3i %8.4g %8.4g\n" i λ norm_r
		push!(eigenvalues, λ)
		push!(residual_norms, norm_r)
		norm_r < tol && break

		# Note: In the Power method this was u = Au
		x = (A - λ*I) \ x
	end
	converged = residual_norms[end] < tol
	!converged && verbose && @warn "RQI not converged."
	
	(; λ, x, converged, eigenvalues, residual_norms)
end

# ╔═╡ 6be8652b-525c-4efa-bcc7-b72285ea2d98
md"""RQI is a very fast algorithm, e.g. notice"""

# ╔═╡ ab7f07f6-20f3-465f-ae46-e694100541b8
Arqi = Diagonal(randn(100));

# ╔═╡ 1f9b91e3-28f4-4e06-8a3d-70c84432007b
x0_rqi = randn(100);

# ╔═╡ ce876450-eea7-4088-9aca-ed3ee462cb09
rayleigh_quotient_iteration(Arqi, x=x0_rqi)

# ╔═╡ 0f7709e7-fec9-4176-bf58-96f2d891e8b4
md"""Notice, how the convergence in the last few steps is *very fast*, meaning that multiple orders of magnitude in the residual norm are gained in each step.

The RQI has fantastic rates of convergence. However, it has two disadvantages:
- The main issue is that *at each RQI iteration* the matrix `A - λ*I` is different and needs to be freshly factorised in order to compute the `\` operation. In other words we pay a full Gaussian elimination each RQI iteration. This in contrast to the *shift-and-invert* power method, where factorisation was computed once and for all for `A - σ*I` *before* the `power_method` was called.
- The second point to note is that it is hard to control *to which eigenpair* the RQI converges. Try re-running the above a few times by updating the initial guess `u0_rqi`. Moreover, it is possible for the RQI to *not* converge --- albeit the set of possible initial guesses that leads to a non-converging RQI is of measure zero. It is thus most useful if one already has a decent guess for the eigenvector.

These observations are summarised in the following
"""

# ╔═╡ 73e61d6c-0762-41d9-a0c1-c2a64926a40a
md"""
**Theorem 5.4.** For all starting vectors `u`, except a set of measure zero, the Rayleigh quotient iteration converges to an eigenpair. When it converges, the convergence rate is cubic: If $λ_k$ is an eigenvalue to $A$ and `u` sufficiently close to the eigenvector $x_k$ corresponding to $λ_k$, then as $i \to \infty$
```math
|λ^{(i+1)}_k - λ_k| = O\big(|λ^{(i)}_k - λ_k|^3\big)
```
and
```math
|x^{(i+1)}_k - x_k| = O\big(|x^{(i)}_k - x_k|^3\big),
```
where $(λ^{(i)}_k, x^{(i)}_k)$ is the approximation to $(λ^{(i)}_k, x^{(i)}_k)$ the algorithm produces in iteration $i$. Here we choose $x_k$ among $\pm x_k$, such that $\langle x^{(i)}_k, x_k \rangle > 0$.

"""

# ╔═╡ c3c91935-60f7-4c4a-b43c-973566d58367
md"""
## Subspace methods

So far we restricted ourselves to the rather basic goal of solving only for a single eigenvector --- essentially by iterating only a single vector. We already saw that some eigenproblems are challenging to solve in this simple setting, namely exactly if the largest-abs eigenvalue has only a small or even a zero gap. 

A natural generalisation to the power-method scheme is to iterate not one, but multiple vectors. Doing this naively is a little flawed as then simply each vector will converge to the eigenvector(s) with the largest-abs eigenvalue. Instead, what is needed is to ensure *in each step* that the obtained eigenvectors are orthogonal to each other.

One way to achieve this is QR-factorisation, i.e. the factorisation $A = Q R$ into an orthogonal matrix $Q$ and an upper triangular matrix $R$. The orthogonalisation thus consists of simply keeping the $Q$ matrix:
"""

# ╔═╡ 41265945-67e6-4e34-8604-20faa6f475b7
ortho_qr(A) = Matrix(qr(A).Q)

# ╔═╡ 8937b497-99f0-4d8f-998a-e2aad0d5aeb9
md"""Built into a power iteration procedure and employing multiple vectors, this gives rise to the **subspace iteration** method (in this example iterating on $k=2$ vectors by default):"""

# ╔═╡ f64ab474-86f6-4048-b013-3c8969d2e125
function subspace_iteration(A; X=randn(eltype(A), size(A, 2), 2), ortho=ortho_qr,
                            tol=1e-6, maxiter=100, verbose=true)
	T = real(eltype(A))

	eigenvalues    = Vector{T}[]
	residual_norms = Vector{T}[]
	λ = T[]
	for i in 1:maxiter
		X = ortho(X)

		AX = A * X
		λ = diag(X'*AX)  # Rayleigh quotient evaluation in all vectors
		push!(eigenvalues, λ)

		residuals = A * X - X * Diagonal(λ)  # compute all residuals
		norm_r = norm.(eachcol(residuals))   # compute all residual norms
		push!(residual_norms, norm_r)

		verbose && @printf "%3i %8.4g %8.4g\n" i λ[end] norm_r[end]
		maximum(norm_r) < tol && break
		
		X = AX
	end
	
	(; λ, X, eigenvalues, residual_norms)
end

# ╔═╡ d2ef7881-bd1b-4bd3-aee4-025c7da5ba71
@bind logsubδ PlutoUI.Slider(-1:0.1:1; show_value=true, default=1)

# ╔═╡ d21b8a45-7129-4038-82ac-250436ff8583
Asubspace = Diagonal([150., 100. + 10^logsubδ, 100., 50., 40., 30.])

# ╔═╡ 6f81e36d-7892-42ec-aeb6-91fca5b6399a
let
	(; residual_norms) = subspace_iteration(Asubspace; X=randn(size(Asubspace, 2), 2),
	                                        verbose=false, tol=1e-6, maxiter=200)
	norms_first  = [n[1] for n in residual_norms]
	norms_second = [n[2] for n in residual_norms]

	p = plot(; yaxis=:log, ylims=[1e-6, 1e2], xlims=[0, 200])
	plot!(p, norms_first,  label="residual first eigenpair",  c=1, lw=2)
	plot!(p, norms_second, label="residual second eigenpair", c=2, lw=2)

	λs = sort(eigvals(Asubspace))
	plot!(p, 1:200, abs(λs[end-1] / λs[end]).^(1:200);
	      label="Rate first", ls=:dash, c=1, lw=2)
	plot!(p, 1:200, abs(λs[end-2] / λs[end-1]).^(1:200);
	      label="Rate second", ls=:dash, c=2, lw=2)
end

# ╔═╡ 137f3b99-d031-4df4-868e-9b5e4d342a83
md"""
We see that each eigenpair $(λ^{(i)}_k, x^{(i)}_k)$ converges with its own speed, which depends on the ratio
```math
\frac{|λ_{k+1}|}{|λ_{k}|},
```
which can also be proven rigorously.
"""

# ╔═╡ 357cd9c6-eb57-48ca-afe5-a752addc6461
md"""
## Projection methods

The above `subspace_iteration` addresses the question of how to compute multiple eigenpairs at once. However, compared to Rayleigh quotient iteration it feels a little like a step back since we are no longer making use of the Rayleigh quotient. Notably, when computing the approximate eigenpairs we are considering each iterated vector separately and allow no "interaction" between them, even though the first iterated vector might well contain a contribution to the eigenvector, which the second iterated vector attempts to approximate. Moreover, we are right now tied to using exactly one iterated vector per eigenvector we want to approximate. We will now see how projection methods allow us to overcome both these restrictions.

The general idea is to somehow construct an $m$-dimensional subspace $\mathcal{S} \subset \mathbb{C}^N$ with $m \ll N$ (details will be discussed later) and fully solve the eigenvalue problem *projected* into this subspace. We consider the basic eigenproblem for $A \in \mathbb{C}^{N \times N}$ Hermitian:
```math
A x_i = λ_i x_i \qquad \text{$(\ddagger)$}
```
An orthogonal projection technique now seeks approximate eigenpairs $(\tilde{\lambda}_i, \tilde{x}_i)$ with $\tilde{\lambda}_i \in \mathbb{R}$, $\tilde{x}_i \in \mathcal{S}$, such that the following **Galerkin condition**
```math
\big(A \tilde{x}_i - \tilde{λ}_i \tilde{x}_i \big) \perp \mathcal{S} \qquad \text{(G)}
```
is satisfied. Equivalently one may enforce
```math
\left\langle v, A \tilde{x}_i - \tilde{λ}_i \tilde{x}_i \right\rangle = 0 \qquad \forall v \in \mathcal{S}. \qquad \text{$(\ast)$}
```
Suppose now that $\{v_1, v_2, \ldots, v_m\}$ is an orthonormal basis for $\mathcal{S}$ and denote by $V = (v_1, v_2, \ldots, v_m) \in \mathbb{C}^{N\times m}$ the matrix with these basis vectors as columns. Expanding $\tilde{x}_i$ in this basis, i.e. introducing a vector $\tilde{y}_i \in \mathbb{C}^m$ with
```math
\tilde{x}_i = V \tilde{y}_i
```
equation $(\ast)$ becomes
```math
\left\langle v_j, A V \tilde{y}_i - \tilde{λ}_i V \tilde{y}_i \right\rangle = 0 \qquad \forall j = 1, \ldots, n.
```
or (since $V^H V = I \in \mathbb{C}^{n \times n}$):
```math
\left(V^H A V\right) \, \tilde{y}_i = \tilde{λ}_i \tilde{y}_i \qquad \text{$(\S)$}
```
Note that this is an $m$-dimensional eigenproblem, which yields an approximate eigenpair $(\tilde{λ}_i, \tilde{y}_i)$ to the original $N$-dimensional eigenpair $(\ddagger)$ above. Thus while $(\ddagger)$ might be too large to be solved by a dense linear algebra routine like `eigen`, the Galerkin approach has allowed us to reduce this to an $m$-dimensional problem, which is easier to treat with `eigen`. Algorithmically this technique is known as the **Rayleigh-Ritz procedure**, which can be summarised as follows:

> 1. Find an orthonormal basis $\{v_i\}_{i=1,\ldots,m}$, spanning a subspace $\mathcal{S}$.
> 2. Compute $A_V = V^H A V$, with $V$ containing the basis vectors of $\mathcal{S}$ as columns.
> 3. Use `eigen` to compute $m$ eigenpairs $(\tilde{λ}_i, \tilde{y}_i)$ of $A_V$.
> 4. Compute the approximate eigenvectors $\tilde{x}_i = V \tilde{y}_i$ of $A$, which yields the approximate eigenpairs of $A$ as $(\tilde{λ}_i, \tilde{x}_i)$.

Note, that the low-dimensional pairs $(\tilde{λ}_i, \tilde{y}_i)$ are exact eigenpairs of the subproblem matrix $A_V$. The corresponding full-dimensional pairs $(\tilde{λ}_i, \tilde{x}_i)$ are usually called **Ritz pair** with $\tilde{λ}_i$ called **Ritz values** and $\tilde{x}_i$ a **Ritz vector**. 

One way to put this into practice is to employ the ideas of the `subspace_iteration` discussed above. Recall that this method generates in each iteration a set of orthogonal vectors `U` as guesses for the eigenpairs to be determined. These vectors span a subspace, which can be employed in a Rayleigh-Ritz procedure, leading to the projected subspace iterations:
"""

# ╔═╡ 7237c2e9-c0f3-4cd8-ad1f-6d5cbff35354
function projected_subspace_iteration(A; tol=1e-6, maxiter=100, verbose=true,
                                      V=randn(eltype(A), size(A, 2), 2),
                                      ortho=ortho_qr)
	T = real(eltype(A))

	eigenvalues    = Vector{T}[]
	residual_norms = Vector{T}[]
	λ = T[]
	Y = nothing
	for i in 1:maxiter
		V = ortho(V)

		AV = A * V
		λ, Y = eigen(Hermitian(V' * AV))  # Notice the change to subspace_iteration
		                                  # This is the Rayleigh-Ritz step
		push!(eigenvalues, λ)
		
		residuals = AV * Y - V * Y * Diagonal(λ)
		norm_r = norm.(eachcol(residuals))
		push!(residual_norms, norm_r)

		verbose && @printf "%3i %8.4g %8.4g\n" i λ[end] norm_r[end]
		maximum(norm_r) < tol && break
		
		V = AV
	end
	X = V * Y
	
	(; λ, X, eigenvalues, residual_norms)
end

# ╔═╡ 2e96d48d-15d9-40cc-8c1d-d8d03e755c2a
Asub = Diagonal(1.0:2.0:100.0)

# ╔═╡ e6a62b4f-9919-43b9-ae74-3cb3c1110829
let
	(; residual_norms, λ) = projected_subspace_iteration(Asub)
	@show λ

	norms_1st = [r[1] for r in residual_norms]
	norms_2nd = [r[2] for r in residual_norms]

	plot( norms_1st; yaxis=:log, label="residual λ₄₉", lw=2)
	plot!(norms_2nd; yaxis=:log, label="residual λ₅₀", lw=2)
end

# ╔═╡ bef6d532-6131-4db7-8b7c-582895a035b1
md"""
An alternative notation for subspace methods, that is sometimes handy, is given in terms of projections.
Making use of our basis for $\mathcal{S}$ it is easy to verify that
```math
P_\mathcal{S} = V V^H \in \mathbb{C}^{N\times N}
```
is indeed a projection onto the subspace $\mathcal{S}$, since the idempotency relation
```math
P_\mathcal{S}^2 = \left( V V^H \right)^2 = V (V^H V) V^H = V I V^H = P_\mathcal{S}
```
is satisfied. Noting further that $\tilde{x}_i = V \tilde{y}_i$ implies $V^H \tilde{x}_i = \tilde{y}_i$ and thus $P_\mathcal{S} \tilde{x}_i = \tilde{x}_i$,
we can rewrite $(\S)$ as
```math
P_\mathcal{S} A \tilde{x}_i = V \left(V^H A V\right) \tilde{y}_i
\stackrel{(\S)}{=} V \tilde{λ}_i \tilde{y}_i
= P_\mathcal{S} \tilde{λ}_i \tilde{x}_i
= \tilde{λ}_i \tilde{x}_i,
```
which in turn implies
```math
P_\mathcal{S} \left( A \tilde{x}_i - \tilde{λ}_i \tilde{x}_i \right) = 0
```
--- a projection version of the Galerkin conditions (G) above.
Yet alternatively we can state
```math
P_\mathcal{S} A P_\mathcal{S} \, \tilde{x}_i
= \tilde{λ}_i P_\mathcal{S} \tilde{x}_i
= \tilde{λ}_i \tilde{x}_i,
```
which is an eigenproblem in the projected operator $P_\mathcal{S} A P_\mathcal{S}$. 
"""

# ╔═╡ 3cea29a7-fd87-4491-81a6-c53f79387e1f
md"""
The appeal of projection methods for computing the eigenpairs of Hermitian matrices results from the connection to the Courant-Fisher min-max principle, which guarantees that the eigenvalues computed by the Rayleigh-Ritz procedure satisfy a strong optimality condition.

First we note that
```math
\langle \tilde{x}, A \tilde{x} \rangle
= \big\langle \tilde{x}, \left(P_\mathcal{S} A P_\mathcal{S} \right)\, \tilde{x} \big\rangle
\qquad \forall \tilde{x} \in \mathcal{S},
```
from which it follows that
if $\tilde{x}$ is inside the subspace employed by an orthogonal projection approach (like the Rayleigh-Ritz procedure), the Rayleigh quotient of the projected matrix $P_\mathcal{S} A P_\mathcal{S}$ and of the original matrix $A$ are identical:
```math
\tilde{λ}_1 = \min_{\tilde{x}\in \mathcal{S}, \tilde{x}\neq 0} \frac{\big\langle \tilde{x}, \left(P_\mathcal{S} A P_\mathcal{S} \right) \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
= \min_{\tilde{x}\in \mathcal{S}, x\neq 0}
\frac{\big\langle P_\mathcal{S} \tilde{x}, A P_\mathcal{S} \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
= \min_{\tilde{x}\in \mathcal{S}, \tilde{x}\neq 0}
\frac{\big\langle \tilde{x}, A \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}.
```
More generally, based on applying Courant-Fisher to the matrix $P_\mathcal{S} A P_\mathcal{S}$, we obtain the following result

**Proposition 5.5.** The $i$-th approximate eigenvalue (counted from small to large including multiplicities) of a Hermitian matrix $A$ obtained from an orthogonal projection method to a subspace $\mathcal{S}$ satisfies
```math
\tilde{λ}_i = \min_{s \subset \mathcal{S},\,\text{dim}(s) = i} \ \max_{\tilde{x}\in s,\,\tilde{x}\neq 0} \frac{\big\langle \tilde{x}, A \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
```

From this we obtain immediately the

**Corollary 5.6.** For $i = 1, 2, \ldots, m$ the inequality
```math
λ_i \leq \tilde{λ}_i
```
holds, i.e. eigenvalues are **approximated from above** in orthogonal projection methods.

**Proof.**
```math
\tilde{λ}_i = \min_{s \subset \textcolor{red}{\mathcal{S}},\,\text{dim}(s) = i} \ \max_{\tilde{x}\in s,\,\tilde{x}\neq 0} \frac{\big\langle \tilde{x}, A \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
\geq
\min_{s \subset \textcolor{red}{\mathbb{C}^N},\,\text{dim}(s) = i} \ \max_{x\in s,\,x\neq 0} \frac{\big\langle x, A x \big\rangle}{\langle x, x \rangle}
= λ_i.
```
This is the desired optimality: The eigenvalue approximation found by a projection method is the best eigenvalue we can possibly find within the chosen subspace, and it provides a rigorous upper bound to the exact one.

Finally, the following characterisation provides a connection between the approximate eigenvalue $\tilde{λ}_i$ and the corresponding approximate eigenvector $\tilde{x}_i$ obtained by an orthogonal projection approach:
```math
\tilde{λ}_i
= \frac{\big\langle \tilde{x}_i, A \tilde{x}_i \big\rangle}{\langle \tilde{x}_i, \tilde{x}_i \rangle}
= \min_{0\neq \tilde{x} \in \mathcal{S},\,\tilde{x} \perp \tilde{\mathcal{X}}_{k-1}} \frac{\big\langle \tilde{x}, A \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
```
where $\tilde{\mathcal{X}}_k = \text{span}(\tilde{x}_1, \ldots, \tilde{x}_k)$,
the space spanned by the $k$ first *approximate* eigenvectors.
"""

# ╔═╡ e1fd4f96-a465-46d8-8c78-cde7c5325e6f
md"""
### Forming a good subspace

Projection methods and the Rayleigh-Ritz procedure are a key ingredient of pretty much all iterative diagonalisation approaches employed nowadays. A full discussion of the standard techniques employed to build the reduced subspace $\mathcal{S}$ is out of the scope of this lecture. An incomplete list of techniques worth mentioning are:
- If a Krylov basis is constructed and employed for $\mathcal{S}$ one obtains diagonalisation methods such as *Lanczos* or *Arnoldi*.
- Employing an orthogonal correction based on the current residual to expand and enrich the subspace leads to the *Davidson* family of methods.
- If the subspace is constructed following the idea of minimising the Rayleigh quotient, we obtain minimisation methods (inverse iterations, PCG, LOPCG). These will be discussed in the following.
"""

# ╔═╡ c56ac72f-643a-4b29-97cb-4f61426d11d5
md"""
## Minimisation methods

In our discussion so far we already employed the variational approach encoded in Courant-Fisher for obtaining eigenvalue problems. E.g. in order to obtain the smallest eigenvalue $λ_1$ we may just
```math
λ_1 = \min_{0 \neq x \in \mathbb{C}^N} R_A(x) \qquad R_A(x) = \frac{x^H A x}{x^H x}
```
The purpose of this subsection is to develop methods, that employ this minimisation property directly for computing the eigenpairs of $A$.

The natural reflex is to employ a steepest descent approach.
We thus compute the gradient of $R_A(x)$, namely
(denoting by $e_i$ the unit vector along coordinate direction $i$):
```math
\begin{aligned}
\frac{\partial R_A(x)}{\partial x_i}
&= \left. \frac{d R_A(x + ε e_i)}{d ε} \right|_{ε=0} \\
&= \frac{d}{d ε} \left. \left(
	\frac{(x + εe_i)^H A (x + εe_i)}{(x + εe_i)^H (x + εe_i)}
\right)\right|_{ε=0}\\
&= \left(
	\frac{\left(2 e_i^H A (x + εe_i) \right) \left( (x + εe_i)^H (x + εe_i) \right) }
	{ \left( (x + εe_i)^H (x + εe_i) \right)^2 } \right. \\
	&\quad\left.\left. -\frac{ \left( 2 e_i^H (x + εe_i)  \right) \left( (x + εe_i)^H A (x + εe_i) \right)}
	{ \left( (x + εe_i)^H (x + εe_i) \right)^2 }
\right) \right|_{ε=0}\\
&= 2\frac{e_i^H A x}{x^Hx} - 2 \frac{R_A(x)}{x^H x} e_i^H x\\
&= \frac{2}{x^H x} \big( (A x)_i - R_A(x) \, x_i \big)
\end{aligned}
```
such that
```math
\nabla R_A(x) = 2 \frac{Ax - R_A(x) x}{\| x \|^2}.
```
We note that the residual $r = Ax - R_A(x) x$ points along the direction of the gradient of $R_A(x)$.

A classic algorithm to solve such minimisation problems in the objective $R_A(x)$ is **preconditioned gradient descent**, in which the iterates are related as
```math
x^{(i+1)} \leftarrow x^{(i)} - \tilde{α} P^{-1} \nabla R_A(x^{(i)})
```
i.e. in which simply a scaled fraction of the gradient is removed from the current iterate. In our case this becomes
```math
x^{(i+1)} \leftarrow x^{(i)} - α P^{-1} r^{(i)} = x^{(i)} - α P^{-1} \left( A x^{(i)} - R_A(x^{(i)}) \, x^{(i)} \right),
```
where $r^{(i)} = Ax^{(i)} - R_A(x^{(i)}) x^{(i)}$.

In practice the step-size $α$ and the preconditioner matrix $P^{-1}$ need to be chosen cleverly to obtain fast convergence. Algorithmically we can write:
"""

# ╔═╡ 1d099205-9d84-4ecc-a2af-81910326a4ba
function preconditioned_gradient_descent(A; x=randn(eltype(A), size(A, 2)),
                                         α=1.0, Pinv=I, tol=1e-6, maxiter=100,
                                         verbose=true)
	T = real(eltype(A))

	eigenvalues    = T[]
	residual_norms = T[]
	λ = NaN
	for i in 1:maxiter
		normalize!(x)
		Ax = A * x
		λ = x'Ax         # Rayleigh quotient
		r = Ax - λ * x   # Residual
		norm_r = norm(r)

		verbose && @printf "%3i %8.4g %8.4g\n" i λ norm_r
		push!(eigenvalues, λ)
		push!(residual_norms, norm_r)
		norm_r < tol && break

		# Precondition the residual and move along the direction:
		x = x - α * Pinv * r
	end
	converged = residual_norms[end] < tol
	!converged && verbose && @warn "PGD not converged."
	
	(; λ, x, converged, eigenvalues, residual_norms)
end

# ╔═╡ f17500cc-6e5f-442d-9b22-2dc406367ad6
Apgd   = Diagonal(abs.(randn(100)))

# ╔═╡ 320a4138-1d2d-4039-adaa-1aeb7f1929cf
let
	Pinv = Diagonal(1 ./ diag(Apgd))  # Diagonal preconditioner for Apgd
	(; residual_norms) = preconditioned_gradient_descent(Apgd; Pinv)
	plot(residual_norms; yaxis=:log, label="prec. grad. desc. residual λ₁", lw=2, mark=:x)
end

# ╔═╡ db913617-2b81-4ea7-a907-033dd31bda5a
md"Convergence is thus great if a good preconditioner is chosen. However, for no preconditioner we obtain:"

# ╔═╡ 68fff309-a3ec-4ebc-bcfe-bafeac6ef3f5
preconditioned_gradient_descent(Apgd)

# ╔═╡ f4bb1d69-1f06-4769-ac6e-081ddaa437d7
md"""
Preconditioned gradient descent is closely related to the inverse power method. Namely, if we chose $P^{-1} = A^{-1}$ and $α = 1$, then we obtain
```math
x^{(i+1)} \leftarrow x^{(i)} - A^{-1} r^{(i)} = x^{(i)} - A^{-1} \left( A x^{(i)} - R_A(x^{(i)}) \, x^{(i)} \right)
= R_A(x^{(i)}) A^{-1} x^{(i)},
```
which is a scaled version of the inverse power method.
If one considered preconditioned gradient ascent one could
similarly draw a connection to the power method.
"""

# ╔═╡ 8481465b-171b-4374-8ec9-d7c19bd23d81
md"""
### Choice of preconditioner

What preconditioner should be chosen is often highly problem dependent. However, there are a number of standard choices, which can often be employed blindly to speed up the calculations. Popular are in particular incomplete factorisations (e.g. Cholesky) or lower-precision version of the original problem. A detailed discussion is out of the scope of these notes. Details can be found, for example, in book of Saad.

"""

# ╔═╡ 179be577-6d0e-4e97-8db1-56e116c502c6
md"""
### Step sizes and LOBPCG

The final aspect we want to briefly address is how to choose the step size $α$. The best step size can be found by a line search,
i.e. by solving the one-dimensional optimisation problem
```math
α_{i,\text{ls}} = \text{argmin}_{α\in \mathbb{R}} R_A(x^{(i)} + α \tilde{r}^{(i)})
```
where we denoted $\tilde{r}^{(i)} = P^{-1} r^{(i)}$.
The next iterate is then formed
as $x^{(i+1)} = x^{(i)} + α_{i,\text{ls}} \tilde{r}^{(i)}$.
Since the Rayleigh quotient contains normalisation in the denominator
we note that
```math
R_A(δ x^{(i)} + ε \tilde{r}^{(i)} ) = R_A(  x^{(i)} + \frac{ε}{δ} \tilde{r}^{(i)} ).
```
Therefore, our full procedure could equivalently be written as
```math
x^{(i+1)} = \text{argmin}_{y \in \text{span} \left\{ x^{(i)}, \tilde{r}^{(i)}  \right\}} R_A(y),
```
i.e. as a minimisation in the subspace spanned by $x^{(i)}$ and $\tilde{r}^{(i)}$.

A slight generalisation of this idea is to make the minimisation subspace even larger by also including the *previous* iterate to the search space, i.e. to employ
```math
x^{(i+1)} = \text{argmin}_{y \in \text{span} \left\{ x^{(i)}, x^{(i-1)}, \tilde{r}^{(i)}  \right\}} R_A(y). \qquad \text{$(\ast)$}
```
This turns out to improve the overall observed convergence rates considerably. This modification has been suggested following a more detailed investigation of minimisation approaches for eigenvalue problems, in particular a comparison with more sophisticated minimisation algorithms such as the conjugate gradient (CG) approach, which is out of the scope of this lecture.

The improved convergence rates of $(\ast)$ unfortunately comes with a downside: As the iterations converge $x^{(i)}$ and $x^{(i-1)}$ *both* converge to the exact eigenvector, such that $x^{(i)}$ and $x^{(i-1)}$ become more and more linearly dependent, which makes the optimisation unstable. This is in practice circumvented by two further modifications:
- Instead of using the basis vectors $\{x^{(i)}, x^{(i-1)}, \tilde{r}^{(i)}\}$ to build the subspace, we compute $p^{(i)} = x^{(i)} - x^{(i-1)}$
  and instead employ $\{ x^{(i)}, p^{(i)}, \tilde{r}^{(i)} \}$,
   which span the same space.
- Instead of first performing the line search to find the optimal $x^{(i+1)}$ and then computing the Rayleigh coefficient $R_A(x^{(i+1)})$, we perform a Rayleigh-Ritz procedure projecting onto the subspace $\mathcal{S} = \text{span} \{ x^{(i)}, p^{(i)}, \tilde{r}^{(i)} \}$. Due to our results on subspace methods (in particular Proposition 5.5.) we know that a Rayleigh-Ritz in the subspace $\mathcal{S}$ is equivalent to finding the minimiser of $R_A(x)$ where $x \in \mathcal{S}$.

These two modifications with our previous discussion lead to the LOPCG (locally optimal preconditioned conjugate gradient) algorithm:
"""

# ╔═╡ 94e6808b-c89c-4431-90d8-5e252e38d834
function lopcg(A; x=randn(eltype(A), size(A, 2)), ortho=ortho_qr,
                  Pinv=I, tol=1e-6, maxiter=100, verbose=true)
	T = real(eltype(A))

	eigenvalues    = T[]
	residual_norms = T[]
	λ = NaN
	p = nothing
	r = nothing
	
	for i in 1:maxiter
		if i > 1
			Z = hcat(x, p, r)
		else
			Z = hcat(x)
		end
		Z = ortho(Z)
		
		# Rayleigh-Ritz step to get smallest eigenvalue
		AZ = A * Z
		λ, Y = eigen(Hermitian(Z' * AZ))
		λ = λ[1]     # Keep only smallest eigenvalue
		y = Y[:, 1]  # Keep only smallest eigenvector
		new_x = Z * y  # x^{(i)} is the eigenvector of the Rayleigh-Ritz,
		               # since this gives the optimal linear combination between
		               # p^{(i-1)}, x^{(i-1)} and ̃r^{(i-1)}

		# Store results and residual
		push!(eigenvalues, λ)
		r = AZ * y - new_x * λ
		norm_r = norm(r)
		push!(residual_norms, norm_r)
		verbose && @printf "%3i %8.4g %8.4g\n" i λ norm_r
		if norm_r < tol
			x = new_x
			break
		end

		# Precondition residual, update x and p
		r = Pinv * r
		p = x - new_x
		x = new_x
	end
	converged = residual_norms[end] < tol
	!converged && verbose && @warn "LOPCG not converged."
	
	(; λ, x, converged, eigenvalues, residual_norms)
end

# ╔═╡ f1398a23-90d9-40ea-8ef2-ebcd0ab67b14
Alopcg = Diagonal(abs.(randn(1000)).^0.1);

# ╔═╡ bba31f5b-7eca-412b-99c7-65adcfa2be1b
let
	Pinv = Diagonal(1 ./ diag(Alopcg))  # Diagonal preconditioner for Apgd
	(; residual_norms) = lopcg(Alopcg; Pinv)

	plot(residual_norms; yaxis=:log, lw=2, label="LOPCG residual norms", mark=:x)
end

# ╔═╡ 9b7bdce7-38cb-4c7c-b1b4-d1ec3c9fa168
md"""
The generalisation to an algorithm capable of computing multiple eigenpairs is straightforward, leading to LOBPCG (locally optimal *block* preconditioned conjugate gradient) is only bookkeeping:
"""

# ╔═╡ 42777c48-0f3a-4c43-974b-7320ad16efb9
function lobpcg(A; X=randn(eltype(A), size(A, 2), 2), ortho=ortho_qr,
                   Pinv=I, tol=1e-6, maxiter=100, verbose=true)
	T = real(eltype(A))
	m = size(X, 2)  # block size

	eigenvalues    = Vector{T}[]
	residual_norms = Vector{T}[]
	λ = NaN
	P = nothing
	R = nothing
	
	for i in 1:maxiter
		if i > 1
			Z = hcat(X, P, R)
		else
			Z = X
		end
		Z = ortho(Z)
		
		# Rayleigh-Ritz step to get smallest eigenvalues
		AZ = A * Z
		λ, Y = eigen(Hermitian(Z' * AZ))
		λ = λ[1:m]
		Y = Y[:, 1:m]
		new_X = Z * Y

		# Store results and residual
		push!(eigenvalues, λ)
		R = AZ * Y - new_X * Diagonal(λ)
		norm_r = norm.(eachcol(R))
		push!(residual_norms, norm_r)
		verbose && @printf "%3i %8.4g %8.4g\n" i λ[end] norm_r[end]
		if maximum(norm_r) < tol
			X = new_X
			break
		end
		
		# Precondition residual, update X and P
		R = Pinv * R
		P = X - new_X
		X = new_X
	end

	(; λ, X, eigenvalues, residual_norms)
end

# ╔═╡ 05175611-34de-44cf-a385-c493a6f3cf19
let
	Pinv = Diagonal(1 ./ diag(Alopcg))  # Diagonal preconditioner for Apgd
	(; residual_norms) = lobpcg(Alopcg; Pinv)

	norms_first  = [n[1] for n in residual_norms]
	norms_second = [n[2] for n in residual_norms]
	
	p = plot(; yaxis=:log)
	plot!(p, norms_first,  label="residual first eigenpair",  c=1, lw=2, mark=:x)
	plot!(p, norms_second, label="residual second eigenpair", c=2, lw=2, mark=:x)
end

# ╔═╡ 15220d0d-2c9f-4daf-a290-672ca439faec
md"""
Note, that this simple implementation of LOBPCG is still not numerically stable and considerable issues remain. In practice one should thus always employ one of the well-tested implementation in standard numerical libraries.

For example the LOBPCG of the DFTK package:
"""

# ╔═╡ 23990899-05d2-4b40-9d52-23caf218c7cf
let
	# Note that DFTK wants the operator P and P⁻¹
	prec = Diagonal(diag(Alopcg))

	X = randn(size(Alopcg, 2), 2)
	(; residual_history) = DFTK.lobpcg_hyper(Alopcg, X;
	                                         prec, tol=1e-6, display_progress=true)

	residual_history .+= eps(Float64)  # to ensure we can use logarithmic axes 
	
	p = plot(; yaxis=:log, ylims=(1e-7, Inf))
	plot!(p, residual_history[1, :],  label="residual first eigenpair",
		  c=1, lw=2, mark=:x)
	plot!(p, residual_history[2, :], label="residual second eigenpair",
		  c=2, lw=2, mark=:x)
end

# ╔═╡ 461162d6-b0e6-41c4-8e60-6293d2fc25ee
md"""
## How do the methods compare ?

To summarise our survey, we want to briefly compare the methods on an example.
We set up a test problem, where we can tune both the size of the matrix and
the gaps between the first eigenvalues. We then try multiple methods to obtain the **smallest** one or two eigenpairs:
"""

# ╔═╡ 92a6feb4-88ad-4681-b6dd-2df2a790360c
md"""
### Single-vector approaches
"""

# ╔═╡ aa9f8c4b-b733-48c6-93f0-4ab8d52b2a3e
md"""
- `logN_s`: Logarithm of the **problem size** $(@bind logN_s PlutoUI.Slider(1.2:0.1:4.0, default=2.0))
- `log1_s`: Log of gap between **1st and 2nd** eigenvalue: $(@bind log1_s PlutoUI.Slider(-6:0.1:3.0; show_value=true, default=1.0))
- `logPrec_s`: **Preconditioner noise** level $(@bind logPrec_s PlutoUI.Slider(-3:0.1:-1.5, default=-2.5, show_value=true))

Plotting options:
- `maxiter_s`: Plotting range is from `0` to `maxiter`: $(@bind maxiter_s PlutoUI.Slider(10:1:100; show_value=true, default=20))
"""

# ╔═╡ 99a01adc-34e2-4284-b918-4bfc020b47d6
begin
	N_s = ceil(Int, 10^logN_s)
	A_s = Diagonal([1.0
	                1.0 + 10^log1_s])
	B_s = Diagonal(A_s[2, 2] .+ 10 .+ randn(N_s-2))
	O_s = zeros(2, N_s-2)
	
	M_s = [A_s  O_s;
	       O_s' B_s]
end

# ╔═╡ 1cf09a9d-2e2b-4a0c-848c-56d3bd1e7f87
x = randn(size(M_s, 2));  # Consistent starting vector for all methods

# ╔═╡ 39215879-264c-4c2a-bafc-a792abf70017
Pnoise_s = 10^logPrec_s * randn(size(M_s, 1));

# ╔═╡ ff50c8b6-ce84-4d74-b5ba-f5abecdd7fb2
let
	function inverse_power_method(A; kwargs...)
		A⁻¹ = InverseMap(factorize(A))
		power_method(A⁻¹; kwargs...)
	end

	p = plot(yaxis=:log, ylims=(1e-7, 10))
	
	non_preconditioned_methods = (
		inverse_power_method, rayleigh_quotient_iteration,
	)
	for method in non_preconditioned_methods
		(; residual_norms) = method(M_s; x, verbose=false, tol=1e-6)
		plot!(p, residual_norms, label=string(method), lw=2, mark=:x)
	end

	preconditioned_method = (
		preconditioned_gradient_descent, lopcg
	)
	for (i, method) in enumerate(preconditioned_method)
		c = i + length(non_preconditioned_methods)
		
		# Perfect preconditioner: The inverse diagonal
		Pinv = Diagonal(1 ./ diag(M_s))
		(; residual_norms) = method(M_s; x, verbose=false, tol=1e-6, Pinv)
		plot!(p, residual_norms; label=string(method) * " (perfect precon)",
			  lw=2, c, mark=:x)

		# Preconditioner plus noise
		Pinv = Diagonal(1 ./ diag(M_s) .+ Pnoise_s)
		(; residual_norms) = method(M_s; x, verbose=false, tol=1e-6, Pinv)
		plot!(p, residual_norms; label=string(method) * " (noisy precon)",
			  lw=2, c, ls=:dash, mark=:x)
	end

	default_lim = xlims(p)
	xlims!(p, 0, min(default_lim[2], maxiter_s))
end

# ╔═╡ f9e4e29c-247d-4896-90d0-9c5f716c203c
md"""By playing with the sliders you should notice:

- Rayleigh quotient iteration generally performs best, but is also the most expensive (Full factorisation required at *every iteration*)
- Inverse power method and preconditioned gradient descent fail when the gap becomes too small.
- Noise in the preconditioner leads to worse convergence in LOPCG and preconditioned gradient descent, but LOPCG is more resistant
- Overall LOPCG offers a pretty good compromise.
"""

# ╔═╡ 2f4000e9-eaca-446d-af43-3b8ea06eb2b1
md"""
### Block and subspace approaches
"""

# ╔═╡ fb887ea1-b4bb-49d5-af84-9a1278791d67
md"""
- `logN_b`: Logarithm of the **problem size** $(@bind logN_b PlutoUI.Slider(1.2:0.1:4.0, default=2.0))
- `log1_b`: Log of gap between **1st and 2nd** eigenvalue: $(@bind log1_b PlutoUI.Slider(-6:0.1:3.0; show_value=true, default=1.0))
- `log2_b`: Log of gap between **2nd and 3rd** eigenvalue: $(@bind log2_b PlutoUI.Slider(-6:0.1:3.0; show_value=true, default=1.0))
- `logPrec_b`: **Preconditioner noise** level $(@bind logPrec_b PlutoUI.Slider(-3:0.1:-1.5, default=-2.5, show_value=true))

Plotting options:
- `maxiter_b`: Plotting range is from `0` to `maxiter`: $(@bind maxiter_b PlutoUI.Slider(10:1:100; show_value=true, default=20))
"""

# ╔═╡ 9d86090f-9680-43fa-b86b-2fc664209faf
begin
	N_b = ceil(Int, 10^logN_b)
	A_b = Diagonal([1.0
	                1.0 + 10^log1_b
	                1.0 + 10^log1_b + 10^log2_b])
	B_b = Diagonal(A_b[3, 3] .+ 10 .+ randn(N_b-3))
	O_b = zeros(3, N_b-3)
	
	M_b = [A_b  O_b;
	       O_b' B_b]
end

# ╔═╡ eb909306-b913-464b-81d6-becc9946e7c1
X = randn(size(M_b, 2), 2);  # Consistent starting vectors for all methods

# ╔═╡ aa4f358a-e8fa-4258-bcb6-b6707ef62b0a
Pnoise_b = 10^logPrec_b * randn(size(M_b, 1));

# ╔═╡ bdcc33c7-eda7-48fe-a73c-84c2e7c64dab
let
	function lobpcg_dftk(A; X, Pinv::Diagonal, tol, verbose)
		prec = Diagonal(1 ./ diag(Pinv))
		(; λ, residual_history) = DFTK.lobpcg_hyper(A, X; prec, tol,
		                                            display_progress=verbose)
		residual_norms = [Array(c) for c in eachcol(residual_history)]
		(; λ, residual_norms)
	end

	function plot_both!(p, λ, residual_norms; label, kwargs...)
		# Extract residual, but make sure that eigenpairs are sorted ascendingly
		sorting = sortperm(λ)
		resid_1st = [resid[sorting[1]] for resid in residual_norms]
		resid_2nd = [resid[sorting[2]] for resid in residual_norms]  

		plot!(p, resid_1st; label=label * " λ₁", m=:x, lw=2, markersize=4, kwargs...)
		plot!(p, resid_2nd; label=label * " λ₂", m=:+, lw=2, markersize=4, kwargs...)
	end
		
	p = plot(yaxis=:log, ylims=(1e-7, 10))
	c = 1

	for method in (subspace_iteration, projected_subspace_iteration)
		M⁻¹ = Diagonal(1 ./ diag(M_b))
		(; λ, residual_norms) = method(M⁻¹; X, verbose=false, tol=1e-6)
		plot_both!(p, 1 ./ λ, residual_norms; label=string(method), c)
		c += 1
	end

	for method in (lobpcg, lobpcg_dftk)
		# Perfect preconditioner: The inverse diagonal
		Pinv = Diagonal(1 ./ diag(M_b))
		(; λ, residual_norms) = method(M_b; X, verbose=false, tol=1e-6, Pinv)
		plot_both!(p, λ, residual_norms; label=string(method) * "(perfect precon)", c)

		# Preconditioner plus noise
		Pinv = Diagonal(1 ./ diag(M_b) .+ Pnoise_b)
		(; λ, residual_norms) = method(M_b; X, verbose=false, tol=1e-6, Pinv)
		plot_both!(p, λ, residual_norms; label=string(method) * "(noisy precon)",
		           c, ls=:dash)

		c += 1
	end

	default_lim = xlims(p)
	xlims!(p, 0, min(default_lim[2], maxiter_b))
end

# ╔═╡ 9759f044-eeba-4263-84e6-0113b7310c0b
md"""
By playing with the sliders you should notice:
- TODO
"""

# ╔═╡ b2cbf4dd-f203-4028-90b1-57d72d21c3c5
begin

toc = 
# sidebar --- DO NOT TOUCH THIS LINE
Markdown.parse( "**Error control in scientific modeling** 
" * read("sidebar.md",String)) 
# sidebar --- DO NOT TOUCH THIS LINE 



Sidebar(elts...; location="upper right") = @htl("""
<aside class="sidebar" style='top: 530px;right: 17px;'>$elts</aside>

<style>
aside.sidebar {
	position: fixed;
	max-width: min(30%, 300px, calc(100vw - 750px));
	padding: 0.4rem;
	border-radius: 10px;
	max-height: calc(100vh - 590px);
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
DFTK = "acf6eb54-70d9-11e9-0013-234b7a5f5337"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearMaps = "7a12625a-238d-50fd-b39a-03d52299707e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
DFTK = "~0.6.11"
HypertextLiteral = "~0.9.5"
LinearMaps = "~3.11.0"
Plots = "~1.39.0"
PlutoTeachingTools = "~0.2.13"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "76166e354d3230512d3d997aa8438523b00d4b08"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AtomsBase]]
deps = ["LinearAlgebra", "PeriodicTable", "Printf", "Requires", "StaticArrays", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "995c2b6b17840cd87b722ce9c6cdd72f47bab545"
uuid = "a963bdd2-2df7-4f54-a1ee-49d51e6be12a"
version = "0.3.5"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bravais]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "08c0613c61dcdd5dd148982cf91b107654e2e811"
uuid = "ada6cbde-b013-4edf-aa94-f6abe8bd6e6b"
version = "0.1.8"

[[deps.Brillouin]]
deps = ["Bravais", "DirectQhull", "DocStringExtensions", "LinearAlgebra", "Reexport", "Requires", "StaticArrays"]
git-tree-sha1 = "d349be305aa8b86593f611409211931a508aeb1f"
uuid = "23470ee3-d0df-4052-8b1a-8cbd6363e7f0"
version = "0.5.11"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CUDA_Driver_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "35a37bb72b35964f2895c12c687ae263b4ac170c"
uuid = "4ee394cb-3365-5eb0-8335-949819d2adfc"
version = "0.6.0+3"

[[deps.CUDA_Runtime_jll]]
deps = ["Artifacts", "CUDA_Driver_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "bfe5a693a11522d58392f742243f2b50dc27afd6"
uuid = "76a88914-d11a-5bdc-97e0-2f5a05c973a2"
version = "0.9.2+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "a1296f0fe01a4c3f9bf0dc2934efbf4416f5db31"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

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
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.ComponentArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "Functors", "LinearAlgebra", "Requires", "StaticArrayInterface"]
git-tree-sha1 = "00380a5de40736c634b867069347b721ca311673"
uuid = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
version = "0.13.15"

    [deps.ComponentArrays.extensions]
    ComponentArraysConstructionBaseExt = "ConstructionBase"
    ComponentArraysGPUArraysExt = "GPUArrays"
    ComponentArraysRecursiveArrayToolsExt = "RecursiveArrayTools"
    ComponentArraysReverseDiffExt = "ReverseDiff"
    ComponentArraysSciMLBaseExt = "SciMLBase"
    ComponentArraysStaticArraysExt = "StaticArrays"

    [deps.ComponentArrays.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SciMLBase = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DFTK]]
deps = ["AbstractFFTs", "Artifacts", "AtomsBase", "Brillouin", "ChainRulesCore", "Dates", "DftFunctionals", "FFTW", "ForwardDiff", "GPUArraysCore", "InteratomicPotentials", "Interpolations", "IterTools", "IterativeSolvers", "LazyArtifacts", "Libxc", "LineSearches", "LinearAlgebra", "LinearMaps", "MPI", "Markdown", "Optim", "OrderedCollections", "PeriodicTable", "PkgVersion", "Polynomials", "PrecompileTools", "Preferences", "Primes", "Printf", "ProgressMeter", "PseudoPotentialIO", "Random", "Requires", "Roots", "SparseArrays", "SpecialFunctions", "Spglib", "StaticArrays", "Statistics", "TimerOutputs", "Unitful", "UnitfulAtomic", "spglib_jll"]
git-tree-sha1 = "135f4f4d8fd9f29ec1c9d67dd009c98eef53621e"
uuid = "acf6eb54-70d9-11e9-0013-234b7a5f5337"
version = "0.6.11"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DftFunctionals]]
deps = ["ComponentArrays", "DiffResults", "ForwardDiff"]
git-tree-sha1 = "171252445375f15afb382fe58ce29a54b3216bf1"
uuid = "6bd331d2-b28d-4fd3-880e-1a1c7f37947f"
version = "0.2.3"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DirectQhull]]
deps = ["Qhull_jll"]
git-tree-sha1 = "5a941ad556ad4d2e310828b0f0b462678887ec2e"
uuid = "c3f9d41a-afcb-471e-bc58-0b8d83bd86f4"
version = "0.2.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "b6def76ffad15143924a2199f72a5cd883a2e8a9"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.9"
weakdeps = ["SparseArrays"]

    [deps.Distances.extensions]
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

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
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "0fa3b52a04a4e210aeb1626def9c90df3ae65268"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.1.0"

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

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "a20eaa3ad64254c61eeb5f230d9306e937405434"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.6.1"
weakdeps = ["SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

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

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9a68d75d466ccc1218d0552a8e1631151c569545"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

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
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad37c091f7d7daf900963171600d7c1c5c3ede32"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InteratomicPotentials]]
deps = ["AtomsBase", "Distances", "LinearAlgebra", "NearestNeighbors", "StaticArrays", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "e52c1cff4fa468972621f0b5dd45ce2ee08dc730"
uuid = "a9efe35a-c65d-452d-b8a8-82646cd5cb04"
version = "0.2.6"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

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
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "81dc6aefcbe7421bd62cb6ca0e700779330acff8"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.25"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.Libxc]]
deps = ["Libxc_GPU_jll", "Libxc_jll", "Requires"]
git-tree-sha1 = "4c3c4e4c0917b3bcfa9bde78b7f1be011a827bac"
uuid = "66e17ffc-8502-11e9-23b5-c9248d0eb96d"
version = "0.3.15"

    [deps.Libxc.extensions]
    LibxcCudaExt = "CUDA"

    [deps.Libxc.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.Libxc_GPU_jll]]
deps = ["Artifacts", "CUDA_Runtime_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg", "TOML"]
git-tree-sha1 = "ee321f68686361802f2ddb978dae441a024e61ea"
uuid = "25af9330-9b41-55d4-a324-1a83c0a0a1ac"
version = "6.1.0+2"

[[deps.Libxc_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c5516f2b1655a103225e69477e3df009347580df"
uuid = "a56a6d9d-ad03-58af-ab61-878bf78270d6"
version = "6.1.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearMaps]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6698ab5e662b47ffc63a82b2f43c1cee015cf80d"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.11.0"
weakdeps = ["ChainRulesCore", "SparseArrays", "Statistics"]

    [deps.LinearMaps.extensions]
    LinearMapsChainRulesCoreExt = "ChainRulesCore"
    LinearMapsSparseArraysExt = "SparseArrays"
    LinearMapsStatisticsExt = "Statistics"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

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
git-tree-sha1 = "0d097476b6c381ab7906460ef1ef1638fbce1d91"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.2"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "eb006abbd7041c28e0d16260e50a24f8f9104913"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.2.0+0"

[[deps.MPI]]
deps = ["Distributed", "DocStringExtensions", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "PkgVersion", "PrecompileTools", "Requires", "Serialization", "Sockets"]
git-tree-sha1 = "df53d0e1e0dbebf2315f4cd35e13e52ad43416c2"
uuid = "da04e1cc-30fd-572f-bb4f-1f8673147195"
version = "0.20.15"

    [deps.MPI.extensions]
    AMDGPUExt = "AMDGPU"
    CUDAExt = "CUDA"

    [deps.MPI.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "8a5b4d2220377d1ece13f49438d71ad20cf1ba83"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.1.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "781916a2ebf2841467cda03b6f1af43e23839d85"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.9"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "6979eccb6a9edbbb62681e158443e79ecc0d056a"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.3.1+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a8027af3d1743b3bfae34e54872359fdebb31422"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.3+4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "f3080f4212a8ba2ceb10a34b938601b862094314"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "4.1.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "963b004d15216f8129f6c0f7d187efa136570be0"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.7"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.PeriodicTable]]
deps = ["Base64", "Test", "Unitful"]
git-tree-sha1 = "9a9731f346797126271405971dfdf4709947718b"
uuid = "7b2266bf-644c-5ea3-82d8-af4bbd25a884"
version = "1.1.4"

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

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

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
git-tree-sha1 = "542de5acb35585afcf202a6d3361b430bc1c3fbd"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.13"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "3aa2bb4982e575acd7583f01531f241af077b163"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.13"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

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

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "4c9f306e5d6603ae203c2000dd460d81a5251489"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.PseudoPotentialIO]]
deps = ["EzXML", "LinearAlgebra"]
git-tree-sha1 = "88cf9598d70015889c99920ff3dacca0eb26ae90"
uuid = "cb339c56-07fa-4cb2-923a-142469552264"
version = "0.1.1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be2449911f4d6cfddacdf7efc895eceda3eee5c1"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1003+0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "7c29f0e8c575428bd84dc3c72ece5178caa67336"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.2+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

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
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "7364d5f608f3492a4352ab1d40b3916955dc6aec"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.5"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "ff42754a57bb0d6dcfe302fd0d4272853190421f"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.19"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Spglib]]
deps = ["StaticArrays", "StructEquality", "spglib_jll"]
git-tree-sha1 = "619f71da7bd903fe7c10408e1ff34d2da17a8ab4"
uuid = "f761d5c5-86db-4880-b97f-9680a7cccfb5"
version = "0.7.0"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f295e0a1da4ca425659c57441bcb59abb035a4bc"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.8"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "03fec6800a986d191f64f5c0996b59ed526eda25"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.4.1"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "d5fb407ec3179063214bc6277712928ba78459e2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.4"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

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
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StructEquality]]
deps = ["Compat"]
git-tree-sha1 = "192a9f1de3cfef80ab1a4ba7b150bb0e11ceedcf"
uuid = "6ec83bb0-ed9f-11e9-3b4c-2b04cb4e219c"
version = "2.1.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

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

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulAtomic]]
deps = ["Unitful"]
git-tree-sha1 = "903be579194534af1c4b4778d1ace676ca042238"
uuid = "a7773ee8-282e-5fa2-be4e-bd808c38a91a"
version = "1.0.0"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "04a51d15436a572301b5abbb9d099713327e9fc4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.4+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

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

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

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

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

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

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.spglib_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "Pkg"]
git-tree-sha1 = "aaaa4deac77ded775242a698a46649811bd7a847"
uuid = "ac4a9f1e-bdb2-5204-990c-47c8b2f70d4e"
version = "1.16.5+0"

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
# ╟─c498a8aa-df5c-4536-b54b-db5912fb411c
# ╟─0f8dc304-6ade-4fbe-8b64-ca368defa2ba
# ╠═80ffd650-d95b-471a-9d2f-e2bc781c322e
# ╟─dddf152e-b1a3-470d-93b7-d213544f0215
# ╟─5729f2a4-7557-4e5c-b230-31faf6286ff3
# ╠═f8827f52-cba1-433f-9c14-9a0b621b2c35
# ╠═23f9fd84-0b58-4ade-a72c-b3e303a79480
# ╠═eb66e66c-278f-465e-a56d-cac7432723f7
# ╠═bfff1b7f-bf26-4373-a499-014b871067ad
# ╟─72264fdd-d151-41b3-abd7-3fd6b0b4fb84
# ╟─fa98f7d8-bc72-4ac8-b46a-2456e37e42de
# ╟─291cb850-5089-480e-b803-1e2a20852b95
# ╠═98c91264-19cc-41e5-951a-a201601fbfc5
# ╟─aea3267b-1164-418c-98a8-70cc37d1c959
# ╠═e33b4164-b980-446f-9398-9ec68cf224f9
# ╟─7661b0fa-839a-47a1-a0ea-a89e8d4abc3a
# ╠═f923e4ba-1833-4b2e-968d-8527fc70720d
# ╟─6a049358-0e6c-4ea9-8950-52f5bdda287e
# ╠═b0e6abc4-7581-4a7b-9160-792d2824a1cf
# ╟─6be8652b-525c-4efa-bcc7-b72285ea2d98
# ╠═ab7f07f6-20f3-465f-ae46-e694100541b8
# ╠═1f9b91e3-28f4-4e06-8a3d-70c84432007b
# ╠═ce876450-eea7-4088-9aca-ed3ee462cb09
# ╟─0f7709e7-fec9-4176-bf58-96f2d891e8b4
# ╟─73e61d6c-0762-41d9-a0c1-c2a64926a40a
# ╟─c3c91935-60f7-4c4a-b43c-973566d58367
# ╠═41265945-67e6-4e34-8604-20faa6f475b7
# ╟─8937b497-99f0-4d8f-998a-e2aad0d5aeb9
# ╠═f64ab474-86f6-4048-b013-3c8969d2e125
# ╠═d21b8a45-7129-4038-82ac-250436ff8583
# ╠═d2ef7881-bd1b-4bd3-aee4-025c7da5ba71
# ╠═6f81e36d-7892-42ec-aeb6-91fca5b6399a
# ╟─137f3b99-d031-4df4-868e-9b5e4d342a83
# ╟─357cd9c6-eb57-48ca-afe5-a752addc6461
# ╠═7237c2e9-c0f3-4cd8-ad1f-6d5cbff35354
# ╠═2e96d48d-15d9-40cc-8c1d-d8d03e755c2a
# ╠═e6a62b4f-9919-43b9-ae74-3cb3c1110829
# ╟─bef6d532-6131-4db7-8b7c-582895a035b1
# ╟─3cea29a7-fd87-4491-81a6-c53f79387e1f
# ╟─e1fd4f96-a465-46d8-8c78-cde7c5325e6f
# ╟─c56ac72f-643a-4b29-97cb-4f61426d11d5
# ╠═1d099205-9d84-4ecc-a2af-81910326a4ba
# ╠═f17500cc-6e5f-442d-9b22-2dc406367ad6
# ╠═320a4138-1d2d-4039-adaa-1aeb7f1929cf
# ╟─db913617-2b81-4ea7-a907-033dd31bda5a
# ╠═68fff309-a3ec-4ebc-bcfe-bafeac6ef3f5
# ╟─f4bb1d69-1f06-4769-ac6e-081ddaa437d7
# ╟─8481465b-171b-4374-8ec9-d7c19bd23d81
# ╟─179be577-6d0e-4e97-8db1-56e116c502c6
# ╠═94e6808b-c89c-4431-90d8-5e252e38d834
# ╠═f1398a23-90d9-40ea-8ef2-ebcd0ab67b14
# ╠═bba31f5b-7eca-412b-99c7-65adcfa2be1b
# ╟─9b7bdce7-38cb-4c7c-b1b4-d1ec3c9fa168
# ╠═42777c48-0f3a-4c43-974b-7320ad16efb9
# ╠═05175611-34de-44cf-a385-c493a6f3cf19
# ╟─15220d0d-2c9f-4daf-a290-672ca439faec
# ╠═2a26de48-c3b6-405b-8c07-bae036dd805d
# ╠═23990899-05d2-4b40-9d52-23caf218c7cf
# ╟─461162d6-b0e6-41c4-8e60-6293d2fc25ee
# ╟─92a6feb4-88ad-4681-b6dd-2df2a790360c
# ╠═99a01adc-34e2-4284-b918-4bfc020b47d6
# ╟─aa9f8c4b-b733-48c6-93f0-4ab8d52b2a3e
# ╠═1cf09a9d-2e2b-4a0c-848c-56d3bd1e7f87
# ╠═39215879-264c-4c2a-bafc-a792abf70017
# ╠═ff50c8b6-ce84-4d74-b5ba-f5abecdd7fb2
# ╟─f9e4e29c-247d-4896-90d0-9c5f716c203c
# ╟─2f4000e9-eaca-446d-af43-3b8ea06eb2b1
# ╠═9d86090f-9680-43fa-b86b-2fc664209faf
# ╟─fb887ea1-b4bb-49d5-af84-9a1278791d67
# ╠═eb909306-b913-464b-81d6-becc9946e7c1
# ╠═aa4f358a-e8fa-4258-bcb6-b6707ef62b0a
# ╠═bdcc33c7-eda7-48fe-a73c-84c2e7c64dab
# ╟─9759f044-eeba-4263-84e6-0113b7310c0b
# ╠═b2cbf4dd-f203-4028-90b1-57d72d21c3c5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
