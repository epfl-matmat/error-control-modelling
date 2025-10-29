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

# ╔═╡ 80ffd650-d95b-471a-9d2f-e2bc781c322e
begin
	using LinearAlgebra
	using LinearMaps
	using PlutoUI
	using PlutoTeachingTools
	using Printf
	using Plots
	using LaTeXStrings
	using HypertextLiteral
end

# ╔═╡ 2a26de48-c3b6-405b-8c07-bae036dd805d
using DFTK

# ╔═╡ e6e2ad31-0d6d-4e59-b80d-beb736759894
md"""
# Overview of diagonalisation algorithms
"""

# ╔═╡ c4d6b47c-6ff7-4e3c-ba16-58a2e08815bc
md"""
This lecture will provide a selective overview of standard diagonalisation algorithms, contrasting their respective scope of applicability and the main ingredients. Our main angle will be to understand questions related to numerical stability in the key ingredients of the algorithms. Instead of taking a proof-guided approach, we will mostly follow a computational approach using the increased-precision and interval techniques we discussed in previous lectures.

For a more gentle introduction see the [chapter on eigenvalue problems in the Numerical Analysis lecture](https://teaching.matmat.org/numerical-analysis/11_Eigenvalue_problems.html).
For a more comprehensive treatment of the topic see the book [Numerical Methods for Large Eigenvalue Problems](https://epubs.siam.org/doi/book/10.1137/1.9781611970739) by Youssef Saad as well as the [Lecture notes on Large Scale Eigenvalue Problems](https://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf) by Peter Arbenz.

In our discussion we will always take $A \in \mathbb{C}^{N\times N}$ to be a Hermitian matrix and we will seek approximations to the eigenpairs $(\lambda_i, x_i) \in \mathbb{R} \times \mathbb{C}^N$, i.e.
```math
A x_i = \lambda_i x_i
```
where we impose an ordering $\lambda_1 \leq \lambda_2 \leq \cdots \leq \lambda_N$.
"""

# ╔═╡ 5729f2a4-7557-4e5c-b230-31faf6286ff3
md"""
## Single vector iteration

### Power method

We already discussed the **power method** in previous lectures and exercises as a way to obtain the eigenvector $x_J$ corresponding to the largest-absolute eigenvalue. Via the Rayleigh quotient the corresponding eigenvalue $\lambda_J$ then becomes accessible as well.

**Remark**. Keep in mind that the **largest-absolute eigenvalue** $\lambda_J$ is in general **not the same as the largest eigenvalue** $\lambda_N$, in case a negative eigenvalue is dominating.

For completeness, our algorithmic recipe to obtain the eigenvector to tolerance $\varepsilon$ was:

> 0. Start from a vector $x^{(0)}$, which has non-zero overlap with $x_n$.
> 1. Compute the matrix-vector product: $\tilde{x}^{(i)} \leftarrow A x^{(i)}$
> 2. Normalise: $x^{(i+1)} \leftarrow \tilde{x}^{(i)} / \| \tilde{x}^{(i)} \|$
> 3. If $\| x^{(i+1)} - x^{(i)} \| < \varepsilon$, exit the procedure, else  $i \leftarrow i+1$ and return to 1.

Below we provide a slightly extended version, which additionally employs the Rayleigh quotient to estimate the eigenvalue and the residual
```math
r^{(i)} = A x^{(i)} - \lambda^{(i)} x^{(i)}
```
to check convergence.
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
From playing with the slider we notice that when the two **largest-absolute eigenvalues get too close**, the **convergence deteriorates**. This can be understood using the following theorem:

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
Part 1 of this theorem is the justification for using the **inverse power method** to approximate the smallest eigenpair of $A$ as done already in one of the previous exercises. Part 3 of this theorem, called **shift and invert** provides us with a recipe to compute an arbitrary isolated eigenvalue.

To illustrate this we consider the matrix
"""

# ╔═╡ 98c91264-19cc-41e5-951a-a201601fbfc5
Ashift = [0.4 -0.6 0.2; -0.3 0.7 -0.4; -0.1 -0.4 0.5]

# ╔═╡ aea3267b-1164-418c-98a8-70cc37d1c959
md"""which has eigenvalues"""

# ╔═╡ e33b4164-b980-446f-9398-9ec68cf224f9
eigvals(Ashift)

# ╔═╡ 7661b0fa-839a-47a1-a0ea-a89e8d4abc3a
md"""Suppose we **don't know the exact eigenvalues**, but we know that the eigenpair $(λ_i, x_i)$ we are interested in has an **eigenvalue around $0.4$**:
- Then we set $σ = 0.4$ and employ the matrix $(A - σ I)^{-1}$ in the `power_method`, which as per our previous discussion will converge to the largest absolute eigenvalue of  $(A - σ I)^{-1}$.
- The eigenvalues of the latter matrix are $\left\{\frac{1}{λ_i - σ} \, \middle| \, i =1,\ldots,3 \right\}$. Clearly, the largest-absolute one of these eigenvalues is attained for the $λ_i$ closest to $σ$.
- From Theorem 5.2 we additionally know that the eigenvector of $(A - σ I)^{-1}$ associated to $\frac{1}{λ_i - σ}$ is exactly the eigenvector of $A$ associated to $λ_i$. Therefore, we have found a way to compute the *exact* eigenpair $(λ_i, x_i)$ based on the approximate eigenvalue $σ = 0.4$.

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

# ╔═╡ aa642130-d82c-4f80-938a-01bfec29ce43
md"""For more details, see the [discussion on spectral transformations](https://teaching.matmat.org/numerical-analysis/11_Eigenvalue_problems.html#Spectral-transformations) in the numerical analysis Bachelor class.
"""

# ╔═╡ b4d035cc-b2f3-495e-9895-d0a3c0790f12
md"""
## Optional: Dynamic shifting

See [discussion in numerical analysis](https://teaching.matmat.org/numerical-analysis/11_Eigenvalue_problems.html#Optional:-Dynamic-shifting) Bachelor course.
"""

# ╔═╡ 6a049358-0e6c-4ea9-8950-52f5bdda287e
md"""
## Optional: Rayleigh quotient iterations

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

With this setup we obtain **Rayleigh-quotient iteration (RQI)**:
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

# ╔═╡ 0f7709e7-fec9-4176-bf58-96f2d891e8b4
md"""Notice, how the convergence in the last few steps is *very fast*, meaning that multiple orders of magnitude in the residual norm are gained in each step.

The RQI has fantastic rates of convergence. However, it has two disadvantages:
- The main issue is that *at each RQI iteration* the matrix `A - λ*I` is different and needs to be **freshly factorised** in order to compute the `\` operation. In other words we pay a full Gaussian elimination each RQI iteration. This in contrast to the *shift-and-invert* power method, where factorisation was computed once and for all for `A - σ*I` *before* the `power_method` was called.
- The second point to note is that it is **hard to control to *which* eigenpair** the RQI converges. Try re-running the above a few times by updating the initial guess `u0_rqi`.
- Moreover, it is possible for the RQI to *not* converge --- albeit the set of possible initial guesses that leads to a non-converging RQI is of measure zero. It is thus most useful if one already has a decent guess for the eigenvector.

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

# ╔═╡ b48af18e-c65b-4ddd-b6d1-6b25e4e4bb80
md"""
To illustrate the point about the **target of RQI's convergence being hard to predict** consider the following example. We initialise the RQI using a random starting vector
"""

# ╔═╡ b8bb7d0b-21c4-4746-9dc7-dc47f76a0736
xstart = randn(eltype(Arqi), size(Arqi, 2))

# ╔═╡ 009c780c-d1cc-4421-86c7-74e5a99646f9
md"""
and allow ourself to select an *initial* shift $σ$,
which we store in the variable `σD`:
- `σD = `$(@bind σD Slider(-0.2:0.01:0.2; default=0.01, show_value=true))

Then we define a version of `rayleigh_quotient_iteration`, which starts the procedure based on applying the matrix $(A - σ I)^{-1}$ to this random vector:
"""

# ╔═╡ 4e208b64-6b88-4430-9bc3-7198c54b00ff
function rayleigh_quotient_iteration(A, σ::Number, x; kwargs...)
	x = (A - σ*I) \ x
	rayleigh_quotient_iteration(A; x, kwargs...)
end

# ╔═╡ ce876450-eea7-4088-9aca-ed3ee462cb09
rayleigh_quotient_iteration(Arqi, x=x0_rqi)

# ╔═╡ 8b2f116a-09e6-4bbf-9a83-b7c38fe34680
let
	# exact eigenvalues of Arqi
	λ_rqi = eigvals(Arqi)
	
	# The eigenvalue the iteration targets
	i = argmin(abs.(λ_rqi .- σD))
	λtarget = λ_rqi[i]

	p = plot(title="Convergence",
		ylabel="residual", xlabel=L"i", legend=:bottomleft,
		yaxis=:log, ylims=(1e-10, 3), xlims=(1, 7))

	q = plot(title="History", xlabel=L"i", ylabel=L"λ_k^{(i)}",
	         legend=:bottomleft, ylims=(-0.2, 0.2), xlims=(1, 7))
	hline!(q, λ_rqi, ls=:dash, label=L"eigenvalues of $A$", c=2, lw=1.5)
	hline!(q, [σD], ls=:dash, label=L"initial shift $σ$", c=3, lw=1.5)

	if abs(σD - λtarget) > 1e-10  # Guard to prevent numerical issues
		results = rayleigh_quotient_iteration(Arqi, σD, xstart; maxiter=7)
		plot!(p, results.residual_norms; mark=:o, label="", lw=1.5)
		plot!(q, results.eigenvalues, mark=:o, c=1, label="", lw=1.5)
	end
	plot(p, q; layout=(1, 2))
end

# ╔═╡ c3c91935-60f7-4c4a-b43c-973566d58367
md"""
## Optional: Subspace methods

So far we restricted ourselves to the rather basic goal of solving only for a **single eigenvector** --- essentially by iterating only a single vector. We already saw that some eigenproblems are challenging to solve in this simple setting, namely exactly if the **largest-absolute eigenvalue has only a small or even a zero gap**. 

A natural generalisation to the power-method scheme is to iterate not one, but **multiple vectors**. Doing this naively is a little flawed as then simply each vector will converge to the eigenvector(s) with the largest-abs eigenvalue. Instead, what is needed is to ensure *in each step* that the obtained eigenvectors are orthogonal to each other.

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
		# Full orthogonalisation
		X = ortho(X)

		# Just normalise ... for testing
		# for c in eachcol(X)
		# 	normalize!(c)
		# end

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

# ╔═╡ 9d5e0e5d-9117-4308-b133-20e3d9e0894e
md"The method yields the following eigenvalues:"

# ╔═╡ bd576cd6-ee70-43b3-ad05-9307d6eeb05d
subspace_iteration(Asubspace; X=randn(size(Asubspace, 2), 2), verbose=false, tol=1e-6, maxiter=200).λ

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

The above `subspace_iteration` addresses the question of how to compute multiple eigenpairs at once. However, compared to Rayleigh quotient iteration it feels a little like a step back since we are no longer making use of the Rayleigh quotient. Notably, when computing the approximate eigenpairs we are considering each iterated vector separately and allow no "interaction" between them, even though the first iterated vector might well contain a contribution to the eigenvector, which the second iterated vector attempts to approximate. Moreover, we are right now tied to using exactly one iterated vector per eigenvector we want to approximate.

In the next section we will see how projection methods allow us to overcome both these restrictions.
"""

# ╔═╡ 1b63714a-1e7f-40f0-9f2f-5a03744a44d5
md"""
## Projection methods

The assumption underlying projection methods is that we are given (by some procedure discussed below) an $m$-dimensional subspace $\mathcal{S} \subset \mathbb{C}^N$ with $m \ll N$. Instead of solving the full eigenproblem, we now restrict ourselves to **solve the eigenvalue problem *projected* into this subspace**. We consider the basic eigenproblem for $A \in \mathbb{C}^{N \times N}$ Hermitian:
```math
A x_i = λ_i x_i \qquad \text{$(\ddagger)$}
```
An orthogonal projection technique seeks approximate eigenpairs $(\tilde{\lambda}_i, \tilde{x}_i)$ with $\tilde{\lambda}_i \in \mathbb{R}$, $\tilde{x}_i \in \mathcal{S}$, such that the following **Galerkin condition**
```math
\big(A \tilde{x}_i - \tilde{λ}_i \tilde{x}_i \big) \perp \mathcal{S} \qquad \text{(G)}
```
is satisfied. Equivalently one may enforce
```math
\left\langle v, A \tilde{x}_i - \tilde{λ}_i \tilde{x}_i \right\rangle = 0 \qquad \forall v \in \mathcal{S}. \qquad \text{$(\ast)$}
```
"""

# ╔═╡ 934ec8b3-a0c0-4129-94f1-f3f04c619365
md"""
Suppose now that $\{v_1, v_2, \ldots, v_m\}$ is an orthonormal basis for $\mathcal{S}$ and denote by $V = (v_1, v_2, \ldots, v_m) \in \mathbb{C}^{N\times m}$ the matrix with these basis vectors as columns. Expanding $\tilde{x}_i$ in this basis, i.e. introducing a vector $\tilde{y}_i \in \mathbb{C}^m$ with
```math
\tilde{x}_i = V \tilde{y}_i
```
equivalent to
```math
(\tilde{x}_i)_α = \sum_k V_{αk} (\tilde{y}_i)_k
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
"""

# ╔═╡ 84e8b9a8-eb55-46e9-81ab-df248d98a800
md"""
> 1. Find an orthonormal basis $\{v_i\}_{i=1,\ldots,m}$, spanning a subspace $\mathcal{S}$.
> 2. Compute $A_V = V^H A V$, with $V$ containing the basis vectors of $\mathcal{S}$ as columns.
> 3. Use `eigen` to compute $m$ eigenpairs $(\tilde{λ}_i, \tilde{y}_i)$ of $A_V$.
> 4. Compute the approximate eigenvectors $\tilde{x}_i = V \tilde{y}_i$ of $A$, which yields the approximate eigenpairs of $A$ as $(\tilde{λ}_i, \tilde{x}_i)$.
"""

# ╔═╡ 0d96be04-f200-43c1-8374-8a1ec85be47f
md"""
Note, that the low-dimensional pairs $(\tilde{λ}_i, \tilde{y}_i)$ are exact eigenpairs of the subproblem matrix $A_V$. The corresponding full-dimensional pairs $(\tilde{λ}_i, \tilde{x}_i)$ are usually called **Ritz pair** with $\tilde{λ}_i$ called **Ritz values** and $\tilde{x}_i$ a **Ritz vector**. 
"""

# ╔═╡ b02d908b-2882-42b4-803c-fa2b9ca3833f
md"""As an example we consider the diagonal matrix"""

# ╔═╡ 63ad906a-0fce-4b3e-a5ef-3001c0a1c315
D  = Diagonal(0:0.1:1)

# ╔═╡ cde89859-7504-4730-afb1-32919d4c5dd2
md"""has smallest eigenvalue zero."""

# ╔═╡ 13bdb670-8235-47f0-a2fe-ff9b6d12afeb
md"""
1. We find an orthonormal basis. Here we just take a random one by orthogonalising $nᵥ random vectors.
"""

# ╔═╡ 5e02082e-0b07-404b-a5f9-8f86155f4381
V = ortho_qr(randn(11, 3)); # Generate 3 random orthogonal vectors

# ╔═╡ 0076b5f3-b038-4d8a-8218-dfac8fb8f731
md"""
2. We compute the projection into the basis
"""


# ╔═╡ c8b15905-8529-4da0-9944-313e271e0665
Dᵥ = V'*D*V

# ╔═╡ 5374fdc7-efac-4044-8757-660120992877
md"""
3. and compute the Ritz pairs
"""

# ╔═╡ 56cfee0e-4eec-49e8-9c8d-6ef3552f7c23
Λᵥ, Yᵥ = eigen(Dᵥ)

# ╔═╡ 5ff2388d-bda2-42bb-957d-6605a1e343ac
md"""
4. The approximation to the eigenvalues are thus
"""

# ╔═╡ 36e7afa4-99ed-4dde-9129-e696759cbf95
Λᵥ

# ╔═╡ 899b256e-e392-4785-8097-6bf7fad79076
md"""
4. and the approximation to the eigenvectors are
"""

# ╔═╡ cd79bd02-2d61-4eaf-9491-796ff041bdba
V * Yᵥ

# ╔═╡ 6c1decec-0f96-4d19-b969-29a4344fdbb7
md"""
### Relation to Courant-Fisher

The **appeal of projection methods** for computing the eigenpairs of Hermitian matrices results from the connection to the **Courant-Fisher min-max principle**, which guarantees that the eigenvalues computed by the Rayleigh-Ritz procedure satisfy a strong optimality condition: Eigenvalues are approximated from above.

To see this, let us introduce a handy **alternative notation** for subspace methods given **in terms of projections**.
"""

# ╔═╡ 4413399f-1a77-44f1-8e71-7fdf50233ee8
md"""
Making use of our basis for $\mathcal{S}$ it is easy to verify that
```math
P_\mathcal{S} = V V^H \in \mathbb{C}^{N\times N}
```
is indeed a **projection onto the subspace** $\mathcal{S}$, since the idempotency relation
```math
P_\mathcal{S}^2 = \left( V V^H \right)^2 = V (V^H V) V^H = V I V^H = P_\mathcal{S}
```
is satisfied. 
"""

# ╔═╡ e4a01539-ef46-4e4c-86c4-5c37f34a0336
Foldable("""Rewriting (§) and (G) in terms of projectors""",
md"""
Using projector notation also the Galerkin and projective eigenvalue problems can be rewritten, which we provide here for completeness.
		 
Since $\tilde{x}_i = V \tilde{y}_i$ implies $V^H \tilde{x}_i = \tilde{y}_i$ and thus $P_\mathcal{S} \tilde{x}_i = \tilde{x}_i$,
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
""")

# ╔═╡ befcfdd1-1f57-45ab-8f96-ae4d01895eea
md"""
Based on the above projectors we note
```math
\langle \tilde{x}, A \tilde{x} \rangle
= \big\langle \tilde{x}, \left(P_\mathcal{S} A P_\mathcal{S} \right)\, \tilde{x} \big\rangle
\qquad \forall \tilde{x} \in \mathcal{S},
```
from which it follows that
if $\tilde{x}$ is inside the subspace employed by an orthogonal projection approach (like the Rayleigh-Ritz procedure), the Rayleigh quotient of the projected matrix $P_\mathcal{S} A P_\mathcal{S}$ and of the original matrix $A$ are identical.
Therefore:
```math
\tilde{λ}_1 = \min_{\tilde{x}\in \mathcal{S}, \tilde{x}\neq 0} \frac{\big\langle \tilde{x}, \left(P_\mathcal{S} A P_\mathcal{S} \right) \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
= \min_{\tilde{x}\in \mathcal{S}, x\neq 0}
\frac{\big\langle P_\mathcal{S} \tilde{x}, A P_\mathcal{S} \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
= \min_{\tilde{x}\in \mathcal{S}, \tilde{x}\neq 0}
\frac{\big\langle \tilde{x}, A \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}.
```



More generally, based on applying Courant-Fisher to the matrix $P_\mathcal{S} A P_\mathcal{S}$, we obtain the following result
"""

# ╔═╡ 0233c113-c3ef-456a-aa50-17eb754608ac
md"""
**Proposition 5.5.** The $i$-th approximate eigenvalue (counted from small to large including multiplicities) of a Hermitian matrix $A$ obtained from an orthogonal projection method to a subspace $\mathcal{S}$ satisfies
```math
\tilde{λ}_i = \min_{s \subset \mathcal{S},\,\text{dim}(s) = i} \ \max_{\tilde{x}\in s,\,\tilde{x}\neq 0} \frac{\big\langle \tilde{x}, A \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
```

From this we obtain immediately:
"""

# ╔═╡ e51a1e89-1675-4b1f-9782-7951b60e1658
md"""
**Corollary 5.6.** For $i = 1, 2, \ldots, m$ the inequality
```math
λ_i \leq \tilde{λ}_i
```
holds, i.e. eigenvalues are **approximated from above** in orthogonal projection methods.
"""

# ╔═╡ 413d5e90-0b44-4e3f-8ab8-15ac1538c1c7
md"""
**Proof.**
```math
\tilde{λ}_i = \min_{s \subset \textcolor{red}{\mathcal{S}},\,\text{dim}(s) = i} \ \max_{\tilde{x}\in s,\,\tilde{x}\neq 0} \frac{\big\langle \tilde{x}, A \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
\geq
\min_{s \subset \textcolor{red}{\mathbb{C}^N},\,\text{dim}(s) = i} \ \max_{x\in s,\,x\neq 0} \frac{\big\langle x, A x \big\rangle}{\langle x, x \rangle}
= λ_i.
```
This is the **desired optimality**: The eigenvalue approximation found by a projection method is the **best eigenvalue approximation we can possibly find** within the chosen subspace, and it provides a rigorous upper bound to the exact one.
"""

# ╔═╡ 262eca5b-d55e-4428-a11a-0614751dd2d6
md"""
	Finally, the following characterisation provides a connection between the approximate eigenvalue $\tilde{λ}_i$ and the corresponding approximate eigenvector $\tilde{x}_i$ obtained by an orthogonal projection approach:
	```math
	\tilde{λ}_i
	= \frac{\big\langle \tilde{x}_i, A \tilde{x}_i \big\rangle}{\langle \tilde{x}_i, \tilde{x}_i \rangle}
	= \min_{0\neq \tilde{x} \in \mathcal{S},\,\tilde{x} \perp \tilde{\mathcal{X}}_{k-1}} \frac{\big\langle \tilde{x}, A \tilde{x} \big\rangle}{\langle \tilde{x}, \tilde{x} \rangle}
	```
	where $\tilde{\mathcal{X}}_k = \text{span}(\tilde{x}_1, \ldots, \tilde{x}_k)$,
	the space spanned by the $k$ first *approximate* eigenvectors.
	"""

# ╔═╡ 24dd589f-3275-4924-ba9e-9ad5c7fd075a
md"""
Sticking with our example from above, we indeed see how increasing the subspace size `nᵥ` makes the approximate eigenvalues approach the true eigenvalues from above.
"""

# ╔═╡ 814837a2-761b-4233-81f7-ded2907cb07f
let
	Vr = ortho_qr(randn(size(D)))
	
	p  = plot(; xlabel="Subspace size nᵥ", ylabel="Approx. eigenvalues")
	Λs = fill(NaN, size(D))
	for nv in 1:size(D, 2)
		V  = Vr[:, 1:nv]
		Λs[1:nv, nv] = eigvals(V'*D*V)
	end
	
	for nv in 1:size(D, 2)
		plot!(p, Λs[nv, :]; mark=:o, markersize=5, lw=2, label="")
	end

	p
end

# ╔═╡ e1fd4f96-a465-46d8-8c78-cde7c5325e6f
md"""
### Forming a good subspace

Projection methods and the Rayleigh-Ritz procedure are a key ingredient of pretty much all iterative diagonalisation approaches employed nowadays. A full discussion of the standard techniques employed to build the reduced subspace $\mathcal{S}$ is out of the scope of this lecture. An incomplete list of techniques worth mentioning are:
- If a Krylov basis is constructed and employed for $\mathcal{S}$ one obtains diagonalisation methods such as *Lanczos* or *Arnoldi*.
- Employing an orthogonal correction based on the current residual to expand and enrich the subspace leads to the *Davidson* family of methods.
- If the subspace is constructed following the idea of minimising the Rayleigh quotient, we obtain minimisation methods (inverse iterations, PCG, LOPCG). These will be discussed in the following.
"""

# ╔═╡ bebb6206-9077-40b4-8bf6-63cf5f23d308
md"""
### Optional: Application to power iterations

One way to put this into practice is to employ the ideas of the `subspace_iteration` discussed above. Recall that this method generates in each iteration a set of orthogonal vectors `U` by applying `A` to the previous set of vectors. In the `subspace_iteration` algorithm discussed above, these vectors are directly taken as the current estimates for the eigenvectors.

However, these vectors also span a subspace, enabling them to be employed within a Rayleigh-Ritz procedure, leading to the projected subspace iterations:
"""

# ╔═╡ 7237c2e9-c0f3-4cd8-ad1f-6d5cbff35354
function projected_subspace_iteration(A, V;
                                      tol=1e-6, maxiter=100, verbose=true,
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

# ╔═╡ a9fdeb7d-e4f8-4320-8d32-792bf405a7aa
# Version to automatically select initiol guess
function projected_subspace_iteration(A;
							X=randn(eltype(A), size(A, 2), 2), kwargs...)
	projected_subspace_iteration(A, X; kwargs...)
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

# ╔═╡ ccd01cf4-e25d-4c4d-99dd-63a12de4e65d
md"""
## Minimisation methods

In our discussion on [Projection methods](#Projection-methods) we already employed the variational ansatz encoded in Courant-Fisher to argue why increasing the subspaces consecutively will eventually lead to a perfect approximation of the eigenpairs.
"""

# ╔═╡ 6b997532-b424-4df6-afdd-ded41bfba8a8
md"""
An orthogonal idea is to directly work on minimising the Reighlay quotient. E.g. since $λ_1$ may just be obtained from
```math
\tag{M}
λ_1 = \min_{0 \neq x \in \mathbb{C}^N} R_A(x) \qquad R_A(x) = \frac{x^H A x}{x^H x}
```
we may think of methods to directly solve this optimisation problem.
This will be the purpose of this section. Amongst others this leads to the LOPCG algorithm, which interestingly yields an iteration of 3-dimensional subspaces to be used within the [Projection methods](#Projection-methods) discussed above.
"""

# ╔═╡ aaf81e2f-5e6e-4f5d-8950-831ae85f2476
md"""
When considering the minimisation problem (M)
the natural reflex is to employ a steepest descent approach.
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
"""

# ╔═╡ 37b95aee-ded0-4698-ae1b-4f236ca4cffb
md"""
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
**Preconditioned gradient descent** is **closely related to the inverse power method**. Namely, if we chose $P^{-1} = A^{-1}$ and $α = 1$, then we obtain
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

# ╔═╡ 652c5513-8263-4211-833f-aa8180711d68
TODO("mention Nystrom approximation")

# ╔═╡ 4b0cd6fc-50ec-4e12-82e0-30fe2ceb98ce
md"""
### Step sizes and LOBPCG

The final aspect we want to briefly address is **how to choose the step size** $α$. The best step size can be found by a **line search**,
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
i.e. as a **minimisation in the subspace spanned** by $x^{(i)}$ and $\tilde{r}^{(i)}$.
"""

# ╔═╡ 9fdc92f6-158f-4e62-82e0-b6c5ce96d9a8
md"""
A **slight generalisation** of this idea is to make the minimisation **subspace even larger** by also including the *previous* iterate to the search space, i.e. to employ
```math
x^{(i+1)} = \text{argmin}_{y \in \text{span} \left\{ x^{(i)}, x^{(i-1)}, \tilde{r}^{(i)}  \right\}} R_A(y). \qquad \text{$(\ast)$}
```
This turns out to improve the overall observed convergence rates considerably. This modification has been suggested following a more detailed investigation of minimisation approaches for eigenvalue problems, in particular a comparison with more sophisticated minimisation algorithms such as the conjugate gradient (CG) approach, which is out of the scope of this lecture.

The improved convergence rates of $(\ast)$ unfortunately comes with a downside: As the iterations converge $x^{(i)}$ and $x^{(i-1)}$ *both* converge to the exact eigenvector, such that $x^{(i)}$ and $x^{(i-1)}$ become more and more **linearly dependent**, which makes the optimisation unstable. This is in practice circumvented by two further modifications:
- Instead of using the basis vectors $\{x^{(i)}, x^{(i-1)}, \tilde{r}^{(i)}\}$ to build the subspace, we compute $p^{(i)} = x^{(i)} - x^{(i-1)}$
  and instead employ $\{ x^{(i)}, p^{(i)}, \tilde{r}^{(i)} \}$,
   which span the same space.
- Instead of first performing the line search to find the optimal $x^{(i+1)}$ and then computing the Rayleigh coefficient $R_A(x^{(i+1)})$, we perform a **Rayleigh-Ritz procedure projecting onto the subspace** $\mathcal{S} = \text{span} \{ x^{(i)}, p^{(i)}, \tilde{r}^{(i)} \}$. Due to our results on [projection methods](#Projection-methods) (in particular Proposition 5.5.) we know that a Rayleigh-Ritz in the subspace $\mathcal{S}$ is equivalent to finding the minimiser of $R_A(x)$ where $x \in \mathcal{S}$.

These two modifications with our previous discussion lead to the **LOPCG (locally optimal preconditioned conjugate gradient)** algorithm:
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
The generalisation to an algorithm capable of computing multiple eigenpairs is straightforward, leading to **LOBPCG (locally optimal *block* preconditioned conjugate gradient)** is only bookkeeping:
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
Note, that this simple implementation of LOBPCG is **still not numerically stable** and **considerable issues remain**. In practice one should thus always employ one of the well-tested implementation in standard numerical libraries.

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

- **Rayleigh quotient iteration** generally performs **best**, but is also the most expensive (Full factorisation required at *every iteration*)
- **Inverse power method** and **preconditioned gradient descent fail** when the gap becomes **too small**.
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

# ╔═╡ 6c1c4615-c98d-4fd4-8f5b-c0fe2702204f
TODO("Fill in the blanks here")

# ╔═╡ dddf152e-b1a3-470d-93b7-d213544f0215
TableOfContents()

# ╔═╡ b2cbf4dd-f203-4028-90b1-57d72d21c3c5
let
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	
	Sidebar(Markdown.parse(read("sidebar.md", String)), 540)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DFTK = "acf6eb54-70d9-11e9-0013-234b7a5f5337"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearMaps = "7a12625a-238d-50fd-b39a-03d52299707e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
DFTK = "~0.7.16"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
LinearMaps = "~3.11.4"
Plots = "~1.40.19"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.71"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "a3b0347c287bf84a5471e2cd40d1b00c4b0e321e"

[[deps.ADTypes]]
git-tree-sha1 = "60665b326b75db6517939d0e1875850bc4a54368"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.17.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

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
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "dbd8c3bbbdbb5c2778f85f4422c39960eac65a42"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.20.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AtomsBase]]
deps = ["LinearAlgebra", "PeriodicTable", "Preferences", "Printf", "Requires", "StaticArrays", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "5f76425eb977584353191c41d739e7783f036b90"
uuid = "a963bdd2-2df7-4f54-a1ee-49d51e6be12a"
version = "0.5.1"

    [deps.AtomsBase.extensions]
    AtomsBaseAtomsViewExt = "AtomsView"

    [deps.AtomsBase.weakdeps]
    AtomsView = "ee286e10-dd2d-4ff2-afcb-0a3cd50c8041"

[[deps.AtomsCalculators]]
deps = ["AtomsBase", "Compat", "StaticArrays", "Test", "Unitful"]
git-tree-sha1 = "99c477979252acdb189f3e4161970e42634e426a"
uuid = "a3e0e189-c65a-42c1-833c-339540406eb1"
version = "0.2.3"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bravais]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a39abce7826834ba1a517e44e9ed673e6621192a"
uuid = "ada6cbde-b013-4edf-aa94-f6abe8bd6e6b"
version = "0.2.5"

[[deps.Brillouin]]
deps = ["Bravais", "DirectQhull", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Reexport", "Requires", "StaticArrays"]
git-tree-sha1 = "caf48ac701636e00df517a24ed4bec2943ba8f58"
uuid = "23470ee3-d0df-4052-8b1a-8cbd6363e7f0"
version = "0.5.29"

    [deps.Brillouin.extensions]
    BrillouinMakieExt = "Makie"
    BrillouinPlotlyJSExt = "PlotlyJS"
    BrillouinSpglibExt = "Spglib"

    [deps.Brillouin.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"
    Spglib = "f761d5c5-86db-4880-b97f-9680a7cccfb5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CUDA_Driver_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e6a1d9f5518122c186fd27786b61d2053cfa1b0c"
uuid = "4ee394cb-3365-5eb0-8335-949819d2adfc"
version = "13.0.1+0"

[[deps.CUDA_Runtime_jll]]
deps = ["Artifacts", "CUDA_Driver_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "e24c6de116c0735c37e83b8bc05ed60d4d359693"
uuid = "76a88914-d11a-5bdc-97e0-2f5a05c973a2"
version = "0.19.1+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

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
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "0037835448781bb46feb39866934e243886d756a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ComponentArrays]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "ConstructionBase", "Functors", "LinearAlgebra", "StaticArrayInterface", "StaticArraysCore"]
git-tree-sha1 = "d8b02e2226568644b6758b2d113fe5b08884eec0"
uuid = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
version = "0.15.29"

    [deps.ComponentArrays.extensions]
    ComponentArraysGPUArraysExt = "GPUArrays"
    ComponentArraysKernelAbstractionsExt = "KernelAbstractions"
    ComponentArraysOptimisersExt = "Optimisers"
    ComponentArraysReactantExt = "Reactant"
    ComponentArraysRecursiveArrayToolsExt = "RecursiveArrayTools"
    ComponentArraysReverseDiffExt = "ReverseDiff"
    ComponentArraysSciMLBaseExt = "SciMLBase"
    ComponentArraysTrackerExt = "Tracker"
    ComponentArraysZygoteExt = "Zygote"

    [deps.ComponentArrays.weakdeps]
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Optimisers = "3bd65402-5787-11e9-1adc-39752487f4e2"
    Reactant = "3c362404-f566-11ee-1572-e11a4b42c853"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SciMLBase = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CrystallographyCore]]
deps = ["LinearAlgebra", "StaticArrays", "StructEquality"]
git-tree-sha1 = "521a8ed44455592672d887632eef056c1ba4678c"
uuid = "80545937-1184-4bc9-b283-396e91386b5c"
version = "0.3.3"

[[deps.DFTK]]
deps = ["AbstractFFTs", "Artifacts", "AtomsBase", "AtomsCalculators", "Brillouin", "Dates", "DftFunctionals", "DifferentiationInterface", "DocStringExtensions", "FFTW", "ForwardDiff", "GPUArraysCore", "Interpolations", "IterTools", "KrylovKit", "Libxc", "LineSearches", "LinearAlgebra", "LinearMaps", "MPI", "Markdown", "Optim", "PeriodicTable", "PkgVersion", "Polynomials", "PrecompileTools", "Preferences", "Printf", "PseudoPotentialData", "PseudoPotentialIO", "Random", "Roots", "SparseArrays", "SpecialFunctions", "Spglib", "StaticArrays", "Statistics", "TimerOutputs", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "a05f869e106b85687be5317f8c0b22c61a08972e"
uuid = "acf6eb54-70d9-11e9-0013-234b7a5f5337"
version = "0.7.16"

    [deps.DFTK.extensions]
    DFTKAMDGPUExt = "AMDGPU"
    DFTKCUDAExt = "CUDA"
    DFTKGenericLinearAlgebraExt = "GenericLinearAlgebra"
    DFTKGeometryOptimizationExt = "GeometryOptimization"
    DFTKIntervalArithmeticExt = "IntervalArithmetic"
    DFTKJLD2Ext = "JLD2"
    DFTKJSON3Ext = "JSON3"
    DFTKPlotsExt = "Plots"
    DFTKWannier90Ext = "wannier90_jll"
    DFTKWannierExt = "Wannier"
    DFTKWriteVTKExt = "WriteVTK"

    [deps.DFTK.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GenericLinearAlgebra = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
    GeometryOptimization = "673bf261-a53d-43b9-876f-d3c1fc8329c2"
    IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
    JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
    Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
    Wannier = "2b19380a-1f7e-4d7d-b1b8-8aa60b3321c9"
    WriteVTK = "64499a7a-5c06-52f2-abe2-ccb03c286192"
    wannier90_jll = "c5400fa0-8d08-52c2-913f-1e3f656c1ce9"

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

[[deps.DftFunctionals]]
deps = ["ComponentArrays", "DiffResults", "ForwardDiff"]
git-tree-sha1 = "b904a2c78f0e768187a13f425719557b2d1db47e"
uuid = "6bd331d2-b28d-4fd3-880e-1a1c7f37947f"
version = "0.3.1"

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

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "16946a4d305607c3a4af54ff35d56f0e9444ed0e"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.7"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DirectQhull]]
deps = ["Qhull_jll"]
git-tree-sha1 = "49951cee00be263d4dc187f38aa96333d74f4d1c"
uuid = "c3f9d41a-afcb-471e-bc58-0b8d83bd86f4"
version = "0.2.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

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

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "680a39c9aadce7c721b68d979e66dc65d2021aa6"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.2.2"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "83dc665d0312b41367b7263e8a4d172eac1897f4"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.4"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3a948313e7a41eb1db7a1e733e6335f17b4ab3c4"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "7.1.1+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "797762812ed063b9b94f6cc7742bc8883bb5e69e"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.9.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "31fd32af86234b6b71add76229d53129aa1b87a9"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.28.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "f0090eb9f8e9d151563dd2300fc3ca3f76b90fe8"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.2.0"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

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

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "1828eb7275491981fa5f1752a5e126e8f26f8741"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.17"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "27299071cc29e409488ada41ec7643e0ab19091f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.17+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "35fbd0cefb04a516104b8e183ce0df11b70a3f1a"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.3+0"

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

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "3d468106a05408f9f7b6f161d9e7715159af247b"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.12.2+0"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

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

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "PackageExtensionCompat", "Printf", "Random", "VectorInterface"]
git-tree-sha1 = "38477816f8db29956ea591feb3086d9edabf6f38"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.9.5"
weakdeps = ["ChainRulesCore"]

    [deps.KrylovKit.extensions]
    KrylovKitChainRulesCoreExt = "ChainRulesCore"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

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

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

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
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d3c8af829abaeba27181db4acb485b18d15d89c6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.1+0"

[[deps.Libxc]]
deps = ["Libxc_GPU_jll", "Libxc_jll"]
git-tree-sha1 = "fe7cc52527ba40c87dd5e955efc59d1bf947317c"
uuid = "66e17ffc-8502-11e9-23b5-c9248d0eb96d"
version = "0.3.19"

    [deps.Libxc.extensions]
    LibxcCudaExt = "CUDA"

    [deps.Libxc.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.Libxc_GPU_jll]]
deps = ["Artifacts", "CUDA_Runtime_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "25bde0bdc076ff720ce0b2383e044f0b62048a03"
uuid = "25af9330-9b41-55d4-a324-1a83c0a0a1ac"
version = "7.0.0+2"

[[deps.Libxc_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "94f7feffa99688628909494a890581fe3b6148dd"
uuid = "a56a6d9d-ad03-58af-ab61-878bf78270d6"
version = "7.0.0+2"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "4adee99b7262ad2a1a4bbbc59d993d24e55ea96f"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.4.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearMaps]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7f6be2e4cdaaf558623d93113d6ddade7b916209"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.11.4"
weakdeps = ["ChainRulesCore", "SparseArrays", "Statistics"]

    [deps.LinearMaps.extensions]
    LinearMapsChainRulesCoreExt = "ChainRulesCore"
    LinearMapsSparseArraysExt = "SparseArrays"
    LinearMapsStatisticsExt = "Statistics"

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

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MPI]]
deps = ["Distributed", "DocStringExtensions", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "PkgVersion", "PrecompileTools", "Requires", "Serialization", "Sockets"]
git-tree-sha1 = "a61ecf714d71064b766d481ef43c094d4c6e3c52"
uuid = "da04e1cc-30fd-572f-bb4f-1f8673147195"
version = "0.20.23"

    [deps.MPI.extensions]
    AMDGPUExt = "AMDGPU"
    CUDAExt = "CUDA"

    [deps.MPI.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "d72d0ecc3f76998aac04e446547259b9ae4c265f"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.3.1+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "e214f2a20bdd64c04cd3e4ff62d3c9be7e969a59"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.4+0"

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

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

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

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "ec764453819f802fc1e144bfe750c454181bd66d"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.8+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2ae7d4ddec2e13ad3bddf5c0796f7547cf682391"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.2+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "61942645c38dd2b5b78e2082c9b51ab315315d10"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.2"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

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

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "275a9a6d85dc86c24d03d1837a0010226a96f540"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.3+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.PeriodicTable]]
deps = ["Base64", "Unitful"]
git-tree-sha1 = "238aa6298007565529f911b734e18addd56985e1"
uuid = "7b2266bf-644c-5ea3-82d8-af4bbd25a884"
version = "1.2.1"

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

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

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
git-tree-sha1 = "0c5a5b7e440c008fe31416a3ac9e0d2057c81106"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.19"

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

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.0"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

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

[[deps.PseudoPotentialData]]
deps = ["Artifacts", "Compat", "DocStringExtensions", "LazyArtifacts", "TOML"]
git-tree-sha1 = "754ae22669bfcb466fd4b0b3cacc1060fc0a7c57"
uuid = "5751a51d-ac76-4487-a056-413ecf6fbe19"
version = "0.2.4"
weakdeps = ["AtomsBase"]

    [deps.PseudoPotentialData.extensions]
    PseudoPotentialDataAtomsBaseExt = "AtomsBase"

[[deps.PseudoPotentialIO]]
deps = ["EzXML", "LinearAlgebra"]
git-tree-sha1 = "88cf9598d70015889c99920ff3dacca0eb26ae90"
uuid = "cb339c56-07fa-4cb2-923a-142469552264"
version = "0.1.1"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6f3ac0623e1173c006cc7377798ec3fb33fa504"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1004+0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "eb38d376097f47316fe089fc62cb7c6d85383a52"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "e1d5e16d0f65762396f9ca4644a5f4ddab8d452b"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

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
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "8a433b1ede5e9be9a7ba5b1cc6698daa8d718f1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

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

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Spglib]]
deps = ["CrystallographyCore", "StaticArrays", "StructEquality", "spglib_jll"]
git-tree-sha1 = "e1f719bd2d19b3014d28c29a4ba51eab126ad605"
uuid = "f761d5c5-86db-4880-b97f-9680a7cccfb5"
version = "0.9.5"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "f737d444cb0ad07e61b3c1bef8eb91203c321eff"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.2.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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

[[deps.StructEquality]]
deps = ["Compat"]
git-tree-sha1 = "192a9f1de3cfef80ab1a4ba7b150bb0e11ceedcf"
uuid = "6ec83bb0-ed9f-11e9-3b4c-2b04cb4e219c"
version = "2.1.0"

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

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

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

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

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
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

[[deps.UnitfulAtomic]]
deps = ["Unitful"]
git-tree-sha1 = "903be579194534af1c4b4778d1ace676ca042238"
uuid = "a7773ee8-282e-5fa2-be4e-bd808c38a91a"
version = "1.0.0"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "af305cc62419f9bd61b6644d19170a4d258c7967"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.7.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9166406dedd38c111a6574e9814be83d267f8aec"
uuid = "409d34a3-91d5-4945-b6ec-7529ddf182d8"
version = "0.5.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "59071150afa35787c1656ba234cf03fdf8e2603f"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.8+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

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

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

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

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "c5bf2dad6a03dfef57ea0a170a1fe493601603f2"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.5+0"

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

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

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
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

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

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.spglib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "bc328924cf4975fe49e6416f7e1622e8ceda55e8"
uuid = "ac4a9f1e-bdb2-5204-990c-47c8b2f70d4e"
version = "2.1.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "fbf139bce07a534df0e699dbb5f5cc9346f95cc1"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.9.2+0"
"""

# ╔═╡ Cell order:
# ╟─e6e2ad31-0d6d-4e59-b80d-beb736759894
# ╠═80ffd650-d95b-471a-9d2f-e2bc781c322e
# ╟─c4d6b47c-6ff7-4e3c-ba16-58a2e08815bc
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
# ╟─aa642130-d82c-4f80-938a-01bfec29ce43
# ╟─b4d035cc-b2f3-495e-9895-d0a3c0790f12
# ╟─6a049358-0e6c-4ea9-8950-52f5bdda287e
# ╠═b0e6abc4-7581-4a7b-9160-792d2824a1cf
# ╟─6be8652b-525c-4efa-bcc7-b72285ea2d98
# ╠═ab7f07f6-20f3-465f-ae46-e694100541b8
# ╠═1f9b91e3-28f4-4e06-8a3d-70c84432007b
# ╠═ce876450-eea7-4088-9aca-ed3ee462cb09
# ╟─0f7709e7-fec9-4176-bf58-96f2d891e8b4
# ╟─73e61d6c-0762-41d9-a0c1-c2a64926a40a
# ╟─b48af18e-c65b-4ddd-b6d1-6b25e4e4bb80
# ╠═b8bb7d0b-21c4-4746-9dc7-dc47f76a0736
# ╟─009c780c-d1cc-4421-86c7-74e5a99646f9
# ╠═4e208b64-6b88-4430-9bc3-7198c54b00ff
# ╠═8b2f116a-09e6-4bbf-9a83-b7c38fe34680
# ╟─c3c91935-60f7-4c4a-b43c-973566d58367
# ╠═41265945-67e6-4e34-8604-20faa6f475b7
# ╟─8937b497-99f0-4d8f-998a-e2aad0d5aeb9
# ╠═f64ab474-86f6-4048-b013-3c8969d2e125
# ╠═d21b8a45-7129-4038-82ac-250436ff8583
# ╠═d2ef7881-bd1b-4bd3-aee4-025c7da5ba71
# ╟─9d5e0e5d-9117-4308-b133-20e3d9e0894e
# ╠═bd576cd6-ee70-43b3-ad05-9307d6eeb05d
# ╠═6f81e36d-7892-42ec-aeb6-91fca5b6399a
# ╟─137f3b99-d031-4df4-868e-9b5e4d342a83
# ╟─1b63714a-1e7f-40f0-9f2f-5a03744a44d5
# ╟─934ec8b3-a0c0-4129-94f1-f3f04c619365
# ╟─84e8b9a8-eb55-46e9-81ab-df248d98a800
# ╟─0d96be04-f200-43c1-8374-8a1ec85be47f
# ╟─b02d908b-2882-42b4-803c-fa2b9ca3833f
# ╠═63ad906a-0fce-4b3e-a5ef-3001c0a1c315
# ╟─cde89859-7504-4730-afb1-32919d4c5dd2
# ╟─13bdb670-8235-47f0-a2fe-ff9b6d12afeb
# ╠═5e02082e-0b07-404b-a5f9-8f86155f4381
# ╟─0076b5f3-b038-4d8a-8218-dfac8fb8f731
# ╠═c8b15905-8529-4da0-9944-313e271e0665
# ╟─5374fdc7-efac-4044-8757-660120992877
# ╠═56cfee0e-4eec-49e8-9c8d-6ef3552f7c23
# ╟─5ff2388d-bda2-42bb-957d-6605a1e343ac
# ╠═36e7afa4-99ed-4dde-9129-e696759cbf95
# ╟─899b256e-e392-4785-8097-6bf7fad79076
# ╠═cd79bd02-2d61-4eaf-9491-796ff041bdba
# ╟─6c1decec-0f96-4d19-b969-29a4344fdbb7
# ╟─4413399f-1a77-44f1-8e71-7fdf50233ee8
# ╟─e4a01539-ef46-4e4c-86c4-5c37f34a0336
# ╟─befcfdd1-1f57-45ab-8f96-ae4d01895eea
# ╟─0233c113-c3ef-456a-aa50-17eb754608ac
# ╟─e51a1e89-1675-4b1f-9782-7951b60e1658
# ╟─413d5e90-0b44-4e3f-8ab8-15ac1538c1c7
# ╟─262eca5b-d55e-4428-a11a-0614751dd2d6
# ╟─24dd589f-3275-4924-ba9e-9ad5c7fd075a
# ╠═814837a2-761b-4233-81f7-ded2907cb07f
# ╟─e1fd4f96-a465-46d8-8c78-cde7c5325e6f
# ╟─bebb6206-9077-40b4-8bf6-63cf5f23d308
# ╠═7237c2e9-c0f3-4cd8-ad1f-6d5cbff35354
# ╠═a9fdeb7d-e4f8-4320-8d32-792bf405a7aa
# ╠═2e96d48d-15d9-40cc-8c1d-d8d03e755c2a
# ╠═e6a62b4f-9919-43b9-ae74-3cb3c1110829
# ╟─ccd01cf4-e25d-4c4d-99dd-63a12de4e65d
# ╟─6b997532-b424-4df6-afdd-ded41bfba8a8
# ╟─aaf81e2f-5e6e-4f5d-8950-831ae85f2476
# ╟─37b95aee-ded0-4698-ae1b-4f236ca4cffb
# ╠═1d099205-9d84-4ecc-a2af-81910326a4ba
# ╠═f17500cc-6e5f-442d-9b22-2dc406367ad6
# ╠═320a4138-1d2d-4039-adaa-1aeb7f1929cf
# ╟─db913617-2b81-4ea7-a907-033dd31bda5a
# ╠═68fff309-a3ec-4ebc-bcfe-bafeac6ef3f5
# ╟─f4bb1d69-1f06-4769-ac6e-081ddaa437d7
# ╟─8481465b-171b-4374-8ec9-d7c19bd23d81
# ╠═652c5513-8263-4211-833f-aa8180711d68
# ╟─4b0cd6fc-50ec-4e12-82e0-30fe2ceb98ce
# ╟─9fdc92f6-158f-4e62-82e0-b6c5ce96d9a8
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
# ╠═6c1c4615-c98d-4fd4-8f5b-c0fe2702204f
# ╟─dddf152e-b1a3-470d-93b7-d213544f0215
# ╟─b2cbf4dd-f203-4028-90b1-57d72d21c3c5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
