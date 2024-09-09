### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ fcf3b117-bb14-4849-818c-a0dfc4de15aa
begin
	using PlutoUI
	using Plots
	using PlutoTeachingTools
	using HypertextLiteral
end

# ╔═╡ 3a59c0c6-5496-4faf-ae89-3c6333d02abb
using BFloat16s

# ╔═╡ baa1625c-7faa-4341-a0b6-bb5ad3fdb67a
using LinearAlgebra

# ╔═╡ 3e0b819d-9dbe-45f7-9e92-c58578b85d21
using AccurateArithmetic

# ╔═╡ a06746c0-c3c3-46cf-985f-2d6873f1f141
using BenchmarkTools

# ╔═╡ 67d5a7a0-9866-41bf-a242-70d9ae44a444
begin
	using DoubleFloats
	using GenericLinearAlgebra
end

# ╔═╡ 8da0c6e6-ba1a-41ae-99f1-060a6cb53bbc
using IntervalArithmetic

# ╔═╡ a3d727fa-fce6-43c7-ac54-51b5bcc0cde0
md"""
# Tackling errors due to floating-point arithmetic
"""

# ╔═╡ b1c739ed-957b-4a57-b1aa-6981bb3f76c4
md"""
In this lecture we will discuss some properties of floating-point numbers and floating-point arithmetic. How are numbers represented on the computer? What can we say about the error introduced by rounding? Are there computational tools to estimate or reduce floating-point error in our computation?

To illustrate our points we will dwell from a number of simple problems, not all related to eigenvalue computations per se.

This material is partly taken from the first 3 chapters of the excellent book *Accuracy and Stability of Numerical Algorithms* by Nicholas J. Higham and a recent Acta Numerica paper titled [Floating-point arithmetic](https://doi.org/10.1017/S0962492922000101) by Boldo, Jeannerod, Melquiond and Muller. 

A good resource is also the [Julia documentation](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/) on this topic.
"""

# ╔═╡ 3c06f892-29d9-4df6-bd85-b22a19d12350
TableOfContents()

# ╔═╡ 77c9a481-b314-42f3-944e-b6422b894dd0
md"""
## Basic concepts
### Some notation and definitions

- Our main measure of error are the **absolute error**
  ```math
  e_\text{abs} = |x - \hat{x}|
  ```
  and the **relative error**
  ```math
  e_\text{rel}(\hat{x}) = \frac{|x - \hat{x}|}{|x|}.
  ```
  For scientific computations the order of magnitude of an expression can vary a lot and thus typically the relative error is of main importance. Notice, that an equivalent way to express the relative error is
  ```math
  e_\text{rel}(\hat{x}) = |δ| \qquad \text{where} \qquad \hat{x} = x (1 + δ)
  ```
- Notice that the relative error is related to the number of correct **significant digits**, but generally considering the relative error instead of the notion of significant digits is clearer. In the following examples, $x$ and $y$ agree with their respective approximation to 3 significant digits, but the relative error varies largely:
"""

# ╔═╡ 0e5d0819-d101-4492-8bb0-73d5302ce461
let
	x     = 1.00000
	xhat  = 1.00499
	E_rel = abs(x - xhat) / abs(x)
end

# ╔═╡ ac79b644-449d-4684-bc14-f06ca24ee160
let
	y     = 9.99999
	yhat  = 9.99899
	E_rel = abs(y - yhat) / abs(y)
end

# ╔═╡ fd0a3bea-2178-4df5-96f4-b24a7e46b512
md"""
- The function $fl(\,\cdot\,)$ denotes evaluation of an expression using floating-point arithmetic. For basic arithmetic operations, such as $+$, $-$, $\ast$, $/$ we will assume
  ```math
  fl(x\, \text{op}\, y) = (x\, \text{op}\, y) (1 + δ) \qquad |δ| < u
  ```
  where $u$ is the **unit roundoff** (or machine precision). In other words we assume that $u$ is an upper bound to the relative error of all basic operations. This is the so-called **standard model** of floating-point arithmetic.
- **Computed quantities** will be denoted by a hat, i.e. $\hat{x}$ is the floating-point approximation to $x$.
- Notice the difference between **precision** and **accuracy**: Whereas *accuracy* refers to the *relative error of an approximate computation* versus the exact result, *precision* refers to the accuracy with which basic operations are performed, which itself is bounded by the unit roundoff $u$ as discussed above. Most importantly precision is not an upper bound to accuracy. We will discuss techniques where computations are done in double precision arithmetic, but more accuracy than double precision is achieved.
"""

# ╔═╡ cba52909-94af-4185-a990-36583976643c
md"""
### Forward and backward error

Assume we have a computation $y = f(x)$, where $f$ is a function $\mathbb{R} \to \mathbb{R}$. Our computation yields a result $\hat{y}$. How should we measure the "quality" of $\hat{y}$? Given that we represent both our computation and our answer in a floating-point number system of precision $u$ we are probably happy with an answer $e_\text{rel}(\hat{y}) \leq u$. However, this is hardly ever feasible.

A slightly different point of view is to take a look on the input data $x$. Given an approximate answer $\hat{y}$, for what range of input data do we actually obtain this solution? Or in other words for which $\Delta x$ do we have $\hat{y} = f(x + \Delta x)$. Notice that in general there are many such answers $\Delta x$, and we are typically interested in obtaining bounds on the largest $|\Delta x|$, which achieves $\hat{y} = f(x + \Delta x)$. This process is called **backwards error analysis** and $|\Delta x|$ is referred to as **backwards error**. To distinguish more clearly the errors in $\hat{y}$, we refer to these as **forward errors**.

An algorithm is called **backwards stable** if for any $x$ it produces a computed $\hat{y}$ with a small backwards error, that is $\hat{y} = f(x + \Delta x)$ is only achieved for a small $\Delta x$. What small means depends on context.

Let us consider an example to illustrate these concepts further. We consider the basic operation
```math
	y = f(x) = x - 1.0
```
"""

# ╔═╡ 21d26883-5dcb-409b-bbee-64c6155b4d53
f(x) = x - 1.0

# ╔═╡ 2f7c1a64-bee2-470f-a2f6-887fad89994f
md"""
Based on our assumption of the standard model, the operation $x-1.0$ will yield the exact answer for any perturbed data $x (1 + δ)$ with $|δ| < u$. By assumption this operation is backwards stable. We can test this numerically for example for the case $x = 1 + 10^{-16}$. In double precision
"""

# ╔═╡ 408f1b1c-dd98-4e74-92ad-754b6122e525
md"""Since the result is $\hat{y} = 0.0$ (instead of the exact $10^{-16}$) $f$ has a forward error of $10^{-16}$ in this example (which is less than $u$). Now"""

# ╔═╡ 4751e59e-98fd-4547-bcd1-43ae740e2e6f
f(1.0)

# ╔═╡ 6003c907-3722-4fd0-9eae-ab96ec88d902
md"thus taking a $\Delta x = -10^{-16}$ yields the same as $f(x)$ for our $x=1+10^{-16}$. 
Further for $\Delta x = 2 \cdot 10^{-16}$"

# ╔═╡ fa2932d3-b09b-4ff2-8d59-295c4fe8cd5d
md"""
Therefore the backward error $|\Delta x|$ is between $10^{-16}$ and $2 \cdot 10^{-16}$ in this example.

Following further along the idea of errors of floating-point addition, we find that the first number larger than $1$, which is representable in `Float64` is
"""

# ╔═╡ 6f7058fe-171a-46a9-9d07-ccd65fb61bcd
md"""
where `eps` gives the **unit in the last place** (ulp) of its input, which is defined as
"""

# ╔═╡ 740d9e34-71c0-46dc-a9a0-0e4537e63e69
epsdef(x) = max(x - prevfloat(x), nextfloat(x) - x)

# ╔═╡ 8d126d65-34b3-478b-a1c8-03142fa4e888
md"""
Where `prevfloat` and `nextfloat` return the previous and next representable floating-point number to a given number, e.g.
"""

# ╔═╡ 1cb60604-7d36-4f11-9683-d5a4ca821c7e
prevfloat(1.2345)

# ╔═╡ 2ec1b321-c63d-40d3-a00b-49b7a80098f8
nextfloat(-1e-2)

# ╔═╡ 293e2d97-5e6f-4d69-99f9-1a9dde7e1c30
md"""
Due to its importance, `eps(one(T))` is also called the **machine epsilon** of the floating-point type `T`. For example:
"""

# ╔═╡ c0d6e714-8c25-45eb-afe8-3fae25a174e6
eps(Float32)

# ╔═╡ 9368e779-9104-44e9-971d-95bc30ee38bc
eps(Float64)

# ╔═╡ 8b261d9b-22e3-4596-971d-30e9f51a6a0f
md"""
### Conditioning and condition number

Forward and backward error can be related to each other in most cases. The relationship between them is usually called the **conditioning**, that is the sensitivity of the solution to the perturbation of the input data for a numerical scheme.

Consider our previous example with $y = f(x)$ and $\hat{y} = f(x + \Delta x)$. Then if $f$ is differentiable we can find an expansion
```math
\hat{y} - y = f(x + \Delta x) - f(x) = f'(x) \Delta x + O\big( (\Delta x)^2 \big).
```
Dividing by $y = f(x)$ then gives
```math
\frac{\hat{y} - y}{y} = \left( \frac{x f'(x)}{f(x)}  \right) \frac{\Delta x}{x} + O\big( (\Delta x)^2 \big),
```
which relates the relative forward error to the relative backward error by means of the factor
```math
\kappa = \left| \frac{x f'(x)}{f(x)}  \right|
```
the so-called **condition number**.

When defined consistently, the following useful relationship holds between forward/backward error and condition number:
```math
\text{forward error} \lessapprox \text{condition number} \times \text{backward error}.
```
Crucially, this implies that for an ill-conditioned problem, a small backward error does *not* imply we have solved the problem to high accuracy, meaning that the reminiscent forward error may still be rather large.
In other words even if we employ a backwards-stable algorithm to compute our approximate solution $\hat{y}$ and $\Delta x$ is small our relative error in $\hat{y}$ may still be large. 

For iterative methods this implies for ill-conditioned problems: Even if successive iterates $x^{(i)}$ and $x^{(i+1)}$ are very close and our algorithm might not even be able to improve any more (e.g. $|x^{(i+1)} - x^{(i)}| < u$), still the error in the solution $\hat{y}$ we might still be large.
"""

# ╔═╡ ab77c84b-c0a3-436b-ac49-cd667983ead8
md"""
**Vector functions.**
This argument can be generalised to functions between vectors by taking vector norms to compute (relative) errors and considering the *maximal* possible ratio between forward and backward error.
For example, we want to solve the **linear system** $Ay = x$, which can be done by the solution "algorithm" $f(x) = A^{-1} x$, which clearly is a function $\mathbb{R}^n \to \mathbb{R}^n$ with $n$ the dimension.
If we are able to perturb $x$ to $x + \Delta x$ to obtain the computed $\hat{y}$,
the relative forward error is
```math
e_\text{fw} = \frac{\| \hat{y} - y \|}{\| y  \|} = \frac{\| A^{-1} (x + \Delta x)  - A^{-1} x \|}{  \| A^{-1} x \|} = \frac{\| A^{-1} \Delta x  \|}{  \| A^{-1} x \|}
```
and the relative backward error clearly
```math
e_\text{bw} = \frac{ \| x + \Delta x - x \| }{\| x \|} = \frac{\| \Delta x \|}{\| x \|},
```
Taking the maximal element of the ratio between the expressions over all possible non-zero values of $\Delta x$ and $x$ we obtain the relative condition number as
```math
\begin{align}
\kappa &= \max_{x\neq 0} \max_{\Delta x \neq 0}
\frac{e_\text{fw}}{e_\text{bw}} \\
&= \max_{x\neq 0} \max_{\Delta x \neq 0}
\frac{\| A^{-1} \Delta x  \| \| x \|  }{\| A^{-1} x \| \| \Delta x \| } \\
&= \max_{x\neq 0}\frac{\| x \|}{\| A^{-1} x \| } \max_{\Delta x \neq 0}
\frac{\| A^{-1} \Delta x  \|   }{ \| \Delta x \| } \\
&= \max_{y\neq 0}\frac{\| A y \|}{\| y \| } \max_{\Delta x \neq 0}
\frac{\| A^{-1} \Delta x  \|   }{ \| \Delta x \| } \\
&= \| A \| \, \| A^{-1} \|,
\end{align}
```
which is the familiar expression for the condition number of linear systems.

Notice: For eigenvalue problems the condition number takes a different expression. A detailed discussion is deferred to our discussion of perturbation theory in eigenvalue problems.
"""

# ╔═╡ 3edfa209-d265-4b90-9461-a2d0a2461a5f
md"""
### Cancellation

Subtractive cancellation happens when nearly equal numbers are subtracted. This is often (but not always) a bad thing. We consider the function
"""

# ╔═╡ bc5e8449-bf43-4c4e-8d12-bf035e1413f7
fᶜᵒˢ(x) = (1 - cos(x))/x^2

# ╔═╡ 102f8541-11bd-4b91-a8bd-21acb683ae20
md"""
By considering the Taylor expansion of `cos` one easily shows that
$0 \leq f^\text{cos}(x) < 1/2$ for all $x \neq 0$. However, we obtain
"""

# ╔═╡ 7f21c4b2-6ac0-4bd4-88ee-5a48fa3be773
fᶜᵒˢ(1.2e-8)

# ╔═╡ e61319a9-2a02-452e-96ac-60bc35f12a32
md"""
which is thus clearly wrong ... not even the first digit is correct. The culprit here is the computation in the numerator of this function, since
"""

# ╔═╡ 360b47ab-992c-43a0-9afa-61865a3352ab
c = cos(1.2e-8)

# ╔═╡ 7400e423-da76-4a61-b32e-bbb08c93609d
md"such that"

# ╔═╡ a7f2a54b-e977-4c6f-a24d-994019d57b4b
1 - c

# ╔═╡ a4d4ced6-b71d-4fd2-8a26-262a443bc6c3
md"""which only has 1 significant figure left.
It is worth emphasising that the subtraction `1 - c` is actually *exact*. The error has been introduced in the computation of `c` itself. In other words subtractive cancellation *brings earlier errors into prominence*.

To analyse this we consider the subtraction of two computed numbers $\hat{a}$ and $\hat{b}$ in *exact arithmetic* to form $\hat{x} = \hat{a} - \hat{b}$. Denoting by $\Delta a$ and $\Delta b$ the relative errors of $\hat{a}$ and $\hat{b}$ from previous computations, i.e.
```math
\hat{a} = a (1 + \Delta a) \qquad \hat{b} = b (1 + \Delta b)
```
then we have (using $x = a - b$):
```math
  \left| \frac{x - \hat{x}}{x} \right|
= \left| \frac{-a \Delta a + b \Delta b }{a - b} \right|
\leq \max(|\Delta a|, |\Delta b|) \frac{|a| + |b|}{|a - b|}.
```
Important to note is that the relative error bound in $\hat{x}$ is thus in particular large when $|a - b| \ll |a| + |b|$, i.e. when we have heavy cancellations in the subtraction. Note that cancellation is not *always* a bad thing. In case $a$ and $b$ are exact (i.e. error-free) the subtraction is actually exact (in the exact arithmetic we have assumed here) or its error bounded by the usual $u$ (if the standard model of FP arithmetic applies).
"""

# ╔═╡ 399a584e-5c21-4502-b89e-ea4a3871a3f8
md"""
### Rounding errors

When performing floating-point operations rounding is inevitable to bring computed values back to the valid range of the employed floating-point datatype. While rounding in the IEEE standard has been designed to be reasonably balanced, such that one may expect errors to cancel or evenly distribute (and they often are), if one is not careful a few rounding errors in an iterative procedure can lead to rather catastrophic results. But in fact rounding errors are not just bad as we will see.

**Accumulation of rounding errors.**
As an illustration consider taking finite $n$ in the definition
```math
e = \lim_{n\to\infty}\left(1 + \frac1n\right)^n
```
defined in the function
"""

# ╔═╡ bfff1a42-d18c-4bdb-bf3f-224ba5261eef
md"""
where `T` varies the employed floating-point type. We plot for a few floating-point types:
"""

# ╔═╡ 18754b1e-f8d5-43ee-a668-a13dcd19f2e0
md"""
Clearly the approximation is poor and moreover degrades completely as $n$ approaches the reciprocal of machine precision.

The reason is that in a binary floating-point system fractions like $1/10^n$ are non-terminating numbers, e.g.
"""

# ╔═╡ 188c06ba-2201-40ec-b14f-042da0c12444
bitstring(1.0 / 10.0)  # Shows the binary representation of a variable as a string

# ╔═╡ c19abeaa-7e9f-4a1f-892a-43ef9c19ee1a
md"""
as a result when $1 + 1/n$ is formed for $n$ being large powers of $10$, only a few significant digits of $1/n$ are retained in the sum. The subsequent exponentiation must thus produce an inaccurate approximation to $e$, even if done exactly.

An important takeaway of this example is that a **single rounding error** was enough to cause this behaviour. No accumulation and successive rounding is involved.
"""

# ╔═╡ e8d93fd9-e32b-476f-9468-e637089a625f
md"""
**Rounding errors can be beneficial.**
It is important to keep in mind that rounding errors (or in general floating-point errors) are not just problematic for numerical algorithms. They can also be beneficial. To see this we return to our good old friend the Power method. An implementation we used a few times is as follows:
"""

# ╔═╡ 39e0ed4f-c812-47dd-b55d-3c4b5e07b2f5
function power_method(A, u=randn(eltype(A), size(A, 2));
                      tol=1e-8, maxiter=100, verbose=true)
	norm_Δu = NaN
	for i in 1:maxiter
		u_prev = u
		u = A * u
		normalize!(u)
		norm_Δu = min(norm(u - u_prev), norm(-u - u_prev))
		norm_Δu < tol && break
		verbose && println("$i   $norm_Δu")
	end
	μ = dot(u, A, u)
	norm_Δu ≥ tol && verbose && @warn "Power not converged $norm_Δu"
	(; μ, u)
end

# ╔═╡ 3f83b876-54a6-4900-be2a-d63609ed1947
md"""We consider the matrix"""

# ╔═╡ b822af3d-36e9-4985-a85f-eb34f96a660c
A = [0.4 -0.6 0.2; -0.3 0.7 -0.4; -0.1 -0.4 0.5]

# ╔═╡ 2361d15d-01c9-4e5d-98c4-ddec8ec2a19a
md"""This matrix has eigenpairs"""

# ╔═╡ 5d571dcf-aef0-448e-8eaf-704d4b02b124
λ, X = eigen(A)

# ╔═╡ ef685429-5b8d-43ea-9480-eb2791345a9a
md"""With in particular the eigenvector (corresponding to the zero eigenvalue)"""

# ╔═╡ fd98ab89-275c-49b4-9a46-17880f2b6540
vstart = [1.0; 1.0; 1.0]

# ╔═╡ 53dec758-3c39-4d2c-b9f4-58dba5b45c68
md"""
being orthogonal to the eigenvector of the largest eigenvalue of $A$. In exact arithmetic the power method would thus converge in one step --- namely to the eigenvalue zero and the normalised starting vector.

However, we observe
"""

# ╔═╡ 3edb0d0f-4b75-4fa5-a5e3-80129abe2220
begin
	res = power_method(A, vstart)
	res.μ
end

# ╔═╡ 1b4bc39a-eed6-451a-b16c-fb5025ff7012
md"""
... in other words a nice convergence to the correct largest eigenvalue and corresponding eigenvector.

Interestingly the reason why this works is exactly the reason why our previous example failed, namely the fact that decimal fractions (like $-0.1$ or $0.2$) are not exactly representable in binary floating-point formats. We actually diagonalise $A+ΔA$, where $ΔA$ is a tiny perturbation. For $A+ΔA$ `vstart` is not an exact eigenvector and thus the overlap with its largest eigenvector is non-zero, and the iteration can proceed. Note, that the final eigenpair we obtain is thus similarly an approximate eigenpair to $A+ΔA$ and not $A$, but in fact due to $ΔA$ being small the eigenpairs of $A+ΔA$ are very good approximations to those of $A$.
"""

# ╔═╡ 5c081511-7391-4d83-a1cf-ff83342b856e
md"""
**Rounding errors are not random.** 
This simple fact is indeed exploited in many stable algorithms, which sometimes even heavily rely on the cancellation of rounding errors. While there are attempts to construct approximate statistical models to characterise rounding errors in computations, there are a number of illustrative examples, where rounding errors follow clear patterns. One example (due to W. Kahan) is the rational function
"""

# ╔═╡ 96e1effa-d45e-4b59-b019-93b08b12dca3
r(x) = (622 - x * (751 - x * (324 - x * (59 - 4x)))) /
       (112 - x * (151 - x * ( 72 - x * (14 -  x))))

# ╔═╡ ec14869d-e277-4628-b1ea-3823db0a2a0f
md"""which is expressed in a form how a ratio of quartic polynomials would be evaluated using Horner's method. Horner's method is a very efficient and thus frequently employed way to evaluate rational functions. For example the evaluation of many standard functions such as `sin`, `cos`, `erf`, etc. are under the hood evaluated using Horner's method.

We evaluate `r(x)` for the numbers"""

# ╔═╡ e19d8e1b-7e75-4e76-baa6-4d2c5bc15496
md"""A striking pattern due to the rounding can be observed."""

# ╔═╡ 8c32c50a-82aa-4750-b91f-e393e2eab11e
md"""
### Increasing the precision to check accuracy

A common technique for estimating the accuracy of a computed result is to redo the computation with increased precision and check if the presumably more accurate answer agrees. Many algorithms --- especially if the main source of error is rounding --- do possess an error bound that decreases when the precision is increased. This implies that at infinite precision, the exact answer is obtained. However, there is generally no guarantee that increasing the precision *always* improves the answer.

To illustrate this point numerically, we exploit that the `BigFloat` type in julia actually allows to dynamically modify the number of bits that are employed for the underlying computation:
"""

# ╔═╡ 31be711b-4250-49e8-b9d5-3cf4db45108b
md"""
We observe that over the full range between 22 and 31 bits the error more or less stagnates, only to decrease later on.

Thus: If **no significant change** in the result is obtained when repeating a calculation **at higher precision**, this is **not a sign we got the right answer**!


"""

# ╔═╡ 97db876f-2aca-45fa-8b05-229ac1139a80
md"""
### Common misconceptions about floating-point errors

In this first part of the lecture we provided a few examples illustrating the difficulties of floating-point operations. We want to pause for a summary of common misconceptions, some of which we disproved by numerical examples shown above:

1. **Cancellation in the subtraction of numbers is always a bad thing.** No: It sometimes leads to convergence to the best possible answer.
2. **Rounding errors can only overwhelm a calculation if a vast number of them accumulate.** No: Sometimes a single rounding error is sufficient.
3. **Increasing the precision always improves the accuracy of the final answer.** No: Stagnation or even deterioration of the answer can occur. 
4. **The final answer of an algorithm is always less accurate than its intermediate quantities.** No: Stable algorithms are often even design around the cancellation of rounding errors, for example.
5. **Rounding errors can only hinder, not aid the success of a computation.** No: See the Power method example.
"""

# ╔═╡ 75bed4a9-54f5-44c1-bcc3-74d5fb11c88a
md"""
## Floating-point number systems
Overall our discussion above (e.g. the aspects of rounding errors) shows, how tricky it is to design a good floating point system, keeping all the beneficial aspects of finite-precision arithmetic, but reducing the impact of rounding, cancellation with respect to the final accuracy obtained.
"""

# ╔═╡ 1d5c2471-1343-4706-9894-293da149ae74
md"""
**Structure of floating-point numbers.**
A floating point number system $F \subset \mathbb{R}$ has elements of the form
```math
y = \pm m \, \beta^{e-t}.
```
It is thus characterised by four integer parameters:
- the **base** (or radix) $\beta$ (typically binary, i.e. $β = 2$),
- the **precision** $t$ and
- the **exponent range** $e_\text{min} \leq e \leq e_\text{max}$.

In the definition of $y$ the integer $m$ is the **significand** (or mantissa), which satisfies $0 \leq m \leq β^t - 1$.

To ensure a *unique* representation for each non-zero $y \in F$, we assume $m \geq β^{t-1}$. Such numbers are called **normalised**. The $0$ is a special case, which does not have a normalised representation.

**Example:** To illustrate some characteristics of floating-point numbers, we consider an example system with $\beta = 2$, $t = 3$, $e_\text{min} = -1$ and $e_\text{max} = 3$. Considering only the positive numbers, the following numbers can be represented:
"""

# ╔═╡ 71e0fa21-4cfc-45a2-b6d4-a703a892c15d
md"Plot limits: $(@bind plot_xlim PlutoUI.Slider([0.5, 1.0, 2.0, 4.0, 8.0]; default=8.0))"

# ╔═╡ a9279bfa-7014-499d-9a64-0f6550a665f8
md"""
We can read off from the plot:
- The **range** of non-zero floating-point numbers is $β^{e_\text{min} - 1} \leq |y| \leq β^{e_\text{max}} (1 - β^{-t})$.
- The spacing of floating-point numbers jumps by a factor of $2$ each power of $2$ and thus is clearly not even. This spacing is characterised by the **machine epsilon** (distance to next larger floating-point number at $1.0$), which is $\varepsilon_M = β^{1-t}$. $\varepsilon_M$ is the spacing between numbers in the range $[1.0, β]$, whereas the spacing in the interval $[1/β, 1]$ is $β^{-t} = \varepsilon_M/β$.

In fact one can prove *(Higham, Lemma 2.1)*:

> The spacing between a normalised floating-point number $x$ and an adjacent normalised floating-point number is at least $β^{-1} \varepsilon_M |x|$ and at most $\varepsilon_M |x|$.

The number system $F$ can be extended by including so-called **subnormal numbers**, that is numbers of the form
```math
y = \pm m \, β^{e_\text{min} - t}, \qquad 0 < m < β^{t-1},
```
i.e. which have the minimum exponent and are *not* normalised. These numbers have fewer digits of precision, but are useful to fill the gap between the smallest positive normalised number $\lambda = β^{e_\text{min} - 1}$ and $0$. Note, that the smallest positive subnormal number is $\lambda \varepsilon_M = β^{e_\text{min} - t}$.
Adding these to our example above we obtain (note the different x-axis scale):
"""

# ╔═╡ 0be06f95-8464-479e-8b2c-3840669354b8
md"""

**Rounding.** Having the structure of the floating-point numbers in place, we now consider **rounding**, i.e. the transformation $x \mapsto fl(x)$, which maps $\mathbb{R}$ to numbers of the form $\pm m \, \beta^{e-t}$ (which could be subnormal). Rounding is monotonic, i.e. $x \leq y$ implies $fl(x) \leq fl(y)$. We say that $fl(x)$ **overflows** if $|fl(x)| > \max\{|y| \, | \, y\in F\}$ and we say it **underflows** if $0 < |fl(x)| < \min\{|y| \, | \, 0\neq y\in F \}$.

In practice different ways of rounding are employed, see the Julia docstring of the `RoundingMode` structure. If the rounding mode does not matter we use $fl$ to denote rounding, if it does we employ the more specific notation:
  - `RoundNearest` in Julia or $RN(\,\cdot\,)$: The IEEE default rounding mode, which rounds to the closest representable number in $F$ with ties being round to the value, which has a `0` bit in the least significant bit of $m$.
  - `RoundUp` in Julia or $RU(\,\cdot\,)$: Always round up to the next largest representable number.
  - `RoundDown` in Julia or $RD(\,\cdot\,)$: Always round down to the next smallest representable number.

With this setup one can prove rigorously that every real number that lies within the range of $F$ can be approximated by $F$ with an error no larger than the **unit roundoff** $u = \frac12 β^{1-t} = \frac12 \varepsilon_M$, namely *(Higham, Theorem 2.2 and Theorem 2.3)*:

> If $x \in \mathbb{R}$ lies within the range of $F$ then
> ```math
> \begin{aligned}
> fl(x) &= x (1 + δ), &&|δ| < u \\
> fl(x) &= \frac{x}{1 + \tilde{δ}}. && |\tilde{δ}| \leq u
> \end{aligned}
> ```

This key result and the assumption of the *standard model* (see Basic notation subsection above) are the main ingredients underlying the numerical analysis of floating-point arithmetic.
"""

# ╔═╡ 1dcdb5e6-a4a0-40a9-9cf0-3072bad872f4
md"""
### IEEE 754 and quasi-standard formats

Until the mid-1980s floating-point numbers where not standardised. This meant that pretty much all major manufacturers adopted different solutions for representing numbers. Details like the base used to represent numbers, the ways of rounding computations to representable numbers or how to handle exceptional cases such as forbidden operations, overflow, underflow were *different* on every machine. Imagine writing portable software for this mess .... Clarity finally brought the IEEE Standard 754, which was originally published in 1985 and has been revised a few times since.

The standard notably defines not only the floating-point formats and some standard rounding modes, but also floating-point exceptions and their handling as well as conversions between integers and floats or floats and other formats.

The main focus of the standard has been on *high-precision* computation. However, in more recent applications such as statistical learning the reduced memory requirement of lower-precision floating-point formats is nowadays a key ingredient to obtain peak performance. This has lead to the introduction of a few **non-standard $16$-bit formats**, which have become widespread and *de facto* standards. The main differences of these formats compared to IEEE is their tradeoff between the number of bits employed for representing the significand and the exponent.

Key parameters of these formats are summarised in the following table:

Type  | Size  | Significand bits $t$  | Exponent bits | $e_\text{min}$ | $e_\text{max}$ | Unit roundoff $u$
:---- | :---- | :--: | :------: | -------------: | -------------: | :-------------
`Bfloat16` | 16bit | 7+1  | 8        | -127           | +126           | $2^{-8} \simeq 0.004$
`Float16`  | 16bit | 10+1 | 5        | -15            | +14            | $2^{-11} \simeq 5\cdot 10^{-4}$
`Float32`  | 32bit | 23+1 | 8        | -127           | +126           | $2^{-24} \simeq 6\cdot 10^{-8}$
`Float64`  | 64bit | 52+1 | 11       | -1023          | +1022          | $2^{-53} \simeq 1\cdot 10^{-16}$

All formats discussed here are binary (i.e. $β = 2$) and further enforce that all numbers are normalised or subnormal. As a result the most significant bit of the significand is always $1$ and thus not stored. This is indicated by the $+1$ in the table above.

Note, that in Julia all aforementioned floating-point types are directly supported. Only `Bfloat16` requires an additional package:
"""

# ╔═╡ f2752b0d-b1d3-4f3d-bb0c-192177a4f997
u_Bfloat16 = eps(BFloat16) / 2  # Compute the u for BFloat16

# ╔═╡ 0fb4ef70-e529-41bc-a00d-3095d21d8ea2
md"""
**IEEE arithmetic is closed:** The standard requires every operation on an IEEE number to produce a result, which is an IEEE number --- whether this makes mathematical sense or not. Exceptional operations raise a signal, which by default is handled by setting a flag (e.g. in form of a special bit pattern). 
For example:
"""

# ╔═╡ a7f9953b-0a14-4d6c-bed6-63a826d53a08
floatmax(Float64) * 10.0   # Overflow

# ╔═╡ 2c625ae1-0275-4548-b8c1-a26a8ba50794
10.0 / -0.0  # Divide by zero

# ╔═╡ 64b720d1-c38d-485e-9243-c501c9155636
0.0 / 0.0    # Invalid operation

# ╔═╡ 2901a2eb-fc06-42e8-90fd-b3932529bc60
Inf / -Inf   # Invalid operation

# ╔═╡ 1c4483d3-37f4-4043-93f6-a3ba09bcac36
md"""Note that `NaN` is **not a number** and thus compares as false to everything, including itself:"""

# ╔═╡ 889f1186-de4e-4060-91b1-bc2511c7828e
NaN == NaN

# ╔═╡ 6ffe44f7-e507-4bd2-bf30-b06ccb0b2027
md"Notably the zero can be signed to keep sign information in case of underflow:"

# ╔═╡ 35466aab-f3ba-4953-8fae-d545bc218ccd
floatmin(Float64) / -floatmax(Float64)

# ╔═╡ a8da2081-9fad-4be6-b486-c5cb6316d152
md"Furthermore sensible operations involving zero and infinity are defined:"

# ╔═╡ 1c2c9d4e-259f-4894-ad5b-87b176716829
-0.0 / -Inf

# ╔═╡ c1113c4b-55aa-4f0c-a602-14279d0702f0
md"""
## A taste of rounding error analysis

With the floating-point number system set up, we now want to perform rounding error analysis on a single standard application, namely the computation of an **inner product** $s_n = x^T y$ between two vectors $x, y \in \mathbb{R}^n$. This is a key building block of linear algebra and underlying many operations, including matrix-vector products and orthogonalisation ... i.e. our ingredients for the Power method.

We thus focus in detail on the evaluation of $s_n = x_1 y_1 + \cdots + x_n y_n$, where we take the evaluation order in the sum from left to right. Using the standard model of FP, we have
```math
\begin{aligned}
\hat{s}_1 &= fl(x_1 y_x) = x_1 y_1 (x + δ_1), \\[0.5em]
\hat{s}_2 &= fl(\hat{s}_1 + x_2 y_2) \\
          &= \big( \hat{s}_1 + x_2 y_2 (1 + δ_2) \big) (1 + δ_3) \\
          &= \big( x_1 y_1 (x + δ_1) + x_2 y_2 (1 + δ_2) \big) (1 + δ_3) \\
          &= x_1 y_1 (x + δ_1) (1 + δ_3) + x_2 y_2 (1 + δ_2) (1 + δ_3), \\[0.5em]
\hat{s}_3 &= fl(\hat{s}_2 + x_3 y_3) \\
          &= \big( \hat{s}_2 + x_3 y_3 (1 + δ_4) \big) (1 + δ_5) \\
          &= \big(   x_1 y_1 (x + δ_1) (1 + δ_3)
                   + x_2 y_2 (1 + δ_2) (1 + δ_3)
                   + x_3 y_3 (1 + δ_4) \big) (1 + δ_5) \\
          &=   x_1 y_1 (x + δ_1) (1 + δ_3) (1 + δ_5)
             + x_2 y_2 (1 + δ_2) (1 + δ_3) (1 + δ_5)
             + x_3 y_3 (1 + δ_4) (1 + δ_5)
\end{aligned}
```
where $|δ_i| \leq u$. The pattern becomes clear: Assuming for simplicity $1 + δ_i = 1 + δ$ then
```math
\hat{s}_n = x_1 y_1 (x \pm δ)^n + x_2 y_2 (1 \pm δ)^n
          + x_3 y_3 (1 \pm δ)^{n-1} + \cdots + x_n y_n (1 \pm δ)^2.
```
"""

# ╔═╡ 474ef698-cb16-44a3-b7fd-5a01005278e9
md"""
To simplify this we employ the following *(Higham, Lemma 3.1)*

> **Lemma.** If $|δ_i| \leq u$ and $ρ_i = \pm 1$ for $1 \leq i \leq n$ and $n u < 1$, then
> ```math
> \prod_{i=1}^n (1 + δ_i)^{ρ_i} = 1 + θ_n,
> ```
> where
> ```math
> |θ_n| \leq \underbrace{\frac{nu}{1 - nu}}_{=\gamma_n}.
> ```

which results in

```math
\hat{s}_n = x_1 y_1 (1 + θ_n) + x_2 y_2 (1 + θ'_n)
          + x_3 y_3 (1 + θ_{n-1}) + \cdots + x_n y_n (1 + θ_2).
```

This result depends on the order of summation (since a smaller backward error is associated to $x_n y_n$, which we sum last). It is easy to see, that a less tight result is (in vector notation)
```math
fl(x^T y) = (x + Δx)^T y = x^T (y + Δy), \qquad |Δx| \leq γ_n |x|, \qquad |Δy| \leq γ_n |y|, \qquad \text{$(\ddagger)$},
```
which is notably independent on the summation order. The operations $|\,\cdot\,|$ and $\leq$ are to be understood element-wise.
"""

# ╔═╡ b85874a4-d9ce-4131-9b1f-beaf79659b82
md"""
This is the *backwards error result*, which can be interpreted as follows: The computed inner product is the exact one for a perturbed set of data $y + Δy$, where each component-wise perturbation is bounded by $γ_n = \frac{nu}{1 - nu}$, which is tiny (e.g. around $10^{-12}$ for `Float64` and $n = 10\,000$). Alternatively one could leave $y$ fixed and consider similar perturbations $x + Δx$.

A *forward error* can be easily obtained from $(\ddagger)$ as
```math
\left|x^T y - fl(x^T y)\right| \leq γ_n \sum_{i=1}^n |x_i y_i| = γ_n |x|^T |y|.
```
"""

# ╔═╡ f02e854d-50bb-43dc-8b4d-b85df11dc739
md"""
**Proof of the Lemma.** We consider the case $ρ_i = 1$, the case $ρ_i = -1$ proceeds similarly. Thus, we are asked to prove that if $|δ_i| \leq u$ for $1 \leq i \leq n$ and $n u < 1$, then
```math
\prod_{i=1}^n (1 + δ_i) = 1 + θ_n,
\qquad \text{where} \quad |θ_n| \leq \frac{nu}{1 - nu}.
```
For this we proceed by induction over $n$. First we note that
```math
|δ_1| \leq u \leq \frac{u}{1 - u} \qquad \text{since $u < 1$},
```
therefore for $n = 1$ a valid $\theta_1$ is simply $δ_1$. Now we consider $n > 1$ and note
```math
\prod_{i=1}^n (1 + δ_i) = (1 + δ_n) (1 + θ_{n-1}) = 1 + θ_n
```
Therefore, $θ_n = δ_n + (1 + δ_n) θ_{n-1}$ and
```math
\begin{aligned}
|θ_n| &\leq u + (1+u) \frac{(n-1) u}{1 - (n-1) u} \\
&= \frac{u \big(1 - (n-1)u\big) + (1+u)(n-1) u}{1 - (n-1) u} \\
&= \frac{nu}{1- (n-1)u}\\
&\leq \frac{nu}{1- nu}.
\end{aligned}
```
"""

# ╔═╡ ffb4717a-5059-4381-9de9-5d9e41a19183
md"""
## Taming floating-point errors

Having discussed a number of things that can go wrong with floating-point computations, we now want to give some techniques that aid in detecting or fixing problems. For more details see Section 4 for of *Boldo et al.*.

### Fused multiply-add (FMA)

One of the more recent additions to IEEE-754 is the requirement for vendors to ship a fused multiply-add (FMA) with correct IEEE-conform rounding. This operation can be defined as
```math
\text{FMA}(a, b, c) = ab + c,
```
but where notably no intermediate rounding takes place after the multiplication of $a$ with $b$. Since this sequence of operations is a key component a number of standard primitives (dot products, multiplication of complex numbers, polynomial evaluations, ...), the ability to avoid intermediate rounding can be exploited to improve accuracy.
"""

# ╔═╡ 8f9e4905-5c39-4ccc-b0b5-af58e86df04c
md"""A standard 2x2 determinant function is:"""

# ╔═╡ 1343eaa4-d899-4ccd-8fbc-0f813e3a7ba6
function determinant_2x2(M)
	(M[1, 1] * M[2, 2]) - (M[1, 2] * M[2, 1])
end

# ╔═╡ da63cabb-5f02-467f-9a30-87c4d4e37288
md"An alternative by Kahan using FMAs:"

# ╔═╡ 3491e075-aac0-418d-8c01-dfee9d19bfd8
@bind logε PlutoUI.Slider(-16:1:0; show_value=true, default=-8)

# ╔═╡ 028b6e6f-e514-43ed-85ca-884b77000fb4
md"""Kahan's FMA-based determinant achieves high relative accuracy throughout the range of $ε$ in `Float64` whereas the naive algorithm deteriorates for small $ε$. Similar for `Float32` the quality is overall better."""

# ╔═╡ 30a193eb-9f22-4509-83b8-ea263884c021
md"""
### Error-free transformations and compensated algorithms

A perhaps surprising result, which was realised remarkably early in the design of floating-point arithmetic is *(Property 2.11 in Boldo et al.)*

> Let $a, b \in F$ and $s = RN(a + b)$. If the FP addition of $a$ and $b$ does not overflow, then $s - (a+b) \in F$.

This result guarantees, that for the `RoundNearest` rounding mode, the rounding error can also be represented as a floating-point number. This motivates the formulation of so-called **error-free transformations**, i.e. methods where we employ some extra floating-point operations in order to not only compute our desired result, but also the associated floating-point error.

For example in the sum case, the above result assures us that we can represent the sum $a + b$ by two numbers $s + t$, where $s$ is the floating-point sum and $t$ is its error. An implementation is:
"""

# ╔═╡ 21e1d043-6b31-4632-8f8a-8148a5b18ab8
function fastTwoSum(a, b)
	s = a + b
	t = b - (s - a)
	(; s, t)
end

# ╔═╡ a4d9ccda-9088-4d52-ab4a-584654ef20e8
fastTwoSum(0.1, 0.2)

# ╔═╡ 7683961e-fa0d-4960-80f6-c1af9d4593a7
md"which indeed gives the floating-point error as the second term"

# ╔═╡ 12bfc098-7c8b-4d7a-b1e4-a2361beb70c3
md"""However, this algorithm makes one assumption, namely that the exponent of $a$ (first argument) is larger than the exponent of $b$ (second argument). This condition is important, and otherwise the error can be *very far* from the error of the FP addition. Consider the example"""

# ╔═╡ 5adb95eb-a275-4d4d-9356-6fee5c387328
fastTwoSum(1.0, 2.0^55)

# ╔═╡ 70616c6f-3f61-494c-a8fa-5de02dbe9a9b
md"where $t=0$ albeit clearly $t=1$. In the general case (without making assumptions on $a$ and $b$) we need to employ:"

# ╔═╡ 7f849f7f-2425-4ef9-8e1e-dd36a0aab126
function twoSum(a, b)
    s = a + b
    v = s - a
    t = (a - (s - v)) + (b - v)
    return (; s, t)
end

# ╔═╡ 6cc3184b-5c44-4aa4-84c9-7647872bdcd3
md"Now the answer from our previous example is correct:"

# ╔═╡ 4891a8a5-7941-4c07-8c20-cd2bf94ab827
twoSum(1.0, 2.0^55)

# ╔═╡ 854efe00-c89f-40c0-b31e-42219acf5c44
md"""
As spelled out here these algorithms might not yet seem extremely useful --- in particular since they need three (respectively six) times the number of operations compared to a standard summation.

However, it is important to note that summations involving data of very unequal sizes are not unusual in scientific applications (e.g. finite differences, Euler-type integration schemes) and can easily lead to catastrophic cancellation (i.e. cancellation where no bit is correct). For example:
"""

# ╔═╡ e6fb6742-543b-4d90-b902-2f620f393b0a
md"We employ this naive summation function on a specially crafted vector of numbers of different size, which are constructed to sum to 1.0."

# ╔═╡ 3c0cdb90-1ff6-488a-a6b4-190bf000feb9
function generate(N::Integer)  # Generate a vector that sums to one by construction
    x = randn(N) .* exp.(10 .* randn(N))
    x = [x; -x; 1.0]
    x[sortperm(rand(length(x)))]
end

# ╔═╡ b7b2d7ff-b5b4-48c2-ad56-660cb813814b
md"""Clearly quite a significant error. Now, based on the idea of `fastTwoSum` one can design so-called **compensated summation** algorithms, i.e. algorithms that try to recover some of the accuracy that can be lost in summation due to catastrophic cancellation.

A simple approach is Kahan's compensated summation method:
"""

# ╔═╡ 2d03b229-ba6e-447c-a13c-72132e5ccdc4
md"""This method is notably better, but still far away from the exact answer. But based on Kahan's idea better methods have been proposed and are available as part of the `AccurateArithmetic` Julia package, for example:
"""

# ╔═╡ 6f507158-ef74-456b-b820-faaf6a4042e1
md"""Both methods are spot-on!

Moreover, the careful and optimised implementation of these algorithms in `AccurateArithmetic` employs modern CPU instructions (such as vectorisation) to effectively hide part of the overhead of the additional summation steps that need to be performed in the compensated implementations. The result is that despite achieving significantly better accuracy, these algorithms are only about a factor of two slower than Julia's builtin `sum` function for larger problem sizes:
"""

# ╔═╡ e5269069-4235-4660-aee1-af5e149ef29d
@bind logN PlutoUI.Slider(1:0.5:6; default=3.0, show_value=true)

# ╔═╡ f29adb10-1992-4e77-ad8c-61a4877e78d5
v = generate(ceil(Int, 10^logN));

# ╔═╡ 9e2fc100-4165-43fd-87e4-65ea2bba6368
md"""
### Double-word arithmetic

We saw above that compensated algorithms like `fastTwoSum` or `twoSum` produce exact summation results $r$ in the form of unevaluated sums $r_h + r_l$, where $r_h = RN(r)$ and $r_l$ is the remainder to obtain the exact answer. If the underlying floating-point type is `Float64` this makes it easily feasible to manipulate numbers with a precision beyond 100 bits! However, since only standard floating-point units are employed, the performance is generally faster than `BigFloat` or related arbitrary-precision libraries.

This motivates the investigation of so-called **double-word arithmetic**, which is indeed based on a double-word (DW) number type that represents the unevaluated sum $x = x_h + x_l$ with $x_h = RN(x)$.

A very sketchy Julia implementation, that supports summation can be obtained in a few lines of code:
"""

# ╔═╡ 3bfbecc2-c89c-4323-aff6-c7d7600767a6
begin
	# Composite type representing the unevaluated sum
	struct DoubleWord
		high::Float64
		low::Float64
	end

	# Cast from ordinary float to DoubleWord
	DoubleWord(x::Float64) = DoubleWord(x, 0.0)

	# Function returning canonical zero in this type
	Base.zero(DoubleWord) = DoubleWord(0.0, 0.0)

	# Addition DoubleWord & DoubleWord
	function Base.:+(xw::DoubleWord, yw::DoubleWord)
		(sh, sl) = twoSum(xw.high, yw.high)
		(th, tl) = twoSum(xw.low,  yw.low)
		c = sl + th
		vh, vl = fastTwoSum(sh, c)
		w = tl + vl
		zh, zl = fastTwoSum(vh, w)
		DoubleWord(zh, zl)
	end
end

# ╔═╡ c2daafe3-672c-4b9a-ad32-e4d773ecaeec
let
	x = 1.0 + 1e-16
	y = f(x)
end

# ╔═╡ 9dee1557-ec64-4146-abd0-aaba2ad3afdc
f(1.0 + 3e-16)

# ╔═╡ 65fc5ab7-364b-4e77-924e-bd0f2d3ad44e
1 + eps(1.0)

# ╔═╡ 605a633d-2362-4ebe-8ca0-7db925e52a2e
compute_e(n::T) where {T} = (one(T) + one(T) / n)^n

# ╔═╡ 562abdb2-0a8a-4bc6-aa2c-0266de924ee6
begin
	p = plot(xaxis=:log, yaxis=:log, ylabel="absolute error",
		     xlabel="n", legend=:bottomleft)
	for T in (Float32, Float64, BigFloat)
		ns = T(10).^(1:18)
		errors = abs.(compute_e.(ns) .- ℯ)
		plot!(p, ns, errors, label=string(T))
	end
	vline!([2/eps(Float32)],  ls=:dash, label="inv. Float32 precision", c=1)
	vline!([2/eps(Float64)], ls=:dash, label="inv. Float64 precision", c=2)
	p
end

# ╔═╡ e8d4565a-9a74-473a-b8b7-87a88f933bcb
let
	ks         = (1:361)
	samples    = 1.606 .+ (ks .- 1) .* 2 .^ -52.0
	references = r.(BigFloat.(samples))

	p = plot(       ks, references, label="BigFloat")
	p = scatter!(p, ks, r.(samples), label="Float64")
end

# ╔═╡ ca3562da-7965-441d-bbaf-0b80c7edfa2d
let
	function eval_with_precision(precision::Int)
		setprecision(BigFloat, precision) do
			x = 1 / BigFloat(7)
			a = 10^(-BigFloat(8))
			b = 2 ^ BigFloat(24)
			x + a * sin(b * x)
		end
	end

	reference = eval_with_precision(256)
	precisions = 10:40
	evals = abs.(eval_with_precision.(precisions) .- reference)
	plot(precisions, evals; yaxis=:log, mark=:x, label="")
end

# ╔═╡ bbbf8e57-2ee0-44c6-9d84-d51c9fe5fd1e
let
	β    = 2.0
	t    = 3
	emin = -1
	emax = 3

	numbers = [0.0]
	append!(numbers, [m * β^(e - t) for m in β^(t-1):(β^t-1) for e in emin:emax])
	sort!(numbers)

	@show numbers

	p = scatter(numbers, zero(numbers); label="Normalised numbers", mark=:x, ylims=(-1, 1), ms=6, legend=:topright)
	vline!(p, [0.25, 0.5, 1, 2, 4], label="Factors of 2", c=3)
	xlims!(p, (-0.01 * plot_xlim, plot_xlim))
end

# ╔═╡ 75a6cd93-4609-4b93-a1f8-8e905cdc6602
let
	β    = 2.0
	t    = 3
	emin = -1
	emax = 3

	numbers = [0.0]
	append!(numbers, [m * β^(e - t) for m in β^(t-1):(β^t-1) for e in emin:emax])
	sort!(numbers)
	
	subnormals = [m * β^(emin - t) for m in 1:(β^(t-1)-1)]
	@show subnormals

	p = scatter(numbers, zero(numbers); label="Normalised numbers", mark=:x, ylims=(-1, 1), ms=6)
	scatter!(subnormals, zero(subnormals); label="Subnormal numbers", mark=:+, ms=6)
	vline!(p, [0.25, 0.5, 1, 2, 4], label="Factors of 2", c=3)
	xlims!(p, -0.1, 4.1)
end

# ╔═╡ 7108724e-5973-43bf-ab75-d4c5e6935cd1
Inf + Inf

# ╔═╡ b5bc55d1-00ca-4733-8f32-08b4528ba297
begin
	ε = 10.0^logε
	M = [π exp(1); 355/113 23225/8544+ε]
end

# ╔═╡ 39f6f711-eab3-4fa1-9061-ec7042010c63
md"""
A classic example is the computation of determinants of 2x2 matrices. We consider the matrix $M$ with $ε =$ $(ε):
"""

# ╔═╡ 39e10226-e527-4699-af33-ea64ccdcb066
md"""Use the slider to change the definition of the matrix 
$M =$ $(latexify_md(M)) and see what happens in single and double precision:"""

# ╔═╡ a8cb9060-9fd1-4006-8f65-d1882b6be082
function det_relative_error(det_value::T) where {T}
	detM_ref = det(BigFloat.(M))
	T(  abs(det_value - detM_ref) / abs(detM_ref) )
end

# ╔═╡ 1929e74c-e6b7-4ffe-b569-61248ba075c3
relerror_naive_fp32 = det_relative_error( determinant_2x2(Float32.(M)) )

# ╔═╡ 08a6f5a4-8c3f-4fa2-b40d-6575ba1597d0
relerror_naive_fp64 = det_relative_error( determinant_2x2(M) )

# ╔═╡ ba7ca52e-e293-4515-a327-af2be18e41b1
function determinant_kahan(M)
	w = M[1, 2] * M[2, 1]
	e = fma(-M[1, 2], M[2, 1],  w)
	f = fma( M[1, 1], M[2, 2], -w)
	e + f
end

# ╔═╡ 540686e3-df39-44f4-86ff-315d5e080ead
relerror_kahan_fp32 = det_relative_error( determinant_kahan(Float32.(M)) )

# ╔═╡ 5162027a-b82f-4690-af0a-a6aabeed8652
relerror_kahan_fp64 = det_relative_error( determinant_kahan(M) )

# ╔═╡ 5435303f-2175-4ddd-80a0-7b069a7f8a8d
(BigFloat(0.1) + BigFloat(0.2)) - (0.1 + 0.2)

# ╔═╡ a4e5b19f-913c-406f-afe7-5d0fb6b3df3f
function sum_naive(x)
	accu = zero(eltype(x))
	for xi in x
		accu += xi
	end
	accu
end

# ╔═╡ 95b1ec78-033e-4f03-8adb-8045848ca123
begin
	sample = generate(100)
	sum_naive(sample)
end

# ╔═╡ 73594ce6-2eb0-444d-bd66-b490d9d72d1e
sum_kbn(sample)   # Improvement over sum_kahan by Kahan, Babuska and Neumaier

# ╔═╡ d2755839-8576-4749-8e2e-34f8b6cbda6b
sum_oro(sample)   # Method by Ogita, Rump and Oishi

# ╔═╡ 39f5ca43-fa61-4f27-b610-016638a54eeb
function sum_kahan(x)
	accu  = zero(eltype(x))
	error = zero(eltype(x))
	for xi in x
		temp = accu
		y = xi + error
		accu = temp + y
		error = (temp - accu) + y
	end
	accu
end

# ╔═╡ e8476598-32e4-4c31-8d21-42835ceb36b6
sum_kahan(sample)

# ╔═╡ 8240ab3a-5b57-4d40-8b04-e8a7bcf5d5cb
let
	p = plot(xlabel="Runtime (ns)", ylabel="error", yaxis=:log, legend=:bottomleft)
	for sumfunction in (sum, sum_kbn, sum_oro)
		error = abs(1 - sumfunction(v)) .+ 1e-20  # Offset to avoid log(0.0)
		
		bench   = @benchmark $sumfunction($v)  # run benchmark
		runtime = BenchmarkTools.median(bench).time

		scatter!(p, [runtime], [error]; label=string(sumfunction))
	end
	p
end

# ╔═╡ bbd6e360-fca7-41b4-83c8-f6b03c5b6243
md"""
Of note a posteriori error bounds are available for these algorithms, which moreover have been formally verified using the proof verification system Coq.

Employing such a `DoubleWord` type directly improves our sums on our sample from the previous section --- where the `sample` vector was constructed to sum to `1`:"""

# ╔═╡ 2f6f1ee8-a85c-4246-ad27-19f31eb1ce0f
# Ordinary Float64:
sum(sample)

# ╔═╡ c300dde5-3a4c-43da-82c3-48b24a31ab04
# DoubleWord-based sum
# (Note how this seamlessly integrated into Julia's own sum function!)
sum(DoubleWord.(sample)).high

# ╔═╡ eb7a1a03-7cf1-4eb1-a7ea-b3c4ca4ad268
md"""
A fully feature-complete implementation of double-word arithmetic in Julia is provided by the `DoubleFloats` package. In combination with the `GenericLinearAlgebra` package, which provides pure-Julia implementations for all key linear-algebra operations (e.g. `eigen`, `qr`, `svd`, `cholesky` ...), this provides an excellent environment to run algorithms in increased precision for testing purposes. Notably the **runtime penalty** of using `DoubleFloats` is **usually a factor of 10 smaller** than using `BigFloat`. This makes `DoubleFloats` applicable for realistic calculations. In MatMat group we use it for example for testing our actual quantum-chemical materials simulations for floating-point stability.
"""

# ╔═╡ 95ed87cb-8a10-4b7e-8641-33277d87dab9
md"""First we try the summation on our constructed vector from above:"""

# ╔═╡ f3909e32-5a43-4ac1-8ac1-ac9feb7d7b8b
sum(Double64.(generate(10000)))

# ╔═╡ 25156ed5-8612-4652-ab20-90afd8416c92
md"Note the excellent agreement with the expected result of 1! 

Next we run the Power method to benchmark the performance. We generate a random dense matrix and starting vector:
"

# ╔═╡ 6f2c5dc3-410e-4011-be4c-bb267de6b539
begin
	Msample = diagm(abs.(randn(1000)))
	v0      = randn(size(Msample, 2))
end;

# ╔═╡ 62916272-9ab8-46b2-afce-c8a5c02c5a1c
md"First `Float64`:"

# ╔═╡ 58809041-6b74-4c91-a589-f581d054b88c
@btime power_method($Msample, $v0; verbose=false, maxiter=50, tol=1e-10)

# ╔═╡ bf73a491-7d4e-48d9-8c93-d6fbef1b7c3f
md"Then `Double64`, which is about 100 times slower:"

# ╔═╡ 01f1ef22-b7b4-4a9d-9378-98547f379902
let
	Md  = Double64.(Msample)
	v0d = Double64.(v0)
	@btime power_method($Md, $v0d; verbose=false, maxiter=50, tol=1e-10)
end

# ╔═╡ 804d1560-83fe-4455-8277-f3ce884c4887
md"In contrast with `BigFloat` we pay a factor of 1000:"

# ╔═╡ 7640338a-daf2-4f83-b1d7-0a3b95c86508
let
	Md  = BigFloat.(Msample)
	v0d = BigFloat.(v0)
	@btime power_method($Md, $v0d; verbose=false, maxiter=50, tol=1e-10)
end

# ╔═╡ 35e54531-a0c6-4419-9223-423d86ced2a2
md"""
### Interval arithmetic

The idea of interval arithmetic is to replace computation employing numbers by **computation employing sets**. To motivate why this might be a good idea we will consider yet another great example by William Kahan.

Consider the function
"""

# ╔═╡ 11b8bda7-a311-47bd-8974-780856ee1e54
begin
	f_interval(x) = (1/80) * log(abs(3*(1 - x) + 1)) + x^2 + 1

	plot(0.5:0.01:2.0, f_interval, leg=false)
end

# ╔═╡ 3cbaf02a-ad01-4b45-bbdf-0f89e685bb10
md"The function looks smooth except for a little blip around 1.3. Let's zoom in a little more closely to see what's going on:"

# ╔═╡ e15ce06b-f435-4f64-aaf2-1285faa1e4db
plot(1.2:0.0001:1.5, f_interval, leg=false)

# ╔═╡ 017d6bfc-7ad5-4f44-b30a-a3f9aec70e97
plot(f_interval, 1.2, 1.5, leg=false)

# ╔═╡ 8a388e43-c16c-432a-839d-9ca405331ee1
plot(f_interval, 1.32, 1.34, leg=false)

# ╔═╡ b389079a-ce96-4111-946d-e98d372e529b
md"Looks like a pretty sharp cusp. In fact if we stare at the original function for a while we see that the function diverges to $-\infty$ at $x = 4/3$.

But $4/3$ is not representable as a binary floating-point number. So what is happening at its neighbours:
"

# ╔═╡ edf7d11c-258a-4666-9286-095a8f63fbbc
begin
	x_minus = Float64(4//3, RoundDown)
	f_interval(x_minus)
end

# ╔═╡ b428c6f0-5812-40fc-8b75-1cd8713162ef
begin
	x_plus = nextfloat(x_minus)
	f_interval(x_plus)
end

# ╔═╡ 617b7c50-9041-445e-babf-cb008790af90
md"""
So in between the *neighbouring* floating-point numbers `x_minus` and `x_plus` the function crashes down to $-\infty$ and back. If we stick to `Float64` we can *never* see this.

Unfortunately even changing to `BigFloat` does not change too much:
"""

# ╔═╡ 8baec096-45c7-4300-bd14-0bd8b34eceb8
let
	x_minus_big = BigFloat(4//3, RoundDown)
	f_interval(x_minus_big)
end

# ╔═╡ c50598ab-0041-4d3d-9f4e-c138b68abe8b
let
	x_plus_big = nextfloat(BigFloat(4//3, RoundDown))
	f_interval(x_plus_big)
end

# ╔═╡ 7d3772a8-f335-4944-950e-252403a1ca74
md"""
One way to make sure we at least detect this behaviour is if we were able to compute $f(x)$ for *all* real $x$ in a certain *interval*, say in this case $\mathbf{x} = [1.2, 1.5]$ (we will denote intervals in bold). In other words we want to compute the range of $f$ over the input set $\mathbf{x}$, which we denote by $f(\mathbf{x})$. If $\mathbf{x}$ is a closed and bounded interval and $f$ is continuous than the range will also be a closed and bounded interval.

The idea to be able to determine the value range of functions given an input interval has lead to the development of **interval arithmetic**. Its main contribution is to provide a cheap computational receipe to calculate an *over-approximation* of the range of $f$, i.e. some interval $\mathbf{y}$ satisfying $f(\mathbf{x}) \subset \mathbf{y}$.

Let us first clarify on the term **interval**. In the context of interval arithmetic this always refers to a closed interval $\mathbf{x} = \{x \in \mathbb{R}\ |\ a \leq x \leq b\} \subseteq \overline{\mathbb{R}}$ where $a, b \in \overline{\mathbb{R}} = \mathbb{R} \cup \{\pm \infty\}$.

Now instead of approximating a quantity $x$ by a single floating-point number $\hat{x}$ we can alternatively approximate it by an interval $\mathbf{x}$ spanning the next smallest and next largest representable floating-point number. For example
"""

# ╔═╡ abe3fbe7-8354-4ea7-a73f-4d0361caa58a
four_thirds = I"4//3"   # Make interval representing the rational number 4/3

# ╔═╡ dd9d43ec-d467-4af9-8910-44ef39c6bfa7
typeof(four_thirds)

# ╔═╡ 2e9ca4cd-ade1-46b0-b57c-cb7fbee1bdcd
md"""This is a Julia interval, which has two fields `lo` and `hi`, representing the interval [`lo`, `hi`]."""

# ╔═╡ 21495bca-2cb2-458b-a7d0-3c4fd1470d1b
four_thirds.lo  # The Float64 number smaller than 4/3

# ╔═╡ 23d35d53-20a4-42a1-828d-823dbfa14aa2
four_thirds.hi  # The Float64 number larger than 4/3

# ╔═╡ 6142c138-d5dc-4e62-8d1a-9df8f8f8c7cd
md"""On top of this we define our operations on these intervals in a way the *containment property* is guaranteed. For example the difference between two intervals $\mathbf{x}$ and $\mathbf{y}$ is computed, such that it satisfies the following property:
```math
\forall x, y \in \overline{\mathbb{R}}: \qquad x \in \mathbf{x} \text{ and } y \in \mathbf{y} \Longrightarrow x - y \in \mathbf{x} - \mathbf{y}
\qquad \text{$(\ast)$}
```
One way to *compute* such a derivative interval $\mathbf{x} - \mathbf{y}$ is to ask : What are the smallest and largest values that $x - y$ can take, given that $x\in\mathbf{x}$ and $y\in\mathbf{y}$. More algorithmically if $\mathbf{x} = [x_\text{lo}, x_\text{hi}]$ and $\mathbf{y} = [y_\text{lo}, y_\text{hi}]$, then
```math
	\mathbf{x} - \mathbf{y} = [RD(x_\text{lo} - y_\text{hi}), RU(x_\text{hi} - y_\text{lo})].
```
Note that this choice of rounding modes along with a proper IEEE implementation of overflow and rounding are crucial in order to be able to provably obtain the containment propery $(\ast)$ above. Packages like `IntervalArithmetic` provide such operations adhering rigorously to the containment propery for pretty much all basic operations (`*`, `+`, `/`, functions like `sqrt`, ...). The full strength of interval arithmetic now results from the

**Fundamental theorem of interval arithmetic.** Given a function $f(x)$, which can be decomposed into a sequence of elementary operations $\left(e_1 \circ e_2 \circ \cdots \circ e_n\right)(x)$ and for each elementary operation an interval-analogue $e_i(\textbf{x})$ satisfying the containment property can be found, then by induction the end result of the calculation $\left(e_1 \circ e_2 \circ \cdots \circ e_n\right)(\textbf{x})$ is a guaranteed enclosure for the range of the original function $f$.

**Strength of interval arithmetic.** A direct consequence of the above theorem is that if we perform our computation using intervals and the resulting interval is narrow, then this is a **mathematical proof** that the exact answer can be found within the resulting interval. We can see this as a kind of automatic run-time error analysis.

For example the computation
"""

# ╔═╡ f64e2ffd-075f-4c68-b7c7-433736c0aea7
sum_naive(interval.(sample))

# ╔═╡ 9309b19e-6eed-41a8-b95e-fa04cda1cc4c
md"""
tells us that albeit the `sum_naive` function has some instabilities if the `sample` triggest too much cancellation, in this particular instance a few digits are *provably* still trustworthy.

**Issues with interval arithmetic.** First and foremostly interval arithmetic doubles the storage cost (two floats for each number) and increases the computational cost by a factor of $4$.
But the more problematic aspect is that the computed enclosures are *usually* too wide, which can cause issues when algorithms contain branching. For example suppose one needs to check $x \leq 0$ in order to determine which code branch to follow and the enclosing interval $\mathbf{x}$ contains the $0$, i.e. the number represented by $\mathbf{x}$ might be positve *or* negative. Which path should we follow? Often enough in complicated calculations one ends up with the interval $\mathbf{x} = [-\infty, \infty]$, which of course trivially contains the exact result, but is pointless to determine its value.

**Dependency problem.** A good example to illustrate the above issue is the dependency problem. Indeed, in many numerical computations intermediate sub-expressions can be strongly correlated, e.g. contain even common computation. Thus, their numerical error is similarly correlated --- which is completely ignored in naive interval arithmetic, which always assumes expressions to be fully independent of each other. The archetypical example is $x - x$. Irrespective of the numerical error in $x$ the result is $0$. Thus, the desired interval would just be $[0, 0]$ as well. However, when evaluating an expression $\textbf{x} - \textbf{x}$ interval arithmetic "does not know" that these are the same objects and gives a much larger interval answer, e.g.
"""

# ╔═╡ 8a1dffc9-a649-4e25-aec5-0c2e0ea8962c
let
	x = interval(0.0, 1.0)
	x - x
end

# ╔═╡ 22f0df08-a997-4b5c-9c52-cbcb3c87792c
md"""**Conclusion.**
The containment property of interval arithmetic **allows for a computational proof** that a floating-point **computation gives the correct answer**, and it provides a way to **compute the floating-point error**. Since intervals tend to be too wide, however, it is only a useful technique if applied to small and carefully chosen code segments. Hardly ever is it insightful to fully run a computation using intervals. Good applications for interval arithmetic in practice are therefore
- to rigorously prove expected invariances or convergence conditions, which should be satisfied if a computation has terminated successfully. For example in eigenpair computations one could use intervals solely in the computation of the residual norm $\|r\| = \|Ax - \lambda x\|$ once the iterative solver has terminated. If the interval is narrow around $0$, the problem has been provably solved exactly despite finite-precision arithmetic.
- as a red flag for numerical instability. For example if run in `Float32` the `sum_naive` returns
"""

# ╔═╡ efa2580e-ade3-411d-a9cb-116d4ae6bd91
sum_naive(Float32.(sample))

# ╔═╡ 434240a9-bdb4-4149-9a06-cd9312557edb
md"""which using intervals results in"""

# ╔═╡ 8f045523-9631-47b0-8d41-fe13dcbb7ca9
sum_naive(interval.(Float32.(sample)))

# ╔═╡ 919fa773-7364-4c0e-b712-3d3cbb805509
md"""a result that clearly should make us suspicious. Similarly, Kahan's example function to an interval containing $4/3$ warns us that the function *might* go down all the way to $-\infty$."""

# ╔═╡ bd6676f4-60b7-414d-b9e9-e96895aa285f
f_interval(interval(1.3, 1.4))

# ╔═╡ e75078d3-4cf9-49a0-977d-59200df1b3b2
TableOfContents()

# ╔═╡ 079db78c-2f64-4807-bfd5-c60b44ef0a0b
let
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	
	Sidebar(Markdown.parse(read("sidebar.md", String)), 530)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AccurateArithmetic = "22286c92-06ac-501d-9306-4abd417d9753"
BFloat16s = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
DoubleFloats = "497a8b3b-efae-58df-a0af-a86822472b78"
GenericLinearAlgebra = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
AccurateArithmetic = "~0.3.8"
BFloat16s = "~0.5.0"
BenchmarkTools = "~1.5.0"
DoubleFloats = "~1.4.0"
GenericLinearAlgebra = "~0.3.11"
HypertextLiteral = "~0.9.5"
IntervalArithmetic = "~0.22.14"
Plots = "~1.40.5"
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "b8751276d2476bd549891c9a0130ec11ad180faf"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AccurateArithmetic]]
deps = ["LinearAlgebra", "Random", "VectorizationBase"]
git-tree-sha1 = "07af26e8d08c211ef85918f3e25d4c0990d20d70"
uuid = "22286c92-06ac-501d-9306-4abd417d9753"
version = "0.3.8"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "5c9b74c973181571deb6442d41e5c902e6b9f38e"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.12.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BFloat16s]]
deps = ["LinearAlgebra", "Printf", "Random", "Test"]
git-tree-sha1 = "2c7cc21e8678eff479978a0a2ef5ce2f51b63dff"
uuid = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"
version = "0.5.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1dff6729bc61f4d49e140da1af55dcd1ac97b2f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.5.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

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

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "b8fe8546d52ca154ac556809e10c75e6e7430ac8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.5"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

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
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d8a9c0b6ac2d9081bf76324b39c78ca3ce4f0c98"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.6"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

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

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DoubleFloats]]
deps = ["GenericLinearAlgebra", "LinearAlgebra", "Polynomials", "Printf", "Quadmath", "Random", "Requires", "SpecialFunctions"]
git-tree-sha1 = "98d485da59c3f9d511429bdcb41b0762bf6ee1d5"
uuid = "497a8b3b-efae-58df-a0af-a86822472b78"
version = "1.4.0"

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
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

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

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "3f74912a156096bd8fdbef211eff66ab446e7297"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random", "libblastrampoline_jll"]
git-tree-sha1 = "02be7066f936af6b04669f7c370a31af9036c440"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.3.11"

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

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "433b0bb201cd76cb087b017e49244f10394ebe9c"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.14"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

    [deps.IntervalArithmetic.weakdeps]
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

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
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "5d3a5a206297af3868151bb4a2cf27ebce46f16d"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.33"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

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

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

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

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

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
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

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

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

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
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

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

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9dd97171646850ee607593965ce1f55063d8d3f9"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.0+0"

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
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "082f0c4b70c202c37784ce4bfbc33c9f437685bf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.5"

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
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "1a9cfb2dc2c2f1bd63f1906d72af39a79b49b736"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.11"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

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

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.Quadmath]]
deps = ["Compat", "Printf", "Random", "Requires"]
git-tree-sha1 = "67fe599f02c3f7be5d97310674cd05429d6f1b42"
uuid = "be4d8f0f-7fa4-5f49-b795-2f01399ab2dd"
version = "0.5.10"

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

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "677b65e17aeb6b4a0be1982e281ec03b0f55155c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.16"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "87d51a3ee9a4b0d2fe054bdd3fc2436258db2603"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.1.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "8963e5a083c837531298fc41599182a759a87a6d"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.1"

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

    [deps.StaticArrayInterface.weakdeps]
    OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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

[[deps.TranscodingStreams]]
git-tree-sha1 = "96612ac5365777520c3c5396314c8cf7408f436a"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.1"
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
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "e7f5b81c65eb858bed630fe006837b935518aca5"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.70"

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
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

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

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

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
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

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
# ╟─a3d727fa-fce6-43c7-ac54-51b5bcc0cde0
# ╠═fcf3b117-bb14-4849-818c-a0dfc4de15aa
# ╟─b1c739ed-957b-4a57-b1aa-6981bb3f76c4
# ╟─3c06f892-29d9-4df6-bd85-b22a19d12350
# ╟─77c9a481-b314-42f3-944e-b6422b894dd0
# ╠═0e5d0819-d101-4492-8bb0-73d5302ce461
# ╠═ac79b644-449d-4684-bc14-f06ca24ee160
# ╟─fd0a3bea-2178-4df5-96f4-b24a7e46b512
# ╟─cba52909-94af-4185-a990-36583976643c
# ╠═21d26883-5dcb-409b-bbee-64c6155b4d53
# ╟─2f7c1a64-bee2-470f-a2f6-887fad89994f
# ╠═c2daafe3-672c-4b9a-ad32-e4d773ecaeec
# ╟─408f1b1c-dd98-4e74-92ad-754b6122e525
# ╠═4751e59e-98fd-4547-bcd1-43ae740e2e6f
# ╟─6003c907-3722-4fd0-9eae-ab96ec88d902
# ╠═9dee1557-ec64-4146-abd0-aaba2ad3afdc
# ╟─fa2932d3-b09b-4ff2-8d59-295c4fe8cd5d
# ╠═65fc5ab7-364b-4e77-924e-bd0f2d3ad44e
# ╟─6f7058fe-171a-46a9-9d07-ccd65fb61bcd
# ╠═740d9e34-71c0-46dc-a9a0-0e4537e63e69
# ╟─8d126d65-34b3-478b-a1c8-03142fa4e888
# ╠═1cb60604-7d36-4f11-9683-d5a4ca821c7e
# ╠═2ec1b321-c63d-40d3-a00b-49b7a80098f8
# ╟─293e2d97-5e6f-4d69-99f9-1a9dde7e1c30
# ╠═c0d6e714-8c25-45eb-afe8-3fae25a174e6
# ╠═9368e779-9104-44e9-971d-95bc30ee38bc
# ╟─8b261d9b-22e3-4596-971d-30e9f51a6a0f
# ╟─ab77c84b-c0a3-436b-ac49-cd667983ead8
# ╟─3edfa209-d265-4b90-9461-a2d0a2461a5f
# ╠═bc5e8449-bf43-4c4e-8d12-bf035e1413f7
# ╟─102f8541-11bd-4b91-a8bd-21acb683ae20
# ╠═7f21c4b2-6ac0-4bd4-88ee-5a48fa3be773
# ╟─e61319a9-2a02-452e-96ac-60bc35f12a32
# ╠═360b47ab-992c-43a0-9afa-61865a3352ab
# ╟─7400e423-da76-4a61-b32e-bbb08c93609d
# ╠═a7f2a54b-e977-4c6f-a24d-994019d57b4b
# ╟─a4d4ced6-b71d-4fd2-8a26-262a443bc6c3
# ╟─399a584e-5c21-4502-b89e-ea4a3871a3f8
# ╠═605a633d-2362-4ebe-8ca0-7db925e52a2e
# ╟─bfff1a42-d18c-4bdb-bf3f-224ba5261eef
# ╠═562abdb2-0a8a-4bc6-aa2c-0266de924ee6
# ╟─18754b1e-f8d5-43ee-a668-a13dcd19f2e0
# ╠═188c06ba-2201-40ec-b14f-042da0c12444
# ╟─c19abeaa-7e9f-4a1f-892a-43ef9c19ee1a
# ╟─e8d93fd9-e32b-476f-9468-e637089a625f
# ╠═39e0ed4f-c812-47dd-b55d-3c4b5e07b2f5
# ╟─3f83b876-54a6-4900-be2a-d63609ed1947
# ╠═b822af3d-36e9-4985-a85f-eb34f96a660c
# ╟─2361d15d-01c9-4e5d-98c4-ddec8ec2a19a
# ╠═5d571dcf-aef0-448e-8eaf-704d4b02b124
# ╟─ef685429-5b8d-43ea-9480-eb2791345a9a
# ╠═fd98ab89-275c-49b4-9a46-17880f2b6540
# ╟─53dec758-3c39-4d2c-b9f4-58dba5b45c68
# ╠═3edb0d0f-4b75-4fa5-a5e3-80129abe2220
# ╟─1b4bc39a-eed6-451a-b16c-fb5025ff7012
# ╟─5c081511-7391-4d83-a1cf-ff83342b856e
# ╠═96e1effa-d45e-4b59-b019-93b08b12dca3
# ╟─ec14869d-e277-4628-b1ea-3823db0a2a0f
# ╠═e8d4565a-9a74-473a-b8b7-87a88f933bcb
# ╟─e19d8e1b-7e75-4e76-baa6-4d2c5bc15496
# ╟─8c32c50a-82aa-4750-b91f-e393e2eab11e
# ╠═ca3562da-7965-441d-bbaf-0b80c7edfa2d
# ╟─31be711b-4250-49e8-b9d5-3cf4db45108b
# ╟─97db876f-2aca-45fa-8b05-229ac1139a80
# ╟─75bed4a9-54f5-44c1-bcc3-74d5fb11c88a
# ╟─1d5c2471-1343-4706-9894-293da149ae74
# ╟─71e0fa21-4cfc-45a2-b6d4-a703a892c15d
# ╠═bbbf8e57-2ee0-44c6-9d84-d51c9fe5fd1e
# ╟─a9279bfa-7014-499d-9a64-0f6550a665f8
# ╠═75a6cd93-4609-4b93-a1f8-8e905cdc6602
# ╟─0be06f95-8464-479e-8b2c-3840669354b8
# ╟─1dcdb5e6-a4a0-40a9-9cf0-3072bad872f4
# ╠═3a59c0c6-5496-4faf-ae89-3c6333d02abb
# ╠═f2752b0d-b1d3-4f3d-bb0c-192177a4f997
# ╟─0fb4ef70-e529-41bc-a00d-3095d21d8ea2
# ╠═a7f9953b-0a14-4d6c-bed6-63a826d53a08
# ╠═2c625ae1-0275-4548-b8c1-a26a8ba50794
# ╠═64b720d1-c38d-485e-9243-c501c9155636
# ╠═2901a2eb-fc06-42e8-90fd-b3932529bc60
# ╟─1c4483d3-37f4-4043-93f6-a3ba09bcac36
# ╠═889f1186-de4e-4060-91b1-bc2511c7828e
# ╟─6ffe44f7-e507-4bd2-bf30-b06ccb0b2027
# ╠═35466aab-f3ba-4953-8fae-d545bc218ccd
# ╟─a8da2081-9fad-4be6-b486-c5cb6316d152
# ╠═1c2c9d4e-259f-4894-ad5b-87b176716829
# ╠═7108724e-5973-43bf-ab75-d4c5e6935cd1
# ╟─c1113c4b-55aa-4f0c-a602-14279d0702f0
# ╟─474ef698-cb16-44a3-b7fd-5a01005278e9
# ╟─b85874a4-d9ce-4131-9b1f-beaf79659b82
# ╟─f02e854d-50bb-43dc-8b4d-b85df11dc739
# ╟─ffb4717a-5059-4381-9de9-5d9e41a19183
# ╠═baa1625c-7faa-4341-a0b6-bb5ad3fdb67a
# ╟─39f6f711-eab3-4fa1-9061-ec7042010c63
# ╠═b5bc55d1-00ca-4733-8f32-08b4528ba297
# ╟─8f9e4905-5c39-4ccc-b0b5-af58e86df04c
# ╠═1343eaa4-d899-4ccd-8fbc-0f813e3a7ba6
# ╟─da63cabb-5f02-467f-9a30-87c4d4e37288
# ╠═ba7ca52e-e293-4515-a327-af2be18e41b1
# ╟─39e10226-e527-4699-af33-ea64ccdcb066
# ╠═a8cb9060-9fd1-4006-8f65-d1882b6be082
# ╠═3491e075-aac0-418d-8c01-dfee9d19bfd8
# ╠═1929e74c-e6b7-4ffe-b569-61248ba075c3
# ╠═540686e3-df39-44f4-86ff-315d5e080ead
# ╠═08a6f5a4-8c3f-4fa2-b40d-6575ba1597d0
# ╠═5162027a-b82f-4690-af0a-a6aabeed8652
# ╟─028b6e6f-e514-43ed-85ca-884b77000fb4
# ╟─30a193eb-9f22-4509-83b8-ea263884c021
# ╠═21e1d043-6b31-4632-8f8a-8148a5b18ab8
# ╠═a4d9ccda-9088-4d52-ab4a-584654ef20e8
# ╟─7683961e-fa0d-4960-80f6-c1af9d4593a7
# ╠═5435303f-2175-4ddd-80a0-7b069a7f8a8d
# ╟─12bfc098-7c8b-4d7a-b1e4-a2361beb70c3
# ╠═5adb95eb-a275-4d4d-9356-6fee5c387328
# ╟─70616c6f-3f61-494c-a8fa-5de02dbe9a9b
# ╠═7f849f7f-2425-4ef9-8e1e-dd36a0aab126
# ╟─6cc3184b-5c44-4aa4-84c9-7647872bdcd3
# ╠═4891a8a5-7941-4c07-8c20-cd2bf94ab827
# ╟─854efe00-c89f-40c0-b31e-42219acf5c44
# ╠═a4e5b19f-913c-406f-afe7-5d0fb6b3df3f
# ╟─e6fb6742-543b-4d90-b902-2f620f393b0a
# ╠═3c0cdb90-1ff6-488a-a6b4-190bf000feb9
# ╠═95b1ec78-033e-4f03-8adb-8045848ca123
# ╟─b7b2d7ff-b5b4-48c2-ad56-660cb813814b
# ╠═39f5ca43-fa61-4f27-b610-016638a54eeb
# ╠═e8476598-32e4-4c31-8d21-42835ceb36b6
# ╟─2d03b229-ba6e-447c-a13c-72132e5ccdc4
# ╠═3e0b819d-9dbe-45f7-9e92-c58578b85d21
# ╠═73594ce6-2eb0-444d-bd66-b490d9d72d1e
# ╠═d2755839-8576-4749-8e2e-34f8b6cbda6b
# ╟─6f507158-ef74-456b-b820-faaf6a4042e1
# ╠═a06746c0-c3c3-46cf-985f-2d6873f1f141
# ╠═e5269069-4235-4660-aee1-af5e149ef29d
# ╠═f29adb10-1992-4e77-ad8c-61a4877e78d5
# ╠═8240ab3a-5b57-4d40-8b04-e8a7bcf5d5cb
# ╟─9e2fc100-4165-43fd-87e4-65ea2bba6368
# ╠═3bfbecc2-c89c-4323-aff6-c7d7600767a6
# ╟─bbd6e360-fca7-41b4-83c8-f6b03c5b6243
# ╠═2f6f1ee8-a85c-4246-ad27-19f31eb1ce0f
# ╠═c300dde5-3a4c-43da-82c3-48b24a31ab04
# ╟─eb7a1a03-7cf1-4eb1-a7ea-b3c4ca4ad268
# ╠═67d5a7a0-9866-41bf-a242-70d9ae44a444
# ╟─95ed87cb-8a10-4b7e-8641-33277d87dab9
# ╠═f3909e32-5a43-4ac1-8ac1-ac9feb7d7b8b
# ╟─25156ed5-8612-4652-ab20-90afd8416c92
# ╠═6f2c5dc3-410e-4011-be4c-bb267de6b539
# ╟─62916272-9ab8-46b2-afce-c8a5c02c5a1c
# ╠═58809041-6b74-4c91-a589-f581d054b88c
# ╟─bf73a491-7d4e-48d9-8c93-d6fbef1b7c3f
# ╠═01f1ef22-b7b4-4a9d-9378-98547f379902
# ╟─804d1560-83fe-4455-8277-f3ce884c4887
# ╠═7640338a-daf2-4f83-b1d7-0a3b95c86508
# ╟─35e54531-a0c6-4419-9223-423d86ced2a2
# ╠═11b8bda7-a311-47bd-8974-780856ee1e54
# ╟─3cbaf02a-ad01-4b45-bbdf-0f89e685bb10
# ╠═e15ce06b-f435-4f64-aaf2-1285faa1e4db
# ╠═017d6bfc-7ad5-4f44-b30a-a3f9aec70e97
# ╠═8a388e43-c16c-432a-839d-9ca405331ee1
# ╟─b389079a-ce96-4111-946d-e98d372e529b
# ╠═edf7d11c-258a-4666-9286-095a8f63fbbc
# ╠═b428c6f0-5812-40fc-8b75-1cd8713162ef
# ╟─617b7c50-9041-445e-babf-cb008790af90
# ╠═8baec096-45c7-4300-bd14-0bd8b34eceb8
# ╠═c50598ab-0041-4d3d-9f4e-c138b68abe8b
# ╟─7d3772a8-f335-4944-950e-252403a1ca74
# ╠═8da0c6e6-ba1a-41ae-99f1-060a6cb53bbc
# ╠═abe3fbe7-8354-4ea7-a73f-4d0361caa58a
# ╠═dd9d43ec-d467-4af9-8910-44ef39c6bfa7
# ╟─2e9ca4cd-ade1-46b0-b57c-cb7fbee1bdcd
# ╠═21495bca-2cb2-458b-a7d0-3c4fd1470d1b
# ╠═23d35d53-20a4-42a1-828d-823dbfa14aa2
# ╟─6142c138-d5dc-4e62-8d1a-9df8f8f8c7cd
# ╠═f64e2ffd-075f-4c68-b7c7-433736c0aea7
# ╟─9309b19e-6eed-41a8-b95e-fa04cda1cc4c
# ╠═8a1dffc9-a649-4e25-aec5-0c2e0ea8962c
# ╟─22f0df08-a997-4b5c-9c52-cbcb3c87792c
# ╠═efa2580e-ade3-411d-a9cb-116d4ae6bd91
# ╟─434240a9-bdb4-4149-9a06-cd9312557edb
# ╠═8f045523-9631-47b0-8d41-fe13dcbb7ca9
# ╟─919fa773-7364-4c0e-b712-3d3cbb805509
# ╠═bd6676f4-60b7-414d-b9e9-e96895aa285f
# ╟─e75078d3-4cf9-49a0-977d-59200df1b3b2
# ╟─079db78c-2f64-4807-bfd5-c60b44ef0a0b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
