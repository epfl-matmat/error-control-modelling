### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ a9d71171-dd24-472a-be5f-12619295c4bf
begin
	using HypertextLiteral
	using PlutoUI
	using PlutoTeachingTools
end

# ╔═╡ 9261420d-3754-4820-a4f0-c5bda988157a
md"""
# MATH-500: Error control in scientific modelling
"""

# ╔═╡ 09579796-7e0e-4990-b727-b20f63242bcb
md"""
This repository contains the teaching material of the course
[MATH-500 Error control in scientific modelling](https://edu.epfl.ch/coursebook/en/error-control-in-scientific-modelling-MATH-500) at EPFL.
"""

# ╔═╡ c0da4104-9a5a-4553-8df0-fdadd0194c14
md"""
## Course outline
"""

# ╔═╡ d448115a-6857-4082-893b-37549fe83b8b
md"""
If you spot an error feel free to make a pull request to the
[github repository](https://github.com/epfl-matmat/error-control-modelling)
generating this website.
"""

# ╔═╡ 854d367b-1609-4c00-9c62-9439874a03d6
md"""
## EPFL resources

If you are an EPFL student, you can find the following course resources:
- **Moodle link:** <https://go.epfl.ch/error-control>

"""

# ╔═╡ ebf43572-95f2-401e-aa2e-d7748d3bf639
md"""
## Summary
Errors are ubiquitous in computational science as neither models nor numerical
techniques are perfect. With respect to eigenvalue problems motivated from
materials science and atomistic modelling we discuss,
implement and apply numerical techniques for estimating simulation error.
"""

# ╔═╡ 7e571edb-6206-4138-863e-5a9d804e07ec
md"""
## Content
* Important eigenvalue problems in materials science
* Motivation for studying errors in eigenvalue problems
* Types of simulation error
* Residual-error relationships for eigenvalue problems
* Perturbation theory and parametrised eigenvalue problems
* Subtleties of infinite-dimensional eigenvalue problems
* Discretisation and discretisation error
* Plane-wave basis sets
* Errors due to uncertain parameters (if time permits)
* Non-linear eigenvalue problems (if time permits)

Algorithm demonstrations and implementations will be based on
the [Julia programming language](https://julialang.org/)
and interactive [Pluto](https://plutojl.org/) notebooks.
"""

# ╔═╡ 0ccf77bd-0c5d-4c3c-82b5-b25c4adbc3ba
md"""
## Prerequisites
* Analysis
* Linear algebra
* Exposure to numerical linear algebra
* Some experience with numerical methods for solving differential equations
  (such as finite-element methods, finite-difference approaches, plane-wave methods)
* Exposure to implementing numerical algorithms (e.g. using Python or Julia)

This course delivers a mathematical viewpoint on materials modelling and it is
explicitly intended for an interdisciplinary student audience. To keep it
accessible, the key mathematical and physical concepts will both be revised as
we go along. However, the learning curve will be steep and an interest to learn
about the respective other discipline is required.

The problem sheets and the
projects require a substantial amount of work and feature both theoretical
(proof-oriented) and applied (programming-based and simulation-based) components.
While there is some freedom for students to select their respective
focus, students are encouraged to team up across the disciplines
for the course work.

In the past participants from materials science found it useful to take this course after they followed the lectures on [Fundamentals of solid state materials](https://edu.epfl.ch/coursebook/en/fundamentals-of-solid-state-materials-MSE-423).
"""

# ╔═╡ 1c4b711d-7275-4b23-82b6-bf0e57ff5d51
md"""
## Resources on the Julia programming language

In this class we employ the [Julia programming language](https://julialang.org/) for code examples and programming exercises. Helpful information and pointers to further resources on Julia can be found in the online material of the [MATH-251(b) Numerical Analysis](https://teaching.matmat.org/numerical-analysis) lecture:
- The [Julia introduction notebook](https://teaching.matmat.org/numerical-analysis/02_Julia.html)
- The [Getting started with Julia](https://teaching.matmat.org/numerical-analysis/exercises/ex0_introduction_julia_pluto_statement.html) exercise sheet
- The [Basic plotting with Julia](https://teaching.matmat.org/numerical-analysis/exercises/ex0_introduction_plots_statement.html) exercise sheet
"""

# ╔═╡ b1c510c9-4130-4f5e-8d15-3bfaf0a61033
md"""
## Literature
There is no single textbook corresponding to the content of the course.
Parts of the lectures have substantial overlap with the following resources,
where further information can be found.

- Youssef Saad. *Numerical Methods for Large Eigenvalue Problems*, SIAM (2011).
  A PDF is available free of charge on
  [Youssef Saad's website](https://www-users.cse.umn.edu/~saad/eig_book_2ndEd.pdf).
- Nicholas J. Higham. *Accuracy and Stability of Numerical Algorithms*, SIAM (2002).
- Peter Arbenz. *Lecture notes on solving large scale eigenvalue problems*, ETHZ.
  A PDF is available from
  [Peter Arbenz' website](https://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf).
- Mathieu Lewin. *Théorie spectrale et mécanique quantique*, Springer (2022).
  A PDF for EPFL students is available from
  [Springer Link]([https://link.springer.com/book/10.1007/978-3-030-93436-1).
"""

# ╔═╡ 4e1fa8a0-ad16-417f-ba58-5cb0ccff45a3
TableOfContents()

# ╔═╡ 8d38c66a-9582-4ef0-96f7-73e0b2ae2e2f
begin
	RobustLocalResource("https://teaching.matmat.org/error-control/sidebar.md", "sidebar.md")
	toc = Markdown.parse(read("sidebar.md", String))
	Sidebar(toc, ypos) = @htl("""<aside class="plutoui-toc aside indent"
		style='top:$(ypos)px; max-height: calc(100vh - $(ypos)px - 55px);' >$toc</aside>""")
	Sidebar(toc, 315)
end

# ╔═╡ fefc5308-7caf-465a-afc0-b1a6f2e68ba2
toc

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
HypertextLiteral = "~0.9.5"
PlutoTeachingTools = "~0.4.5"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "a0d9630b3fb957528e77b2be269fb0ff5bd8a3e3"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "52e1296ebbde0db845b356abbbe67fb82a0a116c"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.9"

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
git-tree-sha1 = "85778cdf2bed372008e6646c64340460764a5b85"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.5"

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
# ╟─9261420d-3754-4820-a4f0-c5bda988157a
# ╟─09579796-7e0e-4990-b727-b20f63242bcb
# ╟─c0da4104-9a5a-4553-8df0-fdadd0194c14
# ╟─a9d71171-dd24-472a-be5f-12619295c4bf
# ╟─fefc5308-7caf-465a-afc0-b1a6f2e68ba2
# ╟─d448115a-6857-4082-893b-37549fe83b8b
# ╟─854d367b-1609-4c00-9c62-9439874a03d6
# ╟─ebf43572-95f2-401e-aa2e-d7748d3bf639
# ╟─7e571edb-6206-4138-863e-5a9d804e07ec
# ╟─0ccf77bd-0c5d-4c3c-82b5-b25c4adbc3ba
# ╟─1c4b711d-7275-4b23-82b6-bf0e57ff5d51
# ╟─b1c510c9-4130-4f5e-8d15-3bfaf0a61033
# ╟─4e1fa8a0-ad16-417f-ba58-5cb0ccff45a3
# ╟─8d38c66a-9582-4ef0-96f7-73e0b2ae2e2f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
