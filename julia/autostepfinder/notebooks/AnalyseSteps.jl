### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ 638e7b43-5c8c-45e1-be75-f1bb0f89e0f1
using Pkg; Pkg.activate("/Users/pabloherrero/sabat/AutoStepFinder/julia/autostepfinder")

# ╔═╡ 18b49836-b0a5-11ee-3610-0f1662ab1bc8
begin
	using DelimitedFiles
	using Plots
	using PlutoUI
	using Images
	using ImageView
	using Colors
	using FixedPointNumbers
	using Statistics
	using Revise
	using StatsBase
	using CSV
	using DataFrames
	#using GR
end

# ╔═╡ 41fc28cb-e5aa-4d54-b00d-daa497cebbcc
PlutoUI.TableOfContents(title="AutoStepfinder Notebook", indent=true)

# ╔═╡ 019f1ff7-cdc3-4a80-89b8-493a8379c2e7
Pkg.instantiate()

# ╔═╡ b100f571-e1af-4cd5-b0c2-df3bb078a798
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ ee024a22-4ead-45ad-8e25-b7b4e4b98d74
asf = ingredients("../src/autostepfinder.jl")

# ╔═╡ cec211dc-1e19-4351-be3b-41c384f52e26
md""" # Test one field"""

# ╔═╡ 69770747-fea0-4c27-878c-6976d9a54893
begin
	pathroot = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/chelated 1e-7M/Photobleaching_1mW_500ms"
	img = asf.autostepfinder.get_image(pathroot*"/Field"*string(1), 23)
	nothing
end

# ╔═╡ dc4f9ba1-4268-4e11-8d5f-5749b23a132a
begin
	path = joinpath(pathroot, "results/full_im")
	pathdf = joinpath(path, "Field1/Field1.csv")
	df = DataFrame(CSV.File(pathdf, skipto=3, header=2), )
	nothing
end

# ╔═╡ df8b8c48-e1de-46be-a07c-afbd5015400c
n = 7

# ╔═╡ bbf808c3-7358-4746-895f-fdcd46525042
begin
	heatmap(img, c=:inferno, size=(800, 600), clim=(1700,2200), title="Field1")
	scatter!(df.x, df.y, c="white", label="single steps")
end

# ╔═╡ 484f88d7-20d2-40d1-8d58-9e01194f5e34
begin
	h1 = histogram(df.height, title="Step height", label="Unchelated",
		color=:green, xlim=(-500, 100))
	h2 = histogram(df.width, title="Step width", yscale=:log10, 
		color=:green,  label="Unchelated")
	n_stepsf = [count(==(i), df.x) for i in unique(df.x)]
	h3 = histogram(n_stepsf, title="Steps per pixel", color=:green,  label="Unchelated")
	plot(h1, h2, h3,layout=(1,3), size=(900, 300))
end

# ╔═╡ c2c3e9f8-543b-4a3f-8ad0-28dc68a7b272
md""" # Analyse all fields"""

# ╔═╡ 4ca47e06-759c-4e27-a4ad-a3597cd1192b
let
	ppp = []
	nfig = n
	for fi in 10:15
		try
			img = asf.autostepfinder.get_image(pathroot*"/Field"*string(fi), 23)	
			pathdf = joinpath(path, "Field$fi/Field$fi.csv")
			df = DataFrame(CSV.File(pathdf, skipto=3, header=2), )
			ht = heatmap(img, c=:inferno, size=(500, 600), clim=(1700,2200), legend = :none)
			#cbar = Colorbar(ht, width=20)
			
			scatter!(df.x, df.y, c="white", label="single steps")
			n_stepsf = [count(==(i), df.x) for i in unique(df.x)]
			h1 = histogram(n_stepsf, title="Steps per pixel", label="Chelated", left_margin = 1Plots.mm)
			hth = plot(ht, h1,layout=grid(1,2, widths=[0.6, 0.4]), title=("Field$fi"), right_margin = 10Plots.mm)#, size=(300, 300))
			
			h2 = histogram(df.height, title="Step height", label="Chelated")
			h3 = histogram(df.width, title="Step width", yscale=:log10, label="Chelated")
			pp = plot(h1, h2,layout=(1,2), size=(900, 300))
	
			pht = plot(hth, pp, layout=(2, 1))
			push!(ppp, pht)
		catch
			@warn "Empty field"
			nfig = nfig - 1
		end
	end
	plot(ppp...; size = (700, 700) .* (1, nfig), layout = (nfig, 1), left_margin = 1Plots.mm)
	#
end

# ╔═╡ bfb88fc6-8157-41d6-a764-4cf9eb2a084d
begin
	#dftot = DataType()
	ppp = []
	for fi in 10:15
		try
			img = asf.autostepfinder.get_image(pathroot*"/Field"*string(fi), 23)
			pathdf = joinpath(path, "Field$fi/Field$fi.csv")
			df = DataFrame(CSV.File(pathdf, skipto=3, header=2), )
			ht = heatmap(img, c=:inferno, size=(400, 300), clim=(1700,2200), legend = :none)
			title!("Field$fi")
			
			scatter!(df.x, df.y, c="white", label="single steps")
			
			push!(ppp, ht)
		catch 
			@warn "Empty field"
		end
	end
	plot(ppp...; size = (400, 400) .* (3, 3), layout = (3, 3), left_margin = 1Plots.mm)
	#
end

# ╔═╡ 11af101c-e27c-441d-807d-cbaef7e66a68
range = vcat([1:6...], [10:15...])

# ╔═╡ 638b2335-ec54-4e92-8439-07764713e2e6
begin
	hh, hw, ns = [], [], []
	for fi in range
		try
			pathdf = joinpath(path, "Field$fi/Field$fi.csv")
			df = DataFrame(CSV.File(pathdf, skipto=3, header=2), )
			n_stepsf = [count(==(i), df.x) for i in unique(df.x)]
			
			append!(hh, df.height)
			append!(hw, df.width)
			append!(ns, n_stepsf)
		catch
			@warn "Empty field"
		end
	end
	
end

# ╔═╡ 0a3650ec-445b-4f15-b8eb-86f4f0e41fb5
let
	h1 = histogram(ns, title="Steps per pixel", c=:blue, label="Chelated")
	h2 = histogram(hh, title="Step height",c=:blue, label="Chelated", xlim=(-500, 200))
	h3 = histogram(hw, title="Step width", yscale=:log10,c=:blue, label="Chelated")
	pp = plot(h1, h2, h3 ,layout=(1,3), size=(900, 400))		
end

# ╔═╡ 8a798fcc-01e1-4ab8-8406-e8a787e13369
histogram(hw, title="Step width", yscale=:log10,c=:blue, label="Chelated")

# ╔═╡ 16ce81e5-14a8-40a1-a375-10e099eb1331
md""" # Store results"""

# ╔═╡ 740a4cf5-eedb-42d0-89c8-76da5e53756b
md""" Check to store results: $(@bind zsave CheckBox())""" 

# ╔═╡ e163e9ac-2aae-43cc-b1fb-60e18d971cad
if zsave
	path_results = joinpath(path, "total_hist.csv")
	dfs = DataFrame(height=hh, width=hw)	
	
	CSV.write(path_results, string.(dfs), header=true, append=true)

	path_results2 = joinpath(path, "total_nstep.csv")
	
	
	CSV.write(path_results2, DataFrame(data=ns), header=true, append=true)
end

# ╔═╡ 4eb40e36-a852-438a-9390-546ae5af3d09
md""" # Compare samples from spin-coating"""

# ╔═╡ c91a3f96-66e8-4c53-ba35-aef8fe68c53d
begin
	pathq = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/Quartz_for_statistics/Sample1/results/full_im"
	dfq = DataFrame(CSV.File(joinpath(pathq, "total_hist.csv"), skipto=2) )
	nsq = DataFrame(CSV.File(joinpath(pathq, "total_nstep.csv"), skipto=2) )

	pathu = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/ZFFL62_free_NAPH3_1e-7M_for_statistics/results/full_im"
	dfu = DataFrame(CSV.File(joinpath(pathu, "total_hist.csv"), skipto=2) )
	nsu = DataFrame(CSV.File(joinpath(pathu, "total_nstep.csv"), skipto=2) )

	pathc = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/chelated 1e-7M/Photobleaching_1mW_500ms/results/full_im"
	dfc = DataFrame(CSV.File(joinpath(pathc, "total_hist.csv"), skipto=2) )
	nsc = DataFrame(CSV.File(joinpath(pathc, "total_nstep.csv"), skipto=2) )
	nothing
end

# ╔═╡ 76379a95-01f7-4900-947e-518c343987b5
histogram(dfc.width, yscale=:log10, bins=170, ylim=(8,3000))

# ╔═╡ c7be13b4-e2ec-45f8-b99a-52955715eee6
dfc.width

# ╔═╡ 6612e64d-75f0-4f60-8cfb-72f416476b06
let
	filtered_nsc = filter(x -> x < 40, nsc.data)
	h1 = histogram(filtered_nsc, c=:blue, label="chelated", bins = 20, xlim=(0,40),
		z=3, fillalpha=0.9, strokecolor=:blue, strokewidth=5, stroke=true)

	filtered_nsu = filter(x -> x < 40, nsu.data)
	histogram!(filtered_nsu, c=:green, label="unchelated", z=1, bins = 10)
	histogram!(nsq.data, title="Steps per pixel", c=:red, bins = 1, label="Quartz")
end

# ╔═╡ e2d17ca2-74c8-428c-b7c1-c30f46a9d0ac
let
	filtered_dfc = filter(x -> x > - 1500, dfc.height)
	h1 = histogram(filtered_dfc, c=:blue, label="chelated", bins = 200, xlim=(-500, 200), z=3, fillalpha=0.9, strokecolor=:blue, strokewidth=5, stroke=true)

	histogram!(dfu.height, c=:green, label="unchelated", z=1, bins = 200)
	histogram!(dfq.height, title="Step height", c=:red, bins = 20, label="Quartz")
end

# ╔═╡ 08cca352-a816-4d3b-8992-7c80a5cb8b60
let
	yax = (:log10, (1,1e4))
	filtered_dfc = filter(x -> x > - 1500, dfc.height)
	h1 = histogram(filtered_dfc, c=:blue, label="chelated", bins = 200, xlim=(-500, 200), z=3, fillalpha=0.9, strokecolor=:blue, yaxis = yax, strokewidth=5, stroke=true)

	filtered_nsu = filter(x -> x < 40, nsu.data)
	histogram!(dfu.height, c=:green, label="unchelated", z=1, yaxis = yax, bins = 200)
	histogram!(dfq.height, title="Step height", c=:red, bins = 30, yaxis = yax,  label="Quartz")
end

# ╔═╡ 6c07b4a8-49c5-4479-abb1-e8724c080522
let
	#filtered_dfc = filter(x -> x > - 500, dfc.width)
	h1 = histogram(dfc.width, c=:blue, label="chelated", bins = 200, 
		xlim=(0, 200),  z=4, fillalpha=0.9, yaxis = (:log10, (1,1e4)),
		strokecolor=:blue, strokewidth=5, stroke=true)

	histogram!(dfu.width, c=:green, label="unchelated", 
		yaxis = (:log10, (1,1e4)), fillalpha=1., z=1, bins = 200)
	histogram!(dfq.width, title="Step width", legend=:topright,
		c=:red, yaxis = (:log10, (1,1e4)), bins = 200, label="quartz")
end

# ╔═╡ b50e258b-9d85-4574-b20f-4b3330329325
let
	#filtered_dfc = filter(x -> x > - 1500, dfc.width)
	h1 = histogram(dfc.width, c=:blue, label="chelated", bins = 200, 
		xlim=(0, 200),  z=4, fillalpha=0.9,
		strokecolor=:blue, strokewidth=5, stroke=true)

	histogram!(dfu.width, c=:green, label="unchelated",  z=1, bins = 200)
	histogram!(dfq.width, title="Step width", legend=:topright,
		c=:red,  bins = 200, label="quartz")
end

# ╔═╡ 965d40a9-9cf6-4fa6-acd0-75d73f6b5667
md""" # Compare samples in vacuum"""

# ╔═╡ a54eb573-a106-4d1f-8f06-2229dc31d8fa
begin
	pathf = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/barium_evaporated+free/48/Photobleaching_1mW_BP450650/results/full_im/Field1"
	dff = DataFrame(CSV.File(joinpath(pathf, "Field1.csv"), skipto=3, header=2) )
	#nsu = DataFrame(CSV.File(joinpath(pathf, "total_nstep.csv"), skipto=2) )

	pathba = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/barium_evaporated+free/41/Photobleaching_1mW_BP450650/results/full_im/Field1"
	dfba = DataFrame(CSV.File(joinpath(pathba, "Field1.csv"), skipto=3, header=2) )
	nothing
end

# ╔═╡ b39489bd-5bde-4cb7-8842-d637bccecf0e
let 
	yax = (:log10, (1,1e4))
	filtered_dff = filter(x -> x > - 500, dff.height)
	filtered_dfba = filter(x -> x > - 500, dfba.height)
	
	h1 = histogram(filtered_dff, c=:blue, label="chelated", bins = 200, xlim=(-500, 200), z=1, fillalpha=0.9, strokecolor=:blue, strokewidth=5, stroke=true)
	histogram!(filtered_dfba, c=:green, label="unchelated", title="Step height", legend=:topright, z=2, bins = 50)

	
	h2 = histogram(dff.width, c=:blue, label="chelated", bins = 100, 
		xlim=(0, 200),  z=4, fillalpha=0.9,  title="Step width", legend=:topright,
		strokecolor=:blue, strokewidth=5, stroke=true)

	histogram!(dfba.width, c=:green, label="unchelated",  z=1, bins = 100)

	plot(h1, h2)
end

# ╔═╡ 879fab73-e4ec-4348-9084-403bcfd03f74


# ╔═╡ Cell order:
# ╠═638e7b43-5c8c-45e1-be75-f1bb0f89e0f1
# ╠═18b49836-b0a5-11ee-3610-0f1662ab1bc8
# ╠═41fc28cb-e5aa-4d54-b00d-daa497cebbcc
# ╟─019f1ff7-cdc3-4a80-89b8-493a8379c2e7
# ╟─b100f571-e1af-4cd5-b0c2-df3bb078a798
# ╠═ee024a22-4ead-45ad-8e25-b7b4e4b98d74
# ╠═cec211dc-1e19-4351-be3b-41c384f52e26
# ╠═69770747-fea0-4c27-878c-6976d9a54893
# ╠═dc4f9ba1-4268-4e11-8d5f-5749b23a132a
# ╠═df8b8c48-e1de-46be-a07c-afbd5015400c
# ╠═bbf808c3-7358-4746-895f-fdcd46525042
# ╠═484f88d7-20d2-40d1-8d58-9e01194f5e34
# ╠═c2c3e9f8-543b-4a3f-8ad0-28dc68a7b272
# ╠═4ca47e06-759c-4e27-a4ad-a3597cd1192b
# ╠═bfb88fc6-8157-41d6-a764-4cf9eb2a084d
# ╠═11af101c-e27c-441d-807d-cbaef7e66a68
# ╠═638b2335-ec54-4e92-8439-07764713e2e6
# ╠═0a3650ec-445b-4f15-b8eb-86f4f0e41fb5
# ╠═8a798fcc-01e1-4ab8-8406-e8a787e13369
# ╠═76379a95-01f7-4900-947e-518c343987b5
# ╠═c7be13b4-e2ec-45f8-b99a-52955715eee6
# ╠═16ce81e5-14a8-40a1-a375-10e099eb1331
# ╠═740a4cf5-eedb-42d0-89c8-76da5e53756b
# ╠═e163e9ac-2aae-43cc-b1fb-60e18d971cad
# ╠═4eb40e36-a852-438a-9390-546ae5af3d09
# ╠═c91a3f96-66e8-4c53-ba35-aef8fe68c53d
# ╠═6612e64d-75f0-4f60-8cfb-72f416476b06
# ╠═e2d17ca2-74c8-428c-b7c1-c30f46a9d0ac
# ╠═08cca352-a816-4d3b-8992-7c80a5cb8b60
# ╠═6c07b4a8-49c5-4479-abb1-e8724c080522
# ╠═b50e258b-9d85-4574-b20f-4b3330329325
# ╠═965d40a9-9cf6-4fa6-acd0-75d73f6b5667
# ╠═a54eb573-a106-4d1f-8f06-2229dc31d8fa
# ╠═b39489bd-5bde-4cb7-8842-d637bccecf0e
# ╠═879fab73-e4ec-4348-9084-403bcfd03f74
