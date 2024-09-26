### A Pluto.jl notebook ###
# v0.19.8

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

# ╔═╡ ac949aea-be3f-404d-8a8b-22c25b26ee38
using Pkg; Pkg.activate("/Users/pabloherrero/sabat/AutoStepFinder/julia/autostepfinder")

# ╔═╡ a03bd74c-dc79-4f88-a142-815a6e94954c
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
begin
	Pkg.add("DelimitedFiles")
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
	Pkg.add("Images")
	Pkg.add("ImageView")
	Pkg.add("Colors")
	Pkg.add("FixedPointNumbers")
	Pkg.add("Statistics")
	Pkg.add("Revise")
	#Pkg.add("GR")
end
  ╠═╡ =#

# ╔═╡ 980dcf83-1e5c-4559-a42b-b02b0db094b3
# ╠═╡ show_logs = false
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

# ╔═╡ 17c55175-2c5e-4691-90df-fdcb1b98d700
default(size=(800, 400))

# ╔═╡ ee654cd8-b63a-4e5e-96f8-3104bea4fa60
Pkg.instantiate()

# ╔═╡ 99ea72eb-b944-4d2a-9347-f4d9255ac790
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

# ╔═╡ c22c5abe-980b-11ee-268c-b16eb65f35f4
# ╠═╡ show_logs = false
asf = ingredients("../src/autostepfinder.jl")

# ╔═╡ d3da2ac5-df51-4d7e-99f5-039595f6f75e
PlutoUI.TableOfContents(title="AutoStepfinder Notebook", indent=true)

# ╔═╡ b2b06a68-dcd5-43e1-9649-8574423e30ec
md"""# Plot image"""

# ╔═╡ 6c52cf61-3acc-414d-bd04-303b031ef3b0
default(size=(800, 400))

# ╔═╡ 74183602-063d-4aa9-884e-aff77e9606d9
path_free = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/barium_evaporated+free/48/Photobleaching_1mW_BP450650/Field1"

# ╔═╡ f94daecf-afda-4b61-a5d1-3ae48c035480
begin
	img = asf.autostepfinder.get_image(path_free, 23)
	heatmap(img, c=:inferno, size=(800, 600), clim=(1700,2200))
end

# ╔═╡ eaf9494d-e8d6-424d-92b0-681562640e8b
roi1 = [1, 512, 1, 512]

# ╔═╡ bd225196-cea7-490d-b8bf-d515f51a31a3
md""" Select Threshold $(@bind tresH NumberField(0.1:10., default=0.1))"""

# ╔═╡ 8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
md"""# Analyse large set of data (free NAPH3)"""

# ╔═╡ d551c578-cf0e-4b9a-82fd-9a8b709a3205
md""" Check to analyze all fields: $(@bind zfields CheckBox())""" 

# ╔═╡ 8b2f646e-e511-495a-8df3-e4d24e916627
begin
	path = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/barium_evaporated+free/48/Photobleaching_1mw_BP450650/Field2"
	#evol = temporal_evolution(path, 400, roi1);
	#phs, res_array = analyze_im(evol, tresH);
end

# ╔═╡ 0db693d2-9e78-4dfd-9ca3-370c8ae0b8fa

function temporal_evolution(folder, N, roi=[1, 512, 1, 512])
    im0 = asf.autostepfinder.get_image(folder, 1)[roi[1]:roi[2], roi[3]:roi[4]]
    data = [im0]
    
    for i in 2:N
        d = get_image(folder, i)[roi[1]:roi[2], roi[3]:roi[4]]
        push!(data, d)
    end
    
    evol = cat(data..., dims=3)
    return evol
end

# ╔═╡ e477b05f-2fbb-4bb8-a37a-a56c286445dc
md"""# Functions """

# ╔═╡ e252b9cb-8a45-4a6e-bfb7-3ed0da607f6c
function analyze_im(evol, tresH)
	hh, hw, n_steps, ijy = asf.autostepfinder.get_hist_roi(evol, tresH, 50); 

	res_array = hh, hw, n_steps, ijy
	p1 = histogram(n_steps, label="Number of steps")
	p2 = histogram(hh, label="Step height")
	p3 = histogram(hw, label="Step width")
	phs = plot(p1, p2, p3, layout=(1,3))
	return phs, res_array
end

# ╔═╡ 0618b3c8-f194-4636-bd42-4923a5eff097
function store_df(ijy, n_steps, hh, hw, path )
	
	ijsave = hcat(inverse_rle(ijy[:,1], n_steps), inverse_rle(ijy[:,2], n_steps))

	dfs = DataFrame(x=ijsave[:,1], y=ijsave[:,2], height=hh, width=hw)	
	
	CSV.write(path, string.(dfs), header=true, append=true)
end

# ╔═╡ 8a55f7f9-7e29-486c-a816-82bb23b13de2
function save_main_results(path, tresH, res_array, phs)
	dirroot, fieldN = splitdir(path)
	
	dirres = joinpath(dirroot, "results", fieldN)
	try
		mkdir(dirres)
	catch
		@warn "Directory exists"
	end

	namefile = fieldN*".csv"
	path_results = joinpath(dirres ,namefile);

	header = " ", " ", "ROI position: full", "Threshold $tresH"
	CSV.write(path_results, [header], header=false)

	hh, hw, n_steps, ijy = res_array
	store_df(ijy, n_steps, hh, hw, path_results)
	dirsave, filename = splitdir(path_results)

	savefig(phs, joinpath(dirsave, fieldN*"hist.png"))
	return dirsave, fieldN
end

# ╔═╡ ad2eb7ed-110d-476e-a0fc-73e6a2420d7c
if zfields
	pathroot, _ = splitdir(path_free)
	HS = []
	for FN in 1:1
		path = joinpath(pathroot, "Field$FN")
		println("Analysing Field $FN")
		evol = asf.autostepfinder.temporal_evolution(path, 400, roi1);
		try
			phs, res_array = analyze_im(evol, tresH);
			title!("Field$FN")		
		
			push!(HS, phs)

			dirsave, fieldN,  = save_main_results(path, tresH, res_array, phs)
		catch
			@warn "Empty field"
		end
	end

end

# ╔═╡ f05f9a6c-ccea-497c-bc0e-6a02842cccc6
function plot_pixel(evol, i, j, img, tresH=0.1, N_iter=50, )
	dataX = evol[ i, j, :]
	S_curve, best_shots, Fits, steptable = asf.autostepfinder.AutoStepMain(dataX, tresH, N_iter)
	p1 = plot(dataX, label = "Data pixel ($i,$j)")
	plot!(Fits, label = "AutoStepFinder fit")

	p2 = plot(S_curve)
	p2 = heatmap(img,colorbar=false, c=:thermal)# xaxis=false, yaxis=false, )

	x_start = i - 2
	x_end = i + 2
	y_start = j - 2
	y_end = j + 2
	
	# Draw a rectangle on top of the heatmap using Shape
	plot!([x_start, x_end, x_end, x_start, x_start], [y_start, y_start, y_end, y_end, y_start],
	    seriestype = :shape, linecolor = :white, fillalpha = 0., linealpha = 1.0, linewidth = 4, label="Pixel ($i, $j)")
		
	ppx = plot(p1, p2)
	return ppx
end

# ╔═╡ Cell order:
# ╠═ac949aea-be3f-404d-8a8b-22c25b26ee38
# ╠═980dcf83-1e5c-4559-a42b-b02b0db094b3
# ╠═17c55175-2c5e-4691-90df-fdcb1b98d700
# ╠═a03bd74c-dc79-4f88-a142-815a6e94954c
# ╠═ee654cd8-b63a-4e5e-96f8-3104bea4fa60
# ╠═c22c5abe-980b-11ee-268c-b16eb65f35f4
# ╟─99ea72eb-b944-4d2a-9347-f4d9255ac790
# ╟─d3da2ac5-df51-4d7e-99f5-039595f6f75e
# ╠═b2b06a68-dcd5-43e1-9649-8574423e30ec
# ╠═6c52cf61-3acc-414d-bd04-303b031ef3b0
# ╠═74183602-063d-4aa9-884e-aff77e9606d9
# ╠═f94daecf-afda-4b61-a5d1-3ae48c035480
# ╠═eaf9494d-e8d6-424d-92b0-681562640e8b
# ╠═bd225196-cea7-490d-b8bf-d515f51a31a3
# ╟─8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
# ╠═d551c578-cf0e-4b9a-82fd-9a8b709a3205
# ╠═ad2eb7ed-110d-476e-a0fc-73e6a2420d7c
# ╠═8b2f646e-e511-495a-8df3-e4d24e916627
# ╠═0db693d2-9e78-4dfd-9ca3-370c8ae0b8fa
# ╠═e477b05f-2fbb-4bb8-a37a-a56c286445dc
# ╠═e252b9cb-8a45-4a6e-bfb7-3ed0da607f6c
# ╠═0618b3c8-f194-4636-bd42-4923a5eff097
# ╠═8a55f7f9-7e29-486c-a816-82bb23b13de2
# ╠═f05f9a6c-ccea-497c-bc0e-6a02842cccc6
