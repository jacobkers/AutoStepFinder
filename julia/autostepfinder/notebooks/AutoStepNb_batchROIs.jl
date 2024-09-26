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
path_free = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/chelated 1e-7M/Photobleaching_1mW_500ms/Field7"

# ╔═╡ f94daecf-afda-4b61-a5d1-3ae48c035480
begin
	img = asf.autostepfinder.get_image(path_free, 23)
	heatmap(img, c=:inferno, size=(800, 600), clim=(1700,2200))
end

# ╔═╡ 86710daa-2869-48dc-a923-054ca666c1ea
global roi2 =  [180 280 100 200
		80 180 100 200;
		180 280 200 300;
		80 180 200 300;
		180 280 300 400;
		80 180 300 400;
		 ];

# ╔═╡ 84fd4f73-1055-4a0b-b1ad-eced4689a791


# ╔═╡ bd225196-cea7-490d-b8bf-d515f51a31a3
md""" Select Threshold $(@bind tresH NumberField(0.1:10., default=0.1))"""

# ╔═╡ 8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
md"""# Analyse large set of data (free NAPH3)"""

# ╔═╡ b5207b13-c840-49c3-91e7-d7e764f70700
md""" Check to analyze all ROI: $(@bind zrois CheckBox())""" 

# ╔═╡ b6834d0a-d907-4f36-8a23-8e92737fde39
md""" # Analyse all regions with same ROI"""

# ╔═╡ d551c578-cf0e-4b9a-82fd-9a8b709a3205
md""" Check to analyze all fields and all ROI: $(@bind zfields CheckBox())""" 

# ╔═╡ dba17286-a139-4f97-9d2b-e8c1eb133eaf


# ╔═╡ e477b05f-2fbb-4bb8-a37a-a56c286445dc
md"""# Functions """

# ╔═╡ 2fceb655-a924-408c-92ea-f8cf03f0f6fb
function plot_roi(path, roi1)
	imroi = asf.autostepfinder.get_image(path, 22)[roi1[3]:roi1[4], roi1[1]:roi1[2]]
	gr()
	#c = cgrad([:red,:yellow,:green], [0.01, 0.9995], categorical = true)

	hmap1 = heatmap(imroi, clim=(1800,2200))
	hmap2 = heatmap(img,colorbar=false,  clim=(1700,2200))# xaxis=false, yaxis=false, )

	x_start = roi1[1]
	x_end = roi1[2]
	y_start = roi1[3]
	y_end = roi1[4]
	
	# Draw a rectangle on top of the heatmap using Shape
	plot!([x_start, x_end, x_end, x_start, x_start], [y_start, y_start, y_end, y_end, y_start],
	    seriestype = :shape, linecolor = :black, fillalpha = 0., linealpha = 1.0, linewidth = 2, label="ROI")
		
	pg = plot(hmap1, hmap2, )#size=(800, 400))
	return pg
end

# ╔═╡ 4ef27ec7-f07b-441b-8fbc-49f4cfba0b47
begin
	n = size(roi2)[1]
	PP = []
	for ir in 1:n
		roi1 = roi2[ir, :]
		pg = plot_roi(path_free, roi1)
		push!(PP, pg)
	end
	plot(PP...; size = (1000, 200) .* (1, n), layout = (n, 1), left_margin = 1Plots.mm)
end

# ╔═╡ e252b9cb-8a45-4a6e-bfb7-3ed0da607f6c
function analyze_roi(evol, tresH)
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
function save_main_results(path, tresH, roiN, res_array, plot_array)
	dirroot, fieldN = splitdir(path)
	
	dirres = joinpath(dirroot, "results", fieldN)
	try
		mkdir(dirres)
	catch
		@warn "Directory exists"
	end

	namefile = fieldN*"_roi"*string(roiN)*".csv"
	path_results = joinpath(dirres ,namefile);

	roi1 = roi2[roiN]
	header = " ", " ", "ROI position: "*string(roi1), "Threshold $tresH"
	CSV.write(path_results, [header], header=false)

	hh, hw, n_steps, ijy = res_array
	store_df(ijy, n_steps, hh, hw, path_results)
	dirsave, filename = splitdir(path_results)
	FieldN, _ = split(filename, "_")

	pg, phs = plot_array
	savefig(pg, joinpath(dirsave, FieldN*"_Roi$roiN.png"))
	savefig(phs, joinpath(dirsave, FieldN*"_Roi$roiN,hist.png"))
	return dirsave, fieldN
end

# ╔═╡ 2cf33115-23b3-4c43-87eb-4ab7aceba942
function analyse_all_rois(path, roi2, tresH)
	HS = []
	n = size(roi2)[1]
	for ir in 1:n
		roi1 = roi2[ir, :]
		
		imroi = asf.autostepfinder.get_image(path, 22)[roi1[3]:roi1[4], roi1[1]:roi1[2]]; 

		pg = plot_roi(path, roi1)
		evol = asf.autostepfinder.temporal_evolution(path, 400, roi1);
		
		phs, res_array = analyze_roi(evol, tresH);
		title!("Roi$ir")
		push!(HS, phs)

		dirsave, fieldN,  = save_main_results(path, tresH, ir, res_array, [pg, phs])

		#save_pixel_steps(dirsave, fieldN, evol, tresH, ir, imroi, res_array)		
	end
	return plot(HS...; size = (700, 200) .* (1, n), layout = (n, 1), left_margin = 1Plots.mm)
end

# ╔═╡ 95da7871-26c1-4cc2-b55c-147a2b53566c
# ╠═╡ show_logs = false
if zrois
	try
		analyse_all_rois(path_free, roi2, tresH, )
	catch
		@warn "Empty ROI"
	end
end

# ╔═╡ ad2eb7ed-110d-476e-a0fc-73e6a2420d7c
if zfields
	pathroot, _ = splitdir(path_free)
	
	for FN in 1:15
		path = joinpath(pathroot, "Field$FN")
		println("Analysing Field $FN")
		try
			analyse_all_rois(path, roi2, tresH)
		catch
			@warn "Empty ROI"
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

# ╔═╡ 7bb25cd5-0f74-48ce-aeee-3020221f6b9b
function save_pixel_steps(dirsave, fieldN, evol, tresH, roiN, imroi, res_array)
	hh, hw, n_steps, ijy = res_array
	for is in 1:size(ijy)[1]
		pxis, pxjs = ijy[is, :]
		ppx = plot_pixel(evol, pxis, pxjs, imroi, tresH)

		savefig(ppx, joinpath(dirsave, fieldN*"_Roi$roiN px($pxis,$pxjs).png"))
	end
end

# ╔═╡ Cell order:
# ╠═ac949aea-be3f-404d-8a8b-22c25b26ee38
# ╠═980dcf83-1e5c-4559-a42b-b02b0db094b3
# ╠═17c55175-2c5e-4691-90df-fdcb1b98d700
# ╠═a03bd74c-dc79-4f88-a142-815a6e94954c
# ╠═ee654cd8-b63a-4e5e-96f8-3104bea4fa60
# ╠═c22c5abe-980b-11ee-268c-b16eb65f35f4
# ╟─99ea72eb-b944-4d2a-9347-f4d9255ac790
# ╠═d3da2ac5-df51-4d7e-99f5-039595f6f75e
# ╠═b2b06a68-dcd5-43e1-9649-8574423e30ec
# ╠═6c52cf61-3acc-414d-bd04-303b031ef3b0
# ╠═74183602-063d-4aa9-884e-aff77e9606d9
# ╠═f94daecf-afda-4b61-a5d1-3ae48c035480
# ╠═86710daa-2869-48dc-a923-054ca666c1ea
# ╠═84fd4f73-1055-4a0b-b1ad-eced4689a791
# ╠═bd225196-cea7-490d-b8bf-d515f51a31a3
# ╠═4ef27ec7-f07b-441b-8fbc-49f4cfba0b47
# ╟─8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
# ╠═b5207b13-c840-49c3-91e7-d7e764f70700
# ╠═95da7871-26c1-4cc2-b55c-147a2b53566c
# ╠═b6834d0a-d907-4f36-8a23-8e92737fde39
# ╠═d551c578-cf0e-4b9a-82fd-9a8b709a3205
# ╠═ad2eb7ed-110d-476e-a0fc-73e6a2420d7c
# ╠═dba17286-a139-4f97-9d2b-e8c1eb133eaf
# ╠═e477b05f-2fbb-4bb8-a37a-a56c286445dc
# ╠═2cf33115-23b3-4c43-87eb-4ab7aceba942
# ╠═2fceb655-a924-408c-92ea-f8cf03f0f6fb
# ╠═e252b9cb-8a45-4a6e-bfb7-3ed0da607f6c
# ╠═0618b3c8-f194-4636-bd42-4923a5eff097
# ╠═8a55f7f9-7e29-486c-a816-82bb23b13de2
# ╠═7bb25cd5-0f74-48ce-aeee-3020221f6b9b
# ╠═f05f9a6c-ccea-497c-bc0e-6a02842cccc6
