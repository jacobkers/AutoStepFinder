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
path_free = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/for_statistics/chelated 1e-7M/Photobleaching_1mW_500ms/Field8"

# ╔═╡ f94daecf-afda-4b61-a5d1-3ae48c035480
begin
	img = asf.autostepfinder.get_image(path_free, 23)
	heatmap(img, c=:inferno, size=(800, 600), clim=(1700,2200))
end

# ╔═╡ 86710daa-2869-48dc-a923-054ca666c1ea
roi2 =  [180 280 100 200
		80 180 100 200;
		180 280 200 300;
		80 180 200 300;
		180 280 300 400;
		80 180 300 400;
		1 512 1 512];

# ╔═╡ cd90501e-5a66-49a7-ad72-341072bdb69e
md""" Select ROI number: $(@bind roiN NumberField(1:10, default=1))"""

# ╔═╡ 12b432f4-5942-4958-9de8-473a08c364a5
begin
	roi1 = roi2[roiN, :]
	evol = asf.autostepfinder.temporal_evolution(path_free, 400, roi1);
	####roi1 = [80, 280, 100, 400]
	nothing
end

# ╔═╡ 10b8263d-8dc5-4ddc-bd21-b310ac7fe11b
begin 
	imroi = asf.autostepfinder.get_image(path_free, 22)[roi1[3]:roi1[4], roi1[1]:roi1[2]]
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
		
	pg = plot(hmap1, hmap2, size=(800, 400))#, aspect_ratio=:equal)
	#display(fig)
end

# ╔═╡ 8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
md"""# Analyse large set of data (free NAPH3)"""

# ╔═╡ bd225196-cea7-490d-b8bf-d515f51a31a3
md""" Select Threshold $(@bind tresH NumberField(0.1:10., default=0.1))"""

# ╔═╡ ccab705b-5193-4618-89f1-ced300100ba3
md""" Check to analyze steps using full ROI: $(@bind zroi CheckBox())""" 

# ╔═╡ ed3cf054-e443-4d4f-87de-280d72044e5a


# ╔═╡ e650b1b7-fa2a-4ad9-8c33-836648987001
md"""### Store results"""

# ╔═╡ d31fb342-d3cc-4fa6-99f5-f3a2ac7c7b42
md""" Check to store histogram data: $(@bind zsave CheckBox())"""

# ╔═╡ 680641ef-16b5-4fca-9a57-149bf8c72eb8
md"""## Plot one evolution"""

# ╔═╡ 630597bd-e8d6-45a8-8e1e-4f2fcc02309c
begin
	dataX = evol[2,52, :]
	S_curve, best_shots, Fits, steptable = asf.autostepfinder.AutoStepMain(dataX, 0.4, 50)
end

# ╔═╡ 7535da1a-4d29-4866-a212-5f68faad75de
cumsum(steptable[1:LST, 5])

# ╔═╡ afd1ef5f-a8ce-4bff-b001-69c1e2731d85
steptable[1:LST, 5]

# ╔═╡ 2ee80e6a-ef83-4c78-8b42-0b864405249f
steptable

# ╔═╡ 5241b7d9-1f56-4b87-bf3a-fdd10cbbb71d


# ╔═╡ 94926618-9161-4ea5-bcd1-346ee41380f0


# ╔═╡ 39a8c1dc-7a7e-4958-a997-039d2eeb20bd
md"""# Analyse Naph3+Ba """

# ╔═╡ 2ba85012-706a-4719-ab87-e0fa91d2f36f
path_ba = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/commercial_molecules/Photobleaching_1mW_100ms/Field3_long_sequence"

# ╔═╡ 5954a7be-70c1-4d23-b4c2-719f23220212
begin
	imgba = asf.autostepfinder.get_image(path_ba, 23)
	heatmap(imgba,  clim=(1800,2500))
end

# ╔═╡ e477b05f-2fbb-4bb8-a37a-a56c286445dc
md"""# Functions """

# ╔═╡ 7c6dcca9-fd13-4977-8bf7-7f0f5d2b5636
function get_hist_roi(evol, tresH = 0.1, N_iter = 50)
	hist_height = Vector{Any}()
	hist_widths = Vector{Any}()
	n_steps = Vector{Int}()
	
	I, J, _ = size(evol)
	ijyes = reshape([], 0, 2)
	
	for i in 1:I
	    for j in 1:J
			dataX = evol[ i, j, :]
			
			
			S_curve, best_shots, Fits, steptable = asf.autostepfinder.AutoStepMain(dataX, tresH, N_iter)
			
	            if best_shots > 0 && S_curve[best_shots] > tresH 

					LST = size(steptable)[1]
					if LST > 2
						push!(n_steps, size(steptable, 1))
						push!(hist_height, steptable[1:LST, 4])
						push!(hist_widths, cumsum(steptable[1:LST, 5]))
						ijyes = [ijyes;  [i j]]
						
					end
	                    
	            end

			end
	    end
	
	hw = vcat(hist_widths...)
	hh = vcat(hist_height...)
	return hh, hw, n_steps, ijyes
end

# ╔═╡ 78f8ff68-ec80-4988-a240-064b1e2bc93c
if zroi
	hh, hw, n_steps, ijy = get_hist_roi(evol, tresH, 50); 
	nothing
end

# ╔═╡ 4cbbadc4-9948-4a21-ab03-18ed7c73573d
begin 
	p1 = histogram(n_steps, label="Number of steps")
	p2 = histogram(hh, label="Step height")
	p3 = histogram(hw, label="Step width")
	phs = plot(p1, p2, p3, layout=(1,3))
end

# ╔═╡ b1076909-7591-41aa-bada-4747616cec67
md""" Select pixel: $(@bind ijview NumberField(1:size(ijy)[1], default=1))"""

# ╔═╡ 0618b3c8-f194-4636-bd42-4923a5eff097
function store_df(ijy, n_steps, hh, hw, path )
	
	ijsave = hcat(inverse_rle(ijy[:,1], n_steps), inverse_rle(ijy[:,2], n_steps))
	dfs = DataFrame(x=ijsave[:,1], y=ijsave[:,2], height=hh, width=hw)	
	
	CSV.write(path, string.(dfs), header=true, append=true)
end

# ╔═╡ d365f70f-16d2-4d08-accb-0f713c6440c6
if zsave

	dirroot, fieldN = splitdir(path_free)

	dirres = joinpath(dirroot, "results", fieldN)
	try
		mkdir(dirres)
	catch
		@warn "Directory exists"
	end
	namefile = fieldN*"_roi"*string(roiN)*".csv"
	path_results = joinpath(dirres ,namefile);
	
	header = " ", " ", "ROI position: "*string(roi1), "Threshold $tresH"
	CSV.write(path_results, [header], header=false)
	
	store_df(ijy, n_steps, hh, hw, path_results)
	dirsave, filename = splitdir(path_results)
	fieldN, _ = split(filename, "_")
	
	savefig(pg, joinpath(dirsave, fieldN*"_Roi$roiN.png"))
	savefig(phs, joinpath(dirsave, fieldN*"_Roi$roiN,hist.png"))

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

# ╔═╡ bf75061a-48f5-437d-900e-aa4e1abf8047
if zsave
	for is in 1:size(ijy)[1]
		pxis, pxjs = ijy[is, :]
		ppx = plot_pixel(evol, pxis, pxjs, imroi, tresH)

		savefig(ppx, joinpath(dirsave, fieldN*"_Roi$roiN px($pxis,$pxjs).png"))
	end
end

# ╔═╡ 9a6e62bb-289e-4ca7-9f95-cd0caaac7700
begin
	pxi, pxj = ijy[ijview, :]
	plot_pixel(evol, pxi, pxj, imroi, tresH)
end

# ╔═╡ Cell order:
# ╠═ac949aea-be3f-404d-8a8b-22c25b26ee38
# ╠═980dcf83-1e5c-4559-a42b-b02b0db094b3
# ╠═17c55175-2c5e-4691-90df-fdcb1b98d700
# ╠═a03bd74c-dc79-4f88-a142-815a6e94954c
# ╠═ee654cd8-b63a-4e5e-96f8-3104bea4fa60
# ╠═c22c5abe-980b-11ee-268c-b16eb65f35f4
# ╠═99ea72eb-b944-4d2a-9347-f4d9255ac790
# ╠═d3da2ac5-df51-4d7e-99f5-039595f6f75e
# ╠═b2b06a68-dcd5-43e1-9649-8574423e30ec
# ╠═6c52cf61-3acc-414d-bd04-303b031ef3b0
# ╠═74183602-063d-4aa9-884e-aff77e9606d9
# ╠═f94daecf-afda-4b61-a5d1-3ae48c035480
# ╠═86710daa-2869-48dc-a923-054ca666c1ea
# ╟─cd90501e-5a66-49a7-ad72-341072bdb69e
# ╠═12b432f4-5942-4958-9de8-473a08c364a5
# ╟─10b8263d-8dc5-4ddc-bd21-b310ac7fe11b
# ╟─8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
# ╟─bd225196-cea7-490d-b8bf-d515f51a31a3
# ╠═ccab705b-5193-4618-89f1-ced300100ba3
# ╟─4cbbadc4-9948-4a21-ab03-18ed7c73573d
# ╠═78f8ff68-ec80-4988-a240-064b1e2bc93c
# ╠═ed3cf054-e443-4d4f-87de-280d72044e5a
# ╟─e650b1b7-fa2a-4ad9-8c33-836648987001
# ╟─d31fb342-d3cc-4fa6-99f5-f3a2ac7c7b42
# ╟─d365f70f-16d2-4d08-accb-0f713c6440c6
# ╠═bf75061a-48f5-437d-900e-aa4e1abf8047
# ╠═680641ef-16b5-4fca-9a57-149bf8c72eb8
# ╠═b1076909-7591-41aa-bada-4747616cec67
# ╠═9a6e62bb-289e-4ca7-9f95-cd0caaac7700
# ╠═630597bd-e8d6-45a8-8e1e-4f2fcc02309c
# ╠═7535da1a-4d29-4866-a212-5f68faad75de
# ╠═afd1ef5f-a8ce-4bff-b001-69c1e2731d85
# ╠═2ee80e6a-ef83-4c78-8b42-0b864405249f
# ╠═5241b7d9-1f56-4b87-bf3a-fdd10cbbb71d
# ╠═94926618-9161-4ea5-bcd1-346ee41380f0
# ╠═39a8c1dc-7a7e-4958-a997-039d2eeb20bd
# ╠═2ba85012-706a-4719-ab87-e0fa91d2f36f
# ╠═5954a7be-70c1-4d23-b4c2-719f23220212
# ╠═e477b05f-2fbb-4bb8-a37a-a56c286445dc
# ╠═7c6dcca9-fd13-4977-8bf7-7f0f5d2b5636
# ╠═0618b3c8-f194-4636-bd42-4923a5eff097
# ╠═f05f9a6c-ccea-497c-bc0e-6a02842cccc6
