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

# ╔═╡ fadf633b-621f-4f18-bbb5-df18273d2e6c
using Pkg; Pkg.activate("C:\\Users\\Alfonso\\MEGA\\PhD\\BOLD\\Data_analysis\\WideField_analysis")

# ╔═╡ 3eac8dee-8709-4ed1-8fc4-ad163f1f5426
begin
	using Images
	using Plots
	using Plots.PlotMeasures
	using LaTeXStrings
	using PlutoUI
	using Base.Filesystem
	using DelimitedFiles
	using CSV
	using Colors
	using FixedPointNumbers
	using Statistics
	using Revise
	using StatsBase
	using DataFrames
	using SpecialFunctions
	using EasyFit
	gr();
end

# ╔═╡ 3afbafe0-1371-11ef-0114-3b191c3f6904
md"""
# Program explanation
"""

# ╔═╡ 40b08f2f-2bfd-4601-9ae5-18b199fce144
md"""
This notebook takes data from a field on a sample and the ROI selected in the spectrum analysis (notebook 1). Then, it processes the data inside the ROI (binning and deconvolve functions, chosen by the user), applies AutoStepFinder in the arrays containing the temporal evolution of every pixel inside the ROI and stores the results, that can be processed together with other fields with notebook 3.
"""

# ╔═╡ 78f03196-450c-4c51-afc3-5a664b81fef3
md"""
## Set the working directory and required parameters
"""

# ╔═╡ 737287f7-f3a5-454b-9b02-0d87ed5d2f8b
md"""
Select the location of the SAMPLE to study (not the field):
"""

# ╔═╡ 854a1826-13de-4a80-883f-675413a630a6
folder_sample = raw"C:\Users\Alfonso\BOLD\WideField\Data\april_2024\57"

# ╔═╡ 5b03b13c-a813-4916-b13e-08413b370670
begin
	filenames_all = readdir(joinpath(folder_sample, "PB"))
	filenames_filtered = filter(filename -> startswith(filename, "Field"), filenames_all)
end;

# ╔═╡ d1a0bd8f-f373-4df5-824a-254460657f79
md"""
Select the FIELD to analyze from the list:
"""

# ╔═╡ 1ba10129-eb26-443e-b436-d72cf28e3793
@bind field_number Select(filenames_filtered)

# ╔═╡ 75ea8ed5-4bc0-41f1-9b18-c3e4ad8f1b83
md"""
Introduce the acquisition rate for the photobleaching curve (time between frame acquisition in SECONDS):
"""

# ╔═╡ 40d8b473-96c1-4952-afe8-583cabcb2e7c
@bind acquisition_rate_string TextField(default = "0.5")

# ╔═╡ 18231a39-2bd9-4bef-b443-05f7261acb88
md"""
Select maximum number of active pixels to store (to avoid overloading):
"""

# ╔═╡ 3891abd0-8802-4944-ae74-f13f42e2f12a
@bind max_pixels_to_store TextField(default = "250")

# ╔═╡ 2e0da7f9-e4ee-43a4-a746-4a8be8d60a5b
md"""
Select the threshold value (TresH) for AutoStepFinder algorithm:
"""

# ╔═╡ ba6ed8fe-b97b-41b8-9e10-ea0e85aa50ab
@bind tresH_string TextField(default = "0.4")

# ╔═╡ 95c930a5-84ba-4670-a0b9-ab145751ad0a
md"""
Select the number of iterations (N_iter) for AutoStepFinder:
"""

# ╔═╡ 43dc348a-3896-4e48-b4eb-44983f200006
@bind N_iter_string TextField(default = "50")

# ╔═╡ 4ed4288f-d077-4fd0-af75-29a4d823f254
md"""
Check to analyze all images of field (time consuming) and store results: $(@bind temp_evol CheckBox())
"""

# ╔═╡ 5481b4eb-bd28-4700-ad34-7c425c08b014
md"""
Check to store the trajectories of all active pixels: $(@bind store_trajectories CheckBox())
"""

# ╔═╡ 3bc06e36-0692-4f85-83f8-a1a38c0d9ace
md"""
## Read ROI from spectrum analysis
"""

# ╔═╡ fd5dc8fa-3a84-4e1e-98ae-b6c2020223ea
begin
	tresH = parse(Float64, tresH_string)
	N_iter = parse(Float64, N_iter_string)
	folder_data = joinpath(folder_sample, "PB", field_number)
end;

# ╔═╡ 6adf724f-73ad-424d-a0cb-2869ce73c3d8
begin
	spectra_folder = joinpath(dirname(dirname(folder_data)), "Spectrum")
	if isdir(joinpath(spectra_folder, "Analysis_AYN", field_number))
		if isfile(joinpath(spectra_folder, "Analysis_AYN", field_number, "ROI.txt"))
			roi_string = read(joinpath(spectra_folder, "Analysis_AYN", field_number, "ROI.txt"), String)
			roi = [parse(Int, val) for val in split(roi_string)]
		else
			println("Spectra analysis folder does not exist, analyze the spectrum before PB")
		end
	else
		println("Spectra analysis folder does not exist, analyze the spectrum before PB")
	end
end;

# ╔═╡ 47c48ee0-ca8e-4f58-8fc8-17489d707ff0
md"""
Check to perform data processing before applying AutoStepFinder algorithm: $(@bind preprocess CheckBox())
"""

# ╔═╡ 6ed429a2-41f6-4fa4-8ecb-b39742aa4547
if preprocess
	println("Set the value of binsize, lambda and aperture diameter")
	binsize = 4
	lambda = 450
	aperture_diameter = 2
end;

# ╔═╡ db148fe2-143a-49f1-a4bc-10cf7599ea01
md"""
Set parameters for AutoStepFinder:
"""

# ╔═╡ 0bef7c1e-c14d-4455-951b-d30e206dbdfb
begin
	data = readdir(folder_data)
	num_files = count(endswith(".tif"), data) - 1
end;

# ╔═╡ 34a8690e-8ccd-4b4a-8547-b6b73ead33a9
md"""
# Sets of frames analysis. Statistics
"""

# ╔═╡ a189b589-cc79-4736-a345-0a3d28f23cca
md"""
# Pixel trajectories (individual and average)
"""

# ╔═╡ 0063fa47-f32f-4302-bd66-55e6e0a05a03
md"""
##### Fit of the photobleaching to an exponential
"""

# ╔═╡ 8754abf0-df80-4c2f-9ad7-84350e5161ee
seconds_to_fit = 50; #temporal limit of the fit in seconds

# ╔═╡ bbbf0b41-746c-4f01-a1e9-dee2bdbeab48
md"""
##### Plots
"""

# ╔═╡ ee9a6aa6-f51b-4295-b2ab-f5bbd360a6ef
md"""
##### Average trajectory in ROI
"""

# ╔═╡ ca58450a-8f5b-4b48-8dff-f370a7d640f9
md"""
##### Creation of the results folder (if not existent) and storage of results
"""

# ╔═╡ 21e6bd2b-ee27-4db4-aa5f-09e95deab2c4
md"""# Functions"""

# ╔═╡ 0bc1224b-5a4f-4f6e-a668-0dcec4a1fc7e
"""
Imports the .jl documents from the given location
"""
function Ingredients(path::String)
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

# ╔═╡ b3fafc65-d893-43e9-9818-021e60f85b77
asf = Ingredients("..\\src\\AutoStepFinder.jl");

# ╔═╡ 2a5ee74c-e17b-4088-b4d2-0c3f51dcb50b
"""
Returns an array with the data from image n in the file folder
"""
function Get_Image(folder, frame)
	nstr = lpad(string(frame), 5, '0')
	files = readdir(folder)
	sort!(files)
	filename = filter(file -> occursin(nstr, file[end-9:end-4]), files)
	full_path = joinpath(folder, filename[end])
	#print(full_path) #debug, prints the resulting path of the image to be acquired
	image = Images.load(full_path)
	image_data = reinterpret(UInt16, channelview(image))
	image_array = Float64.(channelview(image_data))
	return image_array
end	

# ╔═╡ 904eab41-e58e-4625-8090-4fa39956be19
begin
	acquisition_rate = parse(Float64, acquisition_rate_string)
	frame = 10
	figure_heatmap = "heatmap_frame_$(frame).png"
	img = Get_Image(folder_data, frame)
	heatmap(img, c=:inferno, aspect_ratio=:equal, yflip=false, size=(600, 600))
	plot!([roi[1], roi[2], roi[2], roi[1], roi[1]], [roi[3], roi[3], roi[4], roi[4], roi[3]], seriestype=:shape, linecolor=:yellow, fillalpha=0., linealpha=0.9, linewidth=3, title=field_number*" at t = $(frame * acquisition_rate) s", titlefontsize = 12, label=false, widen=false, tickdirection=:out, margin=2mm)
end

# ╔═╡ 137fbbfa-a807-4137-a3bc-f52899ea6e2f
"""
Main function that uses the core processes from AutoStepFinder to return S_curve, best_shots, final fits and the table with steps
"""
function AutoStepMain(data_to_analyze, tresH, N_iter)
	Fit = zeros(size(data_to_analyze))
	newFit, _, _, S_curve, best_shot = asf.AutoStepFinder.stepfindcore(data_to_analyze, tresH, N_iter)
	Fit = asf.AutoStepFinder.AppendFitX(newFit, Fit, data_to_analyze)
	step_table = asf.AutoStepFinder.fit2Steps(data_to_analyze, Fit)
	return S_curve, best_shot, Fit, step_table
end

# ╔═╡ 5211e6b5-b445-4a67-980e-70a2a86cbc31
"""
Processes the outputs from AutoStepMain and builds histograms with the resulting steps for the pixels that pass the algorithm (active_pixels)
"""
function Get_Histograms_ROI(temporal_evolution, tresH, N_iter)
	hist_heights = Vector{Any}()
	hist_widths = Vector{Any}()
	n_steps = Vector{Int}()

	I, J = size(temporal_evolution[1]) #Size of the ROI, taken from the first frame
	active_pixels = reshape([], 0, 2)
	for i in 1:I
		for j in 1:J
			pixel_evolution = [temporal_evolution[k][i, j] for k in 1:length(temporal_evolution)] #We construct an array containing the intensity data of the same pixel over all the frames (temporal evolution of the pixel)
			S_curve, best_shot, Fit, steptable = AutoStepMain(pixel_evolution, tresH, N_iter)

			if best_shot > 0 && S_curve[best_shot] > tresH

				LST = size(steptable)[1]
				if LST > 2
					push!(n_steps, size(steptable, 1))
					push!(hist_heights, steptable[1:LST, 4])
					push!(hist_widths, cumsum(steptable[1:LST, 5]))
					active_pixels = [active_pixels; [j i]]
				end
			end
		end
	end
	hw = vcat(hist_widths...)
	hh = vcat(hist_heights...)
	sorting_indices = sortperm(active_pixels[:, 1])
	sorted_active_pixels = active_pixels[sorting_indices, :]
	return hh, hw, n_steps, sorted_active_pixels
end

# ╔═╡ 627ff63b-2e22-4b04-a81d-3816764159c9
"""
Calls the AutoStepFinder analysis functions, plots the results and returns plots and data
"""
function Analyze_Image(temporal_evolution, tresH)
	hh, hw, n_steps, active_pixels = Get_Histograms_ROI(temporal_evolution, tresH, N_iter);
	res_array = hh, hw, n_steps, active_pixels
	
	#if all(n_steps .> 0)
	#	p1 = histogram(n_steps, label="Number of steps")
	#	p2 = histogram(hh, label="Step height")
	#	p3 = histogram(hw, label="Step width")
	#	phs = plot(p1, p2, p3, layout=(1,3))
	#else
	#	@warn "No steps detected"
	#end
	p1 = histogram(n_steps, color=:blue, label = false, xlabel="Number of steps", xlabelfontsize = 10)
	p2 = histogram(hh, color=:green, label = false, xlabel="Step height", xlabelfontsize = 10)
	p3 = histogram(hw, color=:orange, label = false, xlabel="Step width", xlabelfontsize = 10)
	phs = plot(p1, p2, p3, layout=(1,3))
	return phs, res_array
end

# ╔═╡ bcdd5938-b934-4aa8-a416-40a70acb8540
"""
Stores results in a CSV file located in CSV_path
"""
function Store_Dataframe(active_pixels, n_steps, hh, hw, CSV_path)
	ijsave = hcat(inverse_rle(active_pixels[:,1], n_steps), inverse_rle(active_pixels[:,2], n_steps)) #We are just interested in active pixels, so we create two columns with its coordinates
	
	dfs = DataFrame(x=ijsave[:,1], y=ijsave[:,2], height=hh, width=hw) #We assign to the active pixels their corresponding step data

	CSV.write(CSV_path, string.(dfs), header=true, append=true) #Creation of the CSV file
end

# ╔═╡ 478008c7-5e3c-486b-a8f8-ddc7daaaeecf
"""
Saves all the relevant results (CSV file through Store_Dataframe and graphs)
"""
function Save_Results(path, tresH, res_array, phs)
	dirroot, fieldN = splitdir(path) #As path to store the results, we will use the same as that of the field data
	dirres1 = joinpath(dirroot, "Analysis_AYN") #In case the field analyzed is the first in the folder, we create the results directory to use for all fields
	try
		mkdir(dirres1)
	catch
		@warn "Directory already exists"
	end
	dirres2 = joinpath(dirres1, fieldN) #We create a path to store the results with the number of the field analyzed
	println(dirres2) #debug
	try
		mkdir(dirres2)
	catch
		@warn "Directory already exists. Field already analyzed"
	end
	
	filename = fieldN*".csv" #There will be a CSV file created for every field
	path_results = joinpath(dirres2, filename);
	print(path_results) #debug
	
	header = " ", " ", "ROI position: $roi", "Threshold $tresH"
	CSV.write(path_results, [header], header=false)

	hh, hw, n_steps, active_pixels = res_array
	Store_Dataframe(active_pixels, n_steps, hh, hw, path_results)
	dirsave, filename = splitdir(path_results)

	savefig(phs, joinpath(dirsave, fieldN*"hist.png")) #We also save the histogram figure of the field
	return dirsave
end

# ╔═╡ f9962802-25fb-434e-8eaa-baae842af938
"""
Reads the temporal trajectories of the pixels passing the algorithm and fits them to exponential decays to obtain the photobleaching constant: tau
"""
function PhotoBleaching_Fit(temporal_evolution, active_pixels, seconds_to_fit)
	fit_results_all = []
	for i in 1:size(active_pixels)[1]
		pixel_evolution = [temporal_evolution[k][active_pixels[i, 2], active_pixels[i, 1]] for k in 1:length(temporal_evolution)]

		x = [a * acquisition_rate for a in 1:Int(seconds_to_fit/acquisition_rate)]
		y = pixel_evolution[1:Int(seconds_to_fit/acquisition_rate)]
		fit_results = fitexp(x, y)

		push!(fit_results_all, fit_results)
	end
	return fit_results_all
end

# ╔═╡ 25c3886c-d10f-4918-888c-71755f73e193
md"""
# CURRENTLY IN DEVELOPMENT
"""

# ╔═╡ ff932c8a-ad5c-4748-96d8-d9987ea7c6b3
"""
Takes the data of a certain frame of the selected field and the size of the binning, then bins the data for the posterior analysis
"""
function Binning(unbinned_image, binsize)
	binned_image = zeros(div(size(unbinned_image, 1), binsize), div(size(unbinned_image, 2), binsize))
	for i in 1:binsize:size(unbinned_image, 1)
		for j in 1:binsize:size(unbinned_image, 2)
			print("$j ")
			binned_image[div(i-1, binsize) + 1, div(j-1, binsize) + 1] = sum(unbinned_image[i:min(i + binsize - 1, end), j:min(j + binsize - 1, end)])
		end
	end
	return binned_image
end

# ╔═╡ 1154c4ca-65d6-411a-9865-45d0cf37a902
"""
Creates and returns an array of matrices (temp_evol) containing the N successive temporal images of the roi of the field in the folder given
"""
function Temporal_Evolution(folder, N, roi)
	if preprocess
		println("Data will be processed before building the temporal evolution")
		temp_evol = []
		for i in 1:N
			processed_frame = Binning(Get_Image(folder, i)[roi[1]:roi[2] - 1, roi[3]:roi[4] - 1], binsize)
			#DECONVOLVE NEEDS TO BE REVISED
			#binned_frame = Binning(Get_Image(folder, i)[roi[1]:roi[2] - 1, roi[3]:roi[4] - 1], binsize)
			#processed_frame = Deconvolve_PSF(binned_frame, lambda, aperture_diameter)
			push!(temp_evol, processed_frame)
		end
	else
		println("Data will NOT be processed to build temporal evolution")
		temp_evol = [Get_Image(folder, i)[roi[1]:roi[2] - 1, roi[3]:roi[4] - 1] for i in 1:N]
	end
	return temp_evol
end

# ╔═╡ 4896bcfe-fc29-4343-9433-7729dd70d2d0
evol_try = Temporal_Evolution(folder_data, 30, roi);

# ╔═╡ 073859fe-d83d-40bb-8cdb-ec1850e8ed1b
if temp_evol
	figure = [heatmap(evol_try[i], c=:inferno, aspect_ratio=:equal, yflip=false, size=(600, 600)) for i=1:1];
end

# ╔═╡ a9570315-f705-4c41-8309-52d2ad8e0d9f
plot(figure...)

# ╔═╡ dd338429-8ea0-4404-9af3-74a133ac88d2
if temp_evol
	println("Analyzing data in $folder_data")
	temporal_evolution = Temporal_Evolution(folder_data, num_files, roi)
	plots, results_array = Analyze_Image(temporal_evolution, tresH);
	dirsave = Save_Results(folder_data, tresH, results_array, plots)
	plot(plots)
end

# ╔═╡ 71f63fe1-bb8c-43ae-bb33-526199edc4ca
print("Algorithm passed by $(size(results_array[4], 1)) pixels")

# ╔═╡ 4060e35e-61f0-4121-83b8-df90fdca2408
photobleaching_fits = PhotoBleaching_Fit(temporal_evolution, results_array[4], seconds_to_fit);

# ╔═╡ f3b34417-8138-4948-a158-e38f169cfae0
begin
	tau_values = [fit.b for fit in photobleaching_fits]
	average_tau = mean(tau_values)
	println("Average tau = $average_tau s")
end

# ╔═╡ bddb3753-1f78-4cc7-87a8-f291564382b0
begin
	average_trajectory = []
	for (i, frame) in enumerate(temporal_evolution)
		counts_per_px_frame = sum(frame)/(size(frame)[1]*size(frame)[2])
		push!(average_trajectory, counts_per_px_frame)
	end
end

# ╔═╡ dfe9dbd4-6175-4a54-baf4-bfcb746d127f
begin
	x = [a * acquisition_rate for a in 1:length(average_trajectory)]
	y = Float64.(average_trajectory[1:length(average_trajectory)])
	x_fit = [a * acquisition_rate for a in 1:Int(seconds_to_fit/acquisition_rate)]
	y_fit = Float64.(average_trajectory[1:Int(seconds_to_fit/acquisition_rate)])
	fit_average = fitexp(x_fit, y_fit)
end;

# ╔═╡ d9240721-f5cd-43c8-bb29-8454a660643f
begin
	p_avg = plot(x, y, label = "Average trajectory", xlabel = "t (s)", ylabel = "Pixel intensity value (a.u.)", xlabelfontsize = 10, ylabelfontsize = 10, linewidth=2)
	plot!(x_fit, fit_average.ypred, linewidth=2, label = L"\tau = " * string(round(fit_average.b, digits = 1)) * " s")
end

# ╔═╡ 2f4b2a85-0aba-41ac-8c15-f9fafa26f35b
if store_trajectories && size(results_array[4])[1]< max_pixels_to_store
	dirroot = dirname(folder_data)
	results_folder_name = joinpath(dirroot, "Analysis_AYN")
	if isdir(results_folder_name)
		println("Parent directory already exists")
	else
		try
			mkdir(results_folder_name)
			println("Parent directory created")
		catch e
			println("Could not create directory")
		end	
	end
	subdirectory_field = joinpath(results_folder_name, field_number)
	if isdir(subdirectory_field)
		println("Subdirectory already exists")
	else
		try
			mkdir(subdirectory_field)
			println("Subdirectory created")
		catch e
			println("Could not create subdirectory")
		end
	end
	directory_trajectories = joinpath(subdirectory_field, "Trajectories")
	println(directory_trajectories)
	if isdir(directory_trajectories)
		println("Trajectory directory already exists")
	else
		try
			mkdir(directory_trajectories)
			println("Trajectory directory created")
		catch e
			println("Could not create trajectory directory")
		end
	end
else
	if size(results_array[4])[1] >= max_pixels_to_store
		println("More than $max_pixels_to_store active pixels. Set a higher value for the threshold or a higher value for \"max_pixels_to_store\"")
	end
end

# ╔═╡ e920f086-b0d4-4aa6-9e69-33bf24d137fa
begin
	savefig(p3, joinpath(subdirectory_field, figure_heatmap));
end

# ╔═╡ da4e08b8-b89e-4b70-8ef7-11dcecd18863
begin
	savefig(p_avg, joinpath(subdirectory_field, "Average_trajectory.png"));
	open(joinpath(subdirectory_field, "average_trajectory.txt"), "w") do text_file
		for value in y
		println(text_file, value)
		end
	end
end

# ╔═╡ 14e2baff-6528-4f2f-a355-2b9279650144
"""
Analyzes and shows the desired pixel in the ROI and plots its temporal evolution together with the results of the AutoStepFinder analysis for it
"""
function Plot_Pixel_Trajectory(pixel, i, j, frame, tresH, N_iter, fit_results) #i, j are the coordinates of the pixel to analyze
	pixel_evolution = [temporal_evolution[k][j, i] for k in 1:length(temporal_evolution)]
	
	S_curve, best_shot, Fit, step_table = AutoStepMain(pixel_evolution, tresH, N_iter)
	
	p1 = plot([i * acquisition_rate for i in 1:length(temporal_evolution)], pixel_evolution, label = "Evolution of ($i, $j)", xlabel = "t (s)", ylabel = "Pixel intensity value (a.u.)", xlabelfontsize = 10, ylabelfontsize = 10, linewidth=2)

	if best_shot > 0 && S_curve[best_shot] > tresH
		plot!([a * acquisition_rate for a in 1:length(temporal_evolution)], Fit, label = "ASF Fit", linewidth=2)
		plot!([b * acquisition_rate for b in 1:length(fit_results[pixel].ypred)], fit_results[pixel].ypred, linewidth=2, label = L"\tau = " * string(round(fit_results[pixel].b, digits = 1)) * " s")
	end

	p2 = heatmap(frame, colorbar = false, c=:inferno, aspect_ratio=:equal, yflip=false)

	xmin = i - 2
	xmax = i + 2

	ymin = j - 2
	ymax = j + 2
	plot!([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], seriestype =:shape, linecolor =:white, fillalpha = 0., linealpha = 1.0, linewidth = 2, label="Pixel ($i, $j)")

	final_plot = plot(p1, p2)

	return final_plot
end

# ╔═╡ 7120159c-df13-4140-9add-f52c43530112
begin
	pixel = 1
	i = 44
	j = 77
	Plot_Pixel_Trajectory(pixel, i, j, temporal_evolution[1], tresH, N_iter, photobleaching_fits)
end

# ╔═╡ e3e900d6-3696-4b5d-b5f2-a74644dd0975
begin
	plots_active_pixels = []
	for pixel in 1:size(results_array[4])[1]
		push!(plots_active_pixels, Plot_Pixel_Trajectory(pixel, results_array[4][pixel, 1], results_array[4][pixel, 2], temporal_evolution[1], tresH, N_iter, photobleaching_fits))
	end
end

# ╔═╡ ecdc24c0-b4c4-43b8-8f08-3ac97e035686
plots_active_pixels[1]

# ╔═╡ 9af3cd4e-a16c-4f09-949d-5ac515e2f327
if store_trajectories
	for ii in 1:length(plots_active_pixels)
		subplot = plot(plots_active_pixels[ii])
		savefig(joinpath(directory_trajectories, "Pixel_$ii.png"))
		
	end
	#Pixel trajectory to store in txt file
	for pixel in 1:size(results_array[4])[1]
		ii = results_array[4][pixel, 1]
		jj = results_array[4][pixel, 2]
		pixel_trajectory = [temporal_evolution[kk][jj, ii] for kk in 1:length(temporal_evolution)]
		txt_filename = joinpath(directory_trajectories, "Pixel_$(ii)_$jj.txt")
		open(txt_filename, "w") do file
			for y_value in pixel_trajectory
				println(file, y_value)
			end
		end				
	end
end

# ╔═╡ 0b098090-12ab-4f8b-a9dc-3aaddb4d3b49
#CAREFUL HERE WITH THE ORDER OF I,J. CHECK

# ╔═╡ 971825d2-e532-4778-95dc-5092f2f94e19
"""
Defines the intensity of a pixel at a distance r due to the Airy disk for a given wavelength and aperture diameter
"""
function Airy_Disk(r, lambda, aperture_diameter)
	k = 2 * pi / lambda
	intensity = (2 * besselj(1, k * r) / (k * r))^2
	return intensity
end

# ╔═╡ 9b327272-4ee9-428c-9782-eb9dc1878fac
"""
Takes an image to deconvolve, a wavelength and the diameter of the used aperture and performs the deconvolution using the Point Spread Function of the setup, returning the resulting image
"""
function Deconvolve_PSF(image_to_deconvolve, lambda, aperture_diameter)
	width, height = size(image_to_deconvolve)
	PSF = [Airy_Disk(norm([i, j]), lambda, aperture_diameter) for i in 1:width, j in 1:height] #Calculus of the Point Spread Function of an optical setup
	PSF /= sum(PSF) #Normalization of PSF
	deconvolution = imfilter(image_to_deconvolve, PSF)
	return deconvolution
end

# ╔═╡ 5885fb2d-442f-4a62-800e-68145bf5ef09
md"""
# Other utilities
"""

# ╔═╡ a61eb8ba-4923-45d5-8b9d-8b122f3d7266
PlutoUI.TableOfContents(title="WideField images analysis", indent=true)

# ╔═╡ Cell order:
# ╟─3afbafe0-1371-11ef-0114-3b191c3f6904
# ╟─40b08f2f-2bfd-4601-9ae5-18b199fce144
# ╠═fadf633b-621f-4f18-bbb5-df18273d2e6c
# ╠═3eac8dee-8709-4ed1-8fc4-ad163f1f5426
# ╠═b3fafc65-d893-43e9-9818-021e60f85b77
# ╟─78f03196-450c-4c51-afc3-5a664b81fef3
# ╟─737287f7-f3a5-454b-9b02-0d87ed5d2f8b
# ╠═854a1826-13de-4a80-883f-675413a630a6
# ╟─5b03b13c-a813-4916-b13e-08413b370670
# ╟─d1a0bd8f-f373-4df5-824a-254460657f79
# ╟─1ba10129-eb26-443e-b436-d72cf28e3793
# ╟─75ea8ed5-4bc0-41f1-9b18-c3e4ad8f1b83
# ╟─40d8b473-96c1-4952-afe8-583cabcb2e7c
# ╟─18231a39-2bd9-4bef-b443-05f7261acb88
# ╟─3891abd0-8802-4944-ae74-f13f42e2f12a
# ╟─2e0da7f9-e4ee-43a4-a746-4a8be8d60a5b
# ╟─ba6ed8fe-b97b-41b8-9e10-ea0e85aa50ab
# ╟─95c930a5-84ba-4670-a0b9-ab145751ad0a
# ╟─43dc348a-3896-4e48-b4eb-44983f200006
# ╟─4ed4288f-d077-4fd0-af75-29a4d823f254
# ╟─5481b4eb-bd28-4700-ad34-7c425c08b014
# ╟─3bc06e36-0692-4f85-83f8-a1a38c0d9ace
# ╠═fd5dc8fa-3a84-4e1e-98ae-b6c2020223ea
# ╠═6adf724f-73ad-424d-a0cb-2869ce73c3d8
# ╠═904eab41-e58e-4625-8090-4fa39956be19
# ╠═e920f086-b0d4-4aa6-9e69-33bf24d137fa
# ╟─47c48ee0-ca8e-4f58-8fc8-17489d707ff0
# ╠═6ed429a2-41f6-4fa4-8ecb-b39742aa4547
# ╟─db148fe2-143a-49f1-a4bc-10cf7599ea01
# ╠═0bef7c1e-c14d-4455-951b-d30e206dbdfb
# ╟─34a8690e-8ccd-4b4a-8547-b6b73ead33a9
# ╠═4896bcfe-fc29-4343-9433-7729dd70d2d0
# ╠═073859fe-d83d-40bb-8cdb-ec1850e8ed1b
# ╠═a9570315-f705-4c41-8309-52d2ad8e0d9f
# ╠═dd338429-8ea0-4404-9af3-74a133ac88d2
# ╟─a189b589-cc79-4736-a345-0a3d28f23cca
# ╠═71f63fe1-bb8c-43ae-bb33-526199edc4ca
# ╟─0063fa47-f32f-4302-bd66-55e6e0a05a03
# ╠═8754abf0-df80-4c2f-9ad7-84350e5161ee
# ╠═4060e35e-61f0-4121-83b8-df90fdca2408
# ╠═f3b34417-8138-4948-a158-e38f169cfae0
# ╟─bbbf0b41-746c-4f01-a1e9-dee2bdbeab48
# ╠═7120159c-df13-4140-9add-f52c43530112
# ╠═e3e900d6-3696-4b5d-b5f2-a74644dd0975
# ╠═ecdc24c0-b4c4-43b8-8f08-3ac97e035686
# ╟─ee9a6aa6-f51b-4295-b2ab-f5bbd360a6ef
# ╠═bddb3753-1f78-4cc7-87a8-f291564382b0
# ╠═dfe9dbd4-6175-4a54-baf4-bfcb746d127f
# ╠═d9240721-f5cd-43c8-bb29-8454a660643f
# ╟─ca58450a-8f5b-4b48-8dff-f370a7d640f9
# ╠═2f4b2a85-0aba-41ac-8c15-f9fafa26f35b
# ╠═da4e08b8-b89e-4b70-8ef7-11dcecd18863
# ╠═9af3cd4e-a16c-4f09-949d-5ac515e2f327
# ╟─21e6bd2b-ee27-4db4-aa5f-09e95deab2c4
# ╠═0bc1224b-5a4f-4f6e-a668-0dcec4a1fc7e
# ╠═2a5ee74c-e17b-4088-b4d2-0c3f51dcb50b
# ╠═1154c4ca-65d6-411a-9865-45d0cf37a902
# ╠═627ff63b-2e22-4b04-a81d-3816764159c9
# ╠═5211e6b5-b445-4a67-980e-70a2a86cbc31
# ╠═137fbbfa-a807-4137-a3bc-f52899ea6e2f
# ╠═bcdd5938-b934-4aa8-a416-40a70acb8540
# ╠═478008c7-5e3c-486b-a8f8-ddc7daaaeecf
# ╠═14e2baff-6528-4f2f-a355-2b9279650144
# ╠═f9962802-25fb-434e-8eaa-baae842af938
# ╟─25c3886c-d10f-4918-888c-71755f73e193
# ╠═ff932c8a-ad5c-4748-96d8-d9987ea7c6b3
# ╠═0b098090-12ab-4f8b-a9dc-3aaddb4d3b49
# ╠═971825d2-e532-4778-95dc-5092f2f94e19
# ╠═9b327272-4ee9-428c-9782-eb9dc1878fac
# ╟─5885fb2d-442f-4a62-800e-68145bf5ef09
# ╠═a61eb8ba-4923-45d5-8b9d-8b122f3d7266
