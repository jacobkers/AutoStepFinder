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

# ╔═╡ 1d36be1e-12a7-11ef-0d8b-05789779f611
using Pkg; Pkg.activate("C:\\Users\\Alfonso\\MEGA\\PhD\\BOLD\\Data_analysis\\WideField_analysis")

# ╔═╡ 003ae74c-597c-4792-8c84-847942ea28dd
begin
	using Images
	using Plots
	using Plots.PlotMeasures
	using LaTeXStrings
	using PlutoUI
	using Base.Filesystem
	using DelimitedFiles
	using Statistics
	gr();
end

# ╔═╡ 6e448be7-2eda-4243-b260-4b86101575d4
md"""
# Program explanation
"""

# ╔═╡ 297b76bb-abd1-4fff-848a-3a2a3f39ee38
md"""
This program takes the results from the notebooks "Spectrum\_individual\_fields" and "PhotoBleaching", and does the combined analysis of results extracted from them. In particular, this notebook combines all the spectra of the different fields in the same sample and generates the combined plots for that sample. Same applies for the PB curves
"""

# ╔═╡ 0779c935-6d5c-4f4b-9590-89027b39fe67
md"""
## Set the working directory (sample number)
"""

# ╔═╡ f210f763-03f6-4685-8993-579b3e745799
folder_data = raw"C:\Users\Alfonso\BOLD\WideField\Data\april_2024\55"

# ╔═╡ 54f69ebd-37c5-4c0d-bb14-711d9e83193d
begin
	folder_data_spectra = joinpath(folder_data, "Spectrum")
	folder_data_PB = joinpath(folder_data, "PB")
end;

# ╔═╡ 275b1e4d-5fef-40a7-9e35-3ec8d79a029e
md"""Check to store the results in the folder: $(@bind store_results CheckBox())"""

# ╔═╡ 576baf5a-155d-4457-81ed-2712d99a613b
md"""
##### Creation of the results folder (if not existent)
"""

# ╔═╡ 387cafb5-27b0-4321-8144-f94b60a63729
if store_results
	results_folder_name_spectra = joinpath(folder_data_spectra, "Analysis_AYN", "Common")
	if isdir(results_folder_name_spectra)
		println("Spectra directory already exists")
	else
		try
			mkdir(results_folder_name_spectra)
			println("Spectra directory created")
		catch e
			println("Could not create spectra directory")
		end	
	end
	results_folder_name_PB = joinpath(folder_data_PB, "Analysis_AYN", "Common")
	if isdir(results_folder_name_PB)
		println("PB directory already exists")
	else
		try
			mkdir(results_folder_name_PB)
			println("PB directory created")
		catch e
			println("Could not create PB directory")
		end	
	end
end;

# ╔═╡ 1a81dd3d-4abb-480c-b755-d6cc3202bf3e
md"""
### Construction of the combined spectra and average spectrum for the sample
"""

# ╔═╡ a2c19a14-d58a-4635-b3c1-262b0ad47ba8
md"""
# TO DO
"""

# ╔═╡ 6a1ba57a-859f-402d-b6ac-0919274a9bf9
#Add common analysis for PB curves

# ╔═╡ 61999d4f-30d5-4bb3-92bc-8de8a15b9bc9
md"""
# Functions
"""

# ╔═╡ cc1905c7-7172-4dd8-9237-d9ea2214b56d
"""
Sorts the name of the folders containing the fields data by numerical
"""
function sort_fields(dir1, dir2)
	num1 = parse(Int, match(r"Field(\d+)", dir1).captures[1])
	num2 = parse(Int, match(r"Field(\d+)", dir2).captures[1])
	return num1 < num2
end

# ╔═╡ 46d4dd19-6fab-4407-a9bd-2a5818c38518
begin
	all_files = readdir(folder_data_spectra)
	field_list_unordered = filter(file -> occursin(r"^Field\d+$", file), all_files)
	field_list = sort(field_list_unordered, lt=sort_fields)
	#println("Field directories: ", field_list)
end;

# ╔═╡ 1ce5793c-f58f-4998-b2f1-685c91b393e6
begin
	x_values = []
	y_values = []
	x_errors = []
	for field in field_list
		complete_path = joinpath(folder_data_spectra, "Analysis_AYN", field)
		files = readdir(complete_path)
		txt_spectrum_file = filter(file -> contains(file, "spectrum_ROI.txt"), files)
	
		contents = readlines(joinpath(complete_path, "spectrum_ROI.txt"))
		x_values_ind, y_values_ind, x_errors_ind = [[parse(Float64, val) for val in split(line)] for line in contents]
		push!(x_values, x_values_ind)
		push!(y_values, y_values_ind)
		push!(x_errors, x_errors_ind)
	end
	y_average = [mean([y_values[j][i] for j in 1:length(y_values)]) for i in 1:length(y_values[1])]
end;

# ╔═╡ fb1e75a1-0ef2-4c93-ba59-a0fbcb3c8619
begin
	figure_plot_average_spectrum = "spectrum_ROI_average.png"
	p2 = plot(xlabel=L"$\lambda$ (nm)", ylabel=L"$\gamma$/s/$\Delta \lambda$ per pixel", xlabelfontsize = 9, ylabelfontsize = 9, label=false, linecolor=:black, framestyle=:box, framelinewidth=3, xtickfontsize=9, ytickfontsize=9, xticks=400:50:800, size=(600,300), title="Average spectrum", titlefontsize=12, legend=:outertopright, legendfontsize=9, dpi=300)
	
	plot!(p2, x_values[1], y_average, xerror = x_errors[1], marker=:circle, markersize = 5, line=:2.5, label = "$(size(y_values)[1]) spectra")
	p2
end

# ╔═╡ 39edb00a-1b35-4ffd-a9b4-1a271723c67f
if store_results
	savefig(joinpath(results_folder_name_spectra, figure_plot_average_spectrum))
end

# ╔═╡ b8bcedc8-75b7-4419-abe2-b5375ad20cc0
begin
	figure_plot_spectra = "spectra_ROI_combined.png"
	p = plot(xlabel=L"$\lambda$ (nm)", ylabel=L"$\gamma$/s/$\Delta \lambda$ per pixel", xlabelfontsize = 9, ylabelfontsize = 9, label=false, linecolor=:black, framestyle=:box, framelinewidth=3, xtickfontsize=9, ytickfontsize=9, xticks=400:50:800, size=(600,350), title="Combined spectra", titlefontsize=12, legend=:outertopright, legendfontsize=9, dpi=300)
	for i in 1:length(field_list)
		plot!(p, x_values[i], y_values[i], xerror = x_errors[i], marker=:circle, markersize = 5, line=:2.5, label=field_list[i])
	end
	p
end

# ╔═╡ 4db93882-f2ca-4159-b8d7-33140a5dafd0
if store_results
	savefig(joinpath(results_folder_name_spectra, figure_plot_spectra))
end

# ╔═╡ 5dbcd6b1-47c4-4d68-b52a-2eddf1aa452b
md"""
# Other utilities
"""

# ╔═╡ 710b3535-43d8-4eca-bb32-b6855fc22644
PlutoUI.TableOfContents(title="WideField images analysis", indent=true)

# ╔═╡ Cell order:
# ╟─6e448be7-2eda-4243-b260-4b86101575d4
# ╟─297b76bb-abd1-4fff-848a-3a2a3f39ee38
# ╠═1d36be1e-12a7-11ef-0d8b-05789779f611
# ╠═003ae74c-597c-4792-8c84-847942ea28dd
# ╟─0779c935-6d5c-4f4b-9590-89027b39fe67
# ╠═f210f763-03f6-4685-8993-579b3e745799
# ╠═54f69ebd-37c5-4c0d-bb14-711d9e83193d
# ╟─275b1e4d-5fef-40a7-9e35-3ec8d79a029e
# ╟─576baf5a-155d-4457-81ed-2712d99a613b
# ╠═387cafb5-27b0-4321-8144-f94b60a63729
# ╟─1a81dd3d-4abb-480c-b755-d6cc3202bf3e
# ╠═46d4dd19-6fab-4407-a9bd-2a5818c38518
# ╠═1ce5793c-f58f-4998-b2f1-685c91b393e6
# ╠═b8bcedc8-75b7-4419-abe2-b5375ad20cc0
# ╠═4db93882-f2ca-4159-b8d7-33140a5dafd0
# ╠═fb1e75a1-0ef2-4c93-ba59-a0fbcb3c8619
# ╠═39edb00a-1b35-4ffd-a9b4-1a271723c67f
# ╟─a2c19a14-d58a-4635-b3c1-262b0ad47ba8
# ╠═6a1ba57a-859f-402d-b6ac-0919274a9bf9
# ╟─61999d4f-30d5-4bb3-92bc-8de8a15b9bc9
# ╟─cc1905c7-7172-4dd8-9237-d9ea2214b56d
# ╟─5dbcd6b1-47c4-4d68-b52a-2eddf1aa452b
# ╠═710b3535-43d8-4eca-bb32-b6855fc22644
