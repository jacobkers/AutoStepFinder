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

# ╔═╡ 54dcc270-008c-11ef-1c5a-f33cdd664b5d
using Pkg; Pkg.activate("C:\\Users\\Alfonso\\MEGA\\PhD\\BOLD\\Data_analysis\\WideField_analysis")

# ╔═╡ f7e7aa83-1e2e-420d-9227-1fe1c477766f
begin
	using Images
	using Plots
	using Plots.PlotMeasures
	using LaTeXStrings
	using PlutoUI
	using Base.Filesystem
	using DelimitedFiles
	gr();
end

# ╔═╡ c8cbecab-b14f-452b-a7a5-aad6a2ef6b79
md"""
# Program explanation
"""

# ╔═╡ 7f7ffff4-f895-41b1-84ba-cd1df934394b
md"""
This program calculates the intensity received by the CMOS camera when a sample is illuminated by a laser. This intensity is measured independently in different wavelengths $\lambda$ (nm). The program takes into account the exposure time to the laser, the width of each band pass (BP) filter, the quantum efficiency of the camera at each studied $\lambda$, the number of pixels in the image and the conversion between counts and photons, provided by the camera documentation. This way, for each field acquired we can construct spectra with $\lambda$ in the x-axis and photons per pixel, second and $\Delta \lambda$ in the y-axis.
"""

# ╔═╡ afa580a1-e05d-4b8e-889b-58efd9c9fefb
md"""
## Set the working directory and exposure
"""

# ╔═╡ 43e68012-f4f5-48d8-a38a-9e21d4bfa5d3
md"""
Select the loaction of the SAMPLE to study (not the field):
"""

# ╔═╡ df68d895-b454-4938-8d3f-0d96d5688d7e
begin
	folder_sample = raw"C:\Users\Alfonso\BOLD\WideField\Data\april_2024\55_2"
end;

# ╔═╡ a02b6cc7-eb55-48b6-96f9-22c0d70f1c5c
begin
	filenames_all = readdir(joinpath(folder_sample, "Spectrum"))
	filenames_filtered = filter(filename -> startswith(filename, "Field"), filenames_all)
end;

# ╔═╡ 27ed7df7-4c88-4b68-afee-1a39c1ff4222
md"""
Select the FIELD to analyze from the list:
"""

# ╔═╡ 6b308d30-0d5f-4184-9836-72438aebe1c1
@bind field_number Select(filenames_filtered)

# ╔═╡ a85ef4be-1155-4d77-a160-ebe7375e5fc9
md"""
Introduce the exposure time for the spectra (exposure time per filter in SECONDS):
"""

# ╔═╡ 37208205-e4a9-4b7f-8e43-98a3f02778a5
@bind exposure_time_string TextField(default = "0.5")

# ╔═╡ 2fb14587-7ab5-4179-a2c4-b3d0081229a4
md"""Check to store the results in the folder: $(@bind store_results CheckBox())"""

# ╔═╡ e2958294-7734-4923-93c6-718c9e19327b
md"""
Check to perform independent analysis of background (use only for testing): $(@bind background_analysis CheckBox())
"""

# ╔═╡ da468a8d-9dbd-48eb-a23c-88b0c3aab6db
md"""
### Select the ROI
"""

# ╔═╡ 79a6d6bc-1ced-421e-8a8a-bc10943f7131
roi = [170, 270, 190, 290];

# ╔═╡ 5da7fc25-0962-47c7-aec5-8ecfaae1a2a3
md"""
##### Creation of the results folder (if not existent)
"""

# ╔═╡ a26427cb-8436-420c-9614-b0bd76c900c7
begin
	folder_data = joinpath(folder_sample, "Spectrum", field_number)
	exposure_time = parse(Float64, exposure_time_string)
end;

# ╔═╡ a46802e1-63cf-4471-a3d7-129e995c3c56
if store_results
	results_folder_name = joinpath(dirname(folder_data), "Analysis_AYN")
	if isdir(results_folder_name)
		println("Common directory already exists")
	else
		try
			mkdir(results_folder_name)
			println("Common directory created")
		catch e
			println("Could not create directory")
		end	
	end
	results_folder_name_field = joinpath(results_folder_name, basename(folder_data))
	if isdir(results_folder_name_field)
		println("Directory for field already exists")
	else
		try
			mkdir(results_folder_name_field)
			println("Field directory created")
		catch e
			println("Could not create field directory")
		end
	end
end;

# ╔═╡ 1309ae2a-2941-4f7e-b920-1bf5ed1dd7da
md"""
 ##### Definition of the bandpas filters used in the setup
"""

# ╔═╡ 891d150a-7f46-4cfe-966a-fdb2ddae76cf
begin
	filter_wavelengths = zeros(11, 3)
	filter_wavelengths[1,:]=[426,438,450]
    filter_wavelengths[2,:]=[457,472,487]
    filter_wavelengths[3,:]=[500,510,520]
    filter_wavelengths[4,:]=[524,534,544]
    filter_wavelengths[5,:]=[554,561,568]
    filter_wavelengths[6,:]=[576,586,596]
    filter_wavelengths[7,:]=[605,615,625]
    filter_wavelengths[8,:]=[633,640,647]
    filter_wavelengths[9,:]=[662,676,691]
    filter_wavelengths[10,:]=[696,716,736]
    filter_wavelengths[11,:]=[752,775,798]
	band_widths = filter_wavelengths[:,3] - filter_wavelengths[:,1]
end;

# ╔═╡ 9dca06b6-1b12-4976-a834-f3a895b557de
md"""
##### Definition of transmission factors and efficiencies
"""

# ╔═╡ 35bafaf6-3414-4e2f-9175-6e2e8fdc561f
begin
	cmos_qe=[0.57, 0.71, 0.79, 0.81, 0.82, 0.82, 0.81, 0.80, 0.76, 0.68, 0.56]; #orca hamamatsu c13440 v3. This is the efficiency at the central wavelengths of the 11 BP filters used
	t_bpfilters=[0.97, 0.98, 0.97, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.985, 0.984] #edmund optics 11 filters
	t_tubelens=[0.9467, 0.9512, 0.9568, 0.9627, 0.9674, 0.9694, 0.9691, 0.9701, 0.9761, 0.9848, 0.9720] #Thorlabs TTL200-A
	t_hpfilter=[0.975, 0.977, 0.974, 0.975, 0.975, 0.980, 0.980, 0.981, 0.981, 0.980, 0.980] #Semrock LP02 407RU 25
	t_notch=[0.945, 0.942, 0.950, 0.943, 0.936, 0.9405, 0.948, 0.952, 0.948, 0.945, 0.9445] #Semrock NF03-405E-25
	t_dichroic=[0.978, 0.978, 0.976, 0.981, 0.981, 0.980, 0.972, 0.979, 0.977, 0.984, 0.982] #Semrock FF414-DI01-25x36
end;

# ╔═╡ 178625e6-97fb-4f6a-b658-763cf62be062
transmission_correction = t_bpfilters .* t_tubelens .* t_hpfilter .* t_notch .* t_dichroic;

# ╔═╡ f56c5c52-ce9b-4a12-8f78-e165b7872932
md"""
##### Other relevant parameters
"""

# ╔═╡ b6a31e53-c16d-4a7b-ab2b-c659731f8297
begin
	binning_x = 4 #binning factor applied directly in the camera
	binning_y = 4
end;

# ╔═╡ ff5e992d-fdb8-4897-8a77-a162f57b35db
md"""
### Reading of the images and construction of the spectrum for the field
"""

# ╔═╡ 0bec5a9e-6510-453d-863a-e21341a0645a
begin
	filenames = readdir(folder_data)
	folder_background = joinpath(folder_data, "Background")
	background_filenames = readdir(folder_background)
	
	image_list = filter(file -> occursin(".tiff", file[end-4:end]), filenames)
	image_white_light = image_list[end]
	image_list = image_list[1:end-1]

	background_list = filter(file -> occursin(".tiff", file[end-4:end]), background_filenames)
	background_white_light = background_list[end]
	background_list = background_list[1:end-1]
end;

# ╔═╡ b65cf2a2-b85a-482c-872f-84c894908b18
md"""
##### Creation of spectrum and ROI zoom
"""

# ╔═╡ 1a31138b-f047-4694-a69f-a5cdaf56b2f1
#CONVERSION FROM PX TO LENGTH UNITS WOULD BE GREAT

# ╔═╡ 1c9d6617-47c0-4631-abbc-d7a5cfc90534
if store_results
	open(joinpath(results_folder_name_field,"ROI.txt"), "w") do file
		for value in roi
			println(file, value)
		end
	end
end

# ╔═╡ f2c09158-1187-4239-bc3e-be0ab2ca55ea
#We need to discuss if we should apply the corrections pixelwise or to the average, because it could make a difference in the spectra. Also, if we should apply the corrections to the background values too or just to the signal data
corrections_background = true;

# ╔═╡ a21feaec-6c6e-41b3-8ac7-e60c62993de7
md"""
### Plot of the spectrum
"""

# ╔═╡ 89c22890-73a4-4f41-95be-65b2d4dc9e39
md"""
# Functions
"""

# ╔═╡ 32efd1c5-8d47-4041-a673-4e5419febe62
"""
Returns an array with the data from image n in the file folder
"""
function Get_Image(folder, frame_name)
	#nstr = lpad(string(frame), 2, '0')
	files = readdir(folder)
	sort!(files)
	#filename = filter(file -> occursin(nstr, file[end-7:end-5]), files)
	full_path = joinpath(folder, frame_name)
	#print(full_path) #debug
	image = Images.load(full_path)
	image_data = reinterpret(UInt16, channelview(image))
	image_array = Float64.(channelview(image_data))
	return image_array
end

# ╔═╡ 4391085f-4b8b-4929-9a00-d57990261f37
begin
	figure_heatmap_full = "heatmap_full.png"
	img = Get_Image(folder_data, image_white_light)
	heatmap(img, c=:inferno, aspect_ratio=:equal, yflip=false, size=(600, 600))
	plot!([roi[1], roi[2], roi[2], roi[1], roi[1]], [roi[3], roi[3], roi[4], roi[4], roi[3]], seriestype=:shape, linecolor=:yellow, fillalpha=0., linealpha=0.9, linewidth=3, title=field_number*" (white light)", titlefontsize = 12, label=false, widen=false, tickdirection=:out, margin=2mm)
end

# ╔═╡ 7d1b9dd9-4b1e-43e6-ba1a-9ec4ea76b8cd
if store_results
	savefig(joinpath(results_folder_name_field, figure_heatmap_full))
end;

# ╔═╡ 85c7843b-5ee6-4120-af20-7be2eea8ec2f
begin
	figure_heatmap_roi = "heatmap_ROI.png"
	heatmap(img[roi[3]:roi[4], roi[1]:roi[2]], c=:inferno, aspect_ratio=:equal, yflip=false, size=(350, 350), title=field_number*" ROI", titlefontsize = 12, label=false, widen=false, tickdirection=:out, margin=5mm)
end

# ╔═╡ 76dab4f4-c22f-48e1-9ebb-04ee35669dff
if store_results
	savefig(joinpath(results_folder_name_field, figure_heatmap_roi))
end;

# ╔═╡ 90361bfd-7c3d-4f7c-b413-996f77185447
if background_analysis
	imgbkg = Get_Image(folder_background, background_white_light)
	heatmap(imgbkg, c=:inferno, aspect_ratio=:equal, yflip=false, size=(600, 600))
end

# ╔═╡ a8a5cc7c-74f8-4833-bd6f-f1047f892c5f
if background_analysis
	for field in background_filenames
		imagge = Get_Image(folder_background, field)
		println(sum(imagge)/(size(imagge)[1]*size(imagge)[2]))
	end
end

# ╔═╡ 43dd9792-7c8a-4846-b8f8-8b7f1d8a5224
begin
	data_raw = zeros(Float64, size(image_list)) #to store the analysis results for each field
	data_background = zeros(Float64, size(image_list))
	data_sub_background = zeros(Float64, size(image_list))
	
	for i in 1:length(image_list)
		image = Get_Image(folder_data, image_list[i])
		image_roi = image[roi[3]:roi[4], roi[1]:roi[2]]
		
		background_image = Get_Image(folder_background, background_list[i])
		background_image_roi = background_image[roi[3]:roi[4], roi[1]:roi[2]]
		
		total_counts = sum(image_roi)
		counts_per_pixel = total_counts / (size(image_roi)[1] * size(image_roi)[2])
		data_raw[i] = counts_per_pixel
	
		total_counts_bg = sum(background_image_roi)
		counts_per_pixel_bg = total_counts_bg / (size(background_image_roi)[1] * size(background_image_roi)[2])
		data_background[i] = counts_per_pixel_bg

		if(counts_per_pixel > counts_per_pixel_bg)
			data_sub_background[i] = 0.47 * (counts_per_pixel - counts_per_pixel_bg) / cmos_qe[i] / band_widths[i] / transmission_correction[i] / exposure_time   #not all the BP filters have the same width so we need to take that into account
		else
			data_sub_background[i] = 0.47 * (counts_per_pixel - 100 * binning_x * binning_y) / cmos_qe[i] / band_widths[i] / transmission_correction[i] / exposure_time
		end
		#In this last step, we consider the conversion from counts to photons provided by the documentation: photons = (counts - 100 * binning^2) / QE, but we try to consider the background counts as the threshold instead of a fixed value of 100 * binning ^ 2.
	end
end;

# ╔═╡ b579696e-0a6f-4329-8c75-d9d2af09e52d
begin
	figure_plot_spectrum = "spectrum_ROI.png"
	plot(filter_wavelengths[:,2], data_sub_background, xerror = band_widths/2, xlabel=L"$\lambda$ (nm)", ylabel=L"$\gamma$/s/$\Delta \lambda$ per pixel", xlabelfontsize = 9, ylabelfontsize = 9, label=false, line=:2, linecolor=:black, marker=:circle, markersize = 6, framestyle=:box, framelinewidth=2, xtickfontsize=9, ytickfontsize=9, xticks=400:50:800, size=(375,250), title=field_number, titlefontsize=12, dpi=300)
end

# ╔═╡ 8c232bba-3d2a-4cbd-b0ce-5c0dab83b741
if store_results
	savefig(joinpath(results_folder_name_field, figure_plot_spectrum))
	open(joinpath(results_folder_name_field, "spectrum_ROI.txt"), "w") do file
		println(file, join(filter_wavelengths[:,2], "\t"))
		println(file, join(data_sub_background, "\t"))
		println(file, join(band_widths/2, "\t"))
	end
end;

# ╔═╡ 7f38a4ea-9c72-4b29-9c42-d22d588fd9af
md"""
# Other utilities
"""

# ╔═╡ f6db4370-db7c-4fcb-b872-09ddccf8811a
PlutoUI.TableOfContents(title="WideField images analysis", indent=true)

# ╔═╡ Cell order:
# ╟─c8cbecab-b14f-452b-a7a5-aad6a2ef6b79
# ╟─7f7ffff4-f895-41b1-84ba-cd1df934394b
# ╠═54dcc270-008c-11ef-1c5a-f33cdd664b5d
# ╠═f7e7aa83-1e2e-420d-9227-1fe1c477766f
# ╟─afa580a1-e05d-4b8e-889b-58efd9c9fefb
# ╟─43e68012-f4f5-48d8-a38a-9e21d4bfa5d3
# ╠═df68d895-b454-4938-8d3f-0d96d5688d7e
# ╟─a02b6cc7-eb55-48b6-96f9-22c0d70f1c5c
# ╟─27ed7df7-4c88-4b68-afee-1a39c1ff4222
# ╟─6b308d30-0d5f-4184-9836-72438aebe1c1
# ╟─a85ef4be-1155-4d77-a160-ebe7375e5fc9
# ╟─37208205-e4a9-4b7f-8e43-98a3f02778a5
# ╟─2fb14587-7ab5-4179-a2c4-b3d0081229a4
# ╟─e2958294-7734-4923-93c6-718c9e19327b
# ╟─da468a8d-9dbd-48eb-a23c-88b0c3aab6db
# ╠═79a6d6bc-1ced-421e-8a8a-bc10943f7131
# ╟─4391085f-4b8b-4929-9a00-d57990261f37
# ╟─5da7fc25-0962-47c7-aec5-8ecfaae1a2a3
# ╠═a26427cb-8436-420c-9614-b0bd76c900c7
# ╠═a46802e1-63cf-4471-a3d7-129e995c3c56
# ╟─1309ae2a-2941-4f7e-b920-1bf5ed1dd7da
# ╠═891d150a-7f46-4cfe-966a-fdb2ddae76cf
# ╟─9dca06b6-1b12-4976-a834-f3a895b557de
# ╠═35bafaf6-3414-4e2f-9175-6e2e8fdc561f
# ╠═178625e6-97fb-4f6a-b658-763cf62be062
# ╟─f56c5c52-ce9b-4a12-8f78-e165b7872932
# ╠═b6a31e53-c16d-4a7b-ab2b-c659731f8297
# ╟─ff5e992d-fdb8-4897-8a77-a162f57b35db
# ╠═0bec5a9e-6510-453d-863a-e21341a0645a
# ╠═90361bfd-7c3d-4f7c-b413-996f77185447
# ╟─b65cf2a2-b85a-482c-872f-84c894908b18
# ╠═1a31138b-f047-4694-a69f-a5cdaf56b2f1
# ╠═1c9d6617-47c0-4631-abbc-d7a5cfc90534
# ╠═a8a5cc7c-74f8-4833-bd6f-f1047f892c5f
# ╠═7d1b9dd9-4b1e-43e6-ba1a-9ec4ea76b8cd
# ╠═85c7843b-5ee6-4120-af20-7be2eea8ec2f
# ╠═76dab4f4-c22f-48e1-9ebb-04ee35669dff
# ╠═f2c09158-1187-4239-bc3e-be0ab2ca55ea
# ╠═43dd9792-7c8a-4846-b8f8-8b7f1d8a5224
# ╟─a21feaec-6c6e-41b3-8ac7-e60c62993de7
# ╠═b579696e-0a6f-4329-8c75-d9d2af09e52d
# ╠═8c232bba-3d2a-4cbd-b0ce-5c0dab83b741
# ╟─89c22890-73a4-4f41-95be-65b2d4dc9e39
# ╟─32efd1c5-8d47-4041-a673-4e5419febe62
# ╟─7f38a4ea-9c72-4b29-9c42-d22d588fd9af
# ╠═f6db4370-db7c-4fcb-b872-09ddccf8811a
