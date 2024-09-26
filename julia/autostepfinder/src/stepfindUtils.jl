using DelimitedFiles
using Plots
using PlutoUI
using Images
using ImageView
using Colors
using FixedPointNumbers
using Statistics

function AutoStepMain(dataX, tresH = 0.15, N_iter = 100)
    FitX = zeros(size(dataX))  # Assuming dataX is a 2D array
    
    Fits = zeros(0)  # Initialize empty arrays for storage
    S_curves = zeros(0)
    best_shots = []
	residuX = dataX - FitX
	
	newFitX, _, _, S_curve, best_shot = stepfindcore(
            residuX, tresH, N_iter)

	FitX = AppendFitX(newFitX, FitX, dataX)
       

    steptable = fit2Steps(dataX, FitX)  
    
    return S_curve, best_shot, FitX, steptable
end


function plot_pixel(evol, i, j, img, tresH=0.1, N_iter=50, )
	dataX = evol[ i, j, :]
	S_curve, best_shots, Fits, steptable = asf.autostepfinder.AutoStepMain(dataX, tresH, N_iter)
	p1 = plot(dataX, label = "Data pixel ($i,$j)")
	plot!(Fits, label = "AutoStepFinder fit")

	p2 = plot(S_curve)
	p2 = heatmap(img,colorbar=false, c=:grays)# xaxis=false, yaxis=false, )

	x_start = i - 2
	x_end = i + 2
	y_start = j - 2
	y_end = j + 2
	
	# Draw a rectangle on top of the heatmap using Shape
	plot!([x_start, x_end, x_end, x_start, x_start], [y_start, y_start, y_end, y_end, y_start],
	    seriestype = :shape, linecolor = :blue, fillalpha = 0., linealpha = 1.0, linewidth = 4, label="Pixel")
		
	ppx = plot(p1, p2)
	return ppx
end

function get_hist_roi(evol, tresH = 0.1, N_iter = 50)
	hist_height = Vector{Any}()
	hist_widths = Vector{Any}()
	n_steps = Vector{Int}()
	
	I, J, _ = size(evol)
	ijyes = reshape([], 0, 2)
	
	for i in 1:I
	    for j in 1:J
			dataX = evol[ i, j, :]
			
			
			S_curve, best_shots, Fits, steptable = AutoStepMain(dataX, tresH, N_iter)
			
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


function temporal_evolution(folder, N, roi=[1, 512, 1, 512])
    im0 = get_image(folder, 1)[roi[1]:roi[2], roi[3]:roi[4]]
    data = [im0]
    
    for i in 2:N
        d = get_image(folder, i)[roi[1]:roi[2], roi[3]:roi[4]]
        push!(data, d)
    end
    
    evol = cat(data..., dims=3)
    return evol
end

function get_image(folder, n)
	nstr = lpad(string(n), 3, '0')
	files = readdir(folder)
	sort!(files)
	filename = filter(file -> occursin(nstr, file[end-6:end-4]), files)
	fpath = joinpath(folder, filename[end])
	im = Images.load(fpath)
	im64 = reinterpret(UInt16, channelview(im))
	imarray = Float64.(channelview(im64))
 
	return imarray
end