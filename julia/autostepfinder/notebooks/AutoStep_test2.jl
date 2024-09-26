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
using Pkg; Pkg.activate(ENV["JLaserLab"])

# ╔═╡ 40687fb5-d226-46e5-8a01-125adc668924
# ╠═╡ show_logs = false
Pkg.add("FixedPointNumbers")

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
	#using GR
end

# ╔═╡ cce4b21e-81e7-4e6c-ad7c-014f9f1d72ea
ENV["JLaserLab"]

# ╔═╡ ee654cd8-b63a-4e5e-96f8-3104bea4fa60
#Pkg.instantiate()

# ╔═╡ 933774f1-9a00-4ee0-98e0-980f9e8fb930
#Pkg.precompile()

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
asf = ingredients("../../../src/LaserLab.jl")

# ╔═╡ d3da2ac5-df51-4d7e-99f5-039595f6f75e
PlutoUI.TableOfContents(title="AutoStepfinder Notebook", indent=true)

# ╔═╡ c5133014-5af9-42a0-a04b-4bf1317a2170
md"""# Test with single data file"""

# ╔═╡ 5ad0425b-f976-4fab-873b-cb1949b15d09
begin
	file = "/Users/pabloherrero/sabat/LaserLab/Data-analysis/AutoStepFinder/python/test_data/test_photobleaching.txt"
	dataf = vec(readdlm(file, ' '))
end

# ╔═╡ 23670ee0-8137-4fec-b6b4-86e04b59c324
# ╠═╡ show_logs = false
S_curvef, best_shotsf, Fitsf, steptablef = asf.LaserLab.AutoStepMain(dataf, 0.15, 100)

# ╔═╡ b2b06a68-dcd5-43e1-9649-8574423e30ec
md"""# Test with image"""

# ╔═╡ 9dd0f232-5668-4e57-b7f5-529d23613f49
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

# ╔═╡ 0b04b3af-13cf-48a6-b8b1-800679db9e7d
function temporal_evolution(folder, N, roi=[1, 512, 1, 512])
    im0 = get_image(folder, 1)[roi[1]:roi[2], roi[3]:roi[4]]
    data = [im0]
    
    for i in 2:N
        d = get_image(folder, i)[roi[1]:roi[2], roi[3]:roi[4]]
        push!(data, d)
    end
    
    return cat(data..., dims=3)
end

# ╔═╡ 74183602-063d-4aa9-884e-aff77e9606d9
path = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/ZFFL62_NAPH3_5e-8M/Microscope/Wide field/no_barium/35/Field1"

# ╔═╡ f94daecf-afda-4b61-a5d1-3ae48c035480
begin
	img = get_image(path, 23)
	heatmap(img)
end

# ╔═╡ 3f6a3fe5-14f1-4f7a-ae9b-34c8d76f4ed4
img

# ╔═╡ a0e4327f-8529-4011-9b36-01cdca20a4da
typeof(img)

# ╔═╡ 10b8263d-8dc5-4ddc-bd21-b310ac7fe11b
begin
	roi1 = [100,220,100,220]
	imroi = get_image(path, 23)[roi1[1]:roi1[2],roi1[3]:roi1[4]]
	gr()
	#c = cgrad([:red,:yellow,:green], [0.01, 0.9995], categorical = true)

	heatmap(imroi, clim=(1800,2500))
	#display(fig)
end

# ╔═╡ 12b432f4-5942-4958-9de8-473a08c364a5
evol = temporal_evolution(path, 400, roi1);

# ╔═╡ c1b54e5e-d712-45e3-b752-4cb78c618067


# ╔═╡ 45db938c-4934-47d9-8d02-5e40542258de
hist_height

# ╔═╡ 967fd807-dca0-4f64-bca8-e3cc1fd7c3f6
size(evol)

# ╔═╡ 2e0177fb-47fa-48c7-b81b-3942c2d7794f
typeof(evol[3, 1, :])

# ╔═╡ 0772692d-47fa-460d-8bb6-3d397af3f792
# ╠═╡ show_logs = false
begin
	#float_vector = Float64.(reinterpret(Float16, collect(evol[90, 80, :])))
	S_curve2, best_shots2, Fits2, steptable2 = asf.LaserLab.AutoStepMain(evol[15, 10, :], 0.05, 100)
end

# ╔═╡ d20f7145-b292-4cd8-9129-9ddfda4088f7
begin
	plot(evol[15, 10, :])
	plot!(Fits2)
end

# ╔═╡ 4a902f18-561e-491f-974a-e4e533414d5d
# ╠═╡ disabled = true
#=╠═╡

  ╠═╡ =#

# ╔═╡ 302416ac-ea8b-4e0c-8307-0b5e31f298b9
dataX = evol[12, 13, :]

# ╔═╡ 75ac27d2-372a-470f-b2a1-34a2d3b1253d
typeof(dataX), size(dataX)

# ╔═╡ e80e3660-fa5e-4489-8bbc-4712ca585d9a
plot(dataX)

# ╔═╡ 2e727245-96ea-4a34-991b-e5129f3ef0d4


# ╔═╡ cdea0164-5a87-4b41-ab6b-db4fba83b0c4
function A1210_set_up_splitlogtable(dataX)
    Na = length(dataX)
    i_nxt, avl, avr, rankit, errorcurve = asf.LaserLab.splitFast(dataX)  # Assuming st_splitFast has an equivalent in Julia
    new_row1 = [1 Na i_nxt avl avr rankit]
    split_table = [
        -1 -1 -1 1 1 1
        new_row1
        Na Na Na 1 1 1
    ]
    return split_table
end

# ╔═╡ a319b59c-f3b3-4319-984d-95a3e29acc3b
split_table0 = A1210_set_up_splitlogtable(dataX)

# ╔═╡ 8c4e261f-4e1d-478d-ae9b-d6426c0e86be


# ╔═╡ 95fca961-dc50-4e6a-908e-ab6fe6a568da
function A1211_adapt_fit(FitX, split_table, idx)
    rank = split_table[idx, 6]
	#aFitX = similar(FitX)
    ok2adapt = rank > 0
    if ok2adapt
        leftplateau_istart = Int(split_table[idx, 1])
        leftplateau_istop = Int(split_table[idx, 3])
        rightplateau_istart = leftplateau_istop + 1
        rightplateau_istop = Int(split_table[idx, 2])
        leftplateau_level = split_table[idx, 4]
        rightplateau_level = split_table[idx, 5]
        FitX[leftplateau_istart:leftplateau_istop] .= leftplateau_level
        FitX[rightplateau_istart:rightplateau_istop] .= rightplateau_level
    end
    return FitX
end

# ╔═╡ 66e3ad55-c1aa-4cf5-9fcb-54e2a7199e6a

function A1213_adapt_counterfit(cFitX, dataX, split_table, idx)
    rankings = split_table[idx:idx + 1, 6]
    ok2adapt = all(rankings .!= 0)

    if ok2adapt
        leftplateau_istart = Int(split_table[idx - 1, 3]) + 2
        leftplateau_istop = Int(split_table[idx, 3])
        centerplateau_istart = leftplateau_istop + 1
        centerplateau_istop = Int(split_table[idx + 1, 3])
        rightplateau_istart = centerplateau_istop + 1
        rightplateau_istop = Int(split_table[idx + 2, 3])

        leftplateau_level = mean(dataX[leftplateau_istart:leftplateau_istop])
        centerplateau_level = mean(dataX[centerplateau_istart:centerplateau_istop])
        rightplateau_level = mean(dataX[rightplateau_istart:rightplateau_istop])

        cFitX[leftplateau_istart:leftplateau_istop] .= leftplateau_level
        cFitX[centerplateau_istart:centerplateau_istop] .= centerplateau_level
        cFitX[rightplateau_istart:rightplateau_istop] .= rightplateau_level

		
    end
    return cFitX
end

# ╔═╡ 3d84586d-c4c5-4aff-a480-cb81e827436e


# ╔═╡ 22527abe-a572-421d-a2de-a8aeb0e53fc3
begin
	FitXf = zeros(size(dataX)) .+ mean(dataX)
    cFitXf = similar(FitXf)

    i_startf = Int(split_table0[2,1])
    i_stopf = Int(split_table0[2,2])
    i_nextf = Int(split_table0[2,3])
    avlf = split_table0[2,4]
    avrf = split_table0[2,5]
    cFitXf[i_startf:i_nextf] .= avlf
    cFitXf[i_nextf:i_stopf] .= avrf 
end 

# ╔═╡ 85c62eb4-01b4-47a7-8744-f4da00a75689
cFitXf

# ╔═╡ a9a8ef42-e733-4b25-9b9d-af8d756f5720
avlf, avrf, i_nextf, i_stopf, i_startf

# ╔═╡ dd89c023-7f16-4d05-9610-802ae7534be9
function splitFast(segment)
    invalidFit = 0
    Nmin = 2
    Ns = length(segment)
    if Ns > 2
        var_q = ones(Float64, Ns)

        avl = segment[1]
        avr = sum(segment[2:end]) / (Ns - 1)
        ava = sum(segment) / Ns

        for ii in 2:Ns-1
            n_L = ii
            n_R = Ns - ii
            avl = (avl * (n_L - 1) + segment[ii]) / n_L
            avr = (avr * (n_R + 1) - segment[ii]) / n_R

            delta_l = avl - ava
            delta_r = avr - ava

            varcor_L = 1 + 0 * (n_L + 1) / n_L
            varcor_R = 1 + 0 * (n_R + 1) / n_R

            delta_q = delta_l^2 * n_L * varcor_L + delta_r^2 * n_R * varcor_R
            var_q[ii] = -delta_q
        end

        idx = argmin(var_q)

        if (idx < Nmin - 1) || (Ns - idx < Nmin - 1)
            invalidFit = 1
			println("Entered weird condition")
        else
            avl_fin = mean(segment[1:idx])
            avr_fin = mean(segment[idx + 1:end])
            rankit = (avr_fin - avl_fin)^2 * Ns
            errorcurve = var_q / Ns

        end
    else
        invalidFit = 1
    end

    if invalidFit == 1
		#println("invalidFit!!!!")
		idx = 0
        avl_fin = segment[1]
        avr_fin = segment[1]
		
        rankit = 0.
		
        errorcurve = zeros(Float64, length(segment))
        
    end

	
    return idx, avl_fin, avr_fin, rankit, errorcurve
end

# ╔═╡ b3d5892c-d0db-40ea-a3f4-13a076e60e0a
idx, avl_fin, avr_fin, rankit, error_curve = splitFast(dataX)

# ╔═╡ 7a241f86-d20a-4f56-a87a-94aad6c19ef4
function A1212_adapt_splitlog_table(segm=1, split_table=1, idx=2)
    best_row_entry = split_table[idx, :]
	
	LST = size(split_table)[1]
	last_row_entry = split_table[LST, :]
	
    istart_1 = Int(best_row_entry[1])
    istop_1 = Int(best_row_entry[3])

	
	if length(segm[istart_1:istop_1]) > 2
		inxt_1, avl_1, avr_1, rankit_1, errorcurve_1 = splitFast(segm[istart_1:istop_1])
	else
		inxt_1, avl_1, avr_1, rankit_1, errorcurve_1 = splitFast(segm[istart_1])
	end
    new_row1 = [istart_1 istop_1 inxt_1 + istart_1 avl_1 avr_1 rankit_1]

    istart_2 = Int(best_row_entry[3] + 1)
    istop_2 = Int(best_row_entry[2])

	if istart_2 == Int(last_row_entry[1]) + 1
		istart_2 = Int(best_row_entry[3])
	end
	
	if length(segm[istart_2:istop_2]) > 2
	    inxt_2, avl_2, avr_2, rankit_2, errorcurve_2 = splitFast(segm[istart_2:istop_2])
	
	else 
		inxt_2, avl_2, avr_2, rankit_2, errorcurve_2 = splitFast(segm[istart_2])
	end
	new_row2 = [istart_2 istop_2 inxt_2 + istart_2 avl_2 avr_2 rankit_2]
	
    splitlog_entry = [istop_1 inxt_1 + istart_1 inxt_2 + istart_2]

    block_before = split_table[1:idx-1, :]
    block_after = split_table[idx+1:end, :]
	
    new_split_table = [block_before
						new_row1
						new_row2
						block_after]
    return new_split_table, splitlog_entry
end

# ╔═╡ 1fb74702-1e7a-4011-82b0-d0f5035ee670
function A121_split_until_ready(dataX::Vector{Float64}, N_iter::Int64=50)
    split_table = A1210_set_up_splitlogtable(dataX)
    splitlog = zeros(Int, N_iter, 3)

    FitX = zeros(size(dataX)) .+ mean(dataX)
    cFitX = similar(FitX)

    i_start = Int(split_table[2,1])
    i_stop = Int(split_table[2,2])
    i_next = Int(split_table[2,3])
    avl = split_table[2,4]
    avr = split_table[2,5]
    cFitX[i_start:i_next] .= avl
    cFitX[i_next+1:i_stop] .= avr
		
    S_curve = zeros(Float64, N_iter)
    c = 1
	
	#plot!(p, dataX, label = "Data X", xlabel = "time, a.u.", ylabel = "position")

	
	for ii in 1:N_iter
		
        segment_lengths = vec(split_table[:, 2] .- split_table[:, 1] .+ 1)
        rankings = vec(split_table[:, 6])
		
        ix_valids = findall((segment_lengths .> 2) .& (rankings .> 0))
		
        if !isempty(ix_valids)
            subsel_idx = argmax(rankings[ix_valids])
            best_row_idx = ix_valids[subsel_idx]
            best_row_entry = vec(split_table[best_row_idx, :])
			
            FitX = A1211_adapt_fit(FitX, split_table, best_row_idx)
            split_table, splitlog_entry = A1212_adapt_splitlog_table(dataX, split_table, best_row_idx)

            cFitX = A1213_adapt_counterfit(cFitX, dataX, split_table, best_row_idx)

			splitlog[c, :] = splitlog_entry
            corF = c / (c + 1)
            corF = 1
            S_curve[c] = corF * mean((dataX - cFitX).^2) / mean((dataX - FitX).^2)
            c += 1

			#fig = plot(layout=(2,1))
			
			#plot!(p, FitX, label = "Fit X")
			#plot!(p, cFitX, label = "Counter Fit X")
			#title!("Intermediate Stepfit Result")
			#aspect_ratio!(p1, 1)
			
			#plot(p, S_curve, label = "S-curve", xlabel = "no. of steps in fit", ylabel = "S-value, a.u.")
			#title!("S-curve Evaluation")
			#aspect_ratio!(p2, 0.7)
			
			#plot(p1, p2, layout = (1, 2))  # Arrange plots side by side
			#display(fig)
			
            #end
        end
    end

    return S_curve, splitlog, split_table
end

# ╔═╡ ba55be13-7a28-4dde-b786-3b8b8641a25f


# ╔═╡ 9a91475c-ed86-4453-8268-5a0fa08c4db6
S_curve0, splitlog0, split_tableF = A121_split_until_ready(dataX, 15)
	#best_shot0, S_curve_fin0 = asf.LaserLab.A122_eval_Scurve(S_curve0, 0.001)

# ╔═╡ 072342a2-8544-47a7-b04a-077945375013
A1213_adapt_counterfit(cFitXf, dataX, split_tableF, 2)

# ╔═╡ 7bf4f8b2-1743-43fa-be90-2a66c485de9c
split_tableF

# ╔═╡ db08a7ae-8e9f-4856-ac87-30097b8a2ac9
begin
 	lst = size(split_tableF)[1]
	split_tableF[lst, :]
end

# ╔═╡ f16c1ce5-a2a9-4d48-a97b-031745b0f176
splitlog0

# ╔═╡ 5f0327d8-fe44-4e92-86e2-62ca4e0192e5
plot(S_curve0)

# ╔═╡ 1e7467d4-416b-4d82-b355-a0c85dde0d48
best_shot0, S_curve_finf = asf.LaserLab.A122_eval_Scurve(S_curve0, 0.05)

# ╔═╡ 43258a32-6754-4a0f-8d6c-2a2c97f1f847
splitlog0

# ╔═╡ 474c21d6-279d-4364-b7ba-840de1f3025a
indices_bestfit, indices_counterfit = asf.LaserLab.Splitlog2FinalIndices(splitlog0, best_shot0)  

# ╔═╡ afba4fac-0a23-4f5e-9f39-53be42930675
begin 
	FitXff = asf.LaserLab.Indices2Fit(dataX, indices_bestfit, "mean")  
	cFitXff = asf.LaserLab.Indices2Fit(dataX, indices_counterfit, "mean")
end

# ╔═╡ bf590b23-ce9e-443f-9087-fa41cad24564
begin
	plot(dataX, label="data")
	plot!(FitXff, label="fit")
	#plot!(cFitXff)
end

# ╔═╡ 8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
md"""# Test large set of data"""

# ╔═╡ 29589bee-6a0c-4f8c-ae19-e98c837d827e

function A122_eval_Scurve(S_curve, acceptance=0.15)
    ix_low = findall(S_curve .< 1)
    S1 = copy(S_curve)
    S1[ix_low] .= 1

    LS = length(S_curve)
    baseline = LinRange(1, S_curve[end], LS)
    S_curve_fin = S1 .- baseline

    i_max = argmax(S_curve_fin)
    if S_curve_fin[i_max] > acceptance
        best_shot = i_max
    else
        best_shot = 0
    end

    return best_shot, S_curve_fin
end

# ╔═╡ 2482b515-f159-4823-b4d3-8447972c5324
function Splitlog2FinalIndices(splitlog, best_shot)
	
	#println("best_shot: ", best_shot)
	#println("splitlog shape", size(splitlog))
    indices_bestfit = sort(splitlog[1:best_shot, 1])
    indices_counterfit_all = vcat(splitlog[1:best_shot, 2], splitlog[1:best_shot, 3])
    is_it_used = zeros(length(indices_counterfit_all))

    for thisindex in indices_bestfit
        ix = findall(indices_counterfit_all .== thisindex)
        is_it_used[ix] .= 1
    end

    indices_counterfit = sort(indices_counterfit_all[is_it_used .== 0])
    return indices_bestfit, indices_counterfit
end

# ╔═╡ 9300f184-7b6e-45dc-b5ea-ab0481a03f0b

function fit2Steps(dataX, FitX)
    Lx = length(dataX)
    T = collect(1:Lx)
    globalnoise = std(diff(dataX .- FitX)) / sqrt(2)
    ixes0 = findall(diff(FitX) .!= 0)
    Lix = length(ixes0)
    ixes = vcat([0], ixes0, [Lx])
    steptable = zeros(Float64, Lix, 8)
	
    for ii in 2:Lix-1
        ix_pre = ixes[ii - 1] + 1
        ix = ixes[ii]
        ix_aft = ixes[ii + 1]
        lev_pre = FitX[ix]
        lev_aft = FitX[ix + 1]
        step = FitX[ix + 1] - FitX[ix]
        dwell_pre = ix - ix_pre + 1
        dwell_aft = ix_aft - ix
        error_pred = 2 * (globalnoise^2 / dwell_pre + globalnoise^2 / dwell_aft)^0.5 / sqrt(2)

        rms_pre = std(dataX[ix_pre:ix - 1])
        rms_aft = std(dataX[ix:ix_aft - 1])
		
		error_meas = 2 * ((rms_pre^2 / dwell_pre + rms_aft^2 / dwell_aft)^0.5) / sqrt(2)
        new_row_entry = [
            ix, lev_pre, lev_aft, step, dwell_pre, dwell_aft, error_pred, error_meas
        ]

        steptable[ii - 1, :] = new_row_entry
    end
    return steptable
end

# ╔═╡ c14801a5-ebd8-4873-854e-2fe51091c1bd

function AppendFitX(newFitX, FitX, dataX)
    # combine different fit rounds

    combiFitX = FitX + newFitX
    Lx = length(combiFitX)
    ixes0 = findall(diff(combiFitX) .!= 0)
    if !isempty(ixes0)
        # pad with start and end
        ixes = [0; ixes0; Lx]
        # find indices too close together
        Nmin = 2
        whereblips = findall(diff(ixes) .< Nmin)

        # redo a stepfit over this part to pick just one step location
        for ix in whereblips
			#println(FitX)
			#println("newFitX: ", newFitX)
            lo = ix - 1
            ixlo = ixes[lo]
            ixhi = ixes[lo] + 3

			if ixlo < Lx && ixhi < Lx
				#println("Lx: ", Lx)
				#println("whereblips: ", whereblips)
				#println("ixlo: ", ixlo)
				#println("ixhi: ", ixhi)
	            segment = dataX[ixlo + 2:ixhi]
	            idx, avl, avr, rankit, errorcurve = splitFast(segment)
	            combiFitX_old = copy(combiFitX)
	            combiFitX[ixlo + 2:ixlo + idx + 2] .= avl
	            combiFitX[ixlo + idx + 3:ixhi + 1] .= avr
			end
        end
    end
    return combiFitX
end


# ╔═╡ 2b44eb05-8735-45c0-aebf-631cd4aa26ec


# ╔═╡ 27a53d31-2aa9-4ab8-9c60-7c969c81f0f8

function stepfindcore(dataX::Vector{Float64}, tresH=0.15, N_iter=0, debug = false)
   
    # run an iterative step fit deep into overfitting:
    if N_iter == 0 || N_iter > length(dataX) ÷ 4
        N_iter = length(dataX) ÷ 4
    end
    S_curve, splitlog, split_table = A121_split_until_ready(dataX, N_iter) 
	println(split_table)
    # estimate best step number in this iteration range:
    best_shot, S_curve_fin = A122_eval_Scurve(S_curve, tresH)  

    # re-build best fit from 'split log':
    if best_shot > 0
		
        indices_bestfit, indices_counterfit = Splitlog2FinalIndices(splitlog, best_shot)  
		#println(indices_bestfit, indices_counterfit)
        FitX = asf.LaserLab.Indices2Fit(dataX, indices_bestfit, "mean")  
        cFitX = asf.LaserLab.Indices2Fit(dataX, indices_counterfit, "mean")

		if debug
			p1 = plot(dataX)
			plot!(p1, FitX, label="Fit")
			plot!(p1, cFitX, label="Counterfit")

			p2 = plot(S_curve_fin, label="S curve")
			plot(p1, p2)
		end
    else
        FitX = zeros(length(dataX))
        cFitX = zeros(length(dataX))
        best_shot = 0
    end

    return FitX, cFitX, splitlog, S_curve_fin, best_shot
end

# ╔═╡ f05b1661-54d5-498f-b906-547a1a420e3c
function AutoStepMain(dataX, tresH = 0.15, N_iter = 100, debug = false)

    FitX = zeros(size(dataX))  # Assuming dataX is a 2D array
    
    Fits = zeros(0)  # Initialize empty arrays for storage
    S_curves = zeros(0)
    best_shots = []
	residuX = dataX - FitX
	
	newFitX, _, _, S_curve, best_shot = stepfindcore(
            residuX, tresH, N_iter, debug)

	FitX = AppendFitX(newFitX, FitX, dataX)
       

    steptable = fit2Steps(dataX, FitX)  
    
    return S_curve, best_shot, FitX, steptable
end


# ╔═╡ d258cbc5-2ff0-495c-af5b-72d0ef26e9e2
begin 
	current_fit, current_S_curve = Array{Float64}[], Array{Float64}[]
	counts = Float64[]
	tables = Matrix{Float64}[]
end

# ╔═╡ 1529d78b-5d37-4021-8c09-a5a103f2e26f
@bind advance_epoch PlutoUI.Button("Advance epoch 🖐")

# ╔═╡ 846e5743-1008-466a-8381-69fbfbf7e3b8


# ╔═╡ 194c2282-d507-4226-b333-bc7c97426a67
S_curve, best_shots, Fits, steptable = AutoStepMain(evol[85,10,:], 0.05, 20, true)

# ╔═╡ 5fe17c70-0e67-4bd9-874e-bfe9d9401629
typeof(steptable)

# ╔═╡ 3c1bf9bc-ebf2-43f0-883c-fc94fcdd9a81
steptablet = fit2Steps(evol[12,13,:], Fits)  

# ╔═╡ 1f9d0a19-2e19-4efe-94da-68e60cfe17db
S_curve22, splitlog22, split_table22 = A121_split_until_ready(evol[2, 20, :], 20) 


# ╔═╡ 5c362934-b868-44b1-8879-790339887a8f
best_shot22, S_curve_fin22 = asf.LaserLab.A122_eval_Scurve(S_curve22, 0.05)  

# ╔═╡ ff2827f6-1f7e-4a97-95b6-4d81cffd245c
begin
	plot(dataX, label="data")
	plot!(FitXff, label="fit")
	#plot!(cFitXff)
end

# ╔═╡ 6604220f-7f61-4126-8ef1-f78087a7cc70
size(S_curve22)

# ╔═╡ e6d83e3b-03c6-4061-b39b-888ef5db649a
size(steptablet)[1]

# ╔═╡ df73f9c4-c7bc-4ba2-b34f-adeecba48ede
typeof(evol)

# ╔═╡ 177600ef-163d-4fbf-b178-88e0f628e9e3
function get_hist_roi(evol, tresH = 0.1, N_iter = 50)
	hist_height = Vector{Any}()
	hist_widths = Vector{Any}()
	n_steps = Vector{Int}()
	
	I, J, _ = size(evol)
	
	for i in 1:I
	    for j in 1:J
			dataX = evol[ i, j, :]
			#println("i, j: ", i, ", ", j)
			
			S_curve, best_shots, Fits, steptable = AutoStepMain(dataX, tresH, N_iter)
			
	            if best_shots > 0 && S_curve[best_shots] > tresH 

					LST = size(steptable)[1]
					if LST > 2
						push!(n_steps, size(steptable, 1))
						push!(hist_height, steptable[1:LST-2, 4])
						push!(hist_widths, steptable[1:LST-2, 5])
					end
	                    
	            end

			end
	    end
	
	hw = vcat(hist_widths...)
	hh = vcat(hist_height...)
	return hh, hw, n_steps
end

# ╔═╡ 78f8ff68-ec80-4988-a240-064b1e2bc93c
hh, hw, n_steps = get_hist_roi(evol, 0.05, 40)

# ╔═╡ 4cbbadc4-9948-4a21-ab03-18ed7c73573d
begin 
	p1 = histogram(n_steps, label="Number of steps")
	p2 = histogram(hh, label="Step height")
	p3 = histogram(hw, label="Step width")
	plot(p1, p2, p3, layout=(1,3))
end

# ╔═╡ c8dd0947-b6d3-4123-9a8d-a9c8fa029290
begin
	p11 = plot(dataf, label="Data pixel", xlabel="Nr of frames", ylabel="position (arb. units)")
	plot!(p1, Fitsf, label="AutoStepFinder fit")
	
	p21 = plot(S_curvef, label="S curve")
	plot(p11, p21, layout=(1,2))
end

# ╔═╡ 680641ef-16b5-4fca-9a57-149bf8c72eb8
md"""## Plot one evolution"""

# ╔═╡ 8d9ea0c9-ab17-44f2-8e42-da95f3d6c9ff
function plot_pixel(evol, i, j, tresH=0.1, N_iter=50)
	dataX = evol[ i, j, :]
	S_curve, best_shots, Fits, steptable = AutoStepMain(dataX, tresH, N_iter)
	p1 = plot(dataX, label = "Data pixel ($i,$j)")
	plot!(Fits, label = "AutoStepFinder fit")

	p2 = plot(S_curve)
	plot(p1, p2)
end

# ╔═╡ b1076909-7591-41aa-bada-4747616cec67
md""" Select pixel: $(@bind pxi NumberField(1:512, default=1)) $(@bind pxj NumberField(1:512, default=1))"""

# ╔═╡ bd225196-cea7-490d-b8bf-d515f51a31a3
md""" Select Threshold $(@bind tresH NumberField(0.01:10., default=0.1))"""

# ╔═╡ cc7875e8-e968-496e-a239-bc4e04c57c55
begin
	advance_epoch # reference the variable to make this cell react to it
	let
		epoch = length(current_fit) + 3
		println(epoch)
		FitX = zeros(size(dataX))  # Assuming dataX is a 2D array
		residuX = dataX - FitX
		
		newFitX, _, _, S_curve, best_shot = stepfindcore(
            residuX, tresH, epoch, true)

		steptable = fit2Steps(dataX, newFitX)
		
		LCSC = length(S_curve)
		push!(current_fit, newFitX)
		push!(current_S_curve, S_curve[1:LCSC-1])
		#push!(counts, epoch)
		push!(tables, steptable)
	end
	epoch_done = "🥔"
end

# ╔═╡ a4ebdb67-4037-49cc-9a39-e2d8d98f33ae
let
	epoch_done # reference the variable to make this cell react to it
	p1 = plot(dataX)
	plot!(p1, current_fit, label="Fit", xlabel="Frames", ylabel="Counts")
	p2 = plot(current_S_curve, label="S_curve")
	plot(p1, p2)
end

# ╔═╡ 282ce225-edf1-4c8d-8dd7-d1c37c2a4852
let
	epoch_done
	tables
end

# ╔═╡ 9a6e62bb-289e-4ca7-9f95-cd0caaac7700
plot_pixel(evol, pxi, pxj, tresH)

# ╔═╡ Cell order:
# ╠═cce4b21e-81e7-4e6c-ad7c-014f9f1d72ea
# ╠═ac949aea-be3f-404d-8a8b-22c25b26ee38
# ╠═980dcf83-1e5c-4559-a42b-b02b0db094b3
# ╠═40687fb5-d226-46e5-8a01-125adc668924
# ╟─ee654cd8-b63a-4e5e-96f8-3104bea4fa60
# ╟─933774f1-9a00-4ee0-98e0-980f9e8fb930
# ╟─c22c5abe-980b-11ee-268c-b16eb65f35f4
# ╟─99ea72eb-b944-4d2a-9347-f4d9255ac790
# ╠═d3da2ac5-df51-4d7e-99f5-039595f6f75e
# ╠═c5133014-5af9-42a0-a04b-4bf1317a2170
# ╠═5ad0425b-f976-4fab-873b-cb1949b15d09
# ╠═23670ee0-8137-4fec-b6b4-86e04b59c324
# ╠═c8dd0947-b6d3-4123-9a8d-a9c8fa029290
# ╠═b2b06a68-dcd5-43e1-9649-8574423e30ec
# ╠═0b04b3af-13cf-48a6-b8b1-800679db9e7d
# ╠═9dd0f232-5668-4e57-b7f5-529d23613f49
# ╠═74183602-063d-4aa9-884e-aff77e9606d9
# ╠═3f6a3fe5-14f1-4f7a-ae9b-34c8d76f4ed4
# ╠═f94daecf-afda-4b61-a5d1-3ae48c035480
# ╠═a0e4327f-8529-4011-9b36-01cdca20a4da
# ╠═10b8263d-8dc5-4ddc-bd21-b310ac7fe11b
# ╠═12b432f4-5942-4958-9de8-473a08c364a5
# ╠═c1b54e5e-d712-45e3-b752-4cb78c618067
# ╠═45db938c-4934-47d9-8d02-5e40542258de
# ╠═967fd807-dca0-4f64-bca8-e3cc1fd7c3f6
# ╠═2e0177fb-47fa-48c7-b81b-3942c2d7794f
# ╠═75ac27d2-372a-470f-b2a1-34a2d3b1253d
# ╠═0772692d-47fa-460d-8bb6-3d397af3f792
# ╠═d20f7145-b292-4cd8-9129-9ddfda4088f7
# ╠═4a902f18-561e-491f-974a-e4e533414d5d
# ╠═302416ac-ea8b-4e0c-8307-0b5e31f298b9
# ╠═e80e3660-fa5e-4489-8bbc-4712ca585d9a
# ╠═a319b59c-f3b3-4319-984d-95a3e29acc3b
# ╠═b3d5892c-d0db-40ea-a3f4-13a076e60e0a
# ╠═2e727245-96ea-4a34-991b-e5129f3ef0d4
# ╠═cdea0164-5a87-4b41-ab6b-db4fba83b0c4
# ╠═8c4e261f-4e1d-478d-ae9b-d6426c0e86be
# ╠═95fca961-dc50-4e6a-908e-ab6fe6a568da
# ╠═66e3ad55-c1aa-4cf5-9fcb-54e2a7199e6a
# ╠═3d84586d-c4c5-4aff-a480-cb81e827436e
# ╠═22527abe-a572-421d-a2de-a8aeb0e53fc3
# ╠═85c62eb4-01b4-47a7-8744-f4da00a75689
# ╠═072342a2-8544-47a7-b04a-077945375013
# ╠═a9a8ef42-e733-4b25-9b9d-af8d756f5720
# ╠═7bf4f8b2-1743-43fa-be90-2a66c485de9c
# ╠═db08a7ae-8e9f-4856-ac87-30097b8a2ac9
# ╠═7a241f86-d20a-4f56-a87a-94aad6c19ef4
# ╠═1fb74702-1e7a-4011-82b0-d0f5035ee670
# ╠═dd89c023-7f16-4d05-9610-802ae7534be9
# ╠═ba55be13-7a28-4dde-b786-3b8b8641a25f
# ╠═9a91475c-ed86-4453-8268-5a0fa08c4db6
# ╠═f16c1ce5-a2a9-4d48-a97b-031745b0f176
# ╠═5f0327d8-fe44-4e92-86e2-62ca4e0192e5
# ╠═1e7467d4-416b-4d82-b355-a0c85dde0d48
# ╠═43258a32-6754-4a0f-8d6c-2a2c97f1f847
# ╠═474c21d6-279d-4364-b7ba-840de1f3025a
# ╠═afba4fac-0a23-4f5e-9f39-53be42930675
# ╠═bf590b23-ce9e-443f-9087-fa41cad24564
# ╠═8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
# ╠═29589bee-6a0c-4f8c-ae19-e98c837d827e
# ╠═2482b515-f159-4823-b4d3-8447972c5324
# ╠═9300f184-7b6e-45dc-b5ea-ab0481a03f0b
# ╠═c14801a5-ebd8-4873-854e-2fe51091c1bd
# ╠═2b44eb05-8735-45c0-aebf-631cd4aa26ec
# ╠═f05b1661-54d5-498f-b906-547a1a420e3c
# ╠═27a53d31-2aa9-4ab8-9c60-7c969c81f0f8
# ╠═d258cbc5-2ff0-495c-af5b-72d0ef26e9e2
# ╠═cc7875e8-e968-496e-a239-bc4e04c57c55
# ╠═1529d78b-5d37-4021-8c09-a5a103f2e26f
# ╠═a4ebdb67-4037-49cc-9a39-e2d8d98f33ae
# ╠═282ce225-edf1-4c8d-8dd7-d1c37c2a4852
# ╠═846e5743-1008-466a-8381-69fbfbf7e3b8
# ╠═5fe17c70-0e67-4bd9-874e-bfe9d9401629
# ╠═194c2282-d507-4226-b333-bc7c97426a67
# ╠═3c1bf9bc-ebf2-43f0-883c-fc94fcdd9a81
# ╠═1f9d0a19-2e19-4efe-94da-68e60cfe17db
# ╠═5c362934-b868-44b1-8879-790339887a8f
# ╠═ff2827f6-1f7e-4a97-95b6-4d81cffd245c
# ╠═6604220f-7f61-4126-8ef1-f78087a7cc70
# ╠═e6d83e3b-03c6-4061-b39b-888ef5db649a
# ╠═4cbbadc4-9948-4a21-ab03-18ed7c73573d
# ╠═df73f9c4-c7bc-4ba2-b34f-adeecba48ede
# ╠═177600ef-163d-4fbf-b178-88e0f628e9e3
# ╠═78f8ff68-ec80-4988-a240-064b1e2bc93c
# ╠═680641ef-16b5-4fca-9a57-149bf8c72eb8
# ╠═8d9ea0c9-ab17-44f2-8e42-da95f3d6c9ff
# ╠═9a6e62bb-289e-4ca7-9f95-cd0caaac7700
# ╠═b1076909-7591-41aa-bada-4747616cec67
# ╠═bd225196-cea7-490d-b8bf-d515f51a31a3
