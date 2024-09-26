# -*- coding: utf-8 -*-

using Images
using Plots
using Statistics
using PlutoUI
using Printf
using FixedPointNumbers


function stepfindcore(dataX::Vector{Float64}, tresH=0.15, N_iter=0)
   
    # run an iterative step fit deep into overfitting:
    if N_iter == 0 || N_iter > length(dataX) ÷ 4
        N_iter = length(dataX) ÷ 4
    end
    S_curve, splitlog, split_table = A121_split_until_ready(dataX, N_iter) 

    # estimate best step number in this iteration range:
    best_shot, S_curve_fin = A122_eval_Scurve(S_curve, tresH)  

    # re-build best fit from 'split log':
    if best_shot > 0
		
        indices_bestfit, indices_counterfit = Splitlog2FinalIndices(splitlog, best_shot)  
		#println(indices_bestfit, indices_counterfit)
        FitX = Indices2Fit(dataX, indices_bestfit, "mean")  
        cFitX = Indices2Fit(dataX, indices_counterfit, "mean")
        
    else
        FitX = zeros(length(dataX))
        cFitX = zeros(length(dataX))
        best_shot = 0
    end

    return FitX, cFitX, splitlog, S_curve_fin, best_shot
end

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

function A1210_set_up_splitlogtable(dataX)
    Na = length(dataX)
    i_nxt, avl, avr, rankit, errorcurve = splitFast(dataX)  # Assuming st_splitFast has an equivalent in Julia
    new_row1 = [1 Na i_nxt avl avr rankit]
    split_table = [
        -1 -1 -1 1 1 1
        new_row1
        Na Na Na 1 1 1
    ]
    return split_table
end

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

function Splitlog2FinalIndices(splitlog, best_shot)
	
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