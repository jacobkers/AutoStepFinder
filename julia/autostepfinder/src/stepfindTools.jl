# -*- coding: utf-8 -*-
""" 
    Basic tools for stepfinding (Julia):
    -split_fast
    Pablo Herrero 2023 
    """

using Images
using Plots
using Statistics
using PlutoUI
using Printf

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

function Indices2Fit(dataX::Vector{Float64}, indices::Vector{Int}, how="mean")
    ixlo = 1
    LX = length(dataX)
    FitX = zeros(LX)
    indices_ext = vcat([0], indices, LX)
    Lix = length(indices_ext)
    for ii in 1:Lix-1
        ixlo = indices_ext[ii] + 1
        ixhi = indices_ext[ii + 1] 
		#println(ixlo, ",", ixhi)
        if ixhi >= ixlo
            if how == "mean"
                FitX[ixlo:ixhi] .= mean(dataX[ixlo:ixhi])
            elseif how == "median"
                FitX[ixlo:ixhi] .= median(dataX[ixlo:ixhi])
            end
        else
            FitX[ixlo:ixhi] .= dataX[ixlo:ixhi]
        end
    end
    return FitX
end

function fit2Steps(dataX, FitX)
    Lx = length(dataX)
    T = collect(1:Lx)
    globalnoise = std(diff(dataX .- FitX)) / sqrt(2)
    ixes0 = findall(diff(FitX) .!= 0)
    Lix = length(ixes0)
    ixes = vcat([0], ixes0, [Lx])
    steptable = zeros(Float64, Lix, 8)
	
    for ii in 2:Lix+1
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
            lo = ix - 1
            ixlo = ixes[lo]
            ixhi = ixes[lo] + 3
            if ixlo < Lx && ixhi < Lx
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


