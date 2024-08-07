using StatsPlots, Plots, LaTeXStrings
using ProgressMeter
using Distributed
@everywhere using Statistics

const data=180
const epsilon=10
const NZ=2:1:12
remotecall_fetch.(()->data, workers())
remotecall_fetch.(()->epsilon, workers())
remotecall_fetch.(()->NZ, workers())


opt = (
    size=(666,400),
    legendfontsize=12,
    titlefontsize=14,
    tickfontsize=12,
    guidefontsize=12,
    background_color=RGB(46/255, 52/255, 64/255),
    palette=:Set3_11,
    fontfamily="Helvetica",
)

# Models and Gillespie simulations 
@everywhere module EnzymeExample
    include("GillespieEnzyme.jl")
    include("mf.jl")
end
@everywhere using .EnzymeExample

################
# Generate data to plot

@everywhere begin
    Glo(θ) = GillespieEnd(EnzymeReduced(0.1, θ, 1.0))
    Ghi(θ) = GillespieEnd(Enzyme(100, 100, 0.1, θ))

    MF_uncoupled(θ) = MFPair(Glo(θ), Ghi(θ))
    MF_coupled(θ) = MFPair(Glo(θ), Ghi(θ), (3,1))

    f_coupled(θ, n) = sim(MF_coupled(θ), n)
    f_uncoupled(θ, n) = sim(MF_uncoupled(θ), n)
end

trial = sim(MF_coupled(5),10)

using Random
Random.seed!(1234)
const num_sim = fill(1000, size(NZ))
const csim_data = @showprogress pmap(f_coupled, NZ, num_sim)
remotecall_fetch.(()->csim_data, workers())

const bin_width = 2*epsilon
const Ω_max = data+epsilon
const Ω_min = data-epsilon

bin_min = [minimum(dat[1]) for dat in csim_data]
@. bin_min = data + epsilon + (2*epsilon*round((bin_min - data - epsilon)/(2*epsilon), RoundDown))
bin_max = [maximum(dat[1]) for dat in csim_data]
@. bin_max = data + epsilon + (2*epsilon*round((bin_max - data - epsilon)/(2*epsilon), RoundUp))

const bins = range.(bin_min, bin_max; step=bin_width)

# View posteriors
@everywhere isnear(t; epsilon=epsilon) = abs(t-data)<epsilon
inOmega(; epsilon=epsilon) = [isnear.(dat[1], epsilon=epsilon) for dat in csim_data]
const T = map(dat->sum(selectdim(dat[2],1,2)), csim_data)
const That = map(dat->sum(selectdim(dat[2],1,1)), csim_data)

@everywhere function P(epsilon) 
    P = [mean(isnear.(selectdim(dat[1], 1, 2); epsilon=epsilon)) for dat in csim_data]
    P ./= sum(P)
    P
end
@everywhere function Phat(epsilon) 
    P = [mean(isnear.(selectdim(dat[1], 1, 1); epsilon=epsilon)) for dat in csim_data]
    P ./= sum(P)
    P
end

# View sample building exercise
function show_f(nz; show_data = true, kwargs...)
    fig = plot(; title="θ = $(nz)", xlabel="y", yticks=[], opt..., kwargs...)
    I = findall(==(nz), NZ)
    for i in I
        stephist!(fig, selectdim(csim_data[i][1],1,2), normalize=false, seriescolor=i, bins=bins[i], label="f(y | θ = $(nz))")
    end
    if show_data
        YY = ylims(fig)
        bar!(fig, [data], [YY[2]], bar_width=bin_width, label="", seriescolor=:white, seriesalpha=0.3)
    end
    fig
end
function show_f(fn::String; show_lh = false, kwargs...)
    annie = @animate for i in eachindex(csim_data)
        show_f(i+1; xlims=(minimum(bin_min), maximum(bin_max)), ylims=(0, 500), kwargs...)
    end
    return gif(annie, fn, fps=1)
end
function show_ℓ(fn::String; epsilon=epsilon, kwargs...)
    fig = plot(; xticks=NZ, xlabel="θ", xlims=(1,13), ylims=(0,0.4), yticks=[], opt..., kwargs...)
    PP = P(epsilon)
    annie = @animate for i in eachindex(csim_data)
        if i==1
            plot!(fig, NZ[[i]], PP[[i]]; seriescolor=1, seriestype=:bar, label="ABC Likelihood")
        else
            plot!(fig, NZ[[i]], PP[[i]]; seriescolor=1, seriestype=:bar, label="")
        end
    end
    gif(annie, fn, fps=1)
end

function show_f̂(nz; show_data = true, kwargs...)
    fig = plot(; title="θ = $(nz)", xlabel="y", yticks=[], opt..., kwargs...)
    I = findall(==(nz), NZ)
    for i in I
        stephist!(fig, selectdim(csim_data[i][1],1,1), normalize=false, seriescolor=i, bins=bins[i], label="f(y | θ = $(nz))")
    end
    if show_data
        YY = ylims(fig)
        bar!(fig, [data], [YY[2]], bar_width=bin_width, label="", seriescolor=:white, seriesalpha=0.3)
    end
    fig
end
function show_f̂(fn::String; show_lh = false, kwargs...)
    annie = @animate for i in eachindex(csim_data)
        show_f̂(i+1; xlims=(minimum(bin_min), maximum(bin_max)), ylims=(0, 500), kwargs...)
    end
    return gif(annie, fn, fps=1)
end
function show_ℓ̂(fn::String; epsilon=epsilon, kwargs...)
    fig = plot(; xticks=NZ, xlabel="θ", xlims=(1,13), ylims=(0,0.6), yticks=[], opt..., kwargs...)
    PP = Phat(epsilon)
    annie = @animate for i in eachindex(csim_data)
        if i==1
            plot!(fig, NZ[[i]], PP[[i]]; seriescolor=2, seriestype=:bar, label="ABC Likelihood")
        else
            plot!(fig, NZ[[i]], PP[[i]]; seriescolor=2, seriestype=:bar, label="")
        end
    end
    gif(annie, fn, fps=1)
end

using JLD, LinearAlgebra
function compare_π(; epsilon=epsilon, kwargs...)
    muttley = load("budgeted.jld")
    
    Xlo = muttley["Xlo"]
    Xhi = muttley["Xhi"]

    f!(mat) = for v in eachcol(mat)
        normalize!(v,1)
    end
    f!(Xlo)
    f!(Xhi)
   
    fig = groupedbar(NZ, [mean(Xlo,dims=2) mean(Xhi,dims=2)];
        xticks=NZ,
        xlabel="θ",
        opt...,
        label=["Low-fidelity" "High-fidelity"],
        yerror=[std(Xlo,dims=2) std(Xhi,dims=2)],
        seriescolor=[2 1],
    )
    plot!(fig; title="ABC posteriors: ϵ=$(epsilon)", kwargs...)
end
function compare_π_mf(; epsilon=epsilon, kwargs...)
    muttley = load("budgeted.jld")
    
    Xlo = muttley["Xlo"]
    Xhi = muttley["Xhi"]
    Xmf = muttley["Xmf"]

    f!(mat) = for v in eachcol(mat)
        normalize!(v,1)
    end
    f!(Xlo)
    f!(Xhi)
    f!(Xmf)
   
    fig = groupedbar(NZ, [mean(Xlo,dims=2) mean(Xhi,dims=2) mean(Xmf,dims=2)];
        xticks=NZ,
        xlabel="θ",
        opt...,
        label=["Low-fidelity" "High-fidelity" "Multifidelity"],
        yerror=[std(Xlo,dims=2) std(Xhi,dims=2) std(Xmf,dims=2)],
        seriescolor=[2 1 3],
    )
    plot!(fig; title="ABC posteriors: ϵ=$(epsilon)", kwargs...)
end


# Multifidelity analysis

function show_pairs(nz; kwargs...)
    fig = plot(; title="Multifidelity simulation", ylabel="Low-fidelity y", xlabel="High-fidelity y", opt...)
    I = findall(==(nz), NZ)
    for i in I
        scatter!(fig, selectdim(csim_data[i][1],1,2), selectdim(csim_data[i][1],1,1);
            seriescolor=i,
            bins=bins[i],
            linewidth=1,
            serieslinecolor=:black,
            label="θ = $(nz)"
        )
    end
    XX = xlims(fig)
    YY = ylims(fig)
    bar!(fig, [data], [YY[2]], bar_width=bin_width, label="", seriescolor=:white, seriesalpha=0.3)
    bar!(fig, [data], [XX[2]], bar_width=bin_width, label="", orientation=:horizontal, seriescolor=:white, seriesalpha=0.3)
    plot!(fig; xlims=XX, ylims=YY, ratio=:equal, kwargs...)
end

TP(mat) = mean(selectdim(mat, 1, 1) .& selectdim(mat, 1, 2))
FP(mat) = mean(selectdim(mat, 1, 1) .& (!).(selectdim(mat, 1, 2)))
FN(mat) = mean((!).(selectdim(mat, 1, 1)) .& selectdim(mat, 1, 2))
TN(mat) = mean((!).(selectdim(mat, 1, 1)) .& (!).(selectdim(mat, 1, 2)))

roc_vector = map(mat->(TP=TP(mat), FP=FP(mat), FN=FN(mat), TN=TN(mat)), inOmega())
#αstar = map(ConstantCP, roc_vector, T, That)
#ηstar = map((a,b,c)->PWConstCP(a, b, c; isnear=isnear), roc_vector, T, That)
remotecall_fetch.(()->roc_vector, workers())
#remotecall_fetch.(()->αstar, workers())
#remotecall_fetch.(()->ηstar, workers())

const Tbar = map(dat->mean(selectdim(dat[2],1,2)), csim_data)
const Thatbar = map(dat->mean(selectdim(dat[2],1,1)), csim_data)
remotecall_fetch.(()->Tbar, workers())
remotecall_fetch.(()->Thatbar, workers())


# Choose multifidelity
@everywhere begin
    F_nz(k_of_interest) = i -> Int(i==k_of_interest)
    Fbar(k_of_interest) = P(epsilon)[k_of_interest]
    GGval = mean(Thatbar)
    FF(i; k_of_interest) = (F_nz(k_of_interest)(i) - Fbar(k_of_interest))^2 * (roc_vector[i].TP - roc_vector[i].FP)
    FFval(k_of_interest) = mean(i->FF(i; k_of_interest=k_of_interest), 1:11)

    λ(k_of_interest) = sqrt(GGval / FFval(k_of_interest))
    α(i; k_of_interest, αmin=0.001) = max(αmin, λ(k_of_interest) * abs(F_nz(k_of_interest)(i) - Fbar(k_of_interest)) * sqrt((roc_vector[i].FP + roc_vector[i].FN)/Tbar[i]))
    α(; k_of_interest) = map(i->α(i; k_of_interest=k_of_interest), eachindex(NZ))
    alpha_const(; k_of_interest) = maximum(α(; k_of_interest=k_of_interest))

    Tmf(; k_of_interest) = Thatbar .+ α(; k_of_interest=k_of_interest).*Tbar
    Tmf_const(; k_of_interest) = Thatbar .+ alpha_const(; k_of_interest=k_of_interest).*Tbar

    N(sample_budget; k_of_interest) = floor(sample_budget/sum(Tmf(k_of_interest=k_of_interest)))
    N_const(sample_budget; k_of_interest) = floor(sample_budget/sum(Tmf_const(k_of_interest=k_of_interest)))
end

@everywhere function budget_sims(; k_of_interest, epsilon=epsilon, sample_budget::Float64=10.0)
    wlo = Vector{Vector{Float64}}()
    whi = Vector{Vector{Float64}}()
    wmf = Vector{Vector{Float64}}()
    wmf_const = Vector{Vector{Float64}}()

    Tlo_all = sum(Thatbar)
    Thi_all = sum(Tbar)

    Nlo = Int(floor(sample_budget/Tlo_all))
    Nhi = Int(floor(sample_budget/Thi_all))
    Nmf = Int(N(sample_budget; k_of_interest=k_of_interest))
    Nmf_const = Int(N_const(sample_budget; k_of_interest=k_of_interest))

    α_const = ConstantCP(alpha_const(k_of_interest=k_of_interest))
    for (i,nz) in enumerate(NZ)
        α_i = ConstantCP(α(i; k_of_interest=k_of_interest))
        
        # Models
        glo = Glo(nz)
        ghi = Ghi(nz)
        gmf = MFPair(glo, ghi, 3, 1, α_i)
        gmf_const = MFPair(glo, ghi, 3, 1, α_const)

        # Simulations
        ylo, tlo = sim(glo, Nlo)
        yhi, thi = sim(ghi, Nhi)
        ymf, tmf = sim(gmf, Nmf)
        ymf_const, tmf_const = sim(gmf_const, Nmf_const)

        push!(wlo, isnear.(ylo[1,:]; epsilon=epsilon))
        push!(whi, isnear.(yhi[1,:]; epsilon=epsilon))
        push!(wmf, isnear.(ymf[1,:]; epsilon=epsilon) .+ ((!isnan).(ymf[2,:])).*(isnear.(ymf[2,:]; epsilon=epsilon) .- isnear.(ymf[1,:]; epsilon=epsilon))./(α_i.(ymf[1,:])))
        push!(wmf_const, isnear.(ymf_const[1,:]; epsilon=epsilon) .+ ((!isnan).(ymf_const[2,:])).*(isnear.(ymf_const[2,:]; epsilon=epsilon) .- isnear.(ymf_const[1,:]; epsilon=epsilon))./(α_const.(ymf_const[1,:])))
        
    end

    Llo = mean.(wlo)
    Lhi = mean.(whi)
    Lmf = mean.(wmf)
    Lmf_const = mean.(wmf_const)

    return (
        Llo=Llo[k_of_interest]/sum(Llo),
        Lhi=Lhi[k_of_interest]/sum(Lhi),
        Lmf=Lmf[k_of_interest]/sum(Lmf),
        Lmf_const=Lmf_const[k_of_interest]/sum(Lmf_const)
    )
end

using StructArrays
function compare_methods(N=100; k_of_interest)
    X = @showprogress pmap(1:N) do i
        budget_sims(k_of_interest=k_of_interest)
    end
    return StructArray(X)
end



function show_ratio(; k_of_interest, kwargs...)
    alpha = α(; k_of_interest=k_of_interest)
    fig = bar(NZ, alpha;
        c=NZ.-1,
        xticks=NZ,
        xlabel=L"\theta",
        ylabel=L"\alpha(\theta)",
        opt...,
        palette=:Spectral_11,
        label="",
        background_color=:white,
        fontfamily="Helvetica",
    )
    hline!(fig, [mean(alpha)], label="Constant", linestyle=:dash, seriescolor=:black)
    plot!(fig; kwargs...)
end

function show_estimates(; kwargs...)
    muttley = load("compared.jld")
    p = muttley["prob_z0_equals_7"]

    fig = boxplot([p.Lhi p.Lmf_const p.Lmf];
    palette=:Spectral_11,
    seriescolor=[1 9 11],
    xticks = ([1,2,3],["High fidelity" "Constant" "Optimised"]),
    xtickfontsize=10,
    legend=:none,
    ylabel = "Estimate",
    opt...,
    background_color=:white,
    fontfamily="Helvetica",
    kwargs...
    )
end

function proposal_plot(; kwargs...)
    l = @layout [a{0.6w} b]

    A = show_ratio(k_of_interest=6, title="A. Optimised high:low ratio", size=(480,300))
    B = show_estimates(title="B. Estimator sample", size=(320,300))
    fig = plot(A, B; layout=l, size=(800,300), kwargs...)
end