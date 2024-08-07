export ContinuationProbability
export ConstantCP, PWConstCP
export MFPair

# Continuation probabilities
abstract type ContinuationProbability end
function tocontinue(A::ContinuationProbability, args...)::Bool
    rand()<A(args...)
end

struct ConstantCP<:ContinuationProbability
    α::Float64
    function ConstantCP(x=1.0)
        0<x<=1 || error("$(x) is not a probability")
        new(x)
    end
end
(A::ConstantCP)(args...) = A.α

struct PWConstCP <: ContinuationProbability
    η1::Float64
    η2::Float64
    isnear
    function PWConstCP(η1, η2, isnear)
        (all(>(0), (η1,η2)) & all(<=(1), (η1,η2))) || error("Not probabilities")
        new(η1, η2, isnear)
    end
end
(H::PWConstCP)(ytilde...) = H.isnear(ytilde...) ? H.η1 : H.η2

# Optimal continuation probabilities
export ROC

ROC = NamedTuple{(:TP, :FP, :FN, :TN),NTuple{4,Float64}}
function ConstantCP(roc::ROC, t, that; αmin=0.01)
    p = roc.TP + roc.FN
    p_dis = roc.FP + roc.FN
    αstar = (that/t) * (p_dis/(p - p_dis))
    return 0<=αstar<=1 ? ConstantCP(max(αmin, sqrt(αstar))) : ConstantCP()
end
function PWConstCP(roc::ROC, t, that; isnear, ηmin=(0.01, 0.01))
    
    p = roc.TP + roc.FN
    η1min, η2min = ηmin

    cp = t*(roc.TP + roc.FP)
    cn = t*(roc.TN + roc.FN)

    ϕ((η1, η2)) = (p - (1-(1/η1))*roc.FP - (1-(1/η2))*roc.FN)*(that + η1*cp + η2*cn)
    eta1(x) = max(η1min, min(1.0, sqrt(((that + cn*x)/(p-roc.FP-roc.FN*(1-(1/x))))*(roc.FP/cp))))
    eta2(x) = max(η2min, min(1.0, sqrt(((that + cp*x)/(p-roc.FN-roc.FP*(1-(1/x))))*(roc.FN/cn))))

    R0 = p - roc.FP - roc.FN
    if R0 > 0
        η1bar = sqrt((that/R0)*(roc.FP/cp))
        η2bar = sqrt((that/R0)*(roc.FN/cn))
        if (η1min<=η1bar<=1)&(η2min<=η2bar<=1)
            return PWConstCP(η1bar, η2bar, isnear)
        end
    end
    
    v = [(1.0, eta2(1.0)), (eta1(1.0), 1.0), (η1min, eta2(η1min)), (eta1(η2min), η2min)]
    w = map(ϕ, v)
    ϕmin, j = findmin(w)

    return try PWConstCP(v[j]..., isnear)
    catch
        PWConstCP(1.0, 1.0, isnear)
    end
end

#### MULTIFIDELITY PAIR
struct MFPair{L<:BRNSimulation, H<:BRNSimulation, CP<:ContinuationProbability} <: BRNSimulation
    lo::L
    hi::H
    c::Vector{Tuple{Int,Int}} # (i,j) where lofi reaction j couples to hifi reaction i 
    α::CP
    function MFPair(lo::L, hi::H, c::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[], α::CP=ConstantCP()) where {L,H,CP}
        for (i,j) in c
            ((j in eachindex(lo.S)) && (i in eachindex(hi.S))) || error("Coupling specification is wrong")
        end
        return new{L,H,CP}(lo, hi, c, α)
    end
end
MFPair(lo, hi, (i,j)::Tuple{Int,Int}, α...) = MFPair(lo, hi, [(i,j)], α...)
MFPair(lo, hi, i::Int, j::Int, α...) = MFPair(lo, hi, (i,j), α...)

Base.size(MF::MFPair) = size(MF.lo).+size(MF.hi)

function reset!(MF::MFPair)
    reset!(MF.lo)
    reset!(MF.hi)
    return nothing
end

function _sim!(MF::MFPair, y::AbstractVector{T}, t::AbstractVector{T}) where T<:Real

    y_idx_lo = range(1, length=size(MF.lo, 1))
    t_idx_lo = range(1, length=size(MF.lo, 2))

    y_idx_hi = range(size(MF.lo, 1)+1, length=size(MF.hi,1))
    t_idx_hi = range(size(MF.lo, 2)+1, length=size(MF.hi,2))

    _sim!(MF.lo, view(y, y_idx_lo), view(t, t_idx_lo))
    if tocontinue(MF.α, view(y, y_idx_lo)...)
        for ij in MF.c
            couple!(MF.hi, MF.lo, ij)
        end
        _sim!(MF.hi, view(y, y_idx_hi), view(t, t_idx_hi))
    else
        view(y, y_idx_hi) .= NaN
        view(t, t_idx_hi) .= 0.0
    end
    return nothing
end