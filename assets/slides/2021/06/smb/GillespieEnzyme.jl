using Distributions, StatsBase, Random, Statistics

export Gillespie, GillespieModel, GillespieEnd
export sim

abstract type GillespieModel end
Base.length(::GillespieModel) = error("Define a length for the Gillespie model: i.e. the number of reaction channels.")

mutable struct ReactionChannel
    dx::Vector{Float64}
    n::Int
    ReactionChannel() = ReactionChannel(Float64[])
    ReactionChannel(dx::Vector{Float64}) = new(dx, length(dx))
end
function ReactionChannel(n::Int)
    C = Vector{ReactionChannel}(undef, n)
    for i in eachindex(C)
        C[i] = ReactionChannel()
    end
    C
end
function Base.iterate(C::ReactionChannel, i::Int=1)
    if i>C.n
        dx = Random.randexp()
        push!(C.dx, dx)
        C.n += 1
    else
        dx = C.dx[i]
    end
    i += 1
    dx, i
end
Base.IteratorSize(::Type{ReactionChannel}) = IsInfinite()
Base.eltype(::Type{ReactionChannel}) = Float64

abstract type BRNSimulation end
Base.size(B::BRNSimulation, n::Int) = size(B)[n]
Base.size(::BRNSimulation) = error("Define a size (dim_y, dim_t) for the BRN simulation.")

struct GillespieEnd{P<:GillespieModel} <: BRNSimulation 
    S::Vector{Vector{Int}}
    u0::Vector{Int}
    C::Vector{ReactionChannel}
    propfun::P
    function GillespieEnd(S, u0, C, propfun::P) where P<:GillespieModel
        m = length(S)
        n = length(u0)
        length(propfun) == m || error("Stoichiometric matrix contains $(m) columns, but $(length(P)) reactions.")
        all(==(n)∘length, S) || error("Initial state vector is length $(n) but stoichiometric matrix doesn't match this number of rows.")
        length(C) == m || error("Need $(m) reaction channels.")
        return new{P}(S, u0, C, propfun)
    end
end
function GillespieEnd(S, u0, propfun)
    m = length(S)
    GillespieEnd(S, u0, ReactionChannel(m), propfun)
end
Base.size(::GillespieEnd) = (1,1) # length of y and length of t

function reset!(C::ReactionChannel)
    empty!(C.dx)
    C.n = 0
    return nothing
end
function reset!(G::GillespieEnd) 
    for c in G.C
        reset!(c)
    end
    return nothing
end
function couple!(dst::Vector{ReactionChannel}, src::Vector{ReactionChannel}, dst_idx::Int, src_idx::Int)
    dst[dst_idx].n = src[src_idx].n
    copy!(dst[dst_idx].dx, src[src_idx].dx)
    return nothing
end
couple!(dst::GillespieEnd, src::GillespieEnd, i, j) = couple!(dst.C, src.C, i, j)
couple!(dst, src, (d_o, s_o)::Tuple{Int,Int}) = couple!(dst, src, d_o, s_o)


function __sim!(::BRNSimulation, y::AbstractVector{T}) where T<:Real
    error("Define __sim!(G, y) [i.e. untimed simulation] or _sim!(G, y) [fallback is timing __sim!]")
end
function _sim!(G::BRNSimulation, y::AbstractVector{T}, t::AbstractVector{T}) where T<:Real
    t[1] = @elapsed __sim!(G, y)
    return nothing
end
function sim!(G::BRNSimulation, y::DenseMatrix{T}, t::DenseMatrix{T}) where T<:Real
    for i in 1:size(y,2)
        @inbounds sim!(G, view(y, :, i), view(t, :, i))
    end
end

function sim!(G::BRNSimulation, y::AbstractVector{T}, t::AbstractVector{T}) where T<:Real
    _sim!(G, y, t)
    reset!(G)
    return nothing
end
function _sim(G::BRNSimulation)
    y = Vector{Float64}(undef, size(G, 1))
    t = Vector{Float64}(undef, size(G, 2))
    _sim!(G, y, t)
    return y, t
end
function sim(G::BRNSimulation)
    y, t = _sim(G)
    reset!(G)
    y, t
end
function sim(G::BRNSimulation, n::Int)
    y = Matrix{Float64}(undef, size(G, 1), n)
    t = Matrix{Float64}(undef, size(G, 2), n)
    sim!(G, y, t)
    return y, t
end
function sim(G::BRNSimulation, T::Float64)
    y = Vector{Float64}(undef, size(G, 1))
    t = Vector{Float64}(undef, size(G, 2))
    
    y_out = Vector{Float64}()
    t_out = Vector{Float64}()
    
    t_stop_condition = 0.0
    while t_stop_condition<T
        sim!(G, y, t)
        t_stop_condition += sum(t)
        append!(y_out, y)
        append!(t_out, t)
    end

    return reshape(y_out, size(G,1), :), reshape(t_out, size(G,2), :)
end


############## Specific to GillespieEnd

function __sim!(G::GillespieEnd, y::AbstractVector{T}) where T<:Real
    m = length(G.S)
    v = zeros(m)

    u = copy(G.u0)
    @inbounds y[1] = 0.0

    cactus = Iterators.Stateful.(G.C)
    dx = first.(cactus)
    dt = zero(dx)
    a = zero(dx)

    while true
        G.propfun(a, u)
        if all(iszero,a)
            return nothing
        else
            @. dt = dx / a
            δt, i = findmin(dt)

            @. dx -= δt*a
            @inbounds dx[i] = first(cactus[i])

            @inbounds y[1] += δt
            @inbounds u .+= G.S[i]
        end
    end
end


####################
# MODEL TYPES

export Enzyme, EnzymeReduced

struct Enzyme <: GillespieModel
    k1::Float64
    m1::Float64
    k2::Float64
    nz::Int
end
function (M::Enzyme)(a, u)
    a[1] = M.k1*u[1]*u[2]
    a[2] = M.m1*u[3]
    a[3] = M.k2*u[3]
    return nothing
end
Base.length(::Enzyme)=3

struct EnzymeReduced <: GillespieModel
    k2::Float64
    nz::Int
    κ::Float64
end
function (M::EnzymeReduced)(a, u)
    a[1] = M.k2 * M.nz * u[1] / (M.κ + u[1])
    return nothing
end
Base.length(::EnzymeReduced)=1

EnzymeReduced(M::Enzyme) = EnzymeReduced(M.k2, M.nz, M.m1/M.k1)

# FUDGE FACTOR! Add 20 onto the reduced model...
function _sim!(G::GillespieEnd{EnzymeReduced}, y::AbstractVector{T}, t::AbstractVector{T}) where T<:Real
    t[1] = @elapsed __sim!(G, y)
    y[1] += 20.0
    return nothing
end

GillespieEnd(M::Enzyme, s0::Int64=100) = GillespieEnd([[-1, -1, 1, 0], [1, 1, -1, 0], [0, 1, -1, 1]], [s0, M.nz , 0, 0], M)
GillespieEnd(M::EnzymeReduced, s0::Int64=100) = GillespieEnd([[-1, 1]], [s0, 0], M)

