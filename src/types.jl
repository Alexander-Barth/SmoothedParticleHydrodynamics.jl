abstract type Kernel{N}
end

struct KernelSpiky{N,T} <: Kernel{N} where T <: Number
    h::T
    h²::T
    coeff::T
    coeff_grad::T
end

struct KernelPoly6{N,T} <: Kernel{N} where T <: Number
    h::T
    h²::T
    coeff::T
    coeff_grad::T
end

struct Particle{N,T}
    x::SVector{N,T} # position
    v::SVector{N,T} # velocity
    f::SVector{N,T} # force
    rho::T          # density
    p::T            # pressure
end

"""
    Adapter type for SpatialHashing such that `loc[i]` returns `p.x[i]` where
`loc = Location(p)`.
"""
struct Location{N,T,TC <: AbstractVector{Particle{N,T}}} <: AbstractVector{SVector{N,T}}
    p::TC
end

@inline Base.size(l::Location) = size(l.p)
@inline Base.@propagate_inbounds Base.getindex(l::Location,index) = l.p[index].x
