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

mutable struct Particle{N,T}
    x::SVector{N,T} # position
    v::SVector{N,T} # velocity
    f::SVector{N,T} # force
    rho::T          # density
    p::T            # pressure
end
