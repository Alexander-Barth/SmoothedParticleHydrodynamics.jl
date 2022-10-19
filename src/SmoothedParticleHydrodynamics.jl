module SmoothedParticleHydrodynamics


using StaticArrays
using LinearAlgebra
using Random
using SpecialFunctions

abstract type Kernel{N}
end

surface_hypersphere(N) = 2*π^(N/2) / gamma(N/2)

struct KernelSpiky{N,T} <: Kernel{N} where T <: Number
    h::T
    h²::T
    coeff::T
    coeff_grad::T
end

#KernelSpiky(N,h) = KernelSpiky{N,typeof(h)}(h,h^2,15/(π*h^6),-10 / (π * h^5))

"""
```math
\\int_0^h\\left(h-r\\right)^3r^{n-1}dr
=\\frac{h^{n+3}}{n}-\\frac{3h^{n+3}}{n+1}+\\frac{3h^{n+3}}{n+2}-\\frac{h^{n+3}}{n+3}
```

https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B0%7D%5E%7Bh%7D%5Cleft(h-%20r%5Cright)%5E%7B3%7D%20r%5E%7Bn-1%7Ddr%20?or=input

"""
function KernelSpiky(N,h)
    coeff = 1/((1/N - 3/(N+1) + 3/(N+2) - 1/(N+3)) * surface_hypersphere(N) * h^(N+3))
    coeff_grad = -3 * coeff
    KernelSpiky{N,typeof(h)}(h,h^2,coeff,coeff_grad)
end

@inline W(k::KernelSpiky,r2) = k.coeff * (k.h - sqrt(r2))^3
@inline ∇W(k::KernelSpiky,rij,r = norm(rij)) = k.coeff_grad * (k.h - r)^2 * rij / r


struct KernelPoly6{N,T} <: Kernel{N} where T <: Number
    h::T
    h²::T
    coeff::T
    coeff_grad::T
end


#https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B0%7D%5E%7Bh%7D%5Cleft(h%5E%7B2%7D-%20r%5E%7B2%7D%5Cright)%5E%7B3%7D%20r%5E%7Bn-1%7Ddr%20?or=input
function KernelPoly6(N,h)
    coeff = 1/((1/N-3/(N+2)+3/(N+4)-1/(N+6)) * surface_hypersphere(N) * h^(N+6))
    coeff_grad = -6*coeff
    KernelPoly6{N,typeof(h)}(h,h^2,coeff,coeff_grad)
end

@inline W(k::KernelPoly6,r2) =  k.coeff * (k.h² - r2)^3
@inline ∇W(k::KernelPoly6,rij,r = norm(rij)) = k.coeff_grad * (k.h² - r^2)^2 * rij

mutable struct Particle{N,T}
    x::MVector{N,T} # position
    v::MVector{N,T} # velocity
    f::MVector{N,T} # force
    rho::T          # density
    p::T            # pressure
end


function Particle(x::MVector{N,T}) where {N,T}
    v = @MVector zeros(T,N)
    f = @MVector zeros(T,N)
    return Particle(x,v,f,T(0),T(0))
end

function InitSPH()
    N = 2
    T = Float32
    nparticles = 500;
    particles = Vector{Particle{N,T}}(undef,0)
    h = 16.f0		   # kernel radius
    boundary_epsilon = h
    limits = SVector(T(1200),T(900))

	for y = boundary_epsilon:h:(limits[2] - boundary_epsilon * 2.f0)
		for x = (limits[1] / 4):h:(limits[1] / 2)
			if (length(particles) >= nparticles)
                break
            end

	        jitter = rand(T)
			push!(particles,Particle(MVector(x + jitter, y)))
        end
    end


    # "Particle-Based Fluid Simulation for Interactive Applications" by Müller et al.

    params = (
        g = SVector(T(0),T(-10)),
        rest_density = 300.f0,  # rest density
        gas_const = 2000.f0, # const for equation of state
        h = h,		   # kernel radius
        h² = h * h,		   # radius^2 for optimization
        mass = 2.5f0*2,		   # assume all particles have the same mass
        viscosity = 200.f0,	   # viscosity constant
        Δt = 0.0007f0,	   # integration timestep

        # smoothing kernels defined in Müller and their gradients
        # adapted to 2D per "SPH Based Shallow Water Simulation" by Solenthaler et al.
        POLY6 = 4.f0 / (π * h^8),
        SPIKY_GRAD = -10.f0 / (π * h^5),
        viscosity_lap = 40.f0 / (π * h^5),
        limits = limits,
        boundary_epsilon = boundary_epsilon, # boundary epsilon
        boundary_damping = -0.5f0,
    )

    return params,particles,KernelSpiky(N,h),KernelPoly6(N,h)
end

function step!(params,particles::Vector{Particle{N,T}}) where {N,T}
    Δt = params.Δt

	for p in particles
		# forward Euler integration
		p.v .+= Δt .* p.f ./ p.rho
		p.x .+= Δt * p.v

        for n = 1:N
		    # enforce boundary conditions
		    if (p.x[n] < params.boundary_epsilon)
			    p.v[n] *= params.boundary_damping
			    p.x[n] = params.boundary_epsilon
            end

		    if (p.x[n] > params.limits[n] - params.boundary_epsilon)
			    p.v[n] *= params.boundary_damping
			    p.x[n] = params.limits[n] - params.boundary_epsilon
            end
        end
    end
end

#=
@inline W_rho(params,r2) = 4 / (π * params.h^8) * (params.h² - r2)^3

@inline W_spiky(params,r) = 15/(π*params.h^6) * (params.h-r^3)

@inline function ∇W_spiky(params,rij,r = norm(rij))
    h = params.h
    return -10.f0 / (π * h^5) * (h - r)^3 * rij / r
end
=#

function density_pressure(params,W_rho,particles::Vector{Particle{N,T}}) where {N,T}
	for pi in particles
		pi.rho = 0

		for pj in particles
			rij = pj.x - pi.x
			r2 = norm(rij)^2

			if r2 < params.h²
				# this computation is symmetric
				pi.rho += @fastmath params.mass * W(W_rho,r2)
            end
        end
		pi.p = params.gas_const * (pi.rho - params.rest_density)
    end
end

function forces!(params,W_spiky,particles::Vector{Particle{N,T}}) where {N,T}
	∇pressure = @MArray zeros(T,N)
	fvisc = @MArray zeros(T,N)
    h = params.h
    g = params.g
    mass = params.mass

	for pi in particles
        ∇pressure .= 0
        fvisc .= 0

		for pj in particles
			if pi == pj
				continue
            end

			rij = pj.x - pi.x
			r = norm(rij)

			if (r < h)
				# compute pressure force contribution
                # Particle-Based Fluid Simulation for Interactive Applications
                # Matthias Müller, et al. 2003, Eq 10.
				#∇pressure .+= -params.mass * (pi.p + pj.p) / (2 * pj.rho) * ∇W(W_spiky,rij,r)
				∇pressure .+= - pi.rho * params.mass * (pi.p/pi.rho^2 + pj.p/pj.rho^2) * ∇W(W_spiky,rij,r)

				# compute viscosity force contribution
				fvisc .+= params.viscosity * params.mass * (pj.v - pi.v) / pj.rho * params.viscosity_lap * (h - r)
            end
        end
		fgrav = g * mass / pi.rho
		pi.f .= ∇pressure + fvisc + fgrav
    end
end

function update!(params,W_spiky,W_rho,particles)
	density_pressure(params,W_rho,particles)
	forces!(params,W_spiky,particles)
	step!(params,particles)
end



end # module SmoothedParticleHydrodynamics
