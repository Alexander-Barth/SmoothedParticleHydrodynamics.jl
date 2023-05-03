
function Particle(x::SVector{N,T}) where {N,T}
    v = @SVector zeros(T,N)
    f = @SVector zeros(T,N)
    return Particle(x,v,f,T(0),T(0))
end

function case_dam_break(
    N,T;
    nparticles_max = Inf,
    h = nothing,
    limits = nothing,
    g = nothing,
    Δt = nothing,
)

    particles = Vector{Particle{N,T}}(undef,0)
    boundary_epsilon = T(h)
    ranges = ntuple(i -> 2*boundary_epsilon:h:(limits[i] - 2*boundary_epsilon),N)

    for i in CartesianIndices(length.(ranges))
        if ranges[1][i[1]] < limits[1]/3
			if (length(particles) >= nparticles_max)
                break
            end

            x = SVector{N}([ranges[j][i[j]]+rand(T)/100 for j in 1:N])
			push!(particles,Particle(x))
        end
    end

    # "Particle-Based Fluid Simulation for Interactive Applications" by Müller et al.

    params = (
        g = SVector{N}(T.(g)),
        rest_density = 300.f0,  # rest density
        gas_const = 2000.f0, # const for equation of state
        h = T(h),		   # kernel radius
        h² = T(h * h),		   # radius^2 for optimization
        mass = 2.5f0*2,		   # assume all particles have the same mass
        viscosity = 200.f0,	   # viscosity constant
        Δt = T(Δt),	   # integration timestep

        # smoothing kernels defined in Müller and their gradients
        # adapted to 2D per "SPH Based Shallow Water Simulation" by Solenthaler et al.
        viscosity_lap = 40.f0 / (π * h^5),
        limits = SVector{N}(T.(limits)),
        boundary_epsilon = T(boundary_epsilon), # boundary epsilon
        boundary_damping = -0.5f0,
        #boundary_damping = 0.5f0,
    )

    return params,particles,KernelSpiky(N,T(h)),KernelPoly6(N,T(h))
end


function InitSPH()
    N = 2
    T = Float32

    case_dam_break(N,T;
        h = 16,
        limits = [1200,900],
        g = [0, -10],
        Δt = 0.0007,	   # integration timestep
    )
end

function step!(params,particles::AbstractVector{Particle{N,T}}) where {N,T}
    Δt = params.Δt

	#=Threads.@threads=# @inbounds for p in particles
		# forward Euler integration
		p.v += Δt .* p.f ./ p.rho
		p.x += Δt * p.v

		# damping at boundary
        p.v = SVector(
            ntuple(Val(N)) do n
                @inbounds begin
                    vn = p.v[n]
                    xn = p.x[n]
		            if (xn < params.boundary_epsilon) ||
                        (xn > params.limits[n] - params.boundary_epsilon)
			            vn *= params.boundary_damping
                    end
                    return vn
                end
            end)

		# enforce boundary conditions
        p.x = SVector(
            ntuple(Val(N)) do n
                @inbounds begin
                    clamp(p.x[n],
                          params.boundary_epsilon,
                          params.limits[n] - params.boundary_epsilon)
                end
            end)

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

function density_pressure(params,W_rho,particles::AbstractVector{Particle{N,T}}) where {N,T}
	#=Threads.@threads=# @inbounds for pi in particles
		pi.rho = zero(T)

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

function forces!(params,W_spiky,particles::AbstractVector{Particle{N,T}}) where {N,T}
    h = params.h
    g = params.g
    mass = params.mass

	#=Threads.@threads=# @inbounds for pi in particles
	    ∇pressure = @SArray zeros(T,N)
	    fvisc = @SArray zeros(T,N)

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
				#∇pressure += -params.mass * (pi.p + pj.p) / (2 * pj.rho) * ∇W(W_spiky,rij,r)
				∇pressure += - pi.rho * params.mass * (pi.p/pi.rho^2 + pj.p/pj.rho^2) * ∇W(W_spiky,rij,r)

				# compute viscosity force contribution
				fvisc += params.viscosity * params.mass * (pj.v - pi.v) / pj.rho * params.viscosity_lap * (h - r)
            end
        end
		fgrav = g * mass / pi.rho
		pi.f = ∇pressure + fvisc + fgrav
    end
end

function update!(params,W_spiky,W_rho,particles)
	density_pressure(params,W_rho,particles)
	forces!(params,W_spiky,particles)
	step!(params,particles)
end

