using SpatialHashing: spatial_hash!, spatial_hash, each_near
import SpatialHashing: spatial_hash!

function Particle(x::SVector{N,T}) where {N,T}
    v = @SVector zeros(T,N)
    f = @SVector zeros(T,N)
    return Particle(x,v,f,T(0),T(0))
end

function case_dam_break(N,T; nparticles = 1219, kwargs...)
    particles = Vector{Particle{N,T}}(undef,nparticles)
    case_dam_break!(particles; kwargs...)
end

@inline function case_dam_break!(
    particles::AbstractVector{<:Particle{N,T}};
    h = T(16.f0),
    limits = (1200,900),
    g = (0, -10),
    Δt = 0.0007,	   # integration timestep
    x_noise = 0.01,
    init_particles = true,
    rng = Random.GLOBAL_RNG,
) where {N,T}


    boundary_epsilon = T(h)
    sz = unsafe_trunc.(Int,(limits .- 4*boundary_epsilon) ./ h) .+ 1

    if init_particles
        j = 0
        @inbounds for i in CartesianIndices(sz)
            x = (Tuple(i) .- 1) .* h .+ 2 * boundary_epsilon
            if x[1] < limits[1]/3
    		    if (j >= length(particles))
                    break
                end

                xp = SVector{N}(
                    ntuple(Val(N)) do n
                        @inbounds clamp(x[n] + rand(rng,T) * T(x_noise),
                          boundary_epsilon,
                          limits[n] - boundary_epsilon)
                        #@inbounds x[j]
                    end)

                j += 1
    		    particles[j] = Particle(xp)
            end
        end
    end

    # "Particle-Based Fluid Simulation for Interactive Applications" by Müller et al.

    config = (
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
        #viscosity_lap = 40.f0 / (π * T(h)^5),
        viscosity_lap = 40.f0 / (π * T(h)^2*T(h)^3),
        limits = SVector{N}(T.(limits)),
        boundary_epsilon = T(boundary_epsilon), # boundary epsilon
        boundary_damping = -0.5f0,
        search_range = 2,
    )

    W_spiky = KernelSpiky(N,T(h))
    W_rho = KernelPoly6(N,T(h))
    return config,particles,W_spiky,W_rho
end



function setup_hash(config,particles)
    h = config.h
    limits = config.limits
    sz = unsafe_trunc.(Int,limits ./ h) .+ 1
    table = zeros(Int,prod(sz)+1)
    num_particles = zeros(Int,length(particles))
    limits = Tuple(limits)
    spatial_hash!(particles,h,limits,table,num_particles)
    visited = zeros(Bool,length(num_particles))

    spatial_index = (; table, num_particles, h, sz)

    return spatial_index,visited
end

function step!(config,particles::AbstractVector{Particle{N,T}}) where {N,T}
    Δt = config.Δt

	#=Threads.@threads=# @inbounds for i in 1:length(particles)
        p = particles[i]
		v = p.v + Δt * p.f ./ p.rho
		x = p.x + Δt * v

		# damping at boundary
        vc = SVector(
            ntuple(Val(N)) do n
                @inbounds begin
                    vn = v[n]
                    xn = x[n]
		            if (xn < config.boundary_epsilon) ||
                        (xn > config.limits[n] - config.boundary_epsilon)
			            vn *= config.boundary_damping
                    end
                    return vn
                end
            end)

		# enforce boundary conditions
        xc = SVector(
            ntuple(Val(N)) do n
                @inbounds begin
                    clamp(x[n],
                          config.boundary_epsilon,
                          config.limits[n] - config.boundary_epsilon)
                end
            end)

        particles[i] = Particle(xc,vc,p.f,p.rho,p.p)
    end
end

#=
@inline W_rho(config,r2) = 4 / (π * config.h^8) * (config.h² - r2)^3

@inline W_spiky(config,r) = 15/(π*config.h^6) * (config.h-r^3)

@inline function ∇W_spiky(config,rij,r = norm(rij))
    h = config.h
    return -10.f0 / (π * h^5) * (h - r)^3 * rij / r
end
=#

function density_pressure(config,W_rho,particles::AbstractVector{Particle{N,T}},spatial_index,visited) where {N,T}
    #=Threads.@threads=# @inbounds for i in 1:length(particles)

        pi = particles[i]
		rho = Ref(zero(T))

        @inline each_near(pi.x,config.search_range,spatial_index,visited) do j
            pj = particles[j]
			rij = pj.x - pi.x
			r2 = norm(rij)^2

			if (r2 < config.h²)
				# this computation is symmetric
				rho[] += @fastmath config.mass * W(W_rho,r2)
            end
        end
        #@show rho[]
		p = config.gas_const * (rho[] - config.rest_density)
        particles[i] = Particle(pi.x,pi.v,pi.f,rho[],p)

	    # rho = zero(T)
		# for j in 1:length(particles)
        #     pj = particles[j]
		# 	rij = pj.x - pi.x
		# 	r2 = norm(rij)^2

		# 	if r2 < config.h²
		# 		# this computation is symmetric
		# 		rho += @fastmath config.mass * W(W_rho,r2)
        #     end
        # end
		# p = config.gas_const * (rho - config.rest_density)
        # particles[i] = Particle(pi.x,pi.v,pi.f,rho,p)
    end
end

function forces!(config,W_spiky,particles::AbstractVector{Particle{N,T}},spatial_index,visited) where {N,T}
    h = config.h
    g = config.g
    mass = config.mass

	#=Threads.@threads=# @inbounds for i in 1:length(particles)
        pi = particles[i]

	    ∇pressure = Ref(@SArray zeros(T,N))
	    fvisc = Ref(@SArray zeros(T,N))
        @inline each_near(pi.x,config.search_range,spatial_index,visited) do j
            pj = particles[j]
			rij = pj.x - pi.x
			r = norm(rij)

			if (r < h) &&  (i != j)
				# compute pressure force contribution
                # Particle-Based Fluid Simulation for Interactive Applications
                # Matthias Müller, et al. 2003, Eq 10.
				#∇pressure[] += -config.mass * (pi.p + pj.p) / (2 * pj.rho) * ∇W(W_spiky,rij,r)
				∇pressure[] += - pi.rho * config.mass * (pi.p/pi.rho^2 + pj.p/pj.rho^2) * ∇W(W_spiky,rij,r)

				# compute viscosity force contribution
				fvisc[] += config.viscosity * config.mass * (pj.v - pi.v) / pj.rho * config.viscosity_lap * (h - r)
            end
        end
		fgrav = g * mass / pi.rho
		f = ∇pressure[] + fvisc[] + fgrav

	    # ∇pressure = @SArray zeros(T,N)
	    # fvisc = @SArray zeros(T,N)

		# for j in 1:length(particles)
        #     pj = particles[j]
		# 	rij = pj.x - pi.x
		# 	r = norm(rij)

		# 	if (r < h) &&  (i != j)
		# 		# compute pressure force contribution
        #         # Particle-Based Fluid Simulation for Interactive Applications
        #         # Matthias Müller, et al. 2003, Eq 10.
		# 		#∇pressure += -config.mass * (pi.p + pj.p) / (2 * pj.rho) * ∇W(W_spiky,rij,r)
		# 		∇pressure += - pi.rho * config.mass * (pi.p/pi.rho^2 + pj.p/pj.rho^2) * ∇W(W_spiky,rij,r)

		# 		# compute viscosity force contribution
		# 		fvisc += config.viscosity * config.mass * (pj.v - pi.v) / pj.rho * config.viscosity_lap * (h - r)
        #     end
        # end
		# fgrav = g * mass / pi.rho
		# f = ∇pressure + fvisc + fgrav
        particles[i] = Particle(pi.x,pi.v,f,pi.rho,pi.p)
    end
end



function spatial_hash!(particles::AbstractVector{<:Particle},h,limits,table,num_particles)
    spatial_hash!(Location(particles),h,limits,table,num_particles)
end


function update!(config,W_spiky,W_rho,particles,spatial_index,visited)
    @inline spatial_hash!(particles,config.h,config.limits,
                  spatial_index.table,spatial_index.num_particles)

	density_pressure(config,W_rho,particles,spatial_index,visited)
	forces!(config,W_spiky,particles,spatial_index,visited)
	step!(config,particles)
end
