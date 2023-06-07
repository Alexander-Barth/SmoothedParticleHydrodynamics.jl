import SmoothedParticleHydrodynamics
import SmoothedParticleHydrodynamics: forces!, update!, surface_hypersphere, KernelSpiky,KernelPoly6, Kernel, W, ∇W, Particle, Location
using StableRNGs
using Test
import SmoothedParticleHydrodynamics: forces!, update!
using BenchmarkTools
using LinearAlgebra
using SpatialHashing: spatial_hash!, spatial_hash, each_near

rng = StableRNG(123)

N = 2
T = Float32


config,particles,W_spiky,W_rho = SmoothedParticleHydrodynamics.case_dam_break(
    N,T,
    h = 16,
    #nparticles = 2000000,
    nparticles = 2000,
    limits = (1200,900),
    g = (0, -10),
    Δt = 0.0007,
    rng = rng)

h = config.h
limits = config.limits

@show h,limits

sz = unsafe_trunc.(Int,limits ./ h) .+ 1
table = zeros(Int,prod(sz)+1)
num_particles = zeros(Int,length(particles))
limits = Tuple(limits)
@time spatial_hash!(particles,h,limits,table,num_particles)
@time spatial_hash!(particles,h,limits,table,num_particles)

spatial_index = (; table, num_particles, config.h, sz)
visited = zeros(Int,length(particles))

@time update!(config,W_spiky,W_rho,particles,spatial_index,visited)

#@code_warntype spatial_hash!(particles,h,limits,table,num_particles)



#locations(i) = particles[i].x
using StaticArrays
import Base: size, getindex



l = @time Location(particles);


spatial_index = spatial_hash(particles,h,limits)

#@test table[end] == length(particles)


i = 1
particles_ref = particles[i]
x = particles_ref.x
search_range = 1
radius2 = h^2
offset = first(CartesianIndices(ntuple(i -> -search_range:search_range,N)))
visited = falses(length(particles))

function find_near!(spatial_index,particles,x,search_range,r2max,near_indices,visited)
    nfound = 0

    @inline each_near(x,search_range,spatial_index,visited) do j
        pj = particles[j]
	    rij = pj.x - x
	    r2 = norm(rij)^2

	    if r2 < r2max
            nfound += 1
#            if nfound <= length(near_indices)
                #nfound += 1
                @inbounds near_indices[nfound] = j
#            end
        end
        #@inbounds near_indices[1] = j
    end
    return nfound
end

near_indices = zeros(Int,length(particles))
r2max = config.h²;
nfound = @btime find_near!(spatial_index,particles,x,search_range,r2max,near_indices,visited)

near = near_indices[1:nfound]


near_ref = Int[]

for j = 1:length(particles)
    pj = particles[j]
	rij = pj.x - particles_ref.x
	r2 = norm(rij)^2

	if r2 < config.h²
        push!(near_ref,j)
    end
end

@test Set(near_ref) == Set(near)



