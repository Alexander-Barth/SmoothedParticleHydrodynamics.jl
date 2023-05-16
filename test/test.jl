import SmoothedParticleHydrodynamics
import SmoothedParticleHydrodynamics: forces!, update!, surface_hypersphere, KernelSpiky,KernelPoly6, Kernel, W, ∇W
using StableRNGs
using Test
import SmoothedParticleHydrodynamics: forces!, update!
using BenchmarkTools
using LinearAlgebra

rng = StableRNG(123)

N = 2
T = Float32

config,particles,W_spiky,W_rho = SmoothedParticleHydrodynamics.case_dam_break(
    N,T,
    h = 16,
    limits = (1200,900),
    g = (0, -10),
    Δt = 0.0007,
    rng = rng)


for n = 1:100
    update!(config,W_spiky,W_rho,particles)
end

@test particles[200].x[1] == 274.05756f0

# mutable struct
# 3.433 ms
# 3.432 ms
# 3.441 ms

# immutable struct
# 2.675 ms
# 2.672 ms
#@btime update!(config,W_spiky,W_rho,particles)

@time update!(config,W_spiky,W_rho,particles)

@time SmoothedParticleHydrodynamics.density_pressure(config,W_rho,particles)

@time SmoothedParticleHydrodynamics.forces!(config,W_spiky,particles)
@time SmoothedParticleHydrodynamics.step!(config,particles)

h = config.h
limits = config.limits

lindex(ind,offsets) = sum((ind .- 1) .* offsets) + 1


function hash(ind,max)
    # https://matthias-research.github.io/pages/publications/tetraederCollision.pdf
    large_primes = (73856093, 19349663, 83492791)[1:length(ind)]
    return abs(reduce(xor,ind .* large_primes)) % max + 1
end

indices(x,h) = unsafe_trunc.(Int,x ./ h) .+ 1



function spatial_hash!(particles,h,limits,table,num_particles)
    table .= 0
    num_particles .= 0
    sz = unsafe_trunc.(Int,limits ./ h) .+ 1
    #offsets = @inbounds (1, cumprod(sz[2:end])...)

    # count particles

    for i = 1:length(particles)
        ind = indices(particles[i].x,h)
        #    l = lindex(ind,offsets)
        l = hash(ind,length(table))
        table[l] += 1
    end


    # partial sums

    for i = 2:length(table)
        table[i] += table[i-1]
    end


    # fill-in

    for i = 1:length(particles)
        ind = indices(particles[i].x,h)
        #    l = lindex(ind,offsets)
        l = hash(ind,length(table))
        num_particles[table[l]] = i
        table[l] -= 1
    end

    return nothing
end

table = zeros(Int,prod(sz)+1)
num_particles = zeros(Int,length(particles))
@time spatial_hash!(particles,h,limits,table,num_particles)
@time spatial_hash!(particles,h,limits,table,num_particles)

#@code_warntype spatial_hash!(particles,h,limits,table,num_particles)


#=

function spatial_hash(particles,h,limits)
    table = zeros(Int,prod(sz)+1)
    num_particles = zeros(Int,length(particles))
    @time spatial_hash!(particles,h,limits,table,num_particles)
    return (; table, num_particles, h)
end

spatial_index = spatial_hash(particles,h,limits)

@test table[end] == length(particles)

function each_near(fun,x,search_range,spatial_index)
    table,num_particles,h = spatial_index
    N = length(x)
    ind = indices(x,h)

    for offset in CartesianIndices(ntuple(i -> -search_range:search_range,N))
        #@show ind .+ Tuple(offset)
        l = hash(ind .+ Tuple(offset),length(table))

        for i = table[l]:table[l+1]
            fun(num_particles[i])
        end
        #num_particles[]
        #num_particles[table[l+1]]
    end
end

i = 100
pi = particles[i]
x = pi.x
search_range = 1
radius2 = h^2
offset = first(CartesianIndices(ntuple(i -> -search_range:search_range,N)))


near = Int[]

each_near(x,search_range,table,num_particles,h) do j
    pj = particles[j]
	rij = pj.x - pi.x
	r2 = norm(rij)^2

	if r2 < config.h²
        push!(near,j)
    end
end

near_ref = Int[]

for j = 1:length(particles)
    pj = particles[j]
	rij = pj.x - pi.x
	r2 = norm(rij)^2

	if r2 < config.h²
        push!(near_ref,j)
    end
end

@test Set(near_ref) == Set(near)

=#
