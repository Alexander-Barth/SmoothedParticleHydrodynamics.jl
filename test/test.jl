import SmoothedParticleHydrodynamics
import SmoothedParticleHydrodynamics: InitSPH, forces!, update!, surface_hypersphere, KernelSpiky,KernelPoly6, Kernel, W, ∇W
using StableRNGs
using Test
import SmoothedParticleHydrodynamics: InitSPH, forces!, update!
using BenchmarkTools

rng = StableRNG(123)
params,particles,W_spiky,W_rho = InitSPH(rng = rng)

for n = 1:100
    update!(params,W_spiky,W_rho,particles)
end

@test particles[200].x[1] == 274.05756f0

# mutable struct
# 3.433 ms
# 3.432 ms
# 3.441 ms

# immutable struct
# 2.675 ms
# 2.672 ms
#@btime update!(params,W_spiky,W_rho,particles)

@time update!(params,W_spiky,W_rho,particles)

@time SmoothedParticleHydrodynamics.density_pressure(params,W_rho,particles)

@time SmoothedParticleHydrodynamics.forces!(params,W_spiky,particles)
@time SmoothedParticleHydrodynamics.step!(params,particles)
