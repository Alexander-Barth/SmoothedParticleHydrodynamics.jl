import SmoothedParticleHydrodynamics
import SmoothedParticleHydrodynamics: forces!, update!, surface_hypersphere, KernelSpiky,KernelPoly6, Kernel, W, ∇W
using StableRNGs
using Test
import SmoothedParticleHydrodynamics: forces!, update!
using BenchmarkTools

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
