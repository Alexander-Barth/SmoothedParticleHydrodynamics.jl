import SmoothedParticleHydrodynamics: InitSPH, forces!, update!, surface_hypersphere, KernelSpiky,KernelPoly6, Kernel, W, ∇W
using StableRNGs
using Test
import SmoothedParticleHydrodynamics: InitSPH, forces!, update!

rng = StableRNG(123)
params,particles,W_spiky,W_rho = InitSPH(rng = rng)

for n = 1:100
    update!(params,W_spiky,W_rho,particles)
end

@test particles[200].x[1] == 274.05356f0
