using Test
using LinearAlgebra
import SmoothedParticleHydrodynamics: forces!, update!, surface_hypersphere, KernelSpiky,KernelPoly6, Kernel, W, ∇W, Particle, KernelViscosity


function numerical_integration(k::Kernel{N},dr) where N
    integral = 0.
        for r = dr/2:dr:k.h
            r2 = r^2
            integral = integral + W(k,r2)*surface_hypersphere(N)*r^(N-1)*dr
        end
    return integral
end

function check_grad(k,rij)
    eps_ = 1e-8
    r = norm(rij)
    @test (W(k,(r+eps_)^2) - W(k,(r-eps_)^2))/(2eps_) * rij/r ≈ ∇W(k,rij)
    atol=10*eps_^2
end

@testset "SmoothedParticleHydrodynamics" begin
    dr = 0.000001
    h = 2.

    @test surface_hypersphere(2) ≈ 2*π
    @test surface_hypersphere(3) ≈ 4*π
    @test surface_hypersphere(4) ≈ 2*π^2
    @test surface_hypersphere(5) ≈ 8*π^2/3
    @test surface_hypersphere(6) ≈ π^3

    rij = [2.,3.,4.]

    for kernel in (KernelSpiky,KernelPoly6,KernelViscosity)
        for N = 2:4
            local k
            k = kernel(N,h)
            @test numerical_integration(k,dr) ≈ 1 atol=10*dr^2
            check_grad(k,rij)
        end
    end

    include("dam_break.jl")
end
