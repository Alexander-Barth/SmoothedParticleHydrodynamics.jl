# TODO

# https://en.wikipedia.org/wiki/Lanczos_approximation
function lanczos_gamma(z::Float64)
    g = 7
    n = 9
    p = (
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    )


    if z < 0.5
        y = π / (sin(π * z) * lanczos_gamma(1 - z))  # Reflection formula
    else
        z -= 1
        x = sum(
            ntuple(length(p)) do i
                if i == 1
                    @inbounds p[i]
                else
                    @inbounds p[i] / (z + i - 1)
                end
            end
        )
        t = z + g + 0.5
        y = sqrt(2 * π) * t ^ (z + 0.5) * exp(-t) * x
    end
    return y
end

surface_hypersphere(N) = 2*π^(N/2) / lanczos_gamma(N/2)

#gamma(N) = 1 # 2D
#surface_hypersphere(N) = 2*π^(N/2) / gamma(N/2)


#KernelSpiky(N,h) = KernelSpiky{N,typeof(h)}(h,h^2,15/(π*h^6),-10 / (π * h^5))

"""
    k = KernelSpiky(n,h)


Kernel in `n` dimension with characteristic radius `h` given by the following expression:

```math
W(r,h) = c \\; (h - r)^3 \\mbox{for } r \\le h
```

where \$r\$ is the radial distance and where the normalization coefficient \$c\$ is given by:

```math
\\frac{1}{S_n c} = \\int_0^h\\left(h-r\\right)^3 r^{n-1}dr
=\\frac{h^{n+3}}{n}-\\frac{3h^{n+3}}{n+1}+\\frac{3h^{n+3}}{n+2}-\\frac{h^{n+3}}{n+3}

where \$S_n\$ is the surface of the unit n-sphere.
```
"""
function KernelSpiky(N,h::T) where T
    # https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B0%7D%5E%7Bh%7D%5Cleft(h-%20r%5Cright)%5E%7B3%7D%20r%5E%7Bn-1%7Ddr%20?or=input
    coeff = 1/((1/N - 3/(N+1) + 3/(N+2) - 1/(N+3)) * surface_hypersphere(N) * h^(N+3))
    coeff_grad = -3 * coeff
    KernelSpiky{N,T}(h,h^2,T(coeff),T(coeff_grad))
end

@inline W(k::KernelSpiky,r2) = k.coeff * (k.h - sqrt(r2))^3
@inline ∇W(k::KernelSpiky,rij,r = norm(rij)) = k.coeff_grad * (k.h - r)^2 * rij / r




#https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B0%7D%5E%7Bh%7D%5Cleft(h%5E%7B2%7D-%20r%5E%7B2%7D%5Cright)%5E%7B3%7D%20r%5E%7Bn-1%7Ddr%20?or=input
function KernelPoly6(N,h::T) where T
    coeff = 1/((1/N-3/(N+2)+3/(N+4)-1/(N+6)) * surface_hypersphere(N) * h^(N+6))
    coeff_grad = -6*coeff
    KernelPoly6{N,T}(h,h^2,T(coeff),T(coeff_grad))
end

@inline W(k::KernelPoly6,r2) =  k.coeff * (k.h² - r2)^3
@inline ∇W(k::KernelPoly6,rij,r = norm(rij)) = k.coeff_grad * (k.h² - r^2)^2 * rij
