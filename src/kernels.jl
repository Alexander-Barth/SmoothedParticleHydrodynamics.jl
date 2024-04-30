"""
    S = surface_hypersphere(N)
    Surface of a N-sphere. For example if N=3, its surface is 4π.
"""
function surface_hypersphere(N)
    if N == 1
        return 2.0
    elseif N == 2
        return 2*π
    else
        return 2*π*surface_hypersphere(N-2)/(N-2)
    end
end


"""
    k = KernelSpiky(n,h)


Kernel in `n` dimension with characteristic radius `h` given by the following expression:

```math
W(r,h) = c \\; (h - r)^3
```

for \$r \\le h\$ where \$r\$ is the radial distance and where the normalization coefficient \$c\$ is given by:

```math
\\frac{1}{S_n \\, c} = \\int_0^h\\left(h-r\\right)^3 r^{n-1}dr
=\\frac{h^{n+3}}{n}-\\frac{3h^{n+3}}{n+1}+\\frac{3h^{n+3}}{n+2}-\\frac{h^{n+3}}{n+3}
```

where \$S_n\$ is the surface of the unit n-sphere.

"""
function KernelSpiky(N,h::T) where T
    # https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B0%7D%5E%7Bh%7D%5Cleft(h-%20r%5Cright)%5E%7B3%7D%20r%5E%7Bn-1%7Ddr%20?or=input
    coeff = 1/((1/N - 3/(N+1) + 3/(N+2) - 1/(N+3)) * surface_hypersphere(N) * h^(N+3))
    coeff_grad = -3 * coeff
    KernelSpiky{N,T}(h,h^2,T(coeff),T(coeff_grad))
end

@inline W(k::KernelSpiky,r2) = k.coeff * (k.h - sqrt(r2))^3
@inline ∇W(k::KernelSpiky,rij,r = norm(rij)) = k.coeff_grad * (k.h - r)^2 * rij / r


# ```math
# W(r,h) = c \\; (h - r)^3 \\mbox{for } r \\le h
# ```

# where \$r\$ is the radial distance and where the normalization coefficient \$c\$ is given by:

# ```math
# \\frac{1}{S_n c} = \\int_0^h\\left(h-r\\right)^3 r^{n-1}dr
# =\\frac{h^{n+3}}{n}-\\frac{3h^{n+3}}{n+1}+\\frac{3h^{n+3}}{n+2}-\\frac{h^{n+3}}{n+3}

# where \$S_n\$ is the surface of the unit n-sphere.


#https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B0%7D%5E%7Bh%7D%5Cleft(h%5E%7B2%7D-%20r%5E%7B2%7D%5Cright)%5E%7B3%7D%20r%5E%7Bn-1%7Ddr%20?or=input
function KernelPoly6(N,h::T) where T
    coeff = 1/((1/N-3/(N+2)+3/(N+4)-1/(N+6)) * surface_hypersphere(N) * h^(N+6))
    coeff_grad = -6*coeff
    KernelPoly6{N,T}(h,h^2,T(coeff),T(coeff_grad))
end

@inline W(k::KernelPoly6,r2) =  k.coeff * (k.h² - r2)^3
@inline ∇W(k::KernelPoly6,rij,r = norm(rij)) = k.coeff_grad * (k.h² - r^2)^2 * rij
