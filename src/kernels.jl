# TODO
#gamma(N) = 1 # 2D
surface_hypersphere(N) = 2*π^(N/2) / gamma(N/2)


#KernelSpiky(N,h) = KernelSpiky{N,typeof(h)}(h,h^2,15/(π*h^6),-10 / (π * h^5))

"""
```math
\\int_0^h\\left(h-r\\right)^3r^{n-1}dr
=\\frac{h^{n+3}}{n}-\\frac{3h^{n+3}}{n+1}+\\frac{3h^{n+3}}{n+2}-\\frac{h^{n+3}}{n+3}
```

https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B0%7D%5E%7Bh%7D%5Cleft(h-%20r%5Cright)%5E%7B3%7D%20r%5E%7Bn-1%7Ddr%20?or=input

"""
function KernelSpiky(N,h::T) where T
    coeff = 1/((1/N - 3/(N+1) + 3/(N+2) - 1/(N+3)) * surface_hypersphere(N) * h^(N+3))
    coeff_grad = -3 * coeff
    KernelSpiky{N,T}(h,h^2,coeff,coeff_grad)
end

@inline W(k::KernelSpiky,r2) = k.coeff * (k.h - sqrt(r2))^3
@inline ∇W(k::KernelSpiky,rij,r = norm(rij)) = k.coeff_grad * (k.h - r)^2 * rij / r




#https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B0%7D%5E%7Bh%7D%5Cleft(h%5E%7B2%7D-%20r%5E%7B2%7D%5Cright)%5E%7B3%7D%20r%5E%7Bn-1%7Ddr%20?or=input
function KernelPoly6(N,h::T) where T
    coeff = 1/((1/N-3/(N+2)+3/(N+4)-1/(N+6)) * surface_hypersphere(N) * h^(N+6))
    coeff_grad = -6*coeff
    KernelPoly6{N,T}(h,h^2,coeff,coeff_grad)
end

@inline W(k::KernelPoly6,r2) =  k.coeff * (k.h² - r2)^3
@inline ∇W(k::KernelPoly6,rij,r = norm(rij)) = k.coeff_grad * (k.h² - r^2)^2 * rij

