var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SmoothedParticleHydrodynamics","category":"page"},{"location":"#SmoothedParticleHydrodynamics","page":"Home","title":"SmoothedParticleHydrodynamics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SmoothedParticleHydrodynamics.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The fluid is represented by a discrete particles centred at the location mathbf r_i:","category":"page"},{"location":"","page":"Home","title":"Home","text":"fracd mathbf r_idt = mathbf v_i","category":"page"},{"location":"","page":"Home","title":"Home","text":"where mathbf v_i is the velocity of the i-th particle.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SmoothedParticleHydrodynamics]","category":"page"},{"location":"#SmoothedParticleHydrodynamics.KernelSpiky-Union{Tuple{T}, Tuple{Any, T}} where T","page":"Home","title":"SmoothedParticleHydrodynamics.KernelSpiky","text":"k = KernelSpiky(n,h)\n\nKernel in n dimension with characteristic radius h given by the following expression:\n\nW(rh) = c  (h - r)^3\n\nfor r le h where r is the radial distance and where the normalization coefficient c is given by:\n\nfrac1S_n  c = int_0^hleft(h-rright)^3 r^n-1dr\n=frach^n+3n-frac3h^n+3n+1+frac3h^n+3n+2-frach^n+3n+3\n\nwhere S_n is the surface of the unit n-sphere.\n\n\n\n\n\n","category":"method"},{"location":"#SmoothedParticleHydrodynamics.Location","page":"Home","title":"SmoothedParticleHydrodynamics.Location","text":"loc = Location(p)\n\nAdapter type for SpatialHashing such that loc[i] returns p.x[i] where p is a list of particles.\n\n\n\n\n\n","category":"type"}]
}
