```@meta
CurrentModule = SmoothedParticleHydrodynamics
```

# SmoothedParticleHydrodynamics

Documentation for [SmoothedParticleHydrodynamics](https://github.com/Alexander-Barth/SmoothedParticleHydrodynamics.jl).


The fluid is represented by a set of discrete particles. Each particle is centered at the location $\mathbf r_i$. The position varies in time:

```math
\frac{d \mathbf r_i}{dt} = \mathbf v_i
```

where $\mathbf v_i$ is the velocity of the i-th particle.
Each particle can have a property $A$:

```math
A(\mathbf r) = \int A(\mathbf r') W(|\mathbf r - \mathbf r'|,h) d^n r'
```

The function $W$ is the kernel function which has a characteristic length-scale $h$
and it normalized:

```math
\int W(|\mathbf r - \mathbf r'|,h) d^n r' = 1
```

and converges to the Dirac function $\delta(\mathbf r - \mathbf r')$ if $h$ tends to zero.

```math
\lim_{h \rightarrow 0} \int W(|\mathbf r - \mathbf r'|,h) d^n r' = \delta(\mathbf r - \mathbf r')
```

For small $h$ enought and sufficiently many particles, we can therefore approximate the integral using the discrete sum:

```math
A(\mathbf r) = \sum_i V_i A_i W(|\mathbf r - \mathbf r_i|,h)
```

where $V_i$ is the volume associated to the i-th particle. Knowing the mass $m_j$ of every particle, we can calculate the density as:

```math
\rho(r) = \sum_j m_j W(|\mathbf r - \mathbf r_j|,h)
```

since $m_i = V_i \rho_i$. The density at the particle location $\mathbf r_i$ becomes:


```math
\rho_i = \rho(\mathbf r_i) = \sum_j m_j W(|\mathbf r_i - \mathbf r_j|,h)
```

If we differentiate the density relative to time:

```math
\frac{d\rho_i}{dt} = \sum_j m_j (\mathbf v_i - \mathbf v_j)
\nabla W_{ij}
```

which can be seen as a discretized form of:

```math
\frac{d\rho}{dt} = -\rho \nabla \cdot v
```

The pressure gradient is computed using the symmetric form by analogy of the evolution equation of density:

```math
\nabla p(r_i) = \rho_i \sum_j m_j
\left(\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}\right)
∇W(|\mathbf r - \mathbf r_i|,h)
```

Similarly the viscosity force is given by:

```math
F_i = \mu \nabla^2 \mathbf v(r_i)
\mu
\sum_j m_j \frac{\mathbf v_j - \mathbf v_i}{\rho_j}

\left(\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}\right)
∇²W
```



```@index
```

```@autodocs
Modules = [SmoothedParticleHydrodynamics]
```
