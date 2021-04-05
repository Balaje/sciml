PDE to Gridap
======

Julia code to obtain weak form (using [Gridap.jl](https://github.com/gridap)) for the Poisson equation (defined using [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl), *for a diagonal diffusion tensor*. Run `example.jl` for an example problem.

The method is based on the `@rules` macro from the `SymbolicUtils.jl` package


```Julia
isDiff = T -> T isa Differential

rule1 = @rule (~x::isDiff)((~~w)) => ((~~w))*(~x); # Integration by parts @rule
# To get coefficients
rule2 = @rule (~b)*(~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => (~~a)*(~b) # Rule 1
rule3 = @rule (~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => ~~a   # Rule 2
```

# Details
Setup:

``` julia
using ModelingToolkit
using Gridap
using SymbolicUtils

include("pde2gridap.jl");
```

Define the PDESystem using ModelingToolkit:

``` julia
@parameters a b  
@variables u(..)

Dx= Differential(a);  
Dy= Differential(b);

eq= Dx(a^2*Dx(u(a,b))) + Dy(b^2*Dy(u(a,b)))~ b^2*exp(a)*sin(b) - 2*a*exp(a)*sin(b) - a^2*exp(a)*sin(b) - 2*b*exp(a)*cos(b)
bcs = [u(0,b) ~ sin(b),
       u(1,b) ~ exp(1)*sin(b),
       u(a,0) ~ 0,
       u(a,1) ~ exp(a)*sin(1)]

domains = [a ∈ IntervalDomain(0.0,1.0),
           b ∈ IntervalDomain(0.0,1.0)]

pdesys = PDESystem(eq,bcs,domains,[a,b],[u(a,b)])

```

Solve using `FEMProblem`
``` julia
uh,Ω,operator = sym2gridap.FEMProblem(pdesys,(50,50)) # (50,50) partition
writevtk(Ω,"results",cellfields=["uh"=>uh])  # Visualize using Paraview   
```
