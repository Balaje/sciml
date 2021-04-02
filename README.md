PDE to Gridap
======

Julia code to obtain weak form (using [Gridap.jl](https://github.com/gridap)) from an elliptic PDE (from [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl) **for a diagonal diffusion tensor**. Run `example.jl` for an example problem.

# Details
Setup:

``` julia
using ModelingToolkit
using Gridap
using SymbolicUtils

include("pde2gridap.jl");
```

Define the PDE using ModelingToolkit:

``` julia
@syms x y w vh
DD=Differential(y)(Differential(y)(w)*sin(x*y)) + Differential(x)(Differential(x)(w)*(x^2))
```

Function `sym2gridap.IBP()` does the integration by parts and outputs an array containing the weak form.
``` julia

julia> WF=sym2gridap.IBP(DD,vh) # Outputs the terms of the weak form with vh as the test function

2-element Array{Any,1}:
 Differential(y)(vh)*Differential(y)(w)*sin(x*y)
 Differential(x)(vh)*Differential(x)(w)*(x^2)
```

To get the coefficients and the order in which they are stored, use `sym2gridap.wf2coef()`. The order is used to define the diffusion tensor in Gridap.

``` julia
julia> rr,rr1=sym2gridap.wf2coef(WF)

julia> rr
2-element Array{Any,1}:
 sin(x*y)
 x^2

julia> rr1
2-element Array{Any,1}:
 2
 1

```

Finally, use `sym2gridap.pde2gridapWF()` to get the affine operator associated with the PDE.

``` julia
julia> f(x) = 0 #Griap style
julia> dbc(x) = x[1] #Gridap style

julia> op1,symWF,symCoeff=sym2gridap.pde2gridapWF(DD, f, (0,1,0,1), (4,4), dbc);

julia> op1.op.matrix #Gives the stiffness matrix

```

The output can be compared using the full Gridap implementation and can be found in `example.jl`.
