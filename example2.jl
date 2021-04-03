using ModelingToolkit
using SymbolicUtils
using Gridap
import Gridap: ∇

include("pde2gridap.jl")

@parameters x y
@variables u(..)

#Dxx= Differential(x)^2
#Dyy= Differential(y)^2
#eq = Dxx(u(x,y)) ~ 0

Dx= Differential(x);
Dy= Differential(y);
eq= Dx(sin(x*y)*Dx(u(x,y))) + Dy(x^2*Dy(u(x,y))) ~ 0;

DD=eq.lhs


l = a -> substitute(eq.rhs, Dict([x => a[1], y => a[2]]))
dbc(x) = x[1]
op1,symWF,symCoeff=sym2gridap.pde2gridapWF(DD, l, (0,1,0,1), (4,4), dbc);
