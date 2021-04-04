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
#eq= Dx(Dx(u(x,y))) + Dy(Dy(u(x,y))) ~ 0;
eq= Dx(x^2*Dx(u(x,y))) + Dy(y^2*Dy(u(x,y)))~ -sin(y)*(2*x*exp(x) + x^2*exp(x)) - exp(x)*(y^2*cos(y)+2*y*sin(y))
bcs = [u(0,y) ~ sin(y),
       u(1,y) ~ exp(1)*sin(y),
       u(x,0) ~ 0,
       u(x,1) ~ exp(x)*sin(1)]

domains = [x ∈ IntervalDomain(0.0,1.0),
           y ∈ IntervalDomain(0.0,1.0)]

pdesys = PDESystem(eq,bcs,domains,[x,y],[u(x,y)])

uh,Ω = sym2gridap.FEMProblem(pdesys,(4,4))
writevtk(Ω,"results",cellfields=["uh"=>uh])
