using ModelingToolkit
using SymbolicUtils
using Gridap
import Gridap: ∇

include("pde2gridap.jl")

@parameters a b
@variables u(..)

#Dxx= Differential(x)^2
#Dyy= Differential(y)^2
#eq = Dxx(u(x,y)) ~ 0

Dx= Differential(a);
Dy= Differential(b);
#eq= Dx(Dx(u(a,b))) + Dy(Dy(u(a,b))) ~ 0;
#eq= Dx(a^2*Dx(u(a,b))) + Dy(b^2*Dy(u(a,b)))~ b^2*exp(a)*sin(b) - 2*a*exp(a)*sin(b) - a^2*exp(a)*sin(b) - 2*b*exp(a)*cos(b)
eq = Dx((a^2+b^2*exp(a*b))*Dx(u(a,b))) + Dy(b^2*exp(a*b)*Dy(u(a,b))) ~ b^2*exp(a*b)*exp(a)*sin(b) - exp(a)*sin(b)*(b^2*exp(a*b) + a^2) - exp(a)*sin(b)*(2*a + b^3*exp(a*b)) - 2*b*exp(a*b)*exp(a)*cos(b) - a*b^2*exp(a*b)*exp(a)*cos(b)

bcs = [u(0,b) ~ sin(b),
       u(2,b) ~ exp(2)*sin(b),
       u(a,0) ~ 0,
       u(a,1) ~ exp(a)*sin(1)]

domains = [a ∈ IntervalDomain(0.0,2.0),
           b ∈ IntervalDomain(0.0,1.0)]

pdesys = PDESystem(eq,bcs,domains,[a,b],[u(a,b)])

uh,Ω,operator = sym2gridap.FEMProblem(pdesys,(50,50))
writevtk(Ω,"results",cellfields=["uh"=>uh])
