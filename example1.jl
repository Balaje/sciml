using ModelingToolkit
using SymbolicUtils
using Gridap
import Gridap: ∇

    include("pde2gridap.jl")

@parameters x y z
@variables u(..)

Dx=Differential(x)
Dy=Differential(y)
Dz=Differential(z)
eq=Dx(Dx(u(x,y,z))) + Dy(Dy(u(x,y,z))) + Dz(Dz(u(x,y,z))) ~ 0;

#usol = exp(x+y)*sin(sqrt(2)*z)
bcs = [u(0,y,z) ~ exp(y)*sin(sqrt(2)*z),
       u(2,y,z) ~ exp(y+2)*sin(sqrt(2)*z),
       u(x,0,z) ~ exp(x)*sin(sqrt(2)*z),
       u(x,1,z) ~ exp(x+1)*sin(sqrt(2)*z),
       u(x,y,0) ~ 0,
       u(x,y,1) ~ exp(x+y)*sin(sqrt(2))];

domains = [x ∈ IntervalDomain(0.0,2.0),
           y ∈ IntervalDomain(0.0,1.0),
           z ∈ IntervalDomain(0.0,1.0)]

pdesys = PDESystem(eq, bcs, domains, [x,y,z], [u(x,y,z)]);
prob = sym2gridap.FEMProblem(pdesys,(10,10,10));

uh,Ω,op=sym2gridap.FEMSolve(prob);

writevtk(Ω,"res3d",cellfields=["uh"=>uh])
