using ModelingToolkit
using SymbolicUtils
using Gridap
import Gridap: ∇

    @parameters x y z
@variables u(..)

Dx=Differential(x)
Dy=Differential(y)
Dz=Differential(z)
eq=Dx(Dx(u(x,y,z)))+Dy(Dy(u(x,y,z)))+Dz(Dz(u(x,y,z))) + sin(x)*cos(x)*u(x,y,z)~ 0;

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


## Define some rules
rule1 = @rule +(~~y) => sum(filter(x -> x isa SymbolicUtils.Mul,(~~y)))
rule2 = @rule +(~~y) => sum(filter(x -> x.f isa Differential,(~~y)))
rule3 = @rule +(~~y) => sum(filter(x -> !(x.f isa Differential),(~~y)))

eqs_multerms = rule1(pdesys.eqs.lhs)
@show eqs_multerms

updated_eqs = pdesys.eqs.lhs - eqs_multerms
eqs_diffterms = rule2(updated_eqs)
@show eqs_diffterms

eqs_otherterms = try rule3(updated_eqs) catch e 0 end
@show eqs_otherterms

# Define integration by parts
