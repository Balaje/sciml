using Gridap
using ModelingToolkit
using SymbolicUtils

@syms x y z w
include("pde2gridap.jl")

# Define the LHS of the Poisson Equation
Dx=Differential(x);
Dy=Differential(y);
DxkDx=Dx(x^2*Dx(w));
DykDy=Dy(sin(x*y)*Dy(w));  ##---Bug in DykDy=Dy((x*y)*Dy(w)); ----

# -- Working on 3D.
#Dz=Differential(z);
#DzkDz=Dz(z^2*Dz(w))
#DD=DxkDx + DykDy + DzkDz; #LHS of the PDE

DD=DxkDx + DykDy; #LHS of the PDE

# 2D Gridap domain
domain = (0,1,0,1)
partition = (4,4)

f(x) = 0
dbc(x) = x[1]
@show dbc([π,π/2])

# Do it from PDE and find the Gridap affine operator
op1,symWF,symCoeff=sym2gridap.pde2gridapWF(DD, f, domain, partition, dbc)

# Full Gridap application
model = CartesianDiscreteModel(domain,partition)
order = 1;
reffe = ReferenceFE(lagrangian,Float64,order);
V0 = TestFESpace(model, reffe, conformity=:H1,dirichlet_tags="boundary");
U = TrialFESpace(V0,dbc); # FE Spaces
degree = 2; # Degree and domain
Ω = Triangulation(model);
dΩ = Measure(Ω,degree);
K1(x)=TensorValue(x[1]^2, 0, 0, sin(x[1]*x[2]));
ah(u,v) = ∫( K1⋅∇(v)⊙∇(u) )*dΩ;
bh(v) = ∫( v*f )*dΩ;
op2 = AffineFEOperator(ah,bh,U,V0);

# Verify the matrices
@show op1.op.matrix == op2.op.matrix
