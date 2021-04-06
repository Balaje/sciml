module sym2gridap

using ModelingToolkit
using SymbolicUtils
using SymbolicUtils.Rewriters
using Gridap
import Gridap: ∇

    @syms vh boundary_normal

# Constructor.
struct discreteFEM
    #K=diffusion tensor, f=right hand side
    domain; partition; dbc; K; f;
end

# Function to determine if the input is a Differential
isDiff = T -> T isa Differential

"""
Rules for IBP
"""
## The IBP rule for single PD.
IBP_rule_vol = @rule (~x::isDiff)((~~w)) => ((~~w))*(~x);
# For nbc/rbc [TODO]
IBP_rule_face = @rule (~x::isDiff)((~~w)) => ((~~w))*boundary_normal;

"""
Rules for determining coefficients. [TODO: Add more rules - convection term, reaction term]
"""
# Coefficient K(x,y) -∇(K(x,y)∇u)
r1=@rule (~~b)*(~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => [(~~b);(~~a)]
r2=@rule (~x::isDiff)(~y)*(~~b)*(~w::isDiff)(~z)*(~~a) => [(~~a);(~~b)]

# Rule for non-linearity [To implement].
r_nl = @rule (~y)*(~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => (~~a)*(~y)

# Rule to obtain the order of differentiation in the test function
r1_order=@rule (~~b)*(~w::isDiff)(~z)*(~~a) => (~w).x

"""
IBP(T)
        Input:
            1) T = Full PDE.
                    Eg. Differential(x)(Differential(x)(w)*(x^2)) + Differential(y)(Differential(y)(w)*(y^2 + (x^2)*(y^2)))
            2) testFunc : @syms vh

"""
#Function to do the IBP on all terms
function IBP(T)
    if(length(SymbolicUtils.arguments(T))==1)
        opsum=IBP_rule_vol(T).outer*IBP_rule_vol(T).inner(vh)
        return opsum;
    else
        opsum=Array{Any}(undef, length(SymbolicUtils.arguments(T)));
        m=1
        for X=SymbolicUtils.arguments(T)
            op=IBP_rule_vol(X).outer*IBP_rule_vol(X).inner(vh)
            opsum[m]=op[1];
            m+=1;
        end
        return opsum;
    end
end

"""

wf2coef(T, indvars)

    Input:
        1) T = Symbolic weak form terms
        2) indvars = independent variables

    Output:
        2) op: Symbolic coefficients
           order: Order in with the derivatives are arranged.
            [indvars[1] => 1, indvars[2] => 2]
"""

function wf2coef(T, indvars)
    if((r1(T) != nothing) | (r2(T) != nothing))
        op=(r1(T) == nothing ? r2(T) : r1(T)) # Check rule
        return (op==Term{Number,Nothing}[]) ? 1 : prod(op);
    else
        op=Array{Any}(undef,length(T))
        op_order=Array{Any}(undef,length(T))
        count=1
        DD=Dict(indvars[i]=>i for i=1:length(indvars))
        for term=T
            op_order[count]=DD[r1_order(term)]
            # Checking different rules to determine the coefficients.
            op[count]=(r1(term) == nothing ? r2(term) : r1(term))
            op[count] = (op[count]== Term{Number,Nothing}[]) ? 1 : prod(op[count])
            count+=1;
        end
        return op, op_order;
    end
end


"""

discretize_pde_fem(LHS, RHS, dΩ, domain, partition, dbc, nbc=0)

    Only 2D implementation done
    Input:
        1) LHS = LHS of the PDE, PDESystem.lhs
        2) RHS = RHS load function, PDESystem.rhs
        3) domain = computational domain Eg. domain = (0,1,0,1)
        4) partition = partition along x and y axes Eg. (4,4)
        5) dbc = Dirichlet boundary conditions, pdesys.dbc
        6) nbc = Neumann boundary conditions

    Output:
        1) AffineFEOperator, containing the matrix-vector equation
        2) Symbolic Weak Form.
        3) Coefficients of the diffusion tensor.

"""
function discretize_pde_fem(pdelhs, pderhs, domain, partition, pdebcs, indvars, nbc=0)

    # The weak form of the LHS (In symbolic weak form)
    term1=IBP(pdelhs)
    # Obtain coefficients and order
    coeffs,order = wf2coef(term1, indvars)
    # Arrange according to order
    term1[order]=term1
    coeffs[order]=coeffs

    @show coeffs # Show Coeffs

    # Dictionary for substitution
    DD=a->Dict(indvars[i]=>a[i] for i=1:length(indvars))

    # Convert variables to Gridap function
    coeff_func=Array{Any}(undef, length(coeffs))
    for m=1:length(coeffs)
        coeff_func[m] = a-> substitute(coeffs[m], DD(a))
    end

    ## 2D/3D
    diffusion_func = a->
        if length(indvars)==2
            (coeff_func[1](a), 0, 0, coeff_func[2](a)); # For 2D problems only
        elseif length(indvars) == 3
            (coeff_func[1](a), 0, 0, 0, coeff_func[2](a), 0, 0, 0, coeff_func[3](a)); # For 3D problems only
        end
    K1(a)=TensorValue(diffusion_func(a))
    #@show K1([1,2,3])


    # RHS function. Will contain Neumann
    f = a -> substitute(pderhs, DD(a));

    # Filter Dirichlet boundary conditions. TODO: Neumann BC
    dbc=filter(x -> isDiff(x.lhs)==false, pdebcs)

    #DBC to transfer to gridap
    # TODO: Using labels instead of functions like u(1,x).
    function dirichletbc(a)
        bcvals=ones(Float64,length(pdebcs))
        rep_count=0
        for m=1:length(pdebcs)
            bc=dbc[m];
            bcfunc = a-> substitute(bc.rhs, DD(a))
            bcargs=zeros(Float64,length(indvars))
            for n=1:length(indvars)
                bcargs[n]=substitute(bc.lhs.arguments[n],DD(a))
            end
            bcvals[m]*=convert(Int64,(tuple(a ...) == tuple(bcargs ...)));
            rep_count+=convert(Int64,(tuple(a ...) == tuple(bcargs ...)));
            bcvals[m]*=bcfunc(a);
        end
        return sum(bcvals)/rep_count;
    end

    problem=discreteFEM(domain, partition, dirichletbc, K1, f); # Constructor
    return problem;
end

function FEMSolve(problem)
    # Build Gridap discretization
    print("Building FEM Problem \n")
    model = CartesianDiscreteModel(problem.domain,problem.partition)
    order = 1;
    reffe = ReferenceFE(lagrangian,Float64,order);
    V0 = TestFESpace(model, reffe, conformity=:H1,dirichlet_tags="boundary");
    U = TrialFESpace(V0,problem.dbc); # FE Spaces
    degree = 2; # Degree and domain
    Ω = Triangulation(model);
    dΩ = Measure(Ω,degree);
    ah(u,v) = ∫( problem.K⋅∇(v)⊙∇(u) )*dΩ;
    lh(v) = ∫( v*problem.f )*dΩ;
    op = AffineFEOperator(ah,lh,U,V0);
    uh=solve(op)
    return uh,Ω,op
end

function FEMProblem(pdesys,partition)
    pde=pdesys.eqs;
    pdelhs=pde.lhs;
    pderhs=pde.rhs;
    pdebcs=pdesys.bcs;
    indvars=pdesys.indvars

    # Get Gridap style domain in 2D. (TODO: 3D)
    domain=[(pdesys.domain[1].domain.lower); (pdesys.domain[1].domain.upper)];
    for m=2:length(indvars)
        domain=[ domain; (pdesys.domain[m].domain.lower);
                 (pdesys.domain[m].domain.upper) ];
    end
    Gridapdomain = tuple(domain ...)
    @show Gridapdomain

    prob=discretize_pde_fem(pdelhs,pderhs,Gridapdomain,partition,pdebcs,indvars)

    return prob
end

end
