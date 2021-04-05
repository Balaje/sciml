module sym2gridap

using ModelingToolkit
using SymbolicUtils
using SymbolicUtils.Rewriters
using Gridap
import Gridap: ∇

    @syms w x y z
@syms vh boundary_normal

# Function to determine if the input is a Differential
isDiff = T -> T isa Differential

## --- List of rules for FEM algebra
# The IBP rule for single PD.
IBP_rule_vol = @rule (~x::isDiff)((~~w)) => ((~~w))*(~x);
IBP_rule_face = @rule (~x::isDiff)((~~w)) => ((~~w))*boundary_normal;
# Rules for determining coefficients.
# TODO: Try to test more rules.
r1=@rule (~b)*(~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => (~~a)*(~b)
r2=@rule (~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => ~~a
r3=@rule (~x::isDiff)(~y)*(~b)*(~w::isDiff)(~z)*(~~a) => (~~a)*(~b)
r4=@rule (~x::isDiff)(~y)*(~~a)*(~w::isDiff)(~z) => (~~a)
# Rule to obtain the order of differentiation in the test function
r1_order=@rule (~~b)*(~w::isDiff)(~z)*(~~a) => (~w).x

"""
IBP(T, testfunc)
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

wf2coef(T)

    Input:
        1) T = Symbolic weak form terms

    Output:
        2) op: Symbolic coefficients
           order: Order in with the derivatives are arranged.
            [x => 1, y => 2, z => 3]
"""
function wf2coef(T, indvars)
    # Define some rules

    # Rule for non-linearity [To implement].
    # r3 = @rule (~x)*(~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => (~~a)*(~b)
    if((r1(T) != nothing) | (r2(T) != nothing))
        op=(r1(T) == nothing ? r2(T) : r1(T))
        op=(op == nothing ? r3(T) : op)
        op=(op == nothing ? r4(T) : op)
        return (op==Term{Number,Nothing}[]) ? 1 : prod(op);
    else
        op=Array{Any}(undef,length(T))
        op_order=Array{Any}(undef,length(T))
        count=1
        DD=Dict([indvars[1] => 1, indvars[2] => 2])
        for term=T
            op_order[count]=DD[r1_order(term)]
            op[count]=(r1(term) == nothing ? r2(term) : r1(term))
            op[count]=(op[count] == nothing ? r3(term) : op[count])
            op[count]=(op[count] == nothing ? r4(term) : op[count])

            op[count] = (op[count]== Term{Number,Nothing}[]) ? 1 : prod(op[count])
            count+=1;
        end
        return op, op_order;
    end
end


"""

pde2gridapWF(LHS, RHS, dΩ, domain, partition, dbc)

    Input:
        1) LHS = LHS of the PDE, PDESystem.lhs
        2) RHS = RHS load function, PDESystem.rhs
        3) domain = computational domain Eg. domain = (0,1,0,1)
        4) dbc = Dirichlet boundary condition Eg. dbc = x[1]

    Output:
        1) AffineFEOperator, containing the matrix-vector equation
        2) Symbolic Weak Form.
        3) Coefficients of the diffusion tensor.

"""
function pde2gridapWF(pdelhs, pderhs, domain, partition, pdebcs, indvars)
    term1=IBP(pdelhs) # The weak form of the LHS (In symbolic weak form)
    coeffs,order = wf2coef(term1, indvars)
    # Arrange according to dict order
    term1[order]=term1
    coeffs[order]=coeffs


    # Convert syms to function
    coeff_func=Array{Any}(undef, length(coeffs))
    for m=1:length(coeffs)
        coeff_func[m] = a-> substitute(coeffs[m], Dict([indvars[1] => a[1], indvars[2] => a[2]]))
    end
    K1(x)=TensorValue(coeff_func[1](x), 0, 0, coeff_func[2](x)); # For 2D problems only
    # Construct RHS function
    f = a -> substitute(pderhs, Dict([indvars[1] => a[1], indvars[2] => a[2]]));

    #Function to get the DBC to transfer to gridap
    function dbc(a)
        bcvals=zeros(Float64,length(pdebcs))
        rep_count=0;
        for m=1:length(pdebcs)
            bc=pdebcs[m];

            bcargs_x=a-> substitute(bc.lhs.arguments[1],Dict([indvars[1]=>a[1], indvars[2]=>a[2]]))
            bcargs_y=a-> substitute(bc.lhs.arguments[2],Dict([indvars[1]=>a[1], indvars[2]=>a[2]]))
            # --- TODO: Add condition to check if the values are within domain
            # ---
            bcfunc = a-> substitute(bc.rhs, Dict([indvars[1]=>a[1], indvars[2]=>a[2]]))

            bcvals[m]=convert(Int64,(a[1] == bcargs_x(a)))*
                convert(Int64,(a[2] == bcargs_y(a)))*bcfunc(a);
            if((a[1]==bcargs_x(a)) & (a[2]==bcargs_y(a)))
                rep_count+=1;
            end
        end
        return sum(bcvals)/rep_count;
    end

    # Build gridap discretization
    print("Building FEM Problem")
    model = CartesianDiscreteModel(domain,partition)
    order = 1;
    reffe = ReferenceFE(lagrangian,Float64,order);
    V0 = TestFESpace(model, reffe, conformity=:H1,dirichlet_tags="boundary");
    U = TrialFESpace(V0,dbc); # FE Spaces
    degree = 2; # Degree and domain
    Ω = Triangulation(model);
    dΩ = Measure(Ω,degree);
    ah(u,v) = ∫( K1⋅∇(v)⊙∇(u) )*dΩ;
    lh(v) = ∫( v*f )*dΩ;
    op = AffineFEOperator(ah,lh,U,V0);
    return op, Ω
end

function FEMProblem(pdesys,partition)
    pde=pdesys.eqs;
    pdelhs=pde.lhs;
    pderhs=pde.rhs;
    pdebcs=pdesys.bcs;
    indvars=pdesys.indvars

    # Get Gridap style domain 2D/3D
    if(length(indvars)==2)
        Gridapdomain=((pdesys.domain[1].variables == convert(Sym,indvars[1]))*(pdesys.domain[1].domain.lower),
                      (pdesys.domain[1].variables == convert(Sym,indvars[1]))*(pdesys.domain[1].domain.upper),
                      (pdesys.domain[2].variables == convert(Sym,indvars[2]))*(pdesys.domain[2].domain.lower),
                      (pdesys.domain[2].variables == convert(Sym,indvars[2]))*(pdesys.domain[2].domain.upper));
    else
        Gridapdomain=((pdesys.domain[1].variables == convert(Sym,indvars[1]))*(pdesys.domain[1].domain.lower),
                      (pdesys.domain[1].variables == convert(Sym,indvars[1]))*(pdesys.domain[1].domain.upper),
                      (pdesys.domain[2].variables == convert(Sym,indvars[2]))*(pdesys.domain[2].domain.lower),
                      (pdesys.domain[2].variables == convert(Sym,indvars[2]))*(pdesys.domain[2].domain.upper),
                      (pdesys.domain[3].variables == convert(Sym,indvars[3]))*(pdesys.domain[3].domain.lower),
                      (pdesys.domain[3].variables == convert(Sym,indvars[3]))*(pdesys.domain[3].domain.upper));
    end


    operator,Ω=pde2gridapWF(pdelhs,pderhs,Gridapdomain,partition,pdebcs,indvars);

    uh=solve(operator)
    return uh,Ω
end

end
