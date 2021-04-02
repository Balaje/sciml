module sym2gridap

using ModelingToolkit
using SymbolicUtils
using SymbolicUtils.Rewriters
using Gridap
import Gridap: ∇

    @syms w x y z

# Function to determine if the input is a Differential
isDiff = T -> T isa Differential

# The IBP rule for single PD.
getargDiff = @rule (~x::isDiff)((~~w)) => ((~~w))*(~x);


"""
IBP(T, testfunc)
        Input:
            1) T = Full PDE.
                    Eg. Differential(x)(Differential(x)(w)*(x^2)) + Differential(y)(Differential(y)(w)*(y^2 + (x^2)*(y^2)))
            2) testFunc : @syms vh

"""
#Function to do the IBP on all terms
function IBP(T, testFunc)
    if(length(SymbolicUtils.arguments(T))==1)
        opsum=getargDiff(T).outer*getargDiff(T).inner(testFunc);
        return opsum[1];
    else
        opsum=Array{Any}(undef, length(SymbolicUtils.arguments(T)));
        m=1
        for X=SymbolicUtils.arguments(T)
            op=getargDiff(X).outer*getargDiff(X).inner(testFunc);
            opsum[m]=op[1];
            m+=1;
        end
        return opsum;
    end
end

"""

sym2coef(T)

    Input:
        1) T = Symbolic weak form terms

    Output:
        2) op: Symbolic coefficients
           order: Order in with the derivatives are arranged.
            [x => 1, y => 2, z => 3]
"""
function wf2coef(T)
    # Define some rules
    r1=@rule (~b)*(~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => (~~a)*(~b)
    r2=@rule (~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => ~~a

    r1_order=@rule (~~b)*(~w::isDiff)(~z)*(~~a) => (~w).x

    # Rule for non-linearity [To implement].
    # r3 = @rule (~x)*(~x::isDiff)(~y)*(~w::isDiff)(~z)*(~~a) => (~~a)*(~b)
    if((r1(T) != nothing) | (r2(T) != nothing))
        op=(r1(T) == nothing ? r2(T) : r1(T))
        return (op==Term{Number,Nothing}[]) ? 1 : prod(op);
    else
        op=Array{Any}(undef,length(T))
        op_order=Array{Any}(undef,length(T))
        count=1
        DD=Dict([x => 1, y => 2, z => 3])
        for term=T
            op_order[count]=DD[r1_order(term)]
            op[count]=(r1(term) == nothing ? r2(term) : r1(term))
            op[count] = (op[count]== Term{Number,Nothing}[]) ? 1 : prod(op[count])
            count+=1;
        end
        return op, op_order;
    end
end



"""

pde2gridapWF(LHS, RHS, dΩ, domain, partition, dbc)

    Input:
        1) LHS = Symbolic LHS of the PDE.
            Eg. Differential(x)(Differential(x)(w)*(x^2)) + Differential(y)(Differential(y)(w)*(y^2 + (x^2)*(y^2)))
        2) RHS = RHS load function Eg. f(x) = x[1]*x[2]
        3) domain = computational domain Eg. domain = (0,1,0,1)
        4) dbc = Dirichlet boundary condition Eg. dbc = x[1]

    Output:
        1) AffineFEOperator, containing the matrix-vector equation
        2) Symbolic Weak Form.
        3) Coefficients of the diffusion tensor.

"""
function pde2gridapWF(DD, f, domain, partition, dbc)
    @syms vh
    term1=IBP(DD, vh) # The weak form of the LHS (In symbolic weak form)
    coeffs,order = wf2coef(term1)
    # Arrange according to dict order
    term1[order]=term1
    coeffs[order]=coeffs

    # Convert syms to function
    coeff_func=Array{Any}(undef, length(coeffs))
    for m=1:length(coeffs)
        coeff_func[m] = a-> substitute(coeffs[m], Dict([x => a[1], y => a[2]]))
    end
    K1(x)=TensorValue(coeff_func[1](x), 0, 0, coeff_func[2](x)); # For 2D problems only

    # Build gridap discretization
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
    return op, sum(term1), coeffs
end

end
