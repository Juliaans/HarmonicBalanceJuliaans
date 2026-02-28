using MethodOfLines
using NonlinearSolve
import ModelingToolkit as Model
import Symbolics as Symb
using DomainSets
import ApproxFun as AF
import DifferentialEquations as DE
using Symbolics
using SparseArrays

function create_ansatz(coords::Tuple, t::Symbolics.Num, omega, harmonics::Int, n_fields::Int=1)
    var_names = Symbol[]
    var_exprs = Symbolics.Num[]
    fields = Symbolics.Num[]
    
    letter_idx = 1
    
    for field_idx in 1:n_fields
        u = Num(0)
        j = 1
        base_char = Char('A' + 2*(field_idx - 1))
        for i in 1:(2*harmonics)
            if isodd(i)
                name = Symbol(Char(base_char), div(i + 1, 2))
            else
                name = Symbol(Char(base_char + 1), div(i, 2))
            end
            letter_idx += 1
            
            v = first(@variables $name(..))
            expr = v(coords...)
            
            if isodd(i)
                u += expr * sin(j * omega * t)
            else
                u += expr * cos(j * omega * t)
                j += 1
            end
            
            push!(var_names, name)
            push!(var_exprs, expr)
        end
    
        push!(fields, u)
    end
    
    return var_names, var_exprs, fields
end

function expand_trig_jl(eqn, t, omega)
    y_exp = Symb.expand(Model.expand_derivatives(eqn))
    symbolics_list = Symb.arguments(y_exp, +)
    contains_var(expr, var) = any(v -> isequal(v, var), Symbolics.get_variables(expr))
    finished_terms = Num[]
    
    for (i, term) in enumerate(symbolics_list)
        trig_terms = Num[]
        spatial_terms = Num[]
        
        for mul_term in Symb.arguments(term, *)
            mul_term_num = Num(mul_term)
            if contains_var(mul_term_num, t)
                push!(trig_terms, mul_term_num)
            else
                push!(spatial_terms, mul_term_num)
            end
        end
        
        spatial::Num = isempty(spatial_terms) ? Num(1) : prod(spatial_terms)
        
        if length(trig_terms) == 1
            unwrapped = Symbolics.unwrap(trig_terms[1])
            if SymbolicUtils.operation(unwrapped) !== ^
                push!(finished_terms, Num(term))
                continue
            end
        end
        
        trig::Num = isempty(trig_terms) ? Num(1) : prod(trig_terms)
        
        trig_func = Symbolics.build_function(trig, t, expression=Val{false})
        # trig_func = x -> Symbolics.value(Symbolics.substitute(trig, t => x))
        period = 2π / omega
        F = AF.Fun(trig_func, AF.Fourier(-period/2 .. period/2))
        coeffs = AF.coefficients(F)::Vector{Float64}
        
        ωt = omega * t
        expanded_trig::Num = Num(0)
        
        for (j, c) in enumerate(coeffs)
            abs(c) < 1e-10 && continue
            
            if j == 1
                expanded_trig += c
            elseif iseven(j)
                n = j ÷ 2
                expanded_trig += c * cos(n * ωt)
            else
                n = (j + 1) ÷ 2
                expanded_trig += c * sin(n * ωt)
            end
        end
        
        push!(finished_terms, spatial * -expanded_trig)
    end
    
    return sum(finished_terms)
end

function make_equations(expanded, harmonics, omega, t)
    eqs = []
    for i in 1:harmonics
        sin_coef = Symb.coeff(expanded, sin(i*omega*t))
        cos_coef = Symb.coeff(expanded, cos(i*omega*t))
        push!(eqs, sin_coef)
        push!(eqs, cos_coef)
    end
    return eqs
end

struct HarmonicBalanceProblem
    coords::Tuple
    t::Symbolics.Num
    omega::Float64
    harmonics::Int
    n_fields::Int
    mkpdes::Function
    params::Any
end


function harmonic_balance(hbp::HarmonicBalanceProblem)
    var_names, var_exprs, fields = create_ansatz(hbp.coords, hbp.t, hbp.omega, hbp.harmonics, hbp.n_fields)

    pdes = hbp.mkpdes(hbp.params, fields...)
    
    equations = []

    for pde in pdes
        expanded = expand_trig_jl(pde, hbp.t, hbp.omega)
        eqns = make_equations(expanded, hbp.harmonics, hbp.omega, hbp.t)
        append!(equations, eqns)
    end

    return var_names, var_exprs, equations
end