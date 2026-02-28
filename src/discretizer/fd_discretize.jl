using MethodOfLines
using NonlinearSolve
import ModelingToolkit as Model
import Symbolics as Symb
using DomainSets
import ApproxFun as AF
import DifferentialEquations as DE
using Symbolics
using SparseArrays
using MacroTools: @capture, postwalk, prewalk

function transform_sym(ex::Vector{Any})
    #TODO: Figure out the most suitable supertype for the symbolic expressions
    return prewalk(Meta.parse(string(ex))) do instr
        if @capture(instr, Differential(x)(Differential(x)(s_(x, y))))
            return :(($s[i+1, j] - 2 * $s[i, j] + $s[i-1, j]) / dx^2)
        elseif @capture(instr, Differential(y)(Differential(y)(s_(x, y))))
            return :(($s[i, j+1] - 2 * $s[i, j] + $s[i, j-1]) / dy^2)
        elseif @capture(instr, Differential(x)(s_(x, y)))
            return :(($s[i+1, j] - $s[i-1, j]) / (2 * dx))
        elseif @capture(instr, Differential(y)(s_(x, y)))
            return :(($s[i, j+1] - $s[i, j-1]) / (2 * dy))
        elseif @capture(instr, s_(x, y))
            id = string(s)
            letter = id[1]
            number = parse(Int, id[2:end])
            return :($(Symbol("$(letter)_array"))[i, j, $(number)])
        elseif @capture(instr, s_[i_, j_])
            id = string(s)
            letter = id[1]
            number = parse(Int, id[2:end])
            return :($(Symbol("$(letter)_array"))[$(i), $(j), $(number)])
        elseif @capture(instr, x)
            return :(i * dx)
        elseif @capture(instr, y)
            return :(j * dy)
        end
        return instr
    end
end

function create_residual_function(N, equationsExpr, harmonics)
    equationsExprMapped::Vector{Expr} = []

    println(length(equationsExpr))
    println(harmonics)

    for H in 1:harmonics
        push!(equationsExprMapped, :(F_view[i, j,1,$(H)] = $(equationsExpr[2*H - 1]))) # Related to equationsExpr[2*H-1]
        push!(equationsExprMapped, :(F_view[i, j,2,$(H)] = $(equationsExpr[2*H]))) # Related to equationsExpr[2*H]

    end

    function_code = quote
        function residual!(F, U, p)
            dx, dy = p
            grid_size = $N * $N
            harmonicsNum = $harmonics


            A_array = reshape(@view(U[1:harmonicsNum*grid_size]), $N, $N, harmonicsNum)
            B_array = reshape(@view(U[harmonicsNum*grid_size+1:2*harmonicsNum*grid_size]), $N, $N, harmonicsNum)


            F_view = reshape(@view(F[1:end]), $N, $N, 2, harmonicsNum)


            # Inner points:
            for i in 2:$(N-1)
                for j in 2:$(N-1)
                    $(equationsExprMapped...)
                end
            end

            # BCs (Dirichlet for now):
            for H in 1:harmonicsNum
                F_view[1,:,1,H] .= A_array[1,:,H]; F_view[end,:,1,H] .= A_array[end,:,H]
                F_view[:,1,1,H] .= A_array[:,1,H]; F_view[:,end,1,H] .= A_array[:,end,H]
                F_view[1,:,2,H] .= B_array[1,:,H]; F_view[end,:,2,H] .= B_array[end,:,H]
                F_view[:,1,2,H] .= B_array[:,1,H]; F_view[:,end,2,H] .= B_array[:,end,H]
            end

            return F
        end
    end
    res = eval(function_code)
    return res
end
