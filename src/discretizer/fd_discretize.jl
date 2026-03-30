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

function transform_sym_2D(exprs)
    #TODO - figure out the most suitable supertype for the symbolic expressions
    arr = []
    for exp in exprs
        ne = prewalk(Meta.parse(string(exp))) do instr
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
        push!(arr, ne)
    end
    return arr
end


function transform_sym_1D(exprs)
    #TODO - figure out the most suitable supertype for the symbolic expressions
    arr = []
    for exp in exprs
        ne = prewalk(Meta.parse(string(exp))) do instr
            if @capture(instr, Differential(x)(Differential(x)(s_(x))))
                return :(($s[i+1] - 2 * $s[i] + $s[i-1]) / dx^2)
            elseif @capture(instr, Differential(x)(s_(x)))
                return :(($s[i+1] - $s[i-1]) / (2 * dx))
            elseif @capture(instr, s_(x))
                id = string(s)
                letter = id[1]
                number = parse(Int, id[2:end])
                return :($(Symbol("$(letter)_array"))[i, $(number)])
            elseif @capture(instr, s_[i_])
                id = string(s)
                letter = id[1]
                number = parse(Int, id[2:end])
                return :($(Symbol("$(letter)_array"))[$(i), $(number)])
            elseif @capture(instr, x)
                return :(i * dx)
            end
            return instr
        end
        push!(arr, ne)
    end
    return arr
end
