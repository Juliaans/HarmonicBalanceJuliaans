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
