function create_residual_function(N, equationsExpr, harmonics)
    equationsExprMapped::Vector{Expr} = []

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

        end
    end
    res = eval(function_code)
    return res
end
