include("../src/symbolics/symbolics.jl")
include("../src/discretizer/fd_discretize.jl")
include("../src/residual/residual.jl")
using Plots

const γ = 0.0;
const ω = 1.0;
const γ₃ = 0.0;
const A_forcing = 250;
const λ_forcing = -40;
const harmonics = 1

const x_left = 0.0
const x_right = 1.0
const y_left = 0.0
const y_right = 1.0
const N = 100
step_x = (x_right-x_left)/N;
step_y = (y_right - y_left)/N;
u0 = ones((N) * (N) * harmonics * 2)
par = [step_x, step_y]
Model.@parameters x y t

const Dt = Differential(t);
const Dx = Differential(x);
const Dy = Differential(y);

function mkpdes(params, u)
    (γ, ω, γ₃, A_forcing, λ_forcing, x, y, t, Dt, Dx, Dy) = params
    [Dt(Dt(u)) - 0.25*(Dx(Dx(u)) + Dy(Dy(u))) + γ*Dt(u) + γ₃*Dt(u)*Dt(u)*Dt(u) - A_forcing * exp(λ_forcing*(x/10.0)^2) * sin(ω * t)]
end

params = (γ, ω, γ₃, A_forcing, λ_forcing, x, y, t, Dt, Dx, Dy)

var_names, var_exprs, equations = HarmonicBalanceProblem((x, y), t, ω, harmonics, 1, mkpdes, params) |> harmonic_balance
equations_transformed = transform_sym(equations)

residual = create_residual_function(N, equations_transformed, harmonics)
nonlinear_function = NonlinearFunction(residual)
prob = NonlinearProblem(nonlinear_function, u0, par)

println("Starting nonlinear solver...")
@time sol = NonlinearSolve.solve(prob, NewtonRaphson(), reltol = 1e-5, abstol = 1e-5, maxiters=100)
println("Nonlinear solver finished!")


A_sol = reshape(sol.u[1:N*N], N, N)
B_sol = reshape(sol.u[N*N+1:2*N*N], N, N)


xgrid = range(x_left, x_right, length=N)
ygrid = range(y_left, y_right, length=N)


tgrid = 0.0:0.1:30.0
total_frames = length(tgrid)


anim = @animate for t in tgrid
    u_t = A_sol .* sin(ω * t) .+ B_sol .* cos(ω * t)
    heatmap(xgrid, ygrid, u_t',
            color=:magma,
            xlabel="x", ylabel="y",
            title="u(x,y,t) at t=$t s",
            clims=(-900, 900),
            aspect_ratio=1)
end

gif(anim, "wave_jac.gif", fps=25)
