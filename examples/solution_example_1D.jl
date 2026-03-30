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
u0 = ones((N) * harmonics * 2)
par = [step_x]
Model.@parameters x t

const Dt = Differential(t);
const Dx = Differential(x);

function mkpdes(params, u)
    (γ, ω, γ₃, A_forcing, λ_forcing, x, t, Dt, Dx) = params
    [Dt(Dt(u)) - 0.25*Dx(Dx(u)) + γ*Dt(u) + γ₃*Dt(u)*Dt(u)*Dt(u) - A_forcing * exp(λ_forcing*(x/10.0)^2) * sin(ω * t)]
end

params = (γ, ω, γ₃, A_forcing, λ_forcing, x, t, Dt, Dx)

var_names, var_exprs, equations = HarmonicBalanceProblem((x,), t, ω, harmonics, 1, mkpdes, params) |> harmonic_balance
equations_transformed = transform_sym_1D(equations)

residual = create_residual_function_1D(N, equations_transformed, harmonics)
nonlinear_function = NonlinearFunction(residual)
prob = NonlinearProblem(nonlinear_function, u0, par)

println("Starting nonlinear solver...")
@time sol = NonlinearSolve.solve(prob, NewtonRaphson(), reltol = 1e-5, abstol = 1e-5, maxiters=100)
println("Nonlinear solver finished!")


A_sol = sol.u[1:N]
B_sol = sol.u[N+1:2*N]


xgrid = range(x_left, x_right, length=N)


tgrid = 0.0:0.1:30.0
total_frames = length(tgrid)

anim = @animate for t in tgrid
    u_t = A_sol .* sin(ω * t) .+ B_sol .* cos(ω * t)
    plot(xgrid, u_t,
         xlabel="x", ylabel="u(x,t)",
         title="u(x,t) at t=$(round(t, digits=2)) s",
         ylims=(-900, 900),
         legend=false)
end
gif(anim, "wave_jac.gif", fps=25)
