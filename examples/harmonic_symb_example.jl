using ModelingToolkit 
include("../src/symbolics/symbolics.jl")

const γ = 0.0;
const ω = 1.0;
const γ₃ = 0.0;
const A_forcing = 250;
const λ_forcing = -40;
const harmonics = 3

@parameters x y t

const Dt = Differential(t);
const Dx = Differential(x);
const Dy = Differential(y);

function mkpdes(params, u)
    (γ, ω, γ₃, A_forcing, λ_forcing, x, y, t, Dt, Dx, Dy) = params
    [Dt(Dt(u)) - 0.25*(Dx(Dx(u)) + Dy(Dy(u))) + γ*Dt(u) + γ₃*Dt(u)*Dt(u)*Dt(u) - A_forcing * exp(λ_forcing*(x/10.0)^2) * sin(ω * t)]
end

params = (γ, ω, γ₃, A_forcing, λ_forcing, x, y, t, Dt, Dx, Dy)

HarmonicBalanceProblem((x, y), t, ω, harmonics, 1, mkpdes, params) |> harmonic_balance |> println
