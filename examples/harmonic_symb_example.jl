using ModelingToolkit 
include("../src/symbolics/symbolics.jl")

const gamma = 0.0;
const omega = 1.0;
const gamma3 = 0.0;
const A_forcing = 250;
const lambda_forcing = -40;
const harmonics = 3

@parameters x y t

const Dt = Differential(t);
const Dx = Differential(x);
const Dy = Differential(y);

function mkpdes(params, u)
    (gamma, omega, gamma3, A_forcing, lambda_forcing, x, y, t, Dt, Dx, Dy) = params
    [Dt(Dt(u)) - 0.25*(Dx(Dx(u)) + Dy(Dy(u))) + gamma*Dt(u) + gamma3*Dt(u)*Dt(u)*Dt(u) - A_forcing * exp(lambda_forcing*(x/10.0)^2) * sin(omega * t)]
end

params = (gamma, omega, gamma3, A_forcing, lambda_forcing, x, y, t, Dt, Dx, Dy)

harmonic_balance(HarmonicBalanceProblem((x, y), t, omega, harmonics, 1, mkpdes, params)) |> println