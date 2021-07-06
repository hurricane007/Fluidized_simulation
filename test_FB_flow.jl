using  OrdinaryDiffEq, LinearAlgebra, Sundials, DiffEqOperators
include("opr_cond_FB.jl")

"""
Simulate Fluidized bed reactor with axial disperCarbion model.
"""
function dFdt!(du, u, p, t)
    fbubA = sum(u[:, 1:nsp], dims=2)*Rg*Temp./pL
    feueA = sum(u[:, 1+nsp:2nsp], dims=2)*Rg*Temp./pL
    #=
    if any(x -> x<0, u)
        locr = findfirst(x->x<0, u)
        println("\n t = $t \n at $locr u = $(u[locr])")
        readline()
    end
    =#
    for i = 1:nsp
        du[:,i] = fbubA.*((Δ1*bc[i]*u[:, i]) + Deg*Δ2*bc[i]*(u[:, i]./fbubA))
        du[:,i+nsp] = feueA.*((Δ1*bc[i+nsp]*u[:,i+nsp]) + Deg*Δ2*bc[i+nsp]*
            (u[:,i+nsp]./feueA))
    end
    du[:, 2nsp+1] =  Des*Δ2*bc[2nsp+1]*u[:, 2nsp+1]
    du[:, 2nsp+2] =  Des*Δ2*bc[2nsp+2]*u[:, 2nsp+2]
end

nsp = length(yA_in)
# descretization
dzi = 500E-6 # smallest knot size
    expFactor = 1.0 # expansion factor
    z = [0.]
    dz = dzi
    while dz + z[end] < Lbed
        push!(z, dz+z[end])
        dz *= expFactor
    end
    if z[end] <= Lbed
        push!(z, Lbed)
    end
    u0 = zeros(length(z), 2length(yA_in)+2)
    for i = 1:2length(yA_in)
        u0[:, i] .= FAin[i]*1E-6
    end
    u0[:, length(yA_in)] .= FAin[length(yA_in)]
    u0[:, 2length(yA_in)] .= FAin[2length(yA_in)]
    u0[:, 2length(yA_in)+1] .= 10.
    u0[:, 2length(yA_in)+2] .= 1-1E-4
    # pressure along z direction
    pL = pres .- fes0*g*ρs*z
    zItv = [z; (2z[end]-z[end-1])] - [-dzi; z]
    Δ1 = UpwindDifference(1, 1, zItv, length(z), -1)
    Δ2 = CenteredDifference(2, 2, zItv, length(z))
    bc = [
    [RobinBC((1.0, -1E-4/fbub0, FAin[i]),(0., 1., 0.), zItv, 1) for i = 1:4] ;
    [RobinBC((1.0, -1E-4/umf0, FAin[i+4]),(0., 1., 0.), zItv, 1) for i = 1:4] ;
        [RobinBC((0., 1., 0.), (0., 1., 0.), zItv, 1) for i = 1:2]]


tspan = (0., 1.)
#para= (Δ1up, Δ2ce, bcAll, length(yA_in))
pdeProb = ODEProblem(dFdt!, u0, tspan)

sol = solve(pdeProb, CVODE_BDF())

using Plots
    plotly()
solu = Array(sol)

surface(sol.t, z, solu[:, 1, :]./FAin[1], xlabel="t [s]", ylabel="z [m]")
