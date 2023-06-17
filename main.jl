include("simulation.jl")

V = [[0,0], [0,1], [-1, 1], [-1,-1], [2,-1]]
Δt = 0.1
t = 8

unfold_anim(V, Δt, t)
