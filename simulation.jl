using LinearAlgebra
using Test
using Plots

# Optional; just to make plots nicer to look at.
theme(:dark)

function length_from_arc(V)
    ℓ = []

    for j ∈ 1:(length(V)-1)
        push!(ℓ, norm(V[j+1] - V[j]))
    end

    return ℓ
end

function angle_from_arc(V)
    # Find angle between v and x axis
    function angle(v)
        (x, y) = v
        mod(atan(y, x), 2π)
    end

    Θ = []

    for j ∈ 1:(length(V)-1)
        l = V[j+1] - V[j]
        push!(Θ, angle(l))
    end

    return Θ
end

function arc_from_angle(Θ, ℓ)
    V = [[0.0, 0.0]]
    for j ∈ eachindex(Θ)
        push!(V, V[j] + ℓ[j] * [cos(Θ[j]), sin(Θ[j])])
    end
    return V
end
@test arc_from_angle([1.5707963267948966, 3.9269908169872414], [1.0, 1.4142135623730951]) ≈ [[0.0, 0.0], [0.0, 1.0], [-1.0, 0.0]] atol = 1e-8

function grad_E(V)
    """
    Parameters:
        V (array): Vertices of an arc-and-cycle set
        cycle (bool): true/false value stating whether the vertices are part of a cycle/arc respectively.

    Returns:
        ∇E (vector): The gradient of E.
    """
    # Kronecker delta
    δ(x, y) = ==(x, y)
    n = length(V)
    Θ = angle_from_arc(V)
    ℓ = length_from_arc(V)

    # Calculate the E derivatives for the chain rule
    function ∂E()

        ∂xE = []
        ∂yE = []

        # Calculate ∂E/∂xᵣ and ∂E/∂yᵣ
        for idxᵣ ∈ 1:n
            ∂Exᵣ = 0.
            ∂Eyᵣ = 0.

            # Iterate over all edges
            for edge_idx ∈ 1:n-1
                idxᵥ = edge_idx
                idxw = edge_idx + 1

                v = V[idxᵥ]
                w = V[idxw]

                for (idxᵤ, u) ∈ pairs(V)

                    # Only consider vertices u that are not v or w.
                    if idxᵤ ∈ [idxᵥ, idxw]
                        continue
                    end

                    if idxᵣ ∈ [idxᵤ, idxᵥ, idxw]
                        ∂Exᵣ += -2((((u[1] - v[1]) * (δ(idxᵤ, idxᵣ) - δ(idxᵥ, idxᵣ))) / (norm(u - v))
                                    +
                                    ((u[1] - w[1]) * (δ(idxᵤ, idxᵣ) - δ(idxw, idxᵣ))) / (norm(u - w)))
                                   /
                                   (norm(u - v) + norm(u - w) - norm(v - w))^3)

                        ∂Eyᵣ += -2((((u[2] - v[2]) * (δ(idxᵤ, idxᵣ) - δ(idxᵥ, idxᵣ))) / (norm(u - v))
                                    +
                                    ((u[2] - w[2]) * (δ(idxᵤ, idxᵣ) - δ(idxw, idxᵣ))) / (norm(u - w)))
                                   /
                                   (norm(u - v) + norm(u - w) - norm(v - w))^3)
                    end
                end
            end

            # println(∂Exᵣ, ∂Eyᵣ)
            push!(∂xE, ∂Exᵣ)
            push!(∂yE, ∂Eyᵣ)
        end

        return (∂xE, ∂yE)
    end

    # Calculate the xy derivatives for the chain rule
    function ∂xy()
        # The value of ∂X at index (j,k) is ∂xⱼ/∂θₖ. Similar for ∂Y.
        ∂X = zeros(n, n - 1)
        ∂Y = zeros(n, n - 1)

        # Iterate over xⱼ
        for j ∈ 1:n-1
            # Iterate over θₖ
            for k ∈ 1:n-1
                ∂X[j+1, k] = ∂X[j, k] - ℓ[j] * δ(j, k) * sin(Θ[j])
                ∂Y[j+1, k] = ∂Y[j, k] + ℓ[j] * δ(j, k) * cos(Θ[j])
            end
        end
        return (∂X, ∂Y)
    end

    ∂E = ∂E()
    ∂xy = ∂xy()

    ∇E = zeros(n - 1)

    for i ∈ 1:n-1
        for j ∈ 1:n
            ∇E[i] += ∂E[1][j] * ∂xy[1][j, i] + ∂E[2][j] * ∂xy[2][j, i]
        end
    end

    ∇E
end

function timestep!(Θ, V, Δt)
    ∇E = grad_E(V)
    Θ = Θ - Δt * normalize(∇E)

    return Θ
end

function timestep_cycle!(Θ, V, Δt)

    # function absolute_turn_angle(Θ)
    #     n = length(Θ)

    #     turn_angle = zeros(n)

    #     for j ∈ 2:n
    #         turn_angle[j] = mod(π - (Θ[j] + Θ[j-1]), 2π)
    #     end
    #     turn_angle[1] = mod(π - (Θ[n] + Θ[1]), 2π)

    #     return turn_angle
    # end

    # If first and last entries are the same, remove the last entry
    if first(V) == last(V)
        pop!(V)
    end

    n = length(V)
    ℓₙ = norm(V[end] - V[begin])
    i = findall(x -> x == maximum(Θ), Θ)[1]

    # A = absolute_turn_angle(Θ)

    # # This is the index of the maximum abs. turn angle
    # i = findall(x -> x == maximum(A), A)[1]

    # As described in 5.2, let v_n be the vertex of max. abs. turn angle, and recalculate Θ
    circshift!(V, length(Θ) - i)
    ℓ = length_from_arc(V)
    Θ = angle_from_arc(V)

    # Only compute gradient up to θ_n-2
    Θ = Θ[1:end-1]
    V = V[1:end-1]
    ∇E = grad_E(V)
    Θ = Θ - Δt * normalize(∇E)

    V = arc_from_angle(Θ, ℓ[1:end-1])

    # Find value of Θ_n-1/V_n
    d = V[n-1] - V[1]
    d⁺ = [-d[2], d[1]] # This is d⟂ but the symbol doesn't work lol

    # vₙ has two solutions
    vₙ = (V[1] + d * ((ℓₙ^2 - ℓ[n-1]^2 + norm(d)^2) / (2 * norm(d)^2))
     + d⁺ * sqrt((ℓₙ / norm(d))^2 - ((ℓₙ^2 - ℓ[n-1]^2 + norm(d)^2)^2) / (4 * norm(d)^4)),
        V[1] + d * ((ℓ[n-1]^2 - ℓ[n-2]^2 + norm(d)^2) / (2 * norm(d)^2))
        -
        d⁺ * sqrt((ℓₙ / norm(d))^2 - ((ℓₙ^2 - ℓ[n-1]^2 + norm(d)^2)^2) / (4 * norm(d)^4)))

    push!(V, vₙ[1])
    return V
end

function unfold_anim!(V, Δt, t, snapshot_delay=4, cycle=false)
    """
    Generate an unfolding animation for a given configuration. Exports gif to folder where program was run.

        Parameters:
            V (array): Vertex configuration.
            Δt (float): Timestep. 
            t (float): Total time to run unfolding. 
            snapshot_delay (int): Delay between frames to snapshot the image.
        
        Returns:
            None
    """
    global Θ = angle_from_arc(V)
    ℓ = length_from_arc(V)
    j = 0

    anim = @animate while j * Δt ≤ t
        # Update Θ and V
        if cycle == false
            Θ = timestep!(Θ, V, Δt)
            global V = arc_from_angle(Θ, ℓ)
        else
            global V = timestep_cycle!(Θ, V, Δt)
            push!(V, V[1])
        end

        curr_t = round(Δt * j, digits=2)
        j += 1

        plot(first.(V), last.(V),
            seriestype=:scatter,
            primary=false,
            axis=([], false),
            # xlimits=(-20, 20),
            # ylimits=(-20, 20)
        )
        plot!(first.(V), last.(V), label="Time = $curr_t")
        
        if cycle 
            pop!(V)
        end
    end every snapshot_delay

    gif(anim, "anim.gif", fps=10)
end