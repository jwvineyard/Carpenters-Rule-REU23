using LinearAlgebra
using Test
using Plots

# Optional; just to make plots nicer to look at.
theme(:dark)

function length_from_arc(V)
    ℓ = []

    for j ∈ 1:(length(V) - 1)
        push!(ℓ, norm(V[j + 1] - V[j]))
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

    for j ∈ 1:(length(V) - 1)
        l = V[j + 1] - V[j]
        push!(Θ, angle(l))
    end

    return Θ
end 

function arc_from_angle(Θ, ℓ)
    V = [[0.,0.]]
    for j ∈ eachindex(Θ)
        push!(V, V[j] + ℓ[j]*[cos(Θ[j]), sin(Θ[j])])
    end
    return V
end
@test arc_from_angle([1.5707963267948966, 3.9269908169872414], [1.0, 1.4142135623730951]) ≈ [[0.0, 0.0], [0.0, 1.0], [-1.0, 0.0]] atol = 1e-8

function grad_E(V, cycle=false)
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
    function ∂E(cycle)

        # Code for arc
        if cycle == false
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

        # For cycles (polygons) we may have to modify the code
        else
            return 0
        end
    end

    # Calculate the xy derivatives for the chain rule
    function ∂xy(cycle)

        # Code for cycle
        if cycle == false
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

        # For cycles (polygons) we may have to modify the code
        else
            return 0
        end
    end

    ∂E = ∂E(cycle)
    ∂xy = ∂xy(cycle)

    ∇E = zeros(n-1)

    for i ∈ 1:n-1
        for j ∈ 1:n
            ∇E[i] += ∂E[1][j]*∂xy[1][j,i] + ∂E[2][j]*∂xy[2][j,i]
        end
    end
    
    ∇E
end

function timestep(Θ, V, Δt, cycle=false)
    ∇E = grad_E(V, cycle)
    Θ′ = Θ - Δt * normalize(∇E)
    return Θ′
end

function unfold_anim(V, Δt, t, snapshot_delay=4, cycle=false)
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
        Θ = timestep(Θ, V, Δt, cycle)
        global V = arc_from_angle(Θ, ℓ)

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
    end every snapshot_delay

    gif(anim, "anim.gif", fps=10)
end