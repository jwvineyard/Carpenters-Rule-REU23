using LinearAlgebra
using Test

function length_list(V)
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

    @test V = [[0.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]
    Θ = angle_from_arc(V)
    ℓ = length_list(V)

end
@test arc_from_angle([1.5707963267948966, 3.9269908169872414], [1.0, 1.4142135623730951]) ≈ [[0.0, 0.0], [0.0, 1.0], [-1.0, 0.0]] atol = 1e-8


function ∇E(V, cycle)
    # Kronecker delta
    δ(x, y) = ==(x, y)
    n = length(V)
    Θ = angle_from_arc(V)
    ℓ = length_list(V)

    function ∂E(V, cycle)

        # Code for arc
        if cycle == false

            # Check all vertices in V
            for (idxᵣ, r) ∈ pairs(V)
                ∂xE = 0 
                ∂yE = 0

                # Iterate over all edges
                for edge_idx ∈ 1:(n - 1)
                    idxᵥ = V[edge_idx]
                    idxw = V[edge_idx + 1]

                    v = V[idxᵥ]
                    w = V[idxw]

                    # Iterate over all vertices not equal to edge
                    for (idxᵤ, u) ∈ pairs(V)
                        if idxᵤ == idxᵥ | idxᵤ == idxw
                            continue
                        end

                        if idxᵣ == idxᵤ | idxᵣ == idxᵥ | idxᵣ == idxw
                            ∂xE += (-2((u[0] - v[0]) * (δ(idxᵤ, idxᵣ) - δ(idxᵥ, idxᵣ)) / (norm(u - v)) + (u[0] - w[0]) * (δ(idxᵤ, idxᵣ) - δ(idxw, idxᵣ)) / (norm(u - w)))) / (norm(u-v) + norm(u-w) + norm(v-w))^3

                            ∂yE += (-2((u[1] - v[1]) * (δ(idxᵤ, idxᵣ) - δ(idxᵥ, idxᵣ)) / (norm(u - v)) + (u[1] - w[1]) * (δ(idxᵤ, idxᵣ) - δ(idxw, idxᵣ)) / (norm(u - w)))) / (norm(u - v) + norm(u - w) + norm(v - w))^3
                        end
                    end 
                end
            end

            return (∂xE, ∂yE)

        # For cycles (polygons) we may have to modify the code
        else
            return 0
        end
    end

    function ∂xy(V, cycle)
        # Code for cycle
        if cycle == true
            # The value of ∂X at index (j,k) is ∂xⱼ/∂θₖ. Similar for ∂Y.
            ∂X = zeros(n, n - 1)
            ∂Y = zeros(n, n - 1)

            # Iterate over xⱼ
            for j ∈ 2:n 
                # Iterate over θₖ
                for k ∈ 1:n-1
                    ∂X[j,k] = ∂X[j - 1, k] - ℓ[j - 1] * δ(j,k) * sin(θ[j - 1])
                    ∂Y[j,k] = ∂Y[j - 1, k] + ℓ[j - 1] * δ(j,k) * cos(θ[j - 1])
                end
            end
            return (∂X, ∂Y)

        # For cycles (polygons) we may have to modify the code
        else
            return 0
        end
    end
end