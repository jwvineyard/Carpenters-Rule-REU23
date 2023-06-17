function grad_energy(V, cycle)

    # Kronecker delta
    δ(x, y) = ==(x, y)

    # Code for cycle
    if cycle == false

        # Check all vertices in V
        for (idxᵣ, r) ∈ pairs(V)
            ∂xE = 0 
            ∂yE = 0

            # Iterate over all edges
            for edge_idx ∈ 1:(length(V) - 1)
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
    end
end