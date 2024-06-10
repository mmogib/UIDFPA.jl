function steepest_descent()
    return function findNewDk(Fu0, Fu1, d0, u0, u1, k)
        -Fu1
    end
end



function cruz_raydan()
    return function findNewDk(Fu0, Fu1, d0, u0, u1, k)
        y0 = Fu1 - Fu0
        s0 = u1 - u0
        λ = dot(s0, s0) / dot(s0, y0)
        return -λ * Fu1
    end
end



function ye_zhou(; r=0.1)
    return function findNewDk(Fu0, Fu1, d0, u0, u1, k)
        s0 = u1 - u0
        y0 = Fu1 - Fu0 + r * s0
        λ = dot(s0, s0) / dot(s0, y0)
        return -λ * Fu1
    end
end

function abubaker_mohammad()
    return function findNewDk(Fu0, Fu1, d0, u0, u1, k)
        # println(k)
        k1 = k + 1
        r = 1 / k1^2
        t = 1 / (exp(k1)^k1)
        s0 = u1 - u0
        y0 = Fu1 - Fu0 + r * s0
        sdy = dot(s0, y0)
        ns0 = norm(s0)
        θ1 = (ns0^2) / sdy
        θ2 = ns0 / norm(y0)
        λ = (1 - t) * θ1 + t * θ2
        return -λ * Fu1
    end
end

function compute_search_direction(; r=0.1, ψ=0.2, αmin=1.0e-10, αmax=Inf64)
    return function findNewDk(Fu0, Fu1, d0, u0, u1, k)
        # println("k=$k", " ", norm.([Fu0, Fu1, d0, u0, u1]))

        y0 = Fu1 - Fu0
        s0 = u1 - u0 + r * y0
        # v1 = max(ψ * norm(d0) * norm(y0), dot(d0, y0), norm(Fu0)^2)
        v1 = max(ψ * norm(d0) * norm(y0), norm(Fu0)^2)
        β1 = dot(Fu1, y0) / v1
        α12 = dot(Fu1, d0) / v1
        α11 = min(αmax, max(αmin, dot(s0, y0) / dot(y0, y0)))

        return -α11 * Fu1 + β1 * d0 - α12 * y0
    end
end

export steepest_descent, cruz_raydan, ye_zhou, abubaker_mohammad, compute_search_direction