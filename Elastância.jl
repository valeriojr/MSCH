using Plots

T = 1.0
k_sys = 0.1
t_sh = 0

function eᵥ(t, T_sys0 = 0.54, k_r = 0.83)
    
    T_sys = T_sys0 - (k_sys / T)
    T_r = k_r * T_sys
    
    t %= T
    if ((t_sh + T_sys) <= t < T) || (0.0 <= t < t_sh)    
        return 0.0
    elseif t_sh <= t < (t_sh + T_r)
        num = π * (t - t_sh)
        den = 2T_r
        return sin(num / den)^2
    elseif (t_sh + T_r) <=  t < (t_sh + T_sys)
        num = π * (t - t_sh - T_r)
        den = 2.0 * (T_sys - T_r)
        return cos(num / den)^2
    end
end;

t = 0:0.0001:2T

lim = map(eᵥ, t)
plot(t, lim, label="Lim", size = (600, 200), color="red")

Eₘₐₓ = 2.0
Eₘᵢₙ = 0.06
HR = 60
tc = 60 / HR
Tₘₐₓ = 0.2 + 0.15tc

function Eₙ(tₙ)
    n₁ = (tₙ / 0.7) ^ 1.9
    n₂ = (tₙ / 1.17) ^ 21.9
    return 1.55 * (n₁ / (1 + n₁)) * (1 / (1 + n₂))
end;

function E(t) 
    return (Eₘₐₓ - Eₘᵢₙ) * Eₙ(t/Tₘₐₓ) + Eₘᵢₙ
end;

simaan = map(Eₙ, t .% T ./ Tₘₐₓ)
plot(t, simaan, label="Simaan", size = (600, 200), color="blue")

function compute_error(a, b)
    return sum((eᵥ.(t, a, b) .- simaan) .^ 2)
end

x = 0.0:0.01:1.0
p = Iterators.product(x, x)

grid_values = Array{Float64}(undef, length(p))

i = 1
for t in p
    grid_values[i] = compute_error(t[1], t[2])
    i += 1
end

i = 1
min_i = argmin(grid_values)
min_t = nothing
for t in p
    if i == min_i
        min_t = t
    end
    i += 1
end

println(min_t)

n = Int64(√length(p))
heat = reshape(grid_values, (n, n))
heatmap(x, x, heat, c=cgrad(:viridis, rev=true), title="Erro entre as curvas", xlabel="\$T_{sys0}\$", ylabel="\$k_r\$")
scatter!([min_t[2]],[min_t[1]], legend=false)

plot(t, simaan, label = "E(t) Simaan", color="blue")
lim = eᵥ.(t, min_t[1], min_t[2])
plot!(t, lim, line = :dash, label = "E(t) Lim", size = (600, 200), color="red")
