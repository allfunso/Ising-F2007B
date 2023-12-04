using Random
using Statistics
using Plots


function initialstate(N)
    state = 2 .* rand(Bool, (N, N)) .- 1
    return state
end

function mcmove(config, beta)
    N = size(config, 1)
    for _ in 1:N, _ in 1:N
        a, b = rand(1:N), rand(1:N)
        a, b = Int(a), Int(b)
        s = config[a, b]
        nb = config[mod1(a+1, N), b] + config[a, mod1(b+1, N)] + config[mod1(a-1, N), b] + config[a, mod1(b-1, N)]
        cost = 2s * nb
        if cost < 0 || rand() < exp(-cost * beta)
            s *= -1
        end
        config[a, b] = s
    end
    return config
end

function calcEnergy(config)
    energy = 0
    N = size(config, 1)
    for i in 1:N, j in 1:N
        S = config[i, j]
        nb = config[mod1(i+1, N), j] + config[i, mod1(j+1, N)] + config[mod1(i-1, N), j] + config[i, mod1(j-1, N)]
        energy += -nb * S
    end
    return energy / 4
end

function calcMag(config)
    return sum(config)
end
#N = 8,16,32,64
nt = 88 
N = 16
eqSteps = 1000
mcSteps = 1000

T = range(1.53, stop=3.28, length=nt)
E = zeros(nt)
M = zeros(nt)
C = zeros(nt)
X = zeros(nt)
n1, n2 = 1.0 / (mcSteps * N * N), 1.0 / (mcSteps * mcSteps * N * N)

for tt in 1:nt
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N)
    iT = 1.0 / T[tt]
    iT2 = iT * iT

    for i in 1:eqSteps
        mcmove(config, iT)
    end

    for i in 1:mcSteps
        mcmove(config, iT)
        Ene = calcEnergy(config)
        Mag = calcMag(config)

        E1 += Ene
        M1 += Mag
        M2 += Mag * Mag
        E2 += Ene * Ene
    end

    E[tt] = n1 * E1
    M[tt] = n1 * M1
    C[tt] = (n1 * E2 - n2 * E1 * E1) * iT2
    X[tt] = (n1 * M2 - n2 * M1 * M1) * iT
end

plt_E = scatter(T, E, markersize= 4, color=:IndianRed, xlabel="Temperatura (T)", ylabel="Energía", legend=false,tickfontsize=8,guidefontsize=12)
plt_M = scatter(T, abs.(M), markersize=4, color=:RoyalBlue, xlabel="Temperatura (T)", ylabel="Magnetización", legend=false,tickfontsize=8,guidefontsize=12)
plt_C = scatter(T, C, markersize=4, color=:IndianRed, xlabel="Temperatura (T)", ylabel="Calor específico", legend=false,tickfontsize=8,guidefontsize=12)
plt_X = scatter(T, X, markersize=4, color=:RoyalBlue, xlabel="Temperatura (T)", ylabel="Susceptibilidad", legend=false,tickfontsize=8,guidefontsize=12)
