#= This code is a simple implementation of the paper
"Simulation of the drawdown and its duraton in L\'evy models via stick-breaking Gaussian approximation"
by Gonzalez Cazares and Mijatovic
=#

##############################################################
### Simulations from the Levy measure are assumed to come  ###
### from sampling from the normalised x ↦ ν(dx)/ν([x1,x2)) ###
##############################################################

using Distributions, Statistics, StatsBase, Random, PoissonRandom

# Simulates (X,infX,τX) over [0,t] for X∼BrownianMotion(σ,μ)
function ϕ(t,σ,μ)
    (X,supX,τX) = ϕ_MAXLOCATION(t,σ,-μ)
    return (-X,-supX,τX)
end

function ϕ_MAXLOCATION(t,σ,μ)
    ts = sqrt(t) * σ
    X = rand(Normal()) + μ * sqrt(t) / σ
    supX = (X + sqrt(X^2 - 2 * log(rand()))) / 2
    if supX >= sqrt(2)
        while true
            Z = 1 + (supX - X)^2 / rand(Normal())^2
            if rand() <= Z * exp((1-Z) * supX^2/2)
                τX = 1/Z
                return (X * ts, supX * ts, τX * t)
            end
        end
    elseif supX - X >= sqrt(2)
        while true
            Z = 1 + supX^2 / rand(Normal())^2
            if rand() <= Z * exp((1-Z) * (supX - X)^2/2)
                τX = 1 - 1/Z
                return (X * ts, supX * ts, τX * t)
            end
        end
    else
        while true
            τX = rand(Arcsine())
            if 4 * rand() * τX * (1 - τX) <= supX^2 * (supX - X)^2 * exp(2 - supX^2/(2 * τX) - (supX - X)^2 / (2*(1 - τX)))
                return (X * ts, supX * ts, τX * t)
            end
        end
    end
end

#= Gaussian approximation increment X_t^{(κ)}

Input:
1. b_κ: drift associated to the cutoff function x↦1_{(-κ,κ)}(x),
2. σ: Brownian component σ≥0,
3. sum_sν_p: sampler from the sum of the jumps drawn from the positive tail measure ν(⋅∩[x,∞)), i.e. sum_sν_p(x,t) gives the sum of a Poisson(ν([x,∞))t) number of positive variables with law ν(⋅∩[x,∞))/ν([x,∞)),
4. sum_sν_n: sampler from the sum of the jumps drawn from the negative tail measure ν(⋅∩(-∞,-x]), i.e. sum_sν_n(x,t) gives the sum of a Poisson(ν((-∞,-x])t) number of positive variables with law ν(⋅∩(-∞,-x])/ν((-∞,-x]),
NOTE both functions return positive numbers
5. σ_κ: standard deviation of small jumps √(∫_{(-κ,κ)}x^2 ν(dx)),
6. t: time horizon t > 0,
7. κ: cutoff parameter.

Output:
A real number with the law of X_t^{(κ)}
=#
function rand_G(b_κ::Number,σ::Number,sum_sν_p::Function,sum_sν_n::Function,σ_κ::Number,t::Number,κ::Number)
    if κ<=0 #σ<0 || σ_κ<0 ||
        error("Incorrect domain of parameters: κ>0")
    end

    return b_κ * t + sum_sν_p(κ,t) - sum_sν_n(κ,t) +
        (σ == 0 ? σ_κ : sqrt(σ^2+σ_κ^2)) * rand(Normal()) * sqrt(t)
end

#= Double Gaussian approximation increment Π_t^{(κ1,κ2)}

Input:
1-2. b_κ1, b_κ2: drift associated to the cutoff function x↦1_{(-κ,κ)}(x),
3. σ: Brownian component σ≥0,
3. sum_sν_p: sampler from the sum of the jumps drawn from the positive tail measure ν(⋅∩[y,∞)), i.e. sum_sν_p(x,y,t) gives (a,b), where b is the sum of a Poisson(ν([y,∞))t) number of positive variables with law ν(⋅∩[y,∞))/ν([y,∞)) and a is the variables in b with size at least x,
4. sum_sν_n: sampler from the sum of the jumps drawn from the negative tail measure ν(⋅∩(-∞,-y]), i.e. sum_sν_n(x,y,t) gives (a,b), where b is the sum of a Poisson(ν((-∞,-y])t) number of positive variables with law ν(⋅∩(-∞,-y])/ν((-∞,-y]) and a is the variables in b with size at least x,
NOTE both functions return positive numbers
6-7. σ_κ1, σ_κ2: standard deviation of small jumps √(∫_{(-κ,κ)}x^2 ν(dx)),
8. t: time horizon t > 0,
9-10. κ1, κ2: cutoff parameters.

Output:
A vector with the law of (X_t^{(κ1)},X_t^{(κ2)}) under a certain coupling
=#
function rand_G(b_κ1::Number,b_κ2::Number,σ::Number,sum_sν_p::Function,sum_sν_n::Function,σ_κ1::Number,σ_κ2::Number,t::Number,κ1::Number,κ2::Number)
    if κ1 <= 0 || κ2 <= 0 || κ1 <= κ2#σ<0 || σ_κ1<0 || σ_κ2<0 ||
        error("Incorrect domain of parameters: κ1>κ2>0")
    end

    (J1,J2) = sum_sν_p(κ1,κ2,t) .- sum_sν_n(κ1,κ2,t)

    B = rand(Normal(0,t))
    return (b_κ1 * t + (σ == 0 ? σ_κ1 : sqrt(σ^2+σ_κ1^2)) * B + J1,
        b_κ2 * t + (σ == 0 ? σ_κ2 : sqrt(σ^2+σ_κ2^2)) * B + J2)
end

#= Trivariate Gaussian approximation law (X_t^{(κ)},min X_s^{(κ)},argmin X_s^{(κ)}) (Ver 2)

Input:
1. b_κ: drift associated to the cutoff function x↦1_{(-κ,κ)}(x),
2. σ: Brownian component σ≥0,
3. sν_p: magnitude batch sampler from the positive tail measure ν(⋅∩[x,∞))t, i.e. sν_p(x,t) gives a Poisson(ν([x,∞))t) number of positive variables with law ν(⋅∩[x,∞))/ν([x,∞)),
4. sν_n: magnitude batch sampler from the negative tail measure ν(⋅∩(-∞,-x])t, i.e. sν_n(x,t) gives a Poisson(ν((-∞,-x])t) number of positive variables with law ν((-⋅)∩(-∞,-x])/ν((-∞,-x]),
NOTE both functions return positive numbers
5. σ_κ: standard deviation of small jumps √(∫_{(-κ,κ)}x^2 ν(dx)),
6. t: time horizon t > 0,
7. κ: cutoff parameter.

Output:
A bivariate variable with the law of (X_t^{(κ)},min X_s^{(κ)},argmin X_s^{(κ)})
=#
function rand_Gχ(b_κ::Number,σ::Number,sν_p::Function,sν_n::Function,σ_κ::Number,t::Number,κ::Number)
    if κ<=0 #σ<0 || σ_κ<0 ||
        error("Incorrect domain of parameters: κ>0")
    end

    s = σ == 0 ? σ_κ : sqrt(σ^2 + σ_κ^2)

    J = sν_p(κ,t)
    append!(J, -sν_n(κ,t))
    U = [0.]

    (X,Xmin,τX) = (0.,0.,0.)

    if length(J) > 0
        # Shuffle the jump heights
        J = J[randperm(length(J))]
        append!(U,rand(Exponential(),length(J)))
        U = t .* cumsum(E) ./ (sum(E) - log(rand()))

        for i = 1:length(J)
            (B1,B2,B3) = ϕ(U[i+1]-U[i],s,b_κ)
            if X + B2 < Xmin
                Xmin = X + B2
                τX = U[i] + B3
            end
            X = X + B1 + J[i]
        end
        (B1,B2,B3) = ϕ(t-U[end],s,b_κ)
        if X + B2 < Xmin
            Xmin = X + B2
            τX = U[end] + B3
        end
        X = X + B1
    else
        (X,Xmin,τX) = ϕ(t,s,b_κ)
    end

    return (X,Xmin,τX)
end

#= Double trivariate Gaussian approximation law Π_t^{(κ1,κ2)}

Input:
1-2. b_κ1, b_κ2: drift associated to the cutoff function x↦1_{(-κ,κ)}(x),
3. σ: Brownian component σ≥0,
4. sν_p: magnitude batch sampler from the positive tail measure ν(⋅∩[x,∞)), i.e. sν_p(x) gives a Poisson(ν([x,∞))) number of positive variables with law ν(⋅∩[x,∞))/ν([x,∞)),
5. sν_n: magnitude batch sampler from the negative tail measure ν(⋅∩(-∞,-x]), i.e. sν_n(x) gives a Poisson(ν((-∞,-x])) number of positive variables with law ν(⋅∩(-∞,-x])/ν((-∞,-x]),
NOTE both functions return positive numbers
6-7. σ_κ1, σ_κ2: standard deviation of small jumps √(∫_{(-κ,κ)}x^2 ν(dx)),
8. t: time horizon t > 0,
9-10. κ1, κ2: cutoff parameters.

Output:
A vector with the law of (X_t^{(κ1)},min X_s^{(κ1)},argmin X_s^{(κ1)},X_t^{(κ2)},min X_s^{(κ2)},argmin X_s^{(κ2)}) under a certain coupling
=#
function rand_Gχ(b_κ1::Number,b_κ2::Number,σ::Number,sν_p::Function,sν_n::Function,σ_κ1::Number,σ_κ2::Number,t::Number,κ1::Number,κ2::Number)
    if κ1 <= 0 || κ2 <= 0 || κ1 <= κ2#σ<0 || σ_κ1<0 || σ_κ2<0 ||
        error("Incorrect domain of parameters: κ1>κ2>0")
    end
    #if κ1<κ2
    #    (a1,a2,a3,b1,b2,b3) = rand_Gχ(b_κ2,b_κ2,σ,sν_p,sν_n,σ_κ2,σ_κ1,t,κ2,κ1)
    #    return (b1,b2,b3,a1,a2,a3)
    #end

    s1 = σ == 0 ? σ_κ1 : sqrt(σ^2 + σ_κ1^2)
    s2 = σ == 0 ? σ_κ2 : sqrt(σ^2 + σ_κ2^2)

    J2 = sν_p(κ2,t)
    append!(J2, -sν_n(κ2,t))

    (X1,Xmin1,τX1,X2,Xmin2,τX2) = (0.,0.,0.,0.,0.,0.)

    if length(J2) > 0
        J2 = J2[randperm(length(J2))]
        # Do calculations for X^{κ2}
        U2 = [0.]
        append!(U2,rand(Exponential(),length(J2)))
        U2 = t .* cumsum(U2) ./ (sum(U2) - log(rand()))

        for i = 1:length(J2)
            (B1,B2,B3) = ϕ(U2[i+1]-U2[i],s2,b_κ2)
            if X2 + B2 < Xmin2
                Xmin2 = X2 + B2
                τX2 = U2[i] + B3
            end
            X2 = X2 + B1 + J2[i]
        end
        (B1,B2,B3) = ϕ(t-U2[end],s2,b_κ2)
        if X2 + B2 < Xmin2
            Xmin2 = X2 + B2
            τX2 = U2[end] + B3
        end
        X2 = X2 + B1 + J2[end]

        # Do calculations for X^{κ1}
        ind = abs.(J2) .> κ1

        if sum(ind) > 0
            
            J1 = J2[ind]
            U1 = [0.]
            append!(U1, U2[2:end][ind])
            
            for i = 1:length(J1)
                (B1,B2,B3) = ϕ(U1[i+1]-U1[i],s1,b_κ1)
                if X1 + B2 < Xmin1
                    Xmin1 = X1 + B2
                    τX1 = U1[i] + B3
                end
                X1 = X1 + B1 + J1[i]
            end
            (B1,B2,B3) = ϕ(t-U1[end],s1,b_κ1)
            if X1 + B2 < Xmin1
                Xmin1 = X1 + B2
                τX1 = U1[end] + B3
            end
            X1 = X1 + B1 + J1[end]
        else
            (X1,Xmin1,τX1) = ϕ(t,s1,b_κ1)
        end
    else
        (X2,Xmin2,τX2) = ϕ(t,s2,b_κ2)
        (X1,Xmin1,τX1) = ϕ(t,s1,b_κ1)
    end

    return (X1,Xmin1,τX1,X2,Xmin2,τX2)
end

#= The main algorithm: SBGaussian approximation law Π_{n,t}^{κ}

Input:
1. b_κ: drift associated to the cutoff function x↦1_{(-κ,κ)}(x),
2. σ: Brownian component σ≥0,
3. sν_p: magnitude batch sampler from the positive tail measure ν(⋅∩[x,∞)), i.e. sν_p(x) gives a Poisson(ν([x,∞))) number of positive variables with law ν(⋅∩[x,∞))/ν([x,∞)),
4. sν_n: magnitude batch sampler from the negative tail measure ν(⋅∩(-∞,-x]), i.e. sν_n(x) gives a Poisson(ν((-∞,-x])) number of positive variables with law ν(⋅∩(-∞,-x])/ν((-∞,-x]),
5. sum_sν_p: sampler from the sum of the jumps drawn from the positive tail measure ν(⋅∩[x,∞)), i.e. sum_sν_p(x,t) gives the sum of a Poisson(ν([x,∞))t) number of positive variables with law ν(⋅∩[x,∞))/ν([x,∞)),
6. sum_sν_n: sampler from the sum of the jumps drawn from the negative tail measure ν(⋅∩(-∞,-x]), i.e. sum_sν_n(x,t) gives the sum of a Poisson(ν((-∞,-x])t) number of positive variables with law ν(⋅∩(-∞,-x])/ν((-∞,-x]),
NOTE both functions return positive numbers
7. σ_κ: standard deviation of small jumps √(∫_{(-κ,κ)}x^2 ν(dx)),
8. n: number of steps n ∈ ℕ,
9. t: time horizon t > 0,
10. κ: cutoff parameters.

Output:
A vector with the law of (X_t^{(κ)},min X_s^{(κ)}, argmin X_s^{(κ)})
=#
function rand_SBG(b_κ::Number,σ::Number,sν_p::Function,sν_n::Function,sum_sν_p::Function,sum_sν_n::Function,σ_κ::Number,n::Integer,t::Number,κ::Number)
    if κ<=0 #σ<0 || σ_κ1<0 || σ_κ2<0 ||
        error("Incorrect domain of parameters: κ>0")
    end

    L = Float64[t]
    append!(L,rand(n-1))
    L = cumprod(L)
    ℓ = L[1:(end-1)] .- L[2:end]

    X = 0.
    Xmin = 0.
    τX = 0.

    for i = 1:(n-1)
        ξ = rand_G(b_κ,σ,sum_sν_p,sum_sν_n,σ_κ,ℓ[i],κ)
        X += ξ
        if ξ < 0
            Xmin += ξ
            τX += ℓ[i]
        end
    end

    (ξ,ξmin,τ) = rand_Gχ(b_κ,σ,sν_p,sν_n,σ_κ,L[end],κ)
    X += ξ
    Xmin += ξmin
    τX += τ

    return (X,Xmin,τX)
end

#= The main algorithm: SBGaussian approximation law Π_{n,t}^{κ1,κ2}

Input:
1-2. b_κ1, b_κ2: drift associated to the cutoff function x↦1_{(-κ,κ)}(x),
3. σ: Brownian component σ≥0,
4. sν_p: magnitude batch sampler from the positive tail measure ν(⋅∩[x,∞)), i.e. sν_p(x) gives a Poisson(ν([x,∞))) number of positive variables with law ν(⋅∩[x,∞))/ν([x,∞)),
5. sν_n: magnitude batch sampler from the negative tail measure ν(⋅∩(-∞,-x]), i.e. sν_n(x) gives a Poisson(ν((-∞,-x])) number of positive variables with law ν(⋅∩(-∞,-x])/ν((-∞,-x]),
6. sum_sν_p: sampler from the sum of the jumps drawn from the positive tail measure ν(⋅∩[y,∞)), i.e. sum_sν_p(x,y,t) gives (a,b), where b is the sum of a Poisson(ν([y,∞))t) number of positive variables with law ν(⋅∩[y,∞))/ν([y,∞)) and a is the variables in b with size at least x,
7. sum_sν_n: sampler from the sum of the jumps drawn from the negative tail measure ν(⋅∩(-∞,-y]), i.e. sum_sν_n(x,y,t) gives (a,b), where b is the sum of a Poisson(ν((-∞,-y])t) number of positive variables with law ν(⋅∩(-∞,-y])/ν((-∞,-y]) and a is the variables in b with size at least x,
NOTE both functions return positive numbers
8-9. σ_κ1, σ_κ2: variance of small jumps √(∫_{(-κ_i,κ_i)}x^2 ν(dx)),
10. n: number of steps n ∈ ℕ,
11. t: time horizon t > 0,
12-13. κ1, κ2: cutoff parameters.

Output:
A vector with the law of (X_t^{(κ1)},min X_s^{(κ1)},argmin X_s^{(κ1)},X_t^{(κ2)},min X_s^{(κ2)},argmin X_s^{(κ2)}) under the SBG coupling
=#
function rand_SBG(b_κ1::Number,b_κ2::Number,σ::Number,sν_p::Function,sν_n::Function,sum_sν_p::Function,sum_sν_n::Function,σ_κ1::Number,σ_κ2::Number,n::Integer,t::Number,κ1::Number,κ2::Number)
    if κ1<=0 || κ2<=0 #σ<0 || σ_κ1<0 || σ_κ2<0 ||
        error("Incorrect domain of parameters: κ1,κ2>0")
    end

    L = Float64[t]
    append!(L,rand(n-1))
    L = cumprod(L)
    ℓ = L[1:(end-1)] .- L[2:end]

    X1 = 0.
    Xmin1 = 0.
    τX1 = 0.
    X2 = 0.
    Xmin2 = 0.
    τX2 = 0.

    for i = 1:(n-1)
        (ξ1,ξ2) = rand_G(b_κ1,b_κ2,σ,sum_sν_p,sum_sν_n,σ_κ1,σ_κ2,ℓ[i],κ1,κ2)
        X1 += ξ1
        if ξ1 < 0
            Xmin1 += ξ1
            τX1 += ℓ[i]
        end
        X2 += ξ2
        if ξ2 < 0
            Xmin2 += ξ2
            τX2 += ℓ[i]
        end
    end

    (ξ1,ξmin1,τ1,ξ2,ξmin2,τ2) = rand_Gχ(b_κ1,b_κ2,σ,sν_p,sν_n,σ_κ1,σ_κ2,L[end],κ1,κ2)
    X1 += ξ1
    Xmin1 += ξmin1
    τX1 += τ1
    X2 += ξ2
    Xmin2 += ξmin2
    τX2 += τ2

    return (X1,Xmin1,τX1,X2,Xmin2,τX2)
end
