#= PROCESSES
Examples of Levy measures (functions ν_p, ν_n and σ_κ)
Tempered stable, truncated stable and Watanabe class
(compatible with the second versions of the functions above)
=#

using SpecialFunctions, GSL, Base.Threads, PoissonRandom

#= TemperedStable (tempered stable), which includes the models: CGMY and KoBoL

Levy measure:
ν(dx)/dx = c_±|x|^{-1-α_±} e^{-λ_±|x|}

Input:
1. b: drift
2-3. c_p, c_n: nonnegative factors such that c_p + c_n > 0
4-5. α_p, α_n: nonnegative factors in [0,2)
6-7. λ_p, λ_n: positive factors

Output:
Functions (b_ν,sν_p,sν_n,sum_sν_p,sum_sν_n,σ_ν), where
1. b_ν(x) = b-∫_{(-1,1)∖(-x,x)}y ν(dy),
2-5. sν_p, sν_n, sum_sν_p and sum_sν_n are as above and
6. σ_ν(x) = √(∫_{(-x,x)}y^2 ν(dy)).
=#
function TemperedStable_init(b::Number,c_p::Number,c_n::Number,α_p::Number,α_n::Number,λ_p::Number,λ_n::Number)
    b_ν(x) = b - ( c_p*λ_p^(α_p-1)*(sf_gamma_inc(1-α_p, x*λ_p)-sf_gamma_inc(1-α_p, λ_p)) -
        c_n*λ_n^(α_n-1)*(sf_gamma_inc(1-α_n, x*λ_n)-sf_gamma_inc(1-α_n, λ_n)) )
    σ_ν(x) = sqrt(c_p*λ_p^(α_p-2)*gamma(2-α_p) * gamma_inc(2-α_p, x*λ_p, 0)[1] +
        c_n*λ_n^(α_n-2)*gamma(2-α_n) * gamma_inc(2-α_n, x*λ_n, 0)[1])
    (sν_p,sum_sν_p) = tempered_power_sampler(α_p,λ_p,c_p)
    (sν_n,sum_sν_n) = tempered_power_sampler(α_n,λ_n,c_n)

    return (b_ν,sν_p,sν_n,sum_sν_p,sum_sν_n,σ_ν)
end

# (Rosinski's series representation)
# Samples the jumps of the measure y ↦ c y^{-α-1}e^{-λy}1_{y>x}dy
function tempered_power_sampler(α::Number,λ::Number,c::Number)
    function s(x::Number,t::Number)
        n = pois_rand(c*t*x^(-α)/α)
        X = min.( rand(n) .^ (-1/α) .* x , rand(Exponential(1/λ),n) .* rand(n) .^ (1/α))
        return X[X .> x]
    end
    function sum_s(x::Number,t::Number)
        n = pois_rand(c*t*x^(-α)/α)
        ap = 1/α
        an = -ap
        dist = Exponential(1/λ)

        if n > 64
            X = zeros(nthreads())
            Y = zeros(nthreads())
            @threads for i = 1:n
                Y[threadid()] = min(rand() ^ an * x , rand(dist) * rand() ^ ap)
                if Y[threadid()] > x
                    X[threadid()] += Y[threadid()]
                end
            end
            return sum(X)
        else
            X = min.(rand(n) .^ an .* x , rand(dist,n) .* rand(n) .^ ap)
            return sum(X[X .> x])
        end
    end
    function sum_s(x::Number,y::Number,t::Number)
        n = pois_rand(c*t*y^(-α)/α)
        ap = 1/α
        an = -ap
        dist = Exponential(1/λ)

        if n > 64
            X1 = zeros(nthreads())
            X2 = zeros(nthreads())
            Y = zeros(nthreads())
            @threads for i = 1:n
                Y[threadid()] = min(rand() ^ an * y , rand(dist) * rand() ^ ap)
                if Y[threadid()] > x
                    X1[threadid()] += Y[threadid()]
                    X2[threadid()] += Y[threadid()]
                elseif Y[threadid()] > y
                    X2[threadid()] += Y[threadid()]
                end
            end
            return (sum(X1),sum(X2))
        else
            X = min.(rand(n) .^ an .* y , rand(dist,n) .* rand(n) .^ ap)
            return (sum(X[X .> x]), sum(X[X .> y]))
        end
    end
    return (s,sum_s)
end

#= Stable process

Levy measure:
ν(dx)/dx = c_±|x|^{-1-α_±}

Input:
1. b: drift
2-3. c_p, c_n: nonnegative factors such that c_p + c_n > 0
4-5. α_p, α_n: nonnegative factors in [0,2)
6-7. λ_p, λ_n: positive factors

Output:
Functions (b_ν,sν_p,sν_n,sum_sν_p,sum_sν_n,σ_ν), where
1. b_ν(x) = b-∫_{(-1,1)∖(-x,x)}y ν(dy),
2-5. sν_p, sν_n, sum_sν_p and sum_sν_n are as above and
6. σ_ν(x) = √(∫_{(-x,x)}y^2 ν(dy)).
=#
function Stable_init(b::Number,c_p::Number,c_n::Number,α_p::Number,α_n::Number)
    b_ν(x) = b + c_p * (α_p != 1 ? (x^(1 - α_p) - 1) / (1-α_p) : log(x)) + c_p * (α_n != 1 ? (x^(1 - α_n) - 1) / (1-α_n) : log(x))
    σ_ν(x) = sqrt(c_p * x^(2 - α_p) / (2-α_p) + c_n * x^(2 - α_n) / (2-α_n))
    (sν_p,sum_sν_p) = power_sampler(α_p,c_p)
    (sν_n,sum_sν_n) = power_sampler(α_n,c_n)

    return (b_ν,sν_p,sν_n,sum_sν_p,sum_sν_n,σ_ν)
end

# Inversion of tail Levy measure
# Samples the jumps of the measure y↦cy^{-α-1}dy
function power_sampler(α::Number,c::Number)
    function s(x::Number,t::Number)
        return rand(pois_rand(c*t*x^(-α)/α)) .^ (-1/α) .* x
    end
    function sum_s(x::Number,t::Number)
        n = pois_rand(c*t*x^(-α)/α)
        a = -1/α
        if n > 32
            X = zeros(nthreads())
            @threads for i = 1:n
                X[threadid()] += rand() ^ a
            end
            return sum(X) * x
        else
            return sum(rand(n) .^ a) * x
        end
    end
    function sum_s(x::Number,y::Number,t::Number)
        n = pois_rand(c*t*y^(-α)/α)
        sc = x/y
        a = -1/α

        if n > 32
            X1 = zeros(nthreads())
            X2 = zeros(nthreads())
            Y = zeros(nthreads())
            @threads for i = 1:n
                Y[threadid()] = rand() ^ a
                X2[threadid()] += Y[threadid()]
                if Y[threadid()] > sc
                    X1[threadid()] += Y[threadid()]
                end
            end
            return (sum(X1) * x, sum(X2) * x)
        else
            X = rand(n) .^ a
            return (sum(X[X .> sc]) * x, sum(X) * x)
        end
    end
    return (s,sum_s)
end

#= (subclass of) Watanabe (continuous singular)

General Levy measure: for some m∈ℕ, m≥2,
ν = ∑_{n∈ℤ}c_{n,+}δ_{m^n}+c_{n,-}δ_{-m^n}
With ∑_{n∈ℕ}(c_{-n,+}+c_{-n,-})b^(-n)+∑_{n∈ℕ}c_{n,+}+c_{n,-}<∞
and max_{n∈ℤ}(c_+ + c_-)<∞.

Our class:
ν = ∑_{n∈ℕ}c_+ δ_{m^(-n)} + c_- δ_{-m^(-n)}.

Note:
Infinite activity iff min(r_+,r_-)≤1
Blumenthal-Getoor index = max(0,-log(max(r_+,r_-))/log(b)))

Input:
1. b: drift
2. m: integer ≥2
3-4. c_p, c_n: nonnegative numbers

Output:
Functions (b_ν,sν_p,sν_n,sum_sν_p,sum_sν_n,σ_ν), where
1. b_ν(x) = b-∫_{(-1,1)∖(-x,x)}y ν(dy),
2-5. sν_p, sν_n, sum_sν_p and sum_sν_n are as above and
6. σ_ν(x) = √(∫_{(-x,x)}y^2 ν(dy)).
=#
function Watanabe_init(b::Number,m::Integer,c_p::Number,c_n::Number)
    function b_ν(x::Number)
        return b - (c_p - c_n) * m^(-( x >= 1 ? 0. : floor(-log(x)/log(m))+1 )) / (1-1/m)
    end
    function σ_ν(x::Number) 
        return sqrt((c_p + c_n) * m^(-2*( x >= 1 ? 0. : floor(-log(x)/log(m))+1 )) / (1-1/m^2))
    end
    (sν_p,sum_sν_p) = wat_sampler(c_p,m)
    (sν_n,sum_sν_n) = wat_sampler(c_n,m)

    return (b_ν,sν_p,sν_n,sum_sν_p,sum_sν_n,σ_ν)
end

# Computes the sum c * |{n∈ℕ:m^n≤1/x}|
function wat(c::Number,m::Integer)
    local f(x) = x >= 1 ? 0. :
        c * ceil(-log(x)/log(m))
    return f
end

# Samples from the density p(k) = 1 / |{n∈ℕ:m^n≤1/x}|
function wat_sampler(c::Number,m::Integer)
    function s(x::Number,t::Number)
        if x >= 1 
            return  Float64[] 
        else
            return [1/m^ceil(rand() * floor(-log(x)/log(m))) for i = 1:Int(pois_rand(c * t * ceil(-log(x)/log(m))))]
        end
    end
    function sum_s(x::Number,t::Number) 
        return sum(s(x,t))
    end
    function sum_s(x::Number,y::Number,t::Number)
        aux = s(y,t)
        return (sum(aux .* (aux .> x)), sum(aux))
    end
    return (s,sum_s)
end
