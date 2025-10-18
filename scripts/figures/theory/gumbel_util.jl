using Distributions
using Infiltrator
using Roots

const C_MIN = 1e-6

struct ReverseGumbel{T} <: ContinuousUnivariateDistribution
    μ::T
    θ::T
end

function Distributions.pdf(d::ReverseGumbel, x::Real)
    z = (-x - d.μ)/d.θ
    return exp.(-z - exp(-z))/d.θ
end

function Distributions.cdf(d::ReverseGumbel, x::Real)
    z = (-x - d.μ)/d.θ
    return 1 - exp(-exp(-z))
end

struct WikipediaAsymptotic end

struct WikipediaExact end

# from:
# https://stats.stackexchange.com/questions/105745/extreme-value-theory-show-normal-to-gumbel
struct QuantileAsymptotic end

# from:
# https://www.panix.com/~kts/Thesis/extreme/extreme2.html#1
struct QuantileExact end


# loc and scale parameters for the Gumbel distribution governing the maximum
# of n iid standard normal variables. 
function gumbel_loc_and_scale(n, ::WikipediaAsymptotic)
    c = sqrt(2*log(n))
    return c, 1/c
end

function gumbel_loc_and_scale(n, ::WikipediaExact)
    c = find_zero(x -> n*exp(-x^2/2)/sqrt(2*pi)/x - 1, (C_MIN, n))
    return c, 1/c
end

function gumbel_loc_and_scale(n, ::QuantileExact)
    b = quantile(Normal(), 1 - 1/n)
    #a = quantile(Normal(), 1 - 1/n/exp(1)) - b
    a = quantile(Normal(), exp(-1/exp(1)/n)) - b
    #a = 1/(n*pdf(Normal(), b))
    return b, a
end

function gumbel_loc_and_scale(n, ::QuantileAsymptotic)
    a = sqrt(2*log(n)*(1 - log(4*pi * log(n))/2/log(n)))
    b = sqrt(log(n^2/4/pi/log(n)))
    return b, 1/a
end