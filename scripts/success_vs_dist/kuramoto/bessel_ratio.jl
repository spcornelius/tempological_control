using SpecialFunctions
using Roots
using ArbNumerics

I₀(κ) = besseli(0, κ)
I₁(κ) = besseli(1, κ)
I₂(κ) = besseli(2, κ)

function get_κ₀(r::T) where {T <: Real}
    # Approximate the inverse, κ, of r = I₁(κ)/I₀(κ) via Taylor series,
    # where I₁, I₀ are the modified Bessel functions of the first kind of
    # the given order. The approximated κ can be refined via, e.g., Newton
    # iteration.

    # Follows the approach detailed in Hill (1981), "Evaluation and Inversion
    # of the Ratios of Modified Bessel Functions".
    0 < r < 1 || error("Must have 0 < r < 1.")

    if r <= 0.85
        rem = 2r^11/(360*(1-r))*(4.6494 - 5.0797r + 5.6076r^2 - r^3)
        return (2r - r^3 -r^5/6 - r^7/24 + r^9/360 + 53r^11/2160 + rem)/(1-r^2)
    else
        y = 2/(1-r)
        if r >= 0.95
            x = 32/(y - 131.5 + 120r)
        else
            x = 2001.035224 + 4317.5526r - 2326r^2
        end
        return (y + 1 + 3/(y - 5 - 12/(y - 10 - x)))/4
    end
end

function inverse_bessel_ratio(r)
    # Numerically calculate the value of κ that solves the equation
    # r = I₁(κ)/I₀(κ), where I₀, I₁ are the modified Bessel functions
    # of the first kind.

    # initial guess via series approximation of inverse
    κ₀ = get_κ₀(r)

    # polish initial guess via Newton iteration
    f(κ) = I₁(κ)/I₀(κ) - r
    f′(κ) = 1/2 + I₂(κ)/2/I₀(κ) - (I₁(κ)/I₀(κ))^2
    κ = find_zero((f, f′), κ₀, Roots.Newton())
    return κ
end
