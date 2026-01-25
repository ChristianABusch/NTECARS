# The only point of these structs is a cleaner interface
function Gaussian(σ)
    PowerVoigt(σ, 0, 1)
end


function Voigt(σ, γ)
    PowerVoigt(σ, γ, 1.0)
end

function PowerVoigt(σ, γ, n)
    # find where profile decreases to less than 1% of the maximum
    ν_max = find_zero(x -> voigt_super(x, σ, γ, n) - 0.0001, (0.0, 1000.0))

    # calculate spectrum in the determined range
    ν     = collect(LinRange(-ν_max, ν_max, 100))
    I     = voigt_super.(ν, σ, γ, n)
    return Spectrum(ν, I .- minimum(I), :wavenumber)
end

##################################################
#               Broanening profiles
##################################################
function voigt(x, σ, γ)
    z = (x + 1im*γ) / (√(2)*σ)
    Faddeeva(z) = erfcx(-1im*z)
    V = 1/(σ*2π) * real(Faddeeva(z))

    z0 = (0.0 + 1im*γ) / (√(2)*σ)
    Faddeeva(z0) = erfcx(-1im*z0)
    V0 = 1/(σ*2π) * real(Faddeeva(z0))
    return V/V0
end

function voigt_super(x, σ, γ, n)
    return voigt(x, σ, γ) ^n
end
