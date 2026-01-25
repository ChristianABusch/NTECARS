abstract type AbstractSpectralProfile end

struct DeltaProfile <: AbstractSpectralProfile end
struct FlatProfile <: AbstractSpectralProfile end


"""
LaserConfiguration(;wavelength_1::Float64, wavelength_2::Float64, stokes_spectrum::Spectrum)
    LaserConfiguration(;wavelength_1::Float64, wavelength_2::Float64, stokes_spectrum::Spectrum, laser_1_profile = nothing, laser_2_profile = nothing)

Contains information on the lasers including central wavelengths and profiles.

`wavelenghts_1` and `wavelenghts_2` have to be provided in units of nm. The spectral range of the `stokes_spectrum`
together with `wavelenghts_1` and `wavelenghts_2` determines the region in which the CARS spectrum is calculated.

# Example
```Julia
lasers = LaserConfiguration(
    wavelength_1    = 532e-9,
    laser_1_profile = Gaussian(0.2/2.35),
    wavelength_2    = 532e-9,
    laser_2_profile = Gaussian(0.2/2.35),
    stokes_spectrum = Spectrum([602.5e-9, 610e-9], [1.0,1.0], :wavelength)
)
```
"""
struct LaserConfiguration
    ν_1::Float64
    ν_2::Float64
    profile_1::Union{DeltaProfile, Spectrum} 
    profile_2::Union{DeltaProfile, Spectrum} 
    stokes_profile::Spectrum 
    ν_aS_limits::Tuple{Float64, Float64}
    
    function LaserConfiguration(;
        wavelength_1::Float64, 
        wavelength_2::Float64, 
        stokes_range::Tuple{Float64, Float64}, 
        profile_1::Union{DeltaProfile, Spectrum} = DeltaProfile(), 
        profile_2::Union{DeltaProfile, Spectrum} = DeltaProfile(),
        stokes_profile::Union{FlatProfile, Spectrum} = FlatProfile()
    )
        if stokes_profile isa FlatProfile
            stokes_profile = Spectrum([stokes_range...], [1.0, 1.0], :wavelength)
        end
        
        ν_1 = wavelength_to_wavenumber(wavelength_1)
        ν_2 = wavelength_to_wavenumber(wavelength_2)
        ν_S = wavenumbers(stokes_profile)
        ν_aS = ν_1 .+ ν_2 .- ν_S
        new(ν_1, ν_2, profile_1, profile_2, stokes_profile, extrema(ν_aS))
    end
end

