@doc raw"""
    blackbody(wave::Vector{<:Unitful.Quantity}, T::Unitful.Quantity)
    blackbody(wave::Vector{<:Real}, T::Real)

Create a blackbody spectrum using Planck's law. The curve follows the mathematical form:

```math
    B_\lambda(T) = \frac{2hc^2}{\lambda^5} \frac{1}{e^{hc/\lambda k_B T} - 1}\, .
```

If `wave` and `T` are not `Unitful.Quantity`, they are assumed to be in angstrom and Kelvin, and the returned flux will be in units `W m^-2 Å^-1`.

The physical constants are the CODATA2018 / SI-2019 exact values.

When `wave` and `T` are `Unitful` or `DynamicQuantities` quantities, the corresponding
units extension handles the calculation and the flux carries the appropriate units.

# References

[Planck's Law](https://en.wikipedia.org/wiki/Planck%27s_law)

# Examples

```jldoctest
julia> using SpectrumBase, Unitful, UnitfulAstro

julia> wave = range(1, 3, length=100)u"μm"
(1.0:0.020202020202020204:3.0) μm

julia> bb = blackbody(wave, 2000u"K")
SingleSpectrum(Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(μm,), 𝐋, nothing}}, Unitful.Quantity{Float64, 𝐌 𝐋^-1 𝐓^-3, Unitful.FreeUnits{(μm^-1, m^-2, W), 𝐌 𝐋^-1 𝐓^-3, nothing}})
  spectral axis (100,): 1.0 μm .. 3.0 μm
  flux axis (100,): 89534.30930426194 W μm^-1 m^-2 .. 49010.54557924032 W μm^-1 m^-2
  meta: Dict{Symbol, Any}(:T => 2000 K, :name => "Blackbody")

julia> blackbody(ustrip.(u"angstrom", wave), 6000)
SingleSpectrum(Float64, Float64)
  spectral axis (100,): 10000.0 .. 30000.0
  flux axis (100,): 1190.9562575755401 .. 40.043256909104144
  meta: Dict{Symbol, Any}(:T => 6000, :name => "Blackbody")

julia> spectral_axis(bb)[argmax(bb)]
1.4444444444444444 μm

julia> 2898u"μm*K" / bb.T # See if it matches up with Wien's law
1.449 μm
```
"""
function blackbody(wave::AbstractVector{<:Real}, T::Real)
    flux = _planck_per_angstrom.(wave, T)
    return spectrum(wave, flux, name = "Blackbody", T = T)
end

# Planck's law Bλ(T) for wavelength `λ` in angstrom and temperature `T` in Kelvin,
# returned in W m^-2 Å^-1. Uses the SI constants defined in `SpectrumBase`; the
# unit-aware `blackbody` methods (Unitful/DynamicQuantities `Quantity` arguments)
# live in the units extensions.
function _planck_per_angstrom(λ, T)
    λ_m = λ * 1e-10  # angstrom -> m
    spectral_radiance = 2 * PLANCK_CONSTANT * SPEED_OF_LIGHT^2 / λ_m^5 /
        expm1(PLANCK_CONSTANT * SPEED_OF_LIGHT / (λ_m * BOLTZMANN_CONSTANT * T))
    return spectral_radiance * 1e-10  # W m^-2 m^-1 -> W m^-2 Å^-1
end

#"""
#    equivalent_width(::AbstractSpectrum)
#
#Calculate the equivalent width of the given continuum-normalized spectrum. Return value has units equal to wavelengths.
#"""
#function equivalent_width(spec::AbstractSpectrum)
#    dx = spectral_axis(spec)[end] - spectral_axis(spec)[1]
#    flux = ustrip(line_flux(spec))
#    return dx - flux * unit(dx)
#end
#
#"""
#    line_flux(::AbstractSpectrum)
#
#Calculate the line flux of the given continuum-normalized spectrum. Return value has units equal to flux.
#"""
#function line_flux(spec::AbstractSpectrum)
#    avg_dx = diff(spectral_axis(spec))
#    return sum(flux_axis(spec)[2:end] .* avg_dx)
#end
