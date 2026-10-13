module SpectrumBaseUnitfulExt

using SpectrumBase
using SpectrumBase: AbstractSpectrum, Spectrum, BinnedSpectrum, spectrum, spectral_axis, flux_axis, meta, blackbody
import SpectrumBase: _ustrip, _velocity_ratio, _spectral_axis_factor
using Unitful: Unitful, @u_str, ustrip, unit, dimension, Quantity, Length, Frequency, Energy, NoUnits
using UnitfulEquivalences: Spectral, SpectralDensity

# Unit-ful physical constants, rebuilt from the SI floats that live in
# `SpectrumBase` so the Unitful backend agrees with the numeric core to the last
# digit.
const h = SpectrumBase.PLANCK_CONSTANT * u"J*s"
const c_0 = SpectrumBase.SPEED_OF_LIGHT * u"m/s"
const k_B = SpectrumBase.BOLTZMANN_CONSTANT * u"J/K"

# Recognized spectral-axis dimensions: wavelength, frequency, photon energy.
const SPECTRAL_DIMENSIONS = (u"𝐋", u"𝐓^-1", u"𝐋^2 * 𝐌 * 𝐓^-2")

# --- units-agnostic hooks from the core -------------------------------------

_ustrip(x::Quantity) = ustrip(x)

_velocity_ratio(v::Quantity) = ustrip(NoUnits, v / c_0)

_spectral_axis_factor(::Type{<:Length}, factor) = factor          # Wavelength: stretches
_spectral_axis_factor(::Type{<:Frequency}, factor) = inv(factor)  # Frequency: compresses
_spectral_axis_factor(::Type{<:Energy}, factor) = inv(factor)     # Energy: compresses
_spectral_axis_factor(::Type{S}, factor) where {S <: Quantity} =
    throw(ArgumentError("cannot shift a spectral axis with dimension $(dimension(S)); expected a wavelength, frequency, or energy"))

# --- spectrum constructors for Unitful axes ---------------------------------

function SpectrumBase.spectrum(spectral_axis::AbstractVector{<:Quantity}, flux_axis::AbstractVector{<:Quantity}; kwds...)
    @assert dimension(eltype(spectral_axis)) ∈ SPECTRAL_DIMENSIONS "spectral_axis not recognized as having dimensions of wavelength, frequency, or energy."
    Spectrum(spectral_axis, flux_axis, Dict{Symbol, Any}(kwds))
end

function SpectrumBase.spectrum(spectral_axis::AbstractMatrix{<:Quantity}, flux_axis::AbstractMatrix{<:Quantity}; kwds...)
    @assert dimension(eltype(spectral_axis)) ∈ SPECTRAL_DIMENSIONS "spectral_axis not recognized as having dimensions of wavelength, frequency, or energy."
    Spectrum(spectral_axis, flux_axis, Dict{Symbol, Any}(kwds))
end

function SpectrumBase.spectrum(spectral_axis::AbstractMatrix{<:Quantity}, flux_axis::AbstractVector{<:Quantity}; kwds...)
    @assert dimension(eltype(spectral_axis)) ∈ SPECTRAL_DIMENSIONS "spectral_axis not recognized as having dimensions of wavelength, frequency, or energy."
    Spectrum(spectral_axis, flux_axis, Dict{Symbol, Any}(kwds))
end

# --- ustrip / unit on whole spectra -----------------------------------------

"""
    Unitful.ustrip(::AbstractSpectrum)

Remove the units from a spectrum. Useful for processing spectra in tools that don't play nicely with `Unitful.jl`

# Examples
```jldoctest
julia> using Random

julia> rng = Random.seed!(0)
TaskLocalRNG()

julia> using Unitful, UnitfulAstro

julia> wave = range(1e4, 3e4, length=1000);

julia> flux = wave .* 10 .+ randn(rng, 1000);

julia> spec = spectrum(wave*u"angstrom", flux*u"W/m^2/angstrom")
SingleSpectrum(Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}, Unitful.Quantity{Float64, 𝐌 𝐋^-1 𝐓^-3, Unitful.FreeUnits{(Å^-1, m^-2, W), 𝐌 𝐋^-1 𝐓^-3, nothing}})
  spectral axis (1000,): 10000.0 Å .. 30000.0 Å
  flux axis (1000,): 99999.76809093042 W Å^-1 m^-2 .. 300000.2474309158 W Å^-1 m^-2
  meta: Dict{Symbol, Any}()

julia> ustrip(spec)
SingleSpectrum(Float64, Float64)
  spectral axis (1000,): 10000.0 .. 30000.0
  flux axis (1000,): 99999.76809093042 .. 300000.2474309158
  meta: Dict{Symbol, Any}()
```
"""
Unitful.ustrip(spec::AbstractSpectrum) = spectrum(ustrip.(spectral_axis(spec)), ustrip.(flux_axis(spec)); meta(spec)...)

"""
    Unitful.unit(::AbstractSpectrum)

Get the units of a spectrum. Returns a tuple of the spectral axis units and flux/sigma units

# Examples
```jldoctest
julia> using Random

julia> rng = Random.seed!(0)
TaskLocalRNG()

julia> using Unitful, UnitfulAstro

julia> wave = range(1e4, 3e4, length=1000);

julia> flux = wave .* 10 .+ randn(rng, 1000);

julia> spec = spectrum(wave * u"angstrom", flux * u"W/m^2/angstrom");

julia> w_unit, f_unit = unit(spec)
(Å, W Å^-1 m^-2)
```
"""
Unitful.unit(spec::AbstractSpectrum) = unit(eltype(spectral_axis(spec))), unit(eltype(flux_axis(spec)))

# --- blackbody with Unitful arguments ---------------------------------------

function SpectrumBase.blackbody(wave::AbstractVector{<:Quantity}, T::Quantity)
    out_unit = u"W/m^2" / unit(eltype(wave))
    flux = _blackbody(wave, T) .|> out_unit
    return spectrum(wave, flux, name = "Blackbody", T = T)
end

"""
    blackbody(T::Unitful.Quantity)

Returns a function for calculating blackbody curves.
"""
SpectrumBase.blackbody(T::Quantity) = w -> 2h * c_0^2 / w^5 / (exp(h * c_0 / (w * k_B * T)) - 1)

_blackbody(wave::AbstractVector{<:Quantity}, T::Quantity) = SpectrumBase.blackbody(T).(wave)

# --- spectral axis conversions ----------------------------------------------

# Same-dimension targets are a plain unit conversion, keeping conversions idempotent.
_spectral_value(spectral_unit, x) =
    dimension(spectral_unit) == dimension(x) ? Unitful.uconvert(spectral_unit, x) :
    Unitful.uconvert(spectral_unit, x, Spectral())

_flux_value(flux_unit, x, at) =
    dimension(flux_unit) == dimension(x) ? Unitful.uconvert(flux_unit, x) :
    Unitful.uconvert(flux_unit, x, SpectralDensity(at))

_convert_spectral_axis(spec::Spectrum, spectral_unit) = _spectral_value.(spectral_unit, spectral_axis(spec))

_convert_flux_axis(spec::Spectrum, ::Nothing) = copy(flux_axis(spec))
function _convert_flux_axis(spec::Spectrum{<:Unitful.Quantity}, flux_unit::Unitful.Units)
    spec isa BinnedSpectrum && throw(ArgumentError(
        "flux density conversion is not defined for binned spectra: bin values are integrated quantities, so converting them requires a bin-integration convention"))
    return _flux_value.(flux_unit, flux_axis(spec), spectral_axis(spec))
end

function _convert_spectrum(spec::Spectrum{<:Unitful.Quantity}, spectral_unit::Unitful.Units, flux_unit::Union{Nothing, Unitful.Units})
    dimension(spectral_unit) ∈ SPECTRAL_DIMENSIONS || throw(ArgumentError(
        "spectral axis unit must be a wavelength, frequency, or energy unit; got $spectral_unit"))
    return Spectrum(_convert_spectral_axis(spec, spectral_unit), _convert_flux_axis(spec, flux_unit), deepcopy(meta(spec)))
end

@doc raw"""
    Unitful.uconvert(spectral_unit::Unitful.Units, spec::Spectrum)
    Unitful.uconvert(units::Tuple{Unitful.Units, Unitful.Units}, spec::Spectrum)

Convert the spectral axis of `spec` to `spectral_unit`, which may be any wavelength, frequency, or
photon-energy unit. Coordinates are related by the photon equivalence ``E = hν = hc/λ``.
A target of the same dimension as the current axis is a plain unit conversion.

The conversion is applied elementwise, so a reciprocal conversion (e.g., wavelength -->
frequency) yields an axis with the opposite direction of monotonicity.

Given a single unit, flux values are copied unchanged. Beware that a flux density then
remains "per" its original coordinate (e.g., a per-wavelength density against a frequency
axis). Pass a `(spectral_unit, flux_unit)` tuple, i.e., the shape returned by `unit(spec)`, to
also convert flux density values into the matching convention via [`UnitfulEquivalences.SpectralDensity`](https://github.com/sostock/UnitfulEquivalences.jl),
which preserves integrals: ``F_ν = F_λ λ^2 / c``.

# Examples

```jldoctest
julia> using Unitful, UnitfulAstro

julia> spec = spectrum([1.0, 1.5, 2.0]u"μm", [1.0, 2.0, 3.0]u"W/m^2/μm");

julia> νspec = uconvert((u"THz", u"Jy"), spec);

julia> issorted(spectral_axis(νspec); rev = true)
true

julia> flux_axis(νspec)[1] ≈ uconvert(u"Jy", 1.0u"W/m^2/μm" * (1.0u"μm")^2 / Unitful.c0)
true
```
"""
Unitful.uconvert(spectral_unit::Unitful.Units, spec::Spectrum{<:Unitful.Quantity}) = _convert_spectrum(spec, spectral_unit, nothing)
Unitful.uconvert(units::Tuple{Unitful.Units, Unitful.Units}, spec::Spectrum{<:Unitful.Quantity}) = _convert_spectrum(spec, units[1], units[2])

const _UNITLESS_AXIS_MSG = "`uconvert` requires a spectral axis with Unitful units. This spectrum's axis is unitless. Attach units first, e.g., `spectrum(axis * u\"angstrom\", flux)`"
Unitful.uconvert(::Unitful.Units, ::AbstractSpectrum) = throw(ArgumentError(_UNITLESS_AXIS_MSG))
Unitful.uconvert(::Tuple{Unitful.Units, Unitful.Units}, ::AbstractSpectrum) = throw(ArgumentError(_UNITLESS_AXIS_MSG))

end # module
