using UnitfulEquivalences: Spectral, SpectralDensity

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
