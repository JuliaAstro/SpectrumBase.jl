module SpectrumBaseDynamicQuantitiesExt

using SpectrumBase
using SpectrumBase: AbstractSpectrum, Spectrum, spectrum, spectral_axis, flux_axis, meta, blackbody
import SpectrumBase: _ustrip, _velocity_ratio, _spectral_axis_factor
using DynamicQuantities: DynamicQuantities, @u_str, ustrip, dimension, uexpand,
    UnionAbstractQuantity, Constants

# DynamicQuantities does not encode a quantity's dimension in its type (all
# `Quantity{Float64,Dimensions}` share one concrete type), so — unlike the
# Unitful backend — every check here inspects the *value*. `uexpand` normalizes
# symbolic units (`us"..."`) to base `Dimensions` so the comparisons below work
# regardless of how the axis was constructed.
const _LENGTH_DIM = dimension(u"m")
const _FREQUENCY_DIM = dimension(u"Hz")
const _ENERGY_DIM = dimension(u"J")

_dimension(q) = dimension(uexpand(q))
_is_wavelength(q) = _dimension(q) == _LENGTH_DIM
_is_frequency(q) = _dimension(q) == _FREQUENCY_DIM
_is_energy(q) = _dimension(q) == _ENERGY_DIM

# --- units-agnostic hooks from the core -------------------------------------

_ustrip(x::UnionAbstractQuantity) = ustrip(x)

_velocity_ratio(v::UnionAbstractQuantity) = ustrip(uexpand(v / Constants.c))

# Resolved from the axis value rather than the element type (see note above): a
# wavelength axis stretches, a frequency or energy axis compresses.
function _spectral_axis_factor(spec::AbstractSpectrum{<:UnionAbstractQuantity}, factor)
    q = first(spectral_axis(spec))
    _is_wavelength(q) && return factor
    (_is_frequency(q) || _is_energy(q)) && return inv(factor)
    throw(ArgumentError("cannot shift a spectral axis with dimension $(_dimension(q)); expected a wavelength, frequency, or energy"))
end

# --- spectrum constructors for DynamicQuantities axes -----------------------

function _assert_spectral_dimension(axis)
    q = first(axis)
    @assert (_is_wavelength(q) || _is_frequency(q) || _is_energy(q)) "spectral_axis not recognized as having dimensions of wavelength, frequency, or energy."
end

function SpectrumBase.spectrum(spectral_axis::AbstractVector{<:UnionAbstractQuantity}, flux_axis::AbstractVector{<:UnionAbstractQuantity}; kwds...)
    _assert_spectral_dimension(spectral_axis)
    Spectrum(spectral_axis, flux_axis, Dict{Symbol, Any}(kwds))
end

function SpectrumBase.spectrum(spectral_axis::AbstractMatrix{<:UnionAbstractQuantity}, flux_axis::AbstractMatrix{<:UnionAbstractQuantity}; kwds...)
    _assert_spectral_dimension(spectral_axis)
    Spectrum(spectral_axis, flux_axis, Dict{Symbol, Any}(kwds))
end

function SpectrumBase.spectrum(spectral_axis::AbstractMatrix{<:UnionAbstractQuantity}, flux_axis::AbstractVector{<:UnionAbstractQuantity}; kwds...)
    _assert_spectral_dimension(spectral_axis)
    Spectrum(spectral_axis, flux_axis, Dict{Symbol, Any}(kwds))
end

# --- ustrip / dimension on whole spectra ------------------------------------

"""
    DynamicQuantities.ustrip(::AbstractSpectrum)

Remove the units from a spectrum, returning a plain-number spectrum with the axis
values expressed in SI base units. The DynamicQuantities analog of the Unitful
`ustrip(::AbstractSpectrum)` method.
"""
DynamicQuantities.ustrip(spec::AbstractSpectrum) = spectrum(ustrip.(spectral_axis(spec)), ustrip.(flux_axis(spec)); meta(spec)...)

"""
    DynamicQuantities.dimension(::AbstractSpectrum)

Get the dimensions of a spectrum as a tuple of the spectral-axis and flux-axis
dimensions. The DynamicQuantities analog of the Unitful `unit(::AbstractSpectrum)`
method (DynamicQuantities has no unit objects, only dimensions).
"""
DynamicQuantities.dimension(spec::AbstractSpectrum) = _dimension(first(spectral_axis(spec))), _dimension(first(flux_axis(spec)))

# --- blackbody with DynamicQuantities arguments -----------------------------

function SpectrumBase.blackbody(wave::AbstractVector{<:UnionAbstractQuantity}, T::UnionAbstractQuantity)
    flux = SpectrumBase.blackbody(T).(wave)
    return spectrum(wave, flux, name = "Blackbody", T = T)
end

"""
    blackbody(T::DynamicQuantities.UnionAbstractQuantity)

Returns a function for calculating blackbody curves. The flux is returned in SI
base units (W m⁻² m⁻¹).
"""
SpectrumBase.blackbody(T::UnionAbstractQuantity) =
    w -> 2 * Constants.h * Constants.c^2 / w^5 / (exp(Constants.h * Constants.c / (w * Constants.k_B * T)) - 1)

end # module
