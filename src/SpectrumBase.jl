module SpectrumBase

# Uniform API
export AbstractSpectrum, Spectrum, spectrum, spectral_axis, flux_axis

# AbstractSpectrum types
export SingleSpectrum, IFUSpectrum, EchelleSpectrum, BinnedSpectrum

# Transforms
export SpectrumResampler, redden, redden!, deredden, deredden!
export redshift, redshift!, doppler_shift, doppler_shift!

# Utilities
export blackbody #, line_flux, equivalent_width

# Fitting
#export continuum, continuum!

using RecipesBase: @recipe
using Measurements: Measurements
using Unitful: Unitful, @u_str, ustrip, unit, dimension
using PhysicalConstants.CODATA2018: h, c_0, k_B

"""
    AbstractSpectrum{S <: Number, F <: Number}

An abstract holder for astronomical spectra. All types inheriting from this must have the following fields:

* `spectral_axis::Array{S, M}`
* `flux_axis::Array{F, N}`
* `meta::Dict{Symbol, Any}`

See [`SingleSpectrum`](@ref), [`EchelleSpectrum`](@ref), and [`IFUSpectrum`](@ref) for different subtypes.
"""
abstract type AbstractSpectrum{S, F} end

"""
    Spectrum <: AbstractSpectrum

A spectrum or spectra stored as arrays of real numbers. For UV/VIS/IR spectra, the `spectral_axis` is assumed to be wavelengths (in angstrom). For X-ray spectra, the `spectral_axis` is assumed to be energies (in keV).
"""
mutable struct Spectrum{S<:Number, F<:Number, M, N, A<:AbstractArray{S, M}, B<:AbstractArray{F, N}} <: AbstractSpectrum{S, F}
    spectral_axis::A
    flux_axis::B
    meta::Dict{Symbol, Any}
    function Spectrum{S, F, M, N, A, B}(s, f, meta) where {S<:Number, F<:Number, M, N, A<:AbstractArray{S, M}, B<:AbstractArray{F, N}}
        # TODO: Investigate using Holy Traits to help with validation
        # Dimension compatibility check
        if size(s, 1) != size(f, 1)
            throw(ArgumentError(
            """
            Spectral axis and flux axis are incompatible sizes. Currently supported sizes are:

            - SingleSpectrum: spectral axis (M-length vector), flux axis (M-length vector)
            - EchelleSpectrum: spectral axis (M x N matrix), flux axis (M x N matrix)
            - IFUSpectrum: spectral axis (M-length vector), flux axis (M x N x K matrix)
            - BinnedSpectrum: spectral axis (M x 2 matrix of bin edges), flux axis (M-length vector)

            See the documentation for each spectrum type for more.
            """))
        end

        # Spectral axis monoticity check
        spec_ax = eachcol(s)
        if !(all(issorted, spec_ax) || all(x -> issorted(x; rev = true), spec_ax))
            throw(ArgumentError("Spectral axis must be strictly increasing or decreasing."))
        end

        # BinnedSpectrum invariants: m × 2 bin edges, each row strictly ordered in the axis direction
        if M == 2 && N == 1
            if size(s, 2) != 2
                throw(ArgumentError("BinnedSpectrum spectral axis must be an m × 2 matrix of bin edges; got size $(size(s))"))
            end
            rows = eachrow(s)
            ordered = all(issorted, spec_ax) ? all(r -> first(r) < last(r), rows) :
                all(r -> first(r) > last(r), rows)
            if !ordered
                throw(ArgumentError("BinnedSpectrum bin edges must be strictly ordered within each row, in the direction of the spectral axis."))
            end
        end

        return new{S, F, M, N, A, B}(s, f, meta)
    end
end

function Spectrum(s, f, meta)
    Spectrum{eltype(s), eltype(f), ndims(s), ndims(f), typeof(s), typeof(f)}(s, f, meta)
end

# Doesn't seem to be used atp
#Spectrum(wave, flux, meta::Dict{Symbol, Any}) = Spectrum(collect(wave), collect(flux), meta)

"""
    spectral_axis(spec::AbstractSpectrum)

Return the spectral axis of `spec`.
"""
spectral_axis(spec::AbstractSpectrum) = spec.spectral_axis

"""
    flux_axis(spec::AbstractSpectrum)

Return the flux axis of `spec`.
"""
flux_axis(spec::AbstractSpectrum) = spec.flux_axis

"""
    meta(spec::AbstractSpectrum)

Return the meta data associated with `spec`.
"""
meta(spec::AbstractSpectrum) = spec.meta

# Fields first: `hasfield` is decidable at compile time, keeping field access inferable.
function Base.getproperty(spec::AbstractSpectrum, nm::Symbol)
    if hasfield(typeof(spec), nm)
        return getfield(spec, nm)
    else
        return getfield(spec, :meta)[nm]
    end
end

function Base.propertynames(spec::AbstractSpectrum)
    natural = (:spectral_axis, :flux_axis, :meta)
    m = keys(meta(spec))
    return (natural..., m...)
end

# Collection
Base.argmax(spec::AbstractSpectrum) = argmax(flux_axis(spec))
Base.argmin(spec::AbstractSpectrum) = argmin(flux_axis(spec))
Base.eltype(spec::AbstractSpectrum) = eltype(flux_axis(spec))
Base.findmax(spec::AbstractSpectrum) = findmax(flux_axis(spec))
Base.findmin(spec::AbstractSpectrum) = findmin(flux_axis(spec))
Base.length(spec::AbstractSpectrum) = length(flux_axis(spec))
Base.maximum(spec::AbstractSpectrum) = maximum(flux_axis(spec))
Base.minimum(spec::AbstractSpectrum) = minimum(flux_axis(spec))
Base.size(spec::AbstractSpectrum) = size(flux_axis(spec))
Base.size(spec::AbstractSpectrum, i) = size(flux_axis(spec), i)
function Base.iterate(spec::AbstractSpectrum, state=0)
    state == length(spec) && return nothing
    return spec[begin + state], state + 1
end

# Arithmetic
Base.:(==)(s::AbstractSpectrum, o::AbstractSpectrum) = spectral_axis(s) == spectral_axis(o) && flux_axis(s) == flux_axis(o) && meta(s) == meta(o)
Base.:+(s::T, A) where {T <: AbstractSpectrum} = T(spectral_axis(s), flux_axis(s) .+ A, meta(s))
Base.:*(s::T, A::Union{Real, AbstractVector}) where {T <: AbstractSpectrum} = T(spectral_axis(s), flux_axis(s) .* A, meta(s))
Base.:/(s::T, A) where {T <: AbstractSpectrum} = T(spectral_axis(s), flux_axis(s) ./ A, meta(s))
Base.:-(s::T) where {T <: AbstractSpectrum} = T(spectral_axis(s), -flux_axis(s), meta(s))
Base.:-(s::AbstractSpectrum, A) = s + -A
Base.:-(A, s::AbstractSpectrum) = s - A
Base.:-(s::AbstractSpectrum, o::AbstractSpectrum) = s - o # Satisfy Aqua

# Multi-Spectrum
Base.:+(s::T, o::T) where {T <: AbstractSpectrum} = T(spectral_axis(s), flux_axis(s) .+ flux_axis(o), meta(s))
Base.:*(s::T, o::T) where {T <: AbstractSpectrum} = T(spectral_axis(s), flux_axis(s) .* flux_axis(o), meta(s))
Base.:/(s::T, o::T) where {T <: AbstractSpectrum} = T(spectral_axis(s), flux_axis(s) ./ flux_axis(o) * unit(s)[2], meta(s))
Base.:-(s::T, o::T) where {T <: AbstractSpectrum} = T(spectral_axis(s), flux_axis(s) .- flux_axis(o), meta(s))

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

# Spectrum types and basic arithmetic
include("spectrum_single.jl")
include("spectrum_binned.jl")
include("spectrum_echelle.jl")
include("spectrum_ifu.jl")

"""
    spectrum(spectral_axis, flux_axis, [meta])

Construct a spectrum given the spectral axis and flux axis. This will automatically dispatch the correct spectrum type given the shape and element type of the given flux. Any keyword arguments will be accessible from the spectrum as properties.

# Examples
```jldoctest
julia> wave = range(1e4, 4e4, length=1000);

julia> flux = 100 .* ones(size(wave));

julia> spec = spectrum(wave, flux)
SingleSpectrum(Float64, Float64)
  spectral axis (1000,): 10000.0 .. 40000.0
  flux axis (1000,): 100.0 .. 100.0
  meta: Dict{Symbol, Any}()

julia> spec = spectrum(wave, flux, name="Just Noise")
SingleSpectrum(Float64, Float64)
  spectral axis (1000,): 10000.0 .. 40000.0
  flux axis (1000,): 100.0 .. 100.0
  meta: Dict{Symbol, Any}(:name => "Just Noise")

julia> spec.name
"Just Noise"
```

There is easy integration with [Unitful.jl](https://github.com/JuliaPhysics/Unitful.jl)
and its sub-projects and [Measurements.jl](https://github.com/juliaphysics/measurements.jl)

```jldoctest
julia> using Random

julia> rng = Random.seed!(0)
TaskLocalRNG()

julia> using Unitful, UnitfulAstro, Measurements

julia> wave = range(1, 4, length=1000)u"μm";

julia> sigma = randn(rng, size(wave));

julia> flux = (100 .± sigma)u"erg/cm^2/s/angstrom";

julia> spec = spectrum(wave, flux)
SingleSpectrum(Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(μm,), 𝐋, nothing}}, Unitful.Quantity{Measurements.Measurement{Float64}, 𝐌 𝐋^-1 𝐓^-3, Unitful.FreeUnits{(Å^-1, erg, cm^-2, s^-1), 𝐌 𝐋^-1 𝐓^-3, nothing}})
  spectral axis (1000,): 1.0 μm .. 4.0 μm
  flux axis (1000,): 100.0 ± -0.23 erg Å^-1 cm^-2 s^-1 .. 100.0 ± 0.25 erg Å^-1 cm^-2 s^-1
  meta: Dict{Symbol, Any}()
```

For an echelle spectrum, all orders must have the same length, so be sure to pad any ragged orders with NaN.

```jldoctest
julia> wave = reshape(range(100, 1e4, length=1000), 100, 10);

julia> flux = repeat(1:10.0, 1, 100)';

julia> spec = spectrum(wave, flux)
EchelleSpectrum(Float64, Float64)
  # orders: 10
  spectral axis (100, 10): 100.0 .. 10000.0
  flux axis (100, 10): 1.0 .. 10.0
  meta: Dict{Symbol, Any}()
```
"""
function spectrum(spectral_axis::AbstractVector{<:Real}, flux_axis::AbstractVector{<:Real}; kwds...)
    Spectrum(spectral_axis, flux_axis, Dict{Symbol,Any}(kwds))
end

function spectrum(spectral_axis::AbstractVector{<:Real}, flux_axis::AbstractArray{<:Real, 3}; kwds...)
    Spectrum(spectral_axis, flux_axis, Dict{Symbol,Any}(kwds))
end

function spectrum(spectral_axis::AbstractMatrix{<:Real}, flux_axis::AbstractMatrix{<:Real}; kwds...)
    Spectrum(spectral_axis, flux_axis, Dict{Symbol,Any}(kwds))
end

function spectrum(spectral_axis::AbstractMatrix{<:Real}, flux_axis::AbstractVector{<:Real}; kwds...)
    Spectrum(spectral_axis, flux_axis, Dict{Symbol,Any}(kwds))
end

# Recognized spectral-axis dimensions: wavelength, frequency, photon energy.
const SPECTRAL_DIMENSIONS = (u"𝐋", u"𝐓^-1", u"𝐋^2 * 𝐌 * 𝐓^-2")

function spectrum(spectral_axis::AbstractVector{<:Unitful.Quantity}, flux_axis::AbstractVector{<:Unitful.Quantity}; kwds...)
    @assert dimension(eltype(spectral_axis)) ∈ SPECTRAL_DIMENSIONS "spectral_axis not recognized as having dimensions of wavelength, frequency, or energy."
    Spectrum(spectral_axis, flux_axis, Dict{Symbol,Any}(kwds))
end

function spectrum(spectral_axis::AbstractMatrix{<:Unitful.Quantity}, flux_axis::AbstractMatrix{<:Unitful.Quantity}; kwds...)
    @assert dimension(eltype(spectral_axis)) ∈ SPECTRAL_DIMENSIONS "spectral_axis not recognized as having dimensions of wavelength, frequency, or energy."
    Spectrum(spectral_axis, flux_axis, Dict{Symbol,Any}(kwds))
end

function spectrum(spectral_axis::AbstractMatrix{<:Unitful.Quantity}, flux_axis::AbstractVector{<:Unitful.Quantity}; kwds...)
    @assert dimension(eltype(spectral_axis)) ∈ SPECTRAL_DIMENSIONS "spectral_axis not recognized as having dimensions of wavelength, frequency, or energy."
    Spectrum(spectral_axis, flux_axis, Dict{Symbol,Any}(kwds))
end

include("utils.jl")
include("transforms/resampler.jl")
include("transforms/transforms.jl")
include("plotting.jl")
#include("fitting/fitting.jl")

end # module
