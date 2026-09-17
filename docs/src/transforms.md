# Transformations

## Extinction

By leveraging [DustExtinction.jl](https://github.com/juliaastro/dustextinction.jl) we can apply common reddening laws to our spectra.

```@example
using SpectrumBase, Unitful, Measurements, Random
rng = Random.seed!(0);
wave = (1:0.5:3)u"μm";
sigma = randn(rng, size(wave));
flux = (100 .± sigma)u"W/m^2/μm"
spec = spectrum(wave, flux)
red = redden(spec, 0.3)
flux_axis(red)
deredden!(red, 0.3)
flux_axis(red) ≈ flux_axis(spec)
```

### API/Reference

```@docs
redden
redden!
deredden
deredden!
```

## Redshift

Spectra can be shifted in wavelength via cosmological redshift or Doppler velocity. Both transformations operate only on the spectral axis; flux values are preserved as-is.

**Cosmological redshift** shifts by a dimensionless parameter `z`:

```jldoctest
julia> using SpectrumBase

julia> spec = spectrum(collect(4000.0:1000.0:8000.0), ones(5));

julia> shifted = redshift(spec, 0.5);

julia> spectral_axis(shifted) ≈ spectral_axis(spec) .* 1.5
true
```

**Doppler shift** shifts by a radial velocity `v`. Pass a `Unitful` velocity or a plain number (interpreted as km/s). Set `relativistic=true` for the full relativistic formula:

```jldoctest
julia> using SpectrumBase, Unitful

julia> spec = spectrum(collect(4000.0:1000.0:8000.0), ones(5));

julia> shifted = doppler_shift(spec, 100u"km/s");

julia> shifted_rel = doppler_shift(spec, 100u"km/s"; relativistic=true);
```

Both `redshift` and `doppler_shift` return a new spectrum. In-place variants `redshift!` and `doppler_shift!` are also available.

Binned spectra are supported as well: both bin edges are shifted and the bin values are carried unchanged. An energy axis is compressed rather than stretched:

```jldoctest
julia> using SpectrumBase, Unitful

julia> spec = spectrum([1.0 2.0; 2.0 3.0; 3.0 4.0]u"keV", [10.0, 20.0, 30.0]u"erg/s/cm^2");

julia> shifted = redshift(spec, 0.5);

julia> spectral_axis(shifted) ≈ spectral_axis(spec) ./ 1.5
true

julia> flux_axis(shifted) == flux_axis(spec)
true
```

### API/Reference

```@docs
redshift
redshift!
doppler_shift
doppler_shift!
```

## Axis conversions

Powered by [UnitfulEquivalences.jl](https://github.com/sostock/UnitfulEquivalences.jl), [`uconvert`](@ref) from [Unitful.jl](https://github.com/JuliaPhysics/Unitful.jl) is extended to convert the spectral axis between wavelength, frequency, and photon energy, along with its associated flux density:

```jldoctest
julia> using SpectrumBase, Unitful, UnitfulAstro

julia> spec = spectrum([1.0, 1.5, 2.0]u"μm", [1.0, 2.0, 3.0]u"W/m^2/μm");

julia> uconvert(u"THz", spec);

julia> νspec = uconvert((u"THz", u"Jy"), spec);

julia> issorted(spectral_axis(νspec); rev = true)
true

julia> flux_axis(νspec)[1] ≈ uconvert(u"Jy", 1.0u"W/m^2/μm" * (1.0u"μm")^2 / Unitful.c0)
true
```

### API/Reference

```@docs
Unitful.uconvert
```

## Resampling

External interpolators, e.g., from [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) or [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl), can be used to resample spectra onto a given wavelength grid. Starting with a sample spectrum `spec`, we first create a [`SpectrumResampler`](@ref) object `resampler` which stores the initial spectrum and interpolator `interp` together. We then apply this object to the wavelength grid of our choice to produce the resampled spectrum. We show example usage in the docstring below:

```@docs
SpectrumResampler
```
