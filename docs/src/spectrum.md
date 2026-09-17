# Spectrum

Here we will go over the different spectral types and how we use them.

## Types

SpectrumBase are defined as possible subtypes of `AbstractSpectrum`. You can use these directly for construction, or use the catch-all [`spectrum`](@ref) function, which is preferred.

```@docs
SpectrumBase.AbstractSpectrum
SpectrumBase.Spectrum
SpectrumBase.SingleSpectrum
SpectrumBase.EchelleSpectrum
SpectrumBase.IFUSpectrum
SpectrumBase.BinnedSpectrum
```

## Constructors

```@docs
SpectrumBase.spectrum
```


## Basic operations

For more advanced transformations, see [Transformations](@ref)

### Getters
```@docs
SpectrumBase.spectral_axis(::AbstractSpectrum)
SpectrumBase.flux_axis(::AbstractSpectrum)
SpectrumBase.meta(::AbstractSpectrum)
```

### Array interface

| Function                           |
|:-----------------------------------|
| `Base.argmax(::AbstractSpectrum)`  |
| `Base.argmin(::AbstractSpectrum)`  |
| `Base.eltype(::AbstractSpectrum)`  |
| `Base.findmax(::AbstractSpectrum)` |
| `Base.findmin(::AbstractSpectrum)` |
| `Base.iterate(::AbstractSpectrum)` |
| `Base.length(::AbstractSpectrum)`  |
| `Base.maximum(::AbstractSpectrum)` |
| `Base.minimum(::AbstractSpectrum)` |
| `Base.size(::AbstractSpectrum)`    |

### Arithmetic

| Function                                            |
|:----------------------------------------------------|
| `+(::AbstractSpectrum, A)`                          |
| `-(::AbstractSpectrum, A)`                          |
| `*(::AbstractSpectrum, A)`                          |
| `/(::AbstractSpectrum, A)`                          |
| `Base.(==)(::AbstractSpectrum, ::AbstractSpectrum)` |

## Unitful helpers

```@docs
Unitful.unit
Unitful.ustrip
```

## Utilities

```@docs
SpectrumBase.blackbody
```

## Plotting

We provide simple plotting recipes for spectra using [Plots.jl](https://github.com/juliaplots/plots.jl)

```@example
using Plots, SpectrumBase

wave = range(1e3, 5e4, length=100)
spec = blackbody(wave, 2000)

plot(spec)
```

## Index

```@index
Pages = ["spectrum.md"]
```
