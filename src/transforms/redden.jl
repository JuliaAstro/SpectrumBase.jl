import DustExtinction

# Extinction laws evaluate at a single wavelength per sample, which bin edges cannot provide.
const _BINNED_REDDEN_MSG = "reddening is not defined for binned spectra: applying an extinction law to bin edges requires a bin-integration convention"
redden!(::BinnedSpectrum, Av; kwargs...) = throw(ArgumentError(_BINNED_REDDEN_MSG))
deredden!(::BinnedSpectrum, Av; kwargs...) = throw(ArgumentError(_BINNED_REDDEN_MSG))

"""
    redden!(::AbstractSpectrum, Av; Rv = 3.1, law = DustExtinction.CCM89)

In-place version of [`redden`](@ref)
"""
function redden!(spec::AbstractSpectrum, Av; Rv = 3.1, law = DustExtinction.CCM89)
    law_instance = DustExtinction.lawinstance(law; Rv)
    DustExtinction.redden!(law_instance, spectral_axis(spec), flux_axis(spec); Av)
    return spec
end

"""
    redden(::AbstractSpectrum, Av; Rv = 3.1, law = DustExtinction.CCM89)

Redden a spectrum using common color laws provided by [DustExtinction.jl](https://github.com/juliaastro/dustextinction.jl).
`Av` is the total extinction, `Rv` is the selective extinction (3.1 is a common
value for the Milky Way) and `law` is the color law to use for determining the
extinction.
"""
redden(spec::AbstractSpectrum, Av; kwargs...) = redden!(deepcopy(spec), Av; kwargs...)

"""
    deredden!(::AbstractSpectrum, Av; Rv = 3.1, law = DustExtinction.CCM89)

In-place version of [`deredden`](@ref)
"""
function deredden!(spec::AbstractSpectrum, Av; Rv = 3.1, law = DustExtinction.CCM89)
    law_instance = DustExtinction.lawinstance(law; Rv)
    DustExtinction.deredden!(law_instance, spectral_axis(spec), flux_axis(spec); Av)
    return spec
end

"""
    deredden(::AbstractSpectrum, Av; Rv = 3.1, law = DustExtinction.CCM89)

Deredden a spectrum using common color laws provided by [DustExtinction.jl](https://github.com/juliaastro/dustextinction.jl).
`Av` is the total extinction, `Rv` is the selective extinction (3.1 is a common
value for the Milky Way) and `law` is the color law to use for determining the
extinction.
"""
deredden(spec::AbstractSpectrum, Av; kwargs...) = deredden!(deepcopy(spec), Av; kwargs...)
