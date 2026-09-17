using ParallelTestRunner: runtests, find_tests, parse_args
import SpectrumBase

const init_code = quote
    using SpectrumBase: SpectrumBase, AbstractSpectrum, spectrum, spectral_axis, flux_axis
    using Measurements: Measurements, ±
    using Unitful: @u_str, uconvert, unit, ustrip
    using UnitfulAstro
    import Random

    const C_MPS = 299792458.0
    const C_KMPS = C_MPS / 1000  # speed of light in km/s, the bare-Real velocity unit
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(SpectrumBase, args; testsuite, init_code)
