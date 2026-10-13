using SpectrumBase: SpectrumResampler
using DataInterpolations: LinearInterpolation, ExtrapolationType

# TODO: See if it makes sense to have an exported API for resample/resample! even though we are using external interpolators.

function resample(spec, new_wave)
    interp = LinearInterpolation(flux_axis(spec), spectral_axis(spec); extrapolation = ExtrapolationType.Constant)
    resampler = SpectrumResampler(spec, interp)
    return resampler(new_wave)
end

function resample(spec1::AbstractSpectrum, spec2::AbstractSpectrum)
    new_wave = spectral_axis(spec2)
    return resample(spec1, new_wave)
end

function mock_spectrum(n::Int = Int(1e3); use_units::Bool = false)
    wave = range(1e4, 3e4, length = n)
    sigma = 0.1 .* sin.(wave)
    T = 6700
    flux = @. 1e14 / (wave^5 * (exp(1 / (wave * T)) - 1)) ± sigma
    if use_units
        wave *= u"angstrom"
        flux *= u"erg/s/cm^2/angstrom"
    end
    spectrum(wave, flux; name = "Test Spectrum")
end

@testset "Resampler" begin
    spec = mock_spectrum()
    s, f = spectral_axis(spec), flux_axis(spec)
    interp = LinearInterpolation(f, s; extrapolation = ExtrapolationType.Constant)
    resampler = SpectrumResampler(spec, interp)
    expected = """
    SpectrumResampler(Float64, Measurements.Measurement{Float64})
      spec: $(typeof(spec))
      interpolator: $(typeof(interp))"""

    @test sprint(show, resampler) == expected
    @test spectral_axis(resampler) == s
    @test flux_axis(resampler) == f
end

@testset "Resampling" begin
    spec = mock_spectrum()
    s, f = spectral_axis(spec), flux_axis(spec)
    new_wave = range(minimum(s), maximum(s); length = Integer(length(s) ÷ 2.4))
    res_spec = resample(spec, new_wave)
    expected = """
    SingleSpectrum(Float64, Measurements.Measurement{Float64})
      spectral axis (416,): 10000.0 .. 30000.0
      flux axis (416,): 67.0 ± 0.031 .. 0.827 ± 0.08
      meta: Dict{Symbol, Any}(:name => "Test Spectrum")"""

    @test sprint(show, res_spec) == expected
    @test spectral_axis(res_spec) == new_wave
    @test length(flux_axis(res_spec)) == length(new_wave)

    # Unitful
    spec = mock_spectrum(; use_units = true)
    s, f = spectral_axis(spec), flux_axis(spec)
    new_wave = range(minimum(s), maximum(s); length = Integer(length(s) ÷ 2.4))
    @assert unit(eltype(new_wave)) == unit(eltype(s))
    res_spec = resample(spec, new_wave)

    @test spectral_axis(res_spec) == new_wave
    @test length(flux_axis(res_spec)) == length(new_wave)

    # Test resampling to another Spectrum
    spec1 = mock_spectrum(Integer(1e4))
    spec2 = mock_spectrum(Integer(1e3))
    res_spec = resample(spec1, spec2)
    @test spectral_axis(res_spec) == spectral_axis(spec2)

    # Unitful
    spec1 = mock_spectrum(Integer(1e4); use_units = true)
    spec2 = mock_spectrum(Integer(1e3); use_units = true)
    res_spec = resample(spec1, spec2)

    @test spectral_axis(res_spec) == spectral_axis(spec2)
    @test unit(eltype(flux_axis(spec1))) == unit(eltype(flux_axis(spec2))) == unit(eltype(flux_axis(res_spec)))

    # Test when spectra2 has different units
    #spec1 = mock_spectrum(Integer(1e4), use_units=true)
    #spec2 = mock_spectrum(Integer(1e3), use_units=true)
    #spec1.wave = uconvert.(u"cm", spec1.wave)
    #resample!(spec1, spec2)
    #@test ustrip.(spectral_axis(spec1)) ≈ ustrip.(unit(eltype(spectral_axis(spec1))), spectral_axis(spec2))

    # address bug when doing the same thing but affecting spec2
    #spec1 = mock_spectrum(Integer(1e4), use_units=true)
    #spec2 = mock_spectrum(Integer(1e3), use_units=true)
    #spec2.wave = uconvert.(u"cm", spec2.wave)
    #resample!(spec1, spec2)
    #@test ustrip.(spectral_axis(spec1)) ≈ ustrip.(unit(eltype(spectral_axis(spec1))), spectral_axis(spec2))
end
