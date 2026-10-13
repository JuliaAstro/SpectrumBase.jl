import DynamicQuantities as DQ
using SpectrumBase: SingleSpectrum, BinnedSpectrum, redshift, redshift!, doppler_shift, doppler_shift!, blackbody

# DynamicQuantities normalizes to SI base units, so `angstrom` and `keV` are
# stored as their SI-base values; `DQ.ustrip` returns those base-unit numbers.
Å = DQ.u"angstrom"
keV = DQ.Constants.keV
Wm2 = DQ.u"W/m^2"

@testset "DynamicQuantities extension" begin
    @testset "spectrum constructor" begin
        wave = collect(range(4000.0, 8000.0, length = 10)) .* Å
        flux = ones(10) .* Wm2
        spec = spectrum(wave, flux)
        @test spec isa SingleSpectrum

        # energy and frequency axes are also accepted
        energy = collect(range(1.0, 5.0, length = 10)) .* keV
        @test spectrum(energy, flux) isa SingleSpectrum
        freq = collect(range(1.0, 5.0, length = 10)) .* DQ.u"THz"
        @test spectrum(freq, flux) isa SingleSpectrum

        # bin edges
        @test spectrum([1.0 2.0; 2.0 3.0] .* keV, [1.0, 2.0] .* Wm2) isa BinnedSpectrum

        # a non-wavelength/frequency/energy axis is rejected
        bad = collect(range(1.0, 10.0, length = 10)) .* DQ.u"s"
        @test_throws AssertionError spectrum(bad, flux)
    end

    @testset "Wavelength axis stretches" begin
        wave = collect(range(4000.0, 8000.0, length = 10)) .* Å
        spec = spectrum(wave, ones(10) .* Wm2)

        shifted = @inferred redshift(spec, 0.5)
        @test all(spectral_axis(shifted) .≈ wave .* 1.5)
        @test flux_axis(shifted) == flux_axis(spec)

        # roundtrip
        @test all(spectral_axis(redshift(shifted, -1 / 3)) .≈ wave)
    end

    @testset "Energy axis compresses" begin
        energy = collect(range(1.0, 5.0, length = 10)) .* keV
        spec = spectrum(energy, ones(10) .* Wm2)

        shifted = @inferred redshift(spec, 0.5)
        @test all(spectral_axis(shifted) .≈ energy ./ 1.5)
        @test all(spectral_axis(redshift(shifted, -1 / 3)) .≈ energy)
    end

    @testset "Frequency axis compresses" begin
        freq = collect(range(1.0, 5.0, length = 10)) .* DQ.u"THz"
        spec = spectrum(freq, ones(10) .* Wm2)

        shifted = @inferred redshift(spec, 0.5)
        @test all(spectral_axis(shifted) .≈ freq ./ 1.5)
    end

    @testset "redshift in-place" begin
        wave = collect(range(4000.0, 8000.0, length = 10)) .* Å
        spec = spectrum(wave, ones(10) .* Wm2)
        expected = spectral_axis(spec) .* 1.5
        result = @inferred redshift!(spec, 0.5)
        @test result === spec
        @test all(spectral_axis(spec) .≈ expected)
    end

    @testset "Doppler shift" begin
        wave = collect(range(4000.0, 8000.0, length = 10)) .* Å
        spec = spectrum(wave, ones(10) .* Wm2)

        # bare Real velocity is km/s, as in the numeric core
        v = 100.0
        @test all(spectral_axis(doppler_shift(spec, v)) .≈ wave .* (1 + v / C_KMPS))

        # DynamicQuantities velocity honors its own units
        vq = 100.0 * DQ.u"km/s"
        β = 1e5 / C_MPS
        @test all(spectral_axis(doppler_shift(spec, vq)) .≈ wave .* (1 + β))

        shifted_rel = doppler_shift(spec, vq; relativistic = true)
        @test all(spectral_axis(shifted_rel) .≈ wave .* sqrt((1 + β) / (1 - β)))

        # superluminal is rejected in the relativistic case
        @test_throws DomainError doppler_shift(spec, DQ.Constants.c; relativistic = true)
    end

    @testset "ustrip and dimension" begin
        wave = collect(range(4000.0, 8000.0, length = 10)) .* Å
        spec = spectrum(wave, ones(10) .* Wm2)

        stripped = DQ.ustrip(spec)
        @test eltype(spectral_axis(stripped)) <: Real
        @test spectral_axis(stripped) ≈ DQ.ustrip.(wave)

        w_dim, f_dim = DQ.dimension(spec)
        @test w_dim == DQ.dimension(DQ.u"m")
        @test f_dim == DQ.dimension(Wm2)
    end

    @testset "blackbody" begin
        wave = collect(range(1.0, 3.0, length = 10)) .* DQ.u"μm"
        bb = blackbody(wave, 2000.0 * DQ.u"K")
        @test bb isa SingleSpectrum
        @test spectral_axis(bb) == wave
        @test all(DQ.ustrip.(flux_axis(bb)) .> 0)
        @test bb.T == 2000.0 * DQ.u"K"
    end
end
