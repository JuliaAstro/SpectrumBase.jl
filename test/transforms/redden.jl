using SpectrumBase:
    redden,
    redden!,
    deredden,
    deredden!

using DustExtinction: DustExtinction

function mock_spectrum(n::Int = Int(1e3); use_units::Bool = false)
    wave = range(1e4, 3e4, length = n)
    sigma = 0.1 .* sin.(wave)
    T = 6700
    flux = @. 1e14 / (wave^5 * (exp(1 / (wave * T)) - 1)) ± sigma
    if use_units
        wave *= u"angstrom"
        flux *= u"erg/s/cm^2/angstrom"
    end
    spectrum(wave, flux, name="Test Spectrum")
end

@testset "redden" begin
    # Custom law
    struct CustomLaw <: DustExtinction.ExtinctionLaw
        Rv::Float64
    end
    (law::CustomLaw)(wave) = sin(law.Rv * wave) / law.Rv

    @testset "Reddening Av=$Av" for Av = [0.5, 1.0, 2.0]
        # Regression usage
        spec = mock_spectrum()
        spec2 = deepcopy(spec)
        reddened = @inferred redden(spec, Av)
        @inferred redden!(spec2, Av)
        @test flux_axis(reddened) ≈ flux_axis(spec2)
        @inferred deredden!(spec2, Av)
        @test flux_axis(spec) ≈ flux_axis(spec2)
        dereddened = @inferred deredden(reddened, Av)
        @test flux_axis(dereddened) ≈ flux_axis(spec)

        # Custom law
        s, f = spectral_axis(spec), flux_axis(spec)
        expected = @. f * 10^(-0.4 * Av * CustomLaw(π)(s))
        @test expected ≈ redden(spec, Av; law=CustomLaw, Rv=π) |> flux_axis

        # Bad law
        @test_throws MethodError redden(spec, Av, law = sin)

        # Unitful
        spec = mock_spectrum(use_units=true)
        spec2 = deepcopy(spec)
        reddened = @inferred redden(spec, Av)
        @inferred redden!(spec2, Av)
        @test flux_axis(reddened) ≈ flux_axis(spec2)
        @inferred deredden!(spec2, Av)
        @test flux_axis(spec) ≈ flux_axis(spec2)
        dereddened = @inferred deredden(reddened, Av)
        @test flux_axis(dereddened) ≈ flux_axis(spec)
    end

    @testset "Reddening Av=0" begin
        # Test no change when Av = 0
        Av = 0.0
        spec = mock_spectrum()
        reddened = redden(spec, Av)
        @test flux_axis(reddened) ≈ flux_axis(spec)
        dereddened = deredden(reddened, Av)
        @test flux_axis(dereddened) ≈ flux_axis(spec)
    end
end

@testset "Binned spectra are rejected" begin
    binned = spectrum([4000.0 4010.0; 4010.0 4020.0], [1.0, 2.0])
    @test_throws ArgumentError redden(binned, 0.3)
    @test_throws ArgumentError redden!(binned, 0.3)
    @test_throws ArgumentError deredden(binned, 0.3)
    @test_throws ArgumentError deredden!(binned, 0.3)
end
