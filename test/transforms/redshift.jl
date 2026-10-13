using SpectrumBase:
    BinnedSpectrum,
    redshift,
    redshift!,
    doppler_shift,
    doppler_shift!

function mock_spectrum_redshift(; use_units::Bool = false)
    wave = collect(range(4000.0, 8000.0, length = 100))
    flux = @. 1e14 / (wave^5 * (exp(1 / (wave * 6700)) - 1))
    if use_units
        wave *= u"angstrom"
        flux *= u"erg/s/cm^2/angstrom"
    end
    spectrum(wave, flux, name = "Test Redshift Spectrum")
end

@testset "redshift" begin
    @testset "Cosmological redshift z=$z" for z in [0.0, 0.1, 0.5, 1.0, 2.0]
        spec = mock_spectrum_redshift()
        original_wave = copy(spectral_axis(spec))
        original_flux = copy(flux_axis(spec))
        original_meta = deepcopy(spec.meta)

        shifted = @inferred redshift(spec, z)

        @test spectral_axis(shifted) ≈ original_wave .* (1 + z)
        @test flux_axis(shifted) ≈ original_flux
        @test shifted.meta == original_meta

        # Original unchanged
        @test spectral_axis(spec) ≈ original_wave
        @test flux_axis(spec) ≈ original_flux
    end

    @testset "Cosmological redshift in-place" begin
        spec = mock_spectrum_redshift()
        expected_wave = spectral_axis(spec) .* 1.5

        result = @inferred redshift!(spec, 0.5)
        @test result === spec
        @test spectral_axis(spec) ≈ expected_wave
    end

    @testset "Cosmological redshift roundtrip" begin
        spec = mock_spectrum_redshift()
        original_wave = copy(spectral_axis(spec))

        shifted = redshift(spec, 1.0)
        unshifted = redshift(shifted, -0.5)
        @test spectral_axis(unshifted) ≈ original_wave
    end

    @testset "Cosmological redshift with Unitful" begin
        spec = mock_spectrum_redshift(use_units = true)
        original_wave = copy(spectral_axis(spec))

        shifted = @inferred redshift(spec, 0.5)
        @test spectral_axis(shifted) ≈ original_wave .* 1.5
        @test flux_axis(shifted) ≈ flux_axis(spec)
    end

    @testset "Energy axis compresses" begin
        # A redshift lowers energy: axis ÷ (1 + z), not × (1 + z).
        axis = collect(range(1.0, 5.0, length = 10)) * u"keV"
        spec = spectrum(axis, ones(10) * u"erg/s/cm^2/keV")
        shifted = @inferred redshift(spec, 0.5)
        @test spectral_axis(shifted) ≈ axis ./ 1.5
        # Round-trips back to the original through the inverse shift
        @test spectral_axis(redshift(shifted, -1 / 3)) ≈ axis
    end

    @testset "Wavelength axis stretches" begin
        axis = collect(range(4000.0, 8000.0, length = 10)) * u"angstrom"
        spec = spectrum(axis, ones(10) * u"erg/s/cm^2/angstrom")
        @test spectral_axis(redshift(spec, 0.5)) ≈ axis .* 1.5
    end

    @testset "Unsupported axis dimension rejected" begin
        # Constructed directly, bypassing the dimension guard in `spectrum`
        axis = collect(range(1.0, 5.0, length = 10)) * u"s"
        spec = SpectrumBase.Spectrum(axis, ones(10) * u"erg/s/cm^2/angstrom", Dict{Symbol, Any}())
        @test_throws ArgumentError redshift(spec, 0.5)
    end

    @testset "Unphysical redshift rejected" begin
        spec = mock_spectrum_redshift()
        @test_throws DomainError redshift(spec, -1.0) # Exactly z = -1
        @test_throws DomainError redshift(spec, -1.5)
        @test_throws DomainError redshift!(spec, -1.0)
        # Blueshifts above the bound remain valid
        @test redshift(spec, -0.5) isa typeof(spec)
    end

    @testset "Binned spectra" begin
        edges = [1.0 2.0; 2.0 3.0; 3.0 4.0]
        flux = [10.0, 20.0, 30.0]
        binned = spectrum(edges, flux)

        shifted = @inferred redshift(binned, 0.5)
        @test shifted isa BinnedSpectrum
        @test spectral_axis(shifted) ≈ edges .* 1.5
        @test flux_axis(shifted) == flux
        @test spectral_axis(binned) == edges

        @test redshift!(binned, 0.5) === binned
        @test spectral_axis(binned) ≈ edges .* 1.5

        # Energy axis compresses
        kev = spectrum(edges * u"keV", flux * u"Jy")
        @test spectral_axis(@inferred redshift(kev, 0.5)) ≈ edges * u"keV" ./ 1.5
        @test spectral_axis(doppler_shift(kev, 1000.0)) ≈ edges * u"keV" ./ (1 + 1000 / C_KMPS)
        @test flux_axis(redshift(kev, 0.5)) == flux_axis(kev)
    end
end

@testset "doppler_shift" begin
    @testset "Non-relativistic Doppler" begin
        spec = mock_spectrum_redshift()
        original_wave = copy(spectral_axis(spec))
        original_flux = copy(flux_axis(spec))
        original_meta = deepcopy(spec.meta)

        v = 100.0  # 100 km/s
        shifted = @inferred doppler_shift(spec, v)
        expected = original_wave .* (1 + v / C_KMPS)

        @test spectral_axis(shifted) ≈ expected
        @test flux_axis(shifted) ≈ original_flux
        @test shifted.meta == original_meta

        # Original unchanged
        @test spectral_axis(spec) ≈ original_wave
        @test flux_axis(spec) ≈ original_flux
    end

    @testset "Relativistic Doppler" begin
        spec = mock_spectrum_redshift()
        original_wave = copy(spectral_axis(spec))

        v = 100.0  # 100 km/s
        β = v / C_KMPS
        shifted = @inferred doppler_shift(spec, v; relativistic = true)
        expected = original_wave .* sqrt((1 + β) / (1 - β))

        @test spectral_axis(shifted) ≈ expected
    end

    @testset "Non-relativistic ≈ relativistic at low v" begin
        spec = mock_spectrum_redshift()

        v = 0.1  # 0.1 km/s — very non-relativistic
        nr = doppler_shift(spec, v; relativistic = false)
        rel = doppler_shift(spec, v; relativistic = true)
        @test spectral_axis(nr) ≈ spectral_axis(rel) rtol = 1e-10
    end

    @testset "Doppler in-place" begin
        spec = mock_spectrum_redshift()
        v = 50.0  # 50 km/s
        expected_wave = spectral_axis(spec) .* (1 + v / C_KMPS)

        result = @inferred doppler_shift!(spec, v)
        @test result === spec
        @test spectral_axis(spec) ≈ expected_wave
    end

    @testset "Doppler with Unitful velocity" begin
        spec = mock_spectrum_redshift()
        original_wave = copy(spectral_axis(spec))

        v = 100.0u"km/s"
        shifted = @inferred doppler_shift(spec, v)
        β = 1e5 / C_MPS  # 100 km/s = 1e5 m/s
        expected = original_wave .* (1 + β)
        @test spectral_axis(shifted) ≈ expected

        shifted_rel = @inferred doppler_shift(spec, v; relativistic = true)
        expected_rel = original_wave .* sqrt((1 + β) / (1 - β))
        @test spectral_axis(shifted_rel) ≈ expected_rel
    end

    @testset "Doppler with Unitful spectrum and velocity" begin
        spec = mock_spectrum_redshift(use_units = true)
        original_wave = copy(spectral_axis(spec))

        v = 100.0u"km/s"
        shifted = @inferred doppler_shift(spec, v)
        β = 1e5 / C_MPS
        @test spectral_axis(shifted) ≈ original_wave .* (1 + β)
        @test flux_axis(shifted) ≈ flux_axis(spec)
    end

    @testset "Zero velocity" begin
        spec = mock_spectrum_redshift()
        shifted = doppler_shift(spec, 0.0)
        @test spectral_axis(shifted) ≈ spectral_axis(spec)
        @test flux_axis(shifted) ≈ flux_axis(spec)
    end

    @testset "Superluminal velocity rejected" begin
        spec = mock_spectrum_redshift()
        @test_throws DomainError doppler_shift(spec, 4e8; relativistic = true)
        @test_throws DomainError doppler_shift(spec, C_KMPS; relativistic = true)  # exactly v = c
        @test_throws DomainError doppler_shift!(spec, 4e8; relativistic = true)
        # Non-relativistic optical convention has no such restriction
        @test doppler_shift(spec, 4e8) isa typeof(spec)
    end
end
