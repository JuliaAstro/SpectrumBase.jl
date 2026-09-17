using SpectrumBase: redshift

const H_JS = 6.62607015e-34u"J*s" # Planck constant, exact by SI definition

@testset "Wavelength and frequency conversion" begin
    wavelength = [1.0, 1.5, 2.0]u"μm"
    flux = [1.0, 2.0, 3.0]u"Jy"
    spec = spectrum(wavelength, flux, name = "example")
    c = C_MPS * u"m/s"

    frequency_spec = @inferred uconvert(u"THz", spec)
    expected_frequency = uconvert.(u"THz", c ./ wavelength)

    # The conversion is elementwise, so an ascending wavelength axis becomes descending
    @test issorted(spectral_axis(frequency_spec); rev = true)
    @test spectral_axis(frequency_spec) ≈ expected_frequency
    @test flux_axis(frequency_spec) == flux
    @test frequency_spec.name == "example"
    @test spectral_axis(spec) == wavelength

    # Round trip restores the original axis and ordering
    wavelength_spec = @inferred uconvert(u"μm", frequency_spec)
    @test spectral_axis(wavelength_spec) ≈ wavelength
    @test flux_axis(wavelength_spec) == flux

    shifted = redshift(frequency_spec, 1.0)
    @test spectral_axis(shifted) ≈ expected_frequency ./ 2

    # Same-dimension targets are plain unit conversions, so conversions are idempotent
    @test spectral_axis(uconvert(u"GHz", frequency_spec)) ≈
        uconvert.(u"GHz", spectral_axis(frequency_spec))
    @test spectral_axis(uconvert(u"nm", spec)) ≈ uconvert.(u"nm", wavelength)

    # A `unit(spec)` tuple converts one spectrum to another's units
    matched = uconvert(unit(frequency_spec), spec)
    @test spectral_axis(matched) == spectral_axis(frequency_spec)
    @test flux_axis(matched) == flux_axis(frequency_spec)

    @test_throws MethodError uconvert(1u"Hz", spec) # a quantity, not a unit
    @test_throws ArgumentError uconvert(u"kg", spec)
    @test_throws ArgumentError uconvert(u"THz", spectrum([1.0, 2.0], [3.0, 4.0]))
end

@testset "Energy axis conversion" begin
    energies = [1.0, 2.0, 5.0]u"keV"
    flux = [1.0, 2.0, 3.0]u"Jy"
    xray = spectrum(energies, flux)
    c = C_MPS * u"m/s"

    λspec = @inferred uconvert(u"angstrom", xray)
    @test spectral_axis(λspec) ≈ uconvert.(u"angstrom", H_JS * c ./ energies)
    @test flux_axis(λspec) == flux

    # E → ν is linear, so the axis direction is preserved
    νspec = @inferred uconvert(u"Hz", xray)
    @test spectral_axis(νspec) ≈ uconvert.(u"Hz", energies ./ H_JS)
    @test flux_axis(νspec) == flux

    back = @inferred uconvert(u"keV", λspec)
    @test spectral_axis(back) ≈ energies
    @test flux_axis(back) == flux
end

@testset "Flux density conversion" begin
    wavelength = [1.0, 1.5, 2.0]u"μm"
    F_λ = [1.0, 2.0, 3.0]u"W/m^2/μm"
    spec = spectrum(wavelength, F_λ)
    c = C_MPS * u"m/s"

    converted = @inferred uconvert((u"THz", u"Jy"), spec)
    expected_F_ν = uconvert.(u"Jy", F_λ .* wavelength .^ 2 ./ c)
    @test spectral_axis(converted) ≈ uconvert.(u"THz", c ./ wavelength)
    @test flux_axis(converted) ≈ expected_F_ν

    # Round trip back to F_λ
    back = uconvert((u"μm", u"W/m^2/μm"), converted)
    @test spectral_axis(back) ≈ wavelength
    @test flux_axis(back) ≈ F_λ

    # Per-energy density, the X-ray convention
    F_E_spec = @inferred uconvert((u"keV", u"W/m^2/keV"), spec)
    @test flux_axis(F_E_spec) ≈ uconvert.(u"W/m^2/keV", expected_F_ν ./ H_JS)

    # A same-dimension flux target is a plain unit conversion
    same = uconvert((u"THz", u"W/m^2/μm"), spec)
    @test flux_axis(same) == F_λ
end

@testset "Binned spectrum axis conversion" begin
    edges = [1.0 2.0; 2.0 3.0; 3.0 4.0]u"μm"
    flux = [10.0, 20.0, 30.0]u"Jy"
    binned = spectrum(edges, flux)
    c = C_MPS * u"m/s"

    converted = uconvert(u"THz", binned)
    @test spectral_axis(converted) ≈ uconvert.(u"THz", c ./ edges)
    @test flux_axis(converted) == flux
    @test all(col -> issorted(col; rev = true), eachcol(spectral_axis(converted)))

    # Bin values are integrated quantities; flux density conversion is undefined for them
    @test_throws ArgumentError uconvert((u"THz", u"W/m^2/μm"), binned)
end
