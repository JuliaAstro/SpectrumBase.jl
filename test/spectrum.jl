using SpectrumBase: SingleSpectrum, BinnedSpectrum, EchelleSpectrum, IFUSpectrum

const RNG = Random.seed!(8675309)

@testset "Spectrum types" begin
    wave, flux = collect(1:5), 1:5
    @test spectrum(wave, flux) isa SingleSpectrum
    wave[3] = 10
    @test_throws ArgumentError spectrum(wave, flux)
    @test_throws ArgumentError spectrum(wave[begin:end-1], flux)

    wave, flux = repeat(1:5, 1, 3), rand(RNG, 5, 3)
    @test spectrum(wave, flux) isa EchelleSpectrum
    wave[3, 3] = 10
    @test_throws ArgumentError spectrum(wave, flux)
    @test_throws ArgumentError spectrum(wave[begin:end-1, :], flux)

    wave, flux = collect(1:5), rand(RNG, 5, 4, 3)
    @test spectrum(wave, flux) isa IFUSpectrum
    wave[3] = 10
    @test_throws ArgumentError spectrum(wave, flux)
    @test_throws ArgumentError spectrum(wave[begin:end-1], flux)

    wave, flux = [1 2; 2 3; 3 4; 4 5; 5 6], 1:5
    @test spectrum(wave, flux) isa BinnedSpectrum
end

@testset "Spectrum - Single" begin
    wave = range(1e4, 5e4, length = 1000)
    sigma = randn(RNG, size(wave))
    sigma[7] = 1
    sigma[134] = 0.1
    flux = 100 .± sigma
    flux[7] = 1000 ± 1
    flux[134] = 1 ± 0.1

    spec = spectrum(wave, flux, name = "test spectrum")
    spec_indexed = spec[begin:end]

    @test propertynames(spec) == (:spectral_axis, :flux_axis, :meta, :name)
    @test [s for s in spec] isa Vector{<:SingleSpectrum{Float64, Measurements.Measurement{Float64}}}
    @test eltype(spec) == eltype(flux_axis(spec))
    @test spectral_axis(spec) == wave
    @test spectral_axis(spec_indexed) == wave
    @test size(spec) === (1000,)
    @test length(spec) == 1000
    @test maximum(spec) == 1000 ± 1
    @test minimum(spec) == 1 ± 0.1
    @test argmax(spec) == 7
    @test argmin(spec) == 134
    @test findmax(spec) == (1000 ± 1, 7)
    @test findmin(spec) == (1 ± 0.1, 134)
    @test flux_axis(spec) == flux
    @test flux_axis(spec_indexed) == flux
    @test Measurements.uncertainty.(flux_axis(spec)) ≈ sigma

    flux_trimmed = flux[200:800]
    @test_throws ArgumentError spectrum(wave, flux_trimmed)
    expected = """
    SingleSpectrum(Float64, Measurements.Measurement{Float64})
      spectral axis (1000,): 10000.0 .. 50000.0
      flux axis (1000,): 100.0 ± -2.8 .. 100.0 ± 0.6
      meta: Dict{Symbol, Any}(:name => "test spectrum")"""
    @test sprint(show, spec) == expected
    @test spec.name == "test spectrum"
end

@testset "Spectrum - Binned" begin
    wave = [1.0 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0]
    flux = [10.0, 20.0, 30.0, 40.0]

    spec = spectrum(wave, flux, name = "test binned spectrum")
    spec_i = spec[2]
    spec_I = spec[2:3]

    @test spec isa BinnedSpectrum
    @test spec_i isa BinnedSpectrum
    @test spec_I isa BinnedSpectrum
    @test spectral_axis(spec_i) == [2.0 3.0]
    @test flux_axis(spec_i) == [20.0]
    @test spectral_axis(spec_I) == [2.0 3.0; 3.0 4.0]
    @test flux_axis(spec_I) == [20.0, 30.0]
    @test firstindex(spec) == 1
    @test lastindex(spec) == 4
    @test spec[firstindex(spec)] == spec[1]
    @test spec[lastindex(spec)] == spec[4]

    # Any integer index selects a single bin
    @test spec[Int32(2)] == spec[2] == spec[2:2]

    # m × 2 shape and row-order invariants
    @test_throws ArgumentError spectrum([1.0 2.0 3.0; 2.0 3.0 4.0], [1.0, 2.0]) # not two columns
    @test_throws ArgumentError spectrum(reshape([1.0, 2.0], 2, 1), [1.0, 2.0]) # not two columns
    @test_throws ArgumentError spectrum([2.0 1.0; 3.0 2.0], [1.0, 2.0]) # edges swapped within rows
    @test_throws ArgumentError spectrum([1.0 2.0; 3.0 2.5], [1.0, 2.0]) # inconsistent row order
    @test_throws ArgumentError spectrum([1.0 1.0; 2.0 2.0], [1.0, 2.0]) # zero-width bins
    @test spectrum([4.0 3.0; 3.0 2.0], [1.0, 2.0]) isa BinnedSpectrum # descending edges are valid

    # Empty selections construct and display
    empty_spec = spec[falses(4)]
    @test length(empty_spec) == 0
    @test occursin("empty", sprint(show, empty_spec))

    expected = """
    BinnedSpectrum(Float64, Float64)
      spectral axis (4, 2): 1.0 .. 5.0
      flux axis (4,): 10.0 .. 40.0
      meta: Dict{Symbol, Any}(:name => "test binned spectrum")"""
    @test sprint(show, spec) == expected
    @test spec.name == "test binned spectrum"
end

@testset "Spectrum - Echelle" begin
    n_wavs = 1000
    n_orders = 3
    wave_1 = range(1e4, 5e4, length=n_wavs)
    wave = repeat(wave_1, 1, n_orders)
    sigma = randn(RNG, size(wave_1))
    sigma[7] = 1
    sigma[134] = 0.1
    flux_1 = 100 .± sigma
    flux_1[7] = 1000 ± 1
    flux_1[134] = 1 ± 0.1
    flux = repeat(flux_1, 1, n_orders)

    spec = spectrum(wave, flux, name = "Test Echelle Spectrum")

    i = 1
    I = 1:3
    spec_i = spec[i]
    spec_i_expected = spectrum(wave_1, flux_1, name=spec.name)
    spec_I = spec[I]
    spec_I_expected = spectrum(wave, flux, name=spec.name)

    @test (spec_i.name, spectral_axis(spec_i), flux_axis(spec_i)) == (spec_i_expected.name, spectral_axis(spec_i_expected), flux_axis(spec_i_expected))
    @test (spec_I.name, spectral_axis(spec_I), flux_axis(spec_I)) == (spec_I_expected.name, spectral_axis(spec_I_expected), flux_axis(spec_I_expected))
    @test propertynames(spec) == (:spectral_axis, :flux_axis, :meta, :name)
    @test propertynames(spec_i) == (:spectral_axis, :flux_axis, :meta, :Order, :name)
    @test propertynames(spec_I) == (:spectral_axis, :flux_axis, :meta, :name, :Orders)
    @test eltype(spec) == eltype(flux_axis(spec))
    @test spectral_axis(spec) == wave
    @test size(spec) == (n_wavs, n_orders)
    @test length(spec) == n_wavs * n_orders
    @test maximum(spec) == 1000 ± 1
    @test minimum(spec) == 1 ± 0.1
    @test argmax(spec) == CartesianIndex(7, 1)
    @test argmin(spec) == CartesianIndex(134, 1)
    @test findmax(spec) == (1000 ± 1, CartesianIndex(7, 1))
    @test findmin(spec) == (1 ± 0.1, CartesianIndex(134, 1))
    @test eachcol(Measurements.uncertainty.(flux_axis(spec))) ≈ fill(sigma, n_orders)

    flux_trimmed = flux[200:800, :]
    @test_throws ArgumentError spectrum(wave, flux_trimmed)
    expected = """
    EchelleSpectrum(Float64, Measurements.Measurement{Float64})
      # orders: 3
      spectral axis (1000, 3): 10000.0 .. 50000.0
      flux axis (1000, 3): 100.0 ± -2.8 .. 100.0 ± 0.6
      meta: Dict{Symbol, Any}(:name => "Test Echelle Spectrum")"""
    @test sprint(show, spec) == expected
    @test spec.name == "Test Echelle Spectrum"
end

@testset "Spectrum - IFU" begin
    wave, flux = [20, 40, 120, 160, 200], rand(RNG, 5, 10, 6)

    spec = spectrum(wave, flux, name = "test spectrum")

    expected = """
    IFUSpectrum(Int64, Float64)
      spectral axis (5,): 20 .. 200
      flux axis (5, 10, 6): 0.9210599764489846 .. 0.47778429984485815
      meta: Dict{Symbol, Any}(:name => "test spectrum")"""

    @test sprint(show, spec) == expected
    @test spec.name == "test spectrum"
    @test propertynames(spec) == (:spectral_axis, :flux_axis, :meta, :name)
    @test eltype(spec) == eltype(flux_axis(spec))
    @test spectral_axis(spec) == wave
    @test size(spec) === (5, 10, 6)
    @test length(spec) == 300
    @test flux_axis(spec) == flux
    @test spec[:, 1, 1] isa SingleSpectrum
    @test spec[:, begin:4, begin:3] isa IFUSpectrum
end

@testset "Unitful Spectrum - Single" begin
    wave = range(1e4, 5e4, length = 1000)
    sigma = randn(RNG, size(wave))
    sigma[7] = 1
    sigma[134] = 0.1
    flux = 100 .± sigma
    flux[7] = 1000 ± 1
    flux[134] = 1 ± 0.1

    funit = u"W/m^2/angstrom"
    spec = spectrum(wave * u"angstrom", flux * funit, name = "test")

    @test spectral_axis(spec) ≈ wave * u"angstrom"

    @test propertynames(spec) == (:spectral_axis, :flux_axis, :meta, :name)
    @test eltype(spec) == eltype(flux_axis(spec))
    @test size(spec) === (1000,)
    @test length(spec) == 1000
    @test maximum(spec) == (1000 ± 1) * funit
    @test minimum(spec) == (1 ± 0.1) * funit
    @test argmax(spec) == 7
    @test argmin(spec) == 134
    @test findmax(spec) == ((1000 ± 1) * funit, 7)
    @test findmin(spec) == ((1 ± 0.1) * funit, 134)
    @test spec.name == "test"

    # Test stripping
    w_unit, f_unit = unit(spec)
    @test w_unit == u"angstrom"
    @test f_unit == u"W/m^2/angstrom"

    strip_spec = ustrip(spec)
    @test spectral_axis(strip_spec) == ustrip.(spectral_axis(spec))
    @test flux_axis(strip_spec) == ustrip.(flux_axis(spec))
    @test strip_spec.meta == spec.meta
    expected = """
    SingleSpectrum(Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}, Unitful.Quantity{Measurements.Measurement{Float64}, 𝐌 𝐋^-1 𝐓^-3, Unitful.FreeUnits{(Å^-1, m^-2, W), 𝐌 𝐋^-1 𝐓^-3, nothing}})
      spectral axis (1000,): 10000.0 Å .. 50000.0 Å
      flux axis (1000,): 100.0 ± -2.8 W Å^-1 m^-2 .. 100.0 ± 0.6 W Å^-1 m^-2
      meta: Dict{Symbol, Any}(:name => "test")"""
    @test sprint(show, spec) == expected
end

@testset "Unitful Spectrum - Echelle" begin
    n_wavs = 1000
    n_orders = 3
    wave_1 = range(1e4, 5e4, length=n_wavs)
    wave = repeat(wave_1, 1, n_orders)
    sigma = randn(RNG, size(wave_1))
    sigma[7] = 1
    sigma[134] = 0.1
    flux_1 = 100 .± sigma
    flux_1[7] = 1000 ± 1
    flux_1[134] = 1 ± 0.1
    flux = repeat(flux_1, 1, n_orders)

    wunit = u"angstrom"
    funit = u"W/m^2/angstrom"

    wave *= wunit
    flux *= funit

    spec = spectrum(wave, flux, name = "test echelle")

    @test spectral_axis(spec) ≈ wave

    @test propertynames(spec) == (:spectral_axis, :flux_axis, :meta, :name)
    @test eltype(spec) == eltype(flux_axis(spec))
    @test size(spec) === (n_wavs, n_orders)
    @test length(spec) == n_wavs * n_orders
    @test maximum(spec) == (1000 ± 1) * funit
    @test minimum(spec) == (1 ± 0.1) * funit
    @test argmax(spec) == CartesianIndex(7, 1)
    @test argmin(spec) == CartesianIndex(134, 1)
    @test findmax(spec) == ((1000 ± 1) * funit, CartesianIndex(7, 1))
    @test findmin(spec) == ((1 ± 0.1) * funit, CartesianIndex(134, 1))
    @test spec.name == "test echelle"

    # Test stripping
    w_unit, f_unit = unit(spec)
    @test w_unit == u"angstrom"
    @test f_unit == u"W/m^2/angstrom"

    strip_spec = ustrip(spec)
    @test spectral_axis(strip_spec) == ustrip.(spectral_axis(spec))
    @test flux_axis(strip_spec) == ustrip.(flux_axis(spec))
    @test strip_spec.meta == spec.meta
    expected = """
    EchelleSpectrum(Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}, Unitful.Quantity{Measurements.Measurement{Float64}, 𝐌 𝐋^-1 𝐓^-3, Unitful.FreeUnits{(Å^-1, m^-2, W), 𝐌 𝐋^-1 𝐓^-3, nothing}})
      # orders: 3
      spectral axis (1000, 3): 10000.0 Å .. 50000.0 Å
      flux axis (1000, 3): 100.0 ± -2.8 W Å^-1 m^-2 .. 100.0 ± 0.6 W Å^-1 m^-2
      meta: Dict{Symbol, Any}(:name => "test echelle")"""
    @test sprint(show, spec) == expected
end

@testset "Arithmetic" begin
    wave = range(1e4, 5e4, length = 1000)
    sigma = randn(RNG, size(wave))
    flux = 100 .± sigma
    flux_2 = 200 .± sigma

    spec = spectrum(wave, flux, name = "test spectrum")
    spec_2 = spectrum(wave, flux_2, name = "test spectrum")

    # negation
    @test -flux_axis(spec) == -flux

    # Scalars / vectors
    values = [10, randn(RNG, size(spec))]
    for A in values
        # addition
        s = spec + A
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) .+ A

        # subtraction
        s = spec - A
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) .- A
        s = A - spec
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) .- A
        s = spec_2 - spec
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec_2) .- flux_axis(spec)

        # multiplication
        s = spec * A
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) .* A

        # division
        s = spec / A
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) ./ A
    end

    # Other spectra
    spec2 = deepcopy(spec)

    s = spec + spec2
    @test spectral_axis(s) == spectral_axis(spec)
    @test flux_axis(s) ≈ 2 .* flux_axis(spec)

    s = spec - spec2
    @test spectral_axis(s) == spectral_axis(spec)
    @test flux_axis(s) ≈ zeros(size(spec))

    s = spec * spec2
    @test spectral_axis(s) == spectral_axis(spec)
    @test flux_axis(s) ≈ flux_axis(spec).^2

    s = spec / spec2
    @test spectral_axis(s) == spectral_axis(spec)
    @test flux_axis(s) ≈ ones(size(spec))


    spec = spectrum(spectral_axis(spec) * u"cm", flux_axis(spec) * u"W/m^2/cm", name = "test unitfulspectrum")

    # Scalars / vectors
    for A in [10u"W/m^2/cm", randn(RNG, size(spec))u"W/m^2/cm"]
        # addition
        s = spec + A
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) .+ A

        # subtraction
        s = spec - A
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) .- A
    end

    for A in [10, randn(RNG, size(spec))]
        # multiplication
        s = spec * A
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) .* A

        # division
        s = spec / 10
        @test spectral_axis(s) == spectral_axis(spec)
        @test flux_axis(s) ≈ flux_axis(spec) ./ 10
    end
end

@testset "Unitful Spectrum - Binned" begin
    edges = [1.0 2.0; 2.0 3.0; 3.0 4.0]u"keV"
    flux = [1.0, 2.0, 3.0]u"Jy"

    spec = spectrum(edges, flux, name = "test xray spectrum")

    @test spec isa BinnedSpectrum
    @test unit(spec) == (u"keV", u"Jy")
    @test spectral_axis(spec) == edges
    @test flux_axis(spec) == flux
    @test spec.name == "test xray spectrum"

    @test_throws AssertionError spectrum([1.0 2.0; 2.0 3.0]u"kg", flux)
end
