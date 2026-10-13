include("redden.jl")
include("redshift.jl")

# ## Broadening ops
# include("kernels.jl")
# using Distributions, DSP

# function _broaden(flux, kernel::Kernel)
#     k = evaluate(kernel)
#     broad_flux = conv(flux, k)
# end

# function broaden(spec::Spectrum, kernel::Kernel)
#     w_unit, f_unit = unit(spec)
#     s = ustrip(spec)
#     broad_flux = _broaden(s.flux, kernel)
#     return Spectrum(spec.wave, broad_flux * f_unit, name = spec.name)
# end
