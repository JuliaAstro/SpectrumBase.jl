using Documenter: DocMeta, doctest
import DynamicQuantities

# Docstrings defined in the units extensions belong to those modules, so doctest them too.
const modules = [
    SpectrumBase,
    Base.get_extension(SpectrumBase, :SpectrumBaseUnitfulExt),
    Base.get_extension(SpectrumBase, :SpectrumBaseDynamicQuantitiesExt),
]

for m in modules
    DocMeta.setdocmeta!(m, :DocTestSetup, :(using SpectrumBase); recursive = true)
end

doctest(joinpath(pkgdir(SpectrumBase), "docs", "src"), modules)
