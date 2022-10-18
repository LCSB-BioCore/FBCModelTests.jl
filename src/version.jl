const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))

include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const FBCMT_VERSION =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])
