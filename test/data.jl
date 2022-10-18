
datadir = joinpath(@__DIR__, "testdata")
mkpath(datadir)

model_file = Dict(
    (
        begin
            fn = joinpath(datadir, filename)
            if !isfile(fn)
                Downloads.download(url, fn)
            end

            cksum = open(fn) do io
                bytes2hex(sha256(io))
            end

            if cksum != hash
                @error "Hash mismatch for `$filename'" downloaded = cksum expected = hash
                @info "Remove test/$fn to force re-downloading"
                error("Downloaded hash mismatch for `$filename'")
            end

            filename => fn

        end
    ) for (filename, url, hash) in [
        (
            "e_coli_core.json",
            "http://bigg.ucsd.edu/static/models/e_coli_core.json",
            "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
        ),
        (
            "e_coli_core.xml",
            "http://bigg.ucsd.edu/static/models/e_coli_core.xml",
            "b4db506aeed0e434c1f5f1fdd35feda0dfe5d82badcfda0e9d1342335ab31116",
        ),
        (
            "iJN746.json",
            "http://bigg.ucsd.edu/static/models/iJN746.json",
            "3da875419f4edff3e3e27424742b8b410b6b778eb0955b8939ae933fad5f1219",
        ),
    ]
)

#TODO remove the globals
model = load_model(model_file["e_coli_core.json"])
iJN746 = load_model(model_file["iJN746.json"])
