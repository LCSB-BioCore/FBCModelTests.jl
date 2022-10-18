
# tried and trusted E. coli core
isfile("e_coli_core.json") || Downloads.download(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
)
model = load_model("e_coli_core.json")

# this model has an energy generating cycle
isfile("iJN746.json") ||
    Downloads.download("http://bigg.ucsd.edu/static/models/iJN746.json", "iJN746.json")
iJN746 = load_model("iJN746.json")

