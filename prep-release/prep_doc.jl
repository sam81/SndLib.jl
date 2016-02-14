using SndLib, Lexicon
include("extract_docstrings.jl")

extract_docstrings(["../src/SndLib.jl",
                    "../src/snd_generate.jl",
                    "../src/snd_process.jl",
                    "../src/utils.jl"],
                   "../docs/API.md")

#Lexicon.save("../docs/API.md", SndLib)
cd("../")
run(`mkdocs build`)
cd("prep-release")
