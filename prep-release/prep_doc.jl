using SndLib, Lexicon


Lexicon.save("../docs/API.md", SndLib)
cd("../")
run(`mkdocs build`)
cd("prep-release")
