using Documenter, DocumenterMarkdown, SndLib
#push!(LOAD_PATH,"src/")
if ispath("build/")
    rm("build/", recursive=true)
end

if ispath("site/")
    rm("site/", recursive=true)
end

include("runWeave.jl") #convert weave file to markdown
makedocs(format = Markdown(), modules = [SndLib], source="src", build="build")
run(`mkdocs build`)
