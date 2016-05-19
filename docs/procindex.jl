using Weave

tangle("index.md", out_path=:pwd, informat="markdown")
include("index.jl")
