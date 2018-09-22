using SndLib

fToProcess = ["snd_generate.jl", "snd_process.jl", "utils.jl"]
for fName in fToProcess
##fName = fToProcess[1]
    fIn = open(string("../src/", fName), "r")
    fOut = open(string("test_", fName), "w")
    lns = readlines(fIn)
    idxStart = (Int)[]
    idxStop = (Int)[]

    for i=1:length(lns)
        if lns[i] == "```julia"
            push!(idxStart, i+1)
        elseif lns[i] == "```"
            push!(idxStop, i-1)
        end
    end

    for i=1:length(idxStart)
        for j=idxStart[i]:idxStop[i]
            write(fOut, lns[j]*"\n")
        end
    end
    
    close(fIn); close(fOut)
end

for fName in fToProcess
    fNameTest = string("test_", fName)
    include(fNameTest)
end

