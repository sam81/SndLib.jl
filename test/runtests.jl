t1 = time_ns() #tic()
if ispath("wavDir") == false
    mkdir("wavDir")
end

include("snd_process_tests.jl")
include("test_AMTone.jl")
include("test_AMToneIPD.jl")
include("test_asynchChord.jl")
include("test_complexTone.jl")
include("test_expAMNoise.jl")
include("test_expSinFMTone.jl")
include("test_FMComplexTone.jl")
include("test_dichotic_pitch.jl")
include("test_FMTone.jl")
include("test_IRN.jl")
include("test_noise.jl")
include("test_pureTone.jl")
include("test_pureToneILD.jl")
include("test_pureToneIPD.jl")
include("test_pureToneITD.jl")
include("test_pureToneIPDILD.jl")
include("test_silence.jl")

include("run_default.jl")
include("run_examples.jl")
include("run_getSpectrum.jl")
t2 = time_ns() #toc()
println((float(t2-t1)*1e-9))


