tic()
if ispath("wavDir") == false
    mkdir("wavDir")
end
include("test_default.jl")
include("test_examples.jl")
include("snd_process_tests.jl")
include("test_AMTone.jl")
include("test_AMToneIPD.jl")
include("test_complexTone.jl")
include("test_FMComplexTone.jl")
include("test_dichotic_pitch.jl")
include("test_FMTone.jl")
include("test_noise.jl")
include("test_pureTone.jl")
include("test_silence.jl")
toc()

