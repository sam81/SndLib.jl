using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
phase = 0
level = 60
dur = 1
rampDur = 0.01
maxLevel = 101
channelOpts = ["right", "left", "diotic", "mono", "odd left", "odd right"]
F0 = 440
harmPhaseOpts = ["sine", "cosine", "alternating", "random", "schroeder-", "schroeder+"]
lowHarm = 1
highHarm = 10
#stretch = 0
MF = 5
FMDepth=100
FMStartPhase = 0
FMStartTime = 0
FMDur = 1
levelAdj = true


for channel in channelOpts
    for harmPhase in harmPhaseOpts
        snd = FMComplex2(midF0=F0, harmPhase=harmPhase,
                 lowHarm=lowHarm, highHarm=highHarm,
                 level=level, dur=dur, rampDur=rampDur,
                 MF=5, FMDepth=FMDepth, FMStartPhase=FMStartPhase,
                 FMStartTime=FMStartTime, FMDur=FMDur,
                 levelAdj=levelAdj, channel=channel,
                 sf=sf, maxLevel=maxLevel)
        wavwrite(snd, wavDir*"FMComplex_"*harmPhase*"_"*channel*".wav",
                 Fs=sf, nbits=nbits)
    end
end


## Test error condition
@test_throws(ErrorException, FMComplex2(dur=1, rampDur=0.6))

@test_throws(ErrorException, FMComplex2(dur=1, FMDur=1, FMStartTime=0.2))

@test_throws(ErrorException, FMComplex2(channel="foo"))

@test_throws(ErrorException, FMComplex2(harmPhase="foo"))
