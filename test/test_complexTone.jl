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
stretch = 0

for channel in channelOpts
    for harmPhase in harmPhaseOpts
        ct = complexTone(F0=F0, harmPhase=harmPhase, lowHarm=lowHarm,
                         highHarm=highHarm, stretch=stretch, level=level,
                         dur=dur, rampDur=rampDur, channel=channel, sf=sf,
                         maxLevel=maxLevel)
        wavwrite(ct, wavDir*"ct_"*harmPhase*"_"*channel*".wav",
                 Fs=sf, nbits=nbits)
    end
end


## Test error condition
@test_throws(ErrorException, complexTone(dur=1, rampDur=0.6))

@test_throws(ErrorException, complexTone(channel="foo"))

@test_throws(ErrorException, complexTone(harmPhase="foo"))
