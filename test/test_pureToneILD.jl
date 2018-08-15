using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
frequency = 1000
phase = 0
level = 60
ILD = 15
dur = 1
rampDur = 0.01
maxLevel = 101
channelOpts = ["right", "left"]

for channel in channelOpts

    pt = pureToneILD(frequency=frequency, phase=phase, level=level, ILD=ILD, dur=dur,
                  rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(pt, wavDir*"pt_ILD_"*channel*".wav", Fs=sf, nbits=nbits)
end

@test_throws(ErrorException, pureToneILD(dur=1, rampDur=0.6))

@test_throws(ErrorException, pureToneILD(channel="foo"))
