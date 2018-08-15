using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
frequency = 1000
phase = 0
level = 60
dur = 1
rampDur = 0.01
maxLevel = 101
channelOpts = ["mono", "right", "left", "diotic"]

for channel in channelOpts

    pt = pureTone(frequency=frequency, phase=phase, level=level, dur=dur,
                  rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(pt, wavDir*"pt_"*channel*".wav", Fs=sf, nbits=nbits)
end

@test_throws(ErrorException, pureTone(dur=1, rampDur=0.6))

@test_throws(ErrorException, pureTone(channel="foo"))
