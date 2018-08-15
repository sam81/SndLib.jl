using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
frequency = 1000
phase = 0
level = 60
IPD = pi
dur = 1
rampDur = 0.01
maxLevel = 101
channelOpts = ["right", "left"]

for channel in channelOpts

    pt = pureToneIPD(frequency=frequency, phase=phase, level=level, IPD=IPD, dur=dur,
                  rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(pt, wavDir*"pt_IPD_"*channel*".wav", Fs=sf, nbits=nbits)
end

@test_throws(ErrorException, pureToneIPD(dur=1, rampDur=0.6))

@test_throws(ErrorException, pureToneIPD(channel="foo"))
