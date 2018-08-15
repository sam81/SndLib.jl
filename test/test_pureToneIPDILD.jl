using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
frequency = 1000
phase = 0
level = 60
IPD = pi
ILD = 10
dur = 1
rampDur = 0.01
maxLevel = 101
channelOpts = ["right", "left"]

for channel in channelOpts

    pt = pureToneIPDILD(frequency=frequency, phase=phase, IPD=IPD, level=level, ILD=ILD, dur=dur,
                  rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(pt, wavDir*"pt_IPD_ILD"*channel*".wav", Fs=sf, nbits=nbits)
end

@test_throws(ErrorException, pureToneIPDILD(dur=1, rampDur=0.6))

@test_throws(ErrorException, pureToneIPDILD(channel="foo"))
