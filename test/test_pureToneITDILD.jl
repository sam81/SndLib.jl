using Base.Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
frequency = 1000
phase = 0
level = 60
ITD = 0.005/10
ILD = 10
dur = 1
rampDur = 0.01
maxLevel = 101
channelOpts = ["right", "left"]

for channel in channelOpts

    pt = pureToneITDILD(frequency=frequency, phase=phase, ITD=ITD, level=level, ILD=ILD, dur=dur,
                  rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(pt, wavDir*"pt_ITD_ILD_"*channel*".wav", Fs=sf, nbits=nbits)
end

@test_throws(ErrorException, pureToneITDILD(dur=1, rampDur=0.6))

@test_throws(ErrorException, pureToneITDILD(channel="foo"))
