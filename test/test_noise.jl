using SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
phase = 0
level = 60
dur = 1
rampDur = 0.01
channelOpts = ["mono", "right", "left", "diotic", "dichotic"]
maxLevel = 101

spectrumLevel = 20
for channel in channelOpts
    bn = broadbandNoise(spectrumLevel=spectrumLevel, dur=dur,
                        rampDur=rampDur, channel=channel, sf=sf,
                        maxLevel=maxLevel)
    wavwrite(bn, wavDir*"white_noise_"*channel*".wav", Fs=sf, nbits=nbits)
end

bn = broadbandNoise(spectrumLevel=spectrumLevel, dur=dur,
                          rampDur=0, channel="diotic", sf=sf,
                          maxLevel=maxLevel)
bn_gated = gate!(bn, rampDur=0.01, sf=sf)
wavwrite(bn_gated, wavDir*"white_noise_gated.wav", Fs=sf, nbits=nbits)

f1 = 900
f2 = 1200
level = 20
for channel in channelOpts
    stn = steepNoise(f1=f1, f2=f2, level=level, dur=dur, rampDur=rampDur,
                     channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(stn, wavDir*"steep_noise_"*channel*".wav", Fs=sf, nbits=nbits)
end



## Test error condition
@test_throws(ErrorException, broadbandNoise(dur=1, rampDur=0.6))

@test_throws(ErrorException, broadbandNoise(channel="foo"))

@test_throws(ErrorException, steepNoise(dur=1, rampDur=0.6))

@test_throws(ErrorException, steepNoise(channel="foo"))

