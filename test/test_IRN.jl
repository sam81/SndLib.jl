using Base.Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
delay = 1/300
gain = 0.8
iterations = 10
spectrumLevel = 25
dur = 1
rampDur = 0.01
maxLevel = 101
configurationOpts = ["add same", "add original"]
channelOpts = ["mono", "right", "left", "diotic"]

for configuration in configurationOpts

    for channel in channelOpts

        tn = IRN(delay=delay, gain=gain, iterations=iterations, configuration=configuration,
                 spectrumLevel=spectrumLevel, dur=dur, rampDur=rampDur, channel=channel,
                 sf=sf, maxLevel=maxLevel)
        wavwrite(tn, wavDir*"IRN_"*channel*"_"*configuration*".wav", Fs=sf, nbits=nbits)
        
    end
end

@test_throws(ErrorException, IRN(dur=1, rampDur=0.6))

@test_throws(ErrorException, IRN(channel="foo"))

@test_throws(ErrorException, IRN(configuration="foo"))
