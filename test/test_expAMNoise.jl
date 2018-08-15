using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32

carrierFreq=200
MF=5
deltaCents=600
FMPhase=pi
AMDepth = 1
spectrumLevel=30
dur=0.4
rampDur=0.01
maxLevel = 101
channelOpts = ["mono", "right", "left", "diotic", "dichotic"]

for channel in channelOpts
    ns = expAMNoise(carrierFreq=carrierFreq, MF=MF, deltaCents=deltaCents, FMPhase=FMPhase, AMDepth=AMDepth, spectrumLevel=spectrumLevel, dur=dur, rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)

    wavwrite(ns, wavDir*"expAMNoise_"*channel*".wav", Fs=sf, nbits=nbits)
end

@test_throws(ErrorException, expAMNoise(dur=1, rampDur=0.6))

@test_throws(ErrorException, expAMNoise(channel="foo"))
