using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32

carrierFreq=200
MF=5
deltaCents=600
FMPhase=pi
phase = 0
level = 50
dur=0.2
rampDur=0.01
maxLevel = 101

channelOpts = ["mono", "right", "left", "diotic"]

for channel in channelOpts

    tn = expSinFMTone(carrierFreq=carrierFreq, MF=MF, deltaCents=deltaCents, FMPhase=FMPhase, phase=phase, level=level, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(tn, wavDir*"expSinFMTone_"*channel*".wav", Fs=sf, nbits=nbits)
end

@test_throws(ErrorException, expSinFMTone(dur=1, rampDur=0.6))

@test_throws(ErrorException, expSinFMTone(channel="foo"))
