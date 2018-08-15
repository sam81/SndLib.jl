 
using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
maxLevel = 101

frequency = 1000
AMFreq = 2
AMDepth = 1
AMPhase = 1.5*pi
phase = 0
level = 65
dur = 1
rampDur = 0.01

channelOpts = ["mono", "right", "left", "diotic"]

for channel in channelOpts

    snd = AMTone(carrierFreq=frequency, AMFreq=AMFreq, AMDepth=AMDepth, carrierPhase=phase, AMPhase=AMPhase, level=level, dur=dur, rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(snd, wavDir*"AMTone_"*channel*".wav", Fs=sf, nbits=nbits)
end


## Test error condition
@test_throws(ErrorException, AMTone(dur=1, rampDur=0.6))

@test_throws(ErrorException, AMTone(channel="foo"))
