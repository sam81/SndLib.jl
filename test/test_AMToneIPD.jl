 
using Base.Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
maxLevel = 101

frequency = 1000
AMFreq = 100
AMDepth = 1
carrierPhase = 0
AMPhase = 1.5*pi
carrierIPD = 0
AMIPD = pi/2
level = 65
dur = 1
rampDur = 0.01

channelOpts = ["right", "left"]

for channel in channelOpts

    snd = AMToneIPD(carrierFreq=frequency, AMFreq=AMFreq, AMDepth=AMDepth, carrierPhase=carrierPhase, AMPhase=AMPhase, carrierIPD=carrierIPD, AMIPD=AMIPD, level=level, dur=dur, rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(snd, wavDir*"AMToneIPD_"*channel*".wav", Fs=sf, nbits=nbits)
end


## Test error condition
@test_throws(ErrorException, AMToneIPD(dur=1, rampDur=0.6))

@test_throws(ErrorException, AMToneIPD(channel="foo"))
