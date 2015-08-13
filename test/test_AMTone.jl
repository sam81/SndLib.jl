 
using SndLib, WAV

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

channelOpts = ["right", "left", "diotic", "mono"]

for channel in channelOpts

    snd = AMTone(carrierFreq=frequency, AMFreq=AMFreq, AMDepth=AMDepth, carrierPhase=phase, AMPhase=AMPhase, level=level, dur=duration, rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)
    wavwrite(snd, wavDir*"AMTone_"*channel*".wav", Fs=sf, nbits=nbits)
end


