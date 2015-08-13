using SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
maxLevel = 101

carrierFreq = 1000
MF = 40
MI = 1
phase = 0
level = 65
dur = 1
rampDur = 0.01

channelOpts = ["mono", "right", "left", "diotic"]

for channel in channelOpts

    snd = FMTone(carrierFreq=carrierFreq, MF=MF, MI=MI, phase=phase,
                 level=level, dur=dur, rampDur=rampDur, channel=channel,
                 sf=sf, maxLevel=maxLevel)
    wavwrite(snd, wavDir*"FMTone_"*channel*".wav", Fs=sf, nbits=nbits)
end

