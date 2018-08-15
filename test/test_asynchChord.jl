using Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32
freqs = [100, 200, 400, 800]
levels = [50, 50, 50, 50]
phases = [0, 0, 0, 0]
level = 60
tonesDur = 0.2
tonesRampDur = 0.01
tonesRampDur=0.01
SOA=0.06
maxLevel = 101
channelOpts = ["mono", "right", "left", "diotic"]

for channel in channelOpts

    asChord = asynchChord(freqs=freqs, levels=levels, phases=phases,
                     tonesDur=0.2, tonesRampDur=0.01, tonesChannel=channel,
                     SOA=0.06, sf=48000, maxLevel=100)
    wavwrite(asChord, wavDir*"asynchChord_"*channel*".wav", Fs=sf, nbits=nbits)

end

@test_throws(ErrorException, asynchChord(tonesDur=1, tonesRampDur=0.6))

@test_throws(ErrorException, asynchChord(tonesChannel="foo"))
