using Test, SndLib

#test when snd2 is embedded in snd1
x = pureTone(dur=1)
y = pureTone(dur=0.2)
addSounds(x,y, delay=0.2)

#test error thrown when snd1 and snd2 have different number of channels
y = pureTone(dur=0.2, channel="mono")

@test_throws(ErrorException, addSounds(x,y, delay=0.2))

noise = broadbandNoise(spectrumLevel=40, dur=0.1, rampDur=0.01, channel="diotic", sf=48000, maxLevel=100)
@test_throws(ErrorException, delayAdd!(noise, delay=1/440, gain=1, iterations=6, configuration="add foo", channel=[1,2], sf=48000))
