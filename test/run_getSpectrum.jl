using DSP, SndLib, Test

sf = 48000
pt = pureTone(sf=sf)
p, freqArr = getSpectrum(pt, sf)
p, freqArr = getSpectrum(pt, sf, window=hamming)
p, freqArr = getSpectrum(pt, sf, powerOfTwo=true)
