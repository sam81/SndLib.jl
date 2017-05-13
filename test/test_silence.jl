using Base.Test, SndLib, WAV

wavDir = "wavDir/"
sf = 48000
nbits = 32

channelOpts = ["mono", "diotic"]

for channel in channelOpts
    sl = silence(dur=0.1, channel="mono", sf=48000)
    wavwrite(sl, wavDir*"silence_"*channel*".wav", Fs=sf, nbits=nbits)
end

## Test error condition
@test_throws(ErrorException, silence(channel="foo"))

