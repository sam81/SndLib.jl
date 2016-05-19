 
# Introduction

``SndLib.jl`` is a Julia module to synthetise sounds for psychoacoustics experiments. Detailed documentation for each function is provided in the [Function Reference](http://samcarcagno.altervista.org/SndLib/site/API/). A few examples of the available functions are given below.

The following code generates a 1 kHz pure tone:

```julia
using SndLib, DSP, PyPlot

pt = pureTone(frequency=1000, dur=0.1, rampDur=0.01, sf=48000)
figure()
plot(pt[:,1], color="black")
xlabel("Time (s)"); ylabel("Amplitude (a.u.)")
savefig("pure_tone.png")
```

![Pure Tone](pure_tone.png)

if you're on Linux or OS X, you can listen to it using the `sound` function:

```julia
sound(pt)
```

The following code generates an amplitude-modulated tone:

```julia
amtone = AMTone(carrierFreq=1000, AMFreq=50, AMDepth=1, AMPhase=-pi/2, dur=0.1, rampDur=0.01, sf=48000)
figure()
plot(amtone[:,1], color="black")
xlabel("Time (s)"); ylabel("Amplitude (a.u.)")
savefig("AM_tone.png")
```
![Amplitude Modulated Tone](AM_tone.png)

we can use the `DSP` module to plot the spectrum of the tone:

```julia
spec = periodogram(amtone[:,1], onesided=true, fs=48000, window=hamming)
idx = find(spec.freq .<= 2000)
figure()
plot(spec.freq[idx], 10*log10(spec.power)[idx], color="black")
xlabel("Frequency (Hz)"); ylabel("Level (dB)")
savefig("AM_tone_spectrum.png")
```
![Amplitude Modulated Tone Spectrum](AM_tone_spectrum.png)
