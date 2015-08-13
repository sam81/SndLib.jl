# SndLib

## Exported

---

<a id="method__amtone.1" class="lexicon_definition"></a>
#### AMTone() [¶](#method__amtone.1)
Generate an amplitude modulated tone.

##### Parameters:

* `frequency`: Carrier frequency in hertz.
* `AMFreq`:  Amplitude modulation frequency in Hz.
* `AMDepth`:  Amplitude modulation depth (a value of 1 corresponds to 100% modulation). 
* `carrierPhase`: Starting phase in radians.
* `AMPhase`: Starting AM phase in radians.
* `level`: Tone level in dB SPL. 
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone will be generated (`mono`, `right`, `left`, or diotic`).
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: array 
       
##### Examples:

```julia
snd = AMTone(carrierFreq=1000, AMFreq=20, AMDepth=1, carrierPhase=0,
AMPhase=0, level=65, dur=1, rampDur=0.01, channel="diotic", sf=48000,
maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:54](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__amtoneipd.1" class="lexicon_definition"></a>
#### AMToneIPD() [¶](#method__amtoneipd.1)
Generate an amplitude modulated tone with an IPD in the carrier and/or
modulation phase.

##### Parameters:

* `frequency`: Carrier frequency in hertz.
* `AMFreq`:  Amplitude modulation frequency in Hz.
* `AMDepth`:  Amplitude modulation depth (a value of 1 corresponds to 100% modulation). 
* `carrierPhase`: Starting phase in radians.
* `AMPhase`: Starting AM phase in radians.
* `carrierIPD`: IPD to apply to the carrier phase in radians.
* `AMIPD`: IPD to apply to the AM phase in radians.
* `level`: Tone level in dB SPL. 
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: channel in which the IPD(s) will be applied (`right`, `left`).
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: array 
       
##### Examples:

```julia
snd = AMToneIPD(carrierFreq=1000, AMFreq=20, AMDepth=1, carrierPhase=0,
AMPhase=0, carrierIPD=0, AMIPD=pi/2, level=65, dur=1, rampDur=0.01,
channel="right", sf=48000, maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:133](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__erbdistance.1" class="lexicon_definition"></a>
#### ERBDistance(f1::Real, f2::Real) [¶](#method__erbdistance.1)
Compute the distance in equivalent rectangular bandwiths (ERBs)
between the frequencies `f1` and `f2`.

##### Parameters

* `f1`: frequency 1 in Hz
* `f2`: frequency 2 in Hz

##### Returns

* `deltaERB`: distance between f1 and f2 in ERBs.

##### References

.. [GM] Glasberg, B. R., & Moore, B. C. J. (1990). Derivation of auditory filter shapes from notched-noise data. Hear. Res., 47(1-2), 103–38.

##### Examples

```julia
ERBDistance(1000, 1200)
```


*source:*
[SndLib/src/utils.jl:147](file:///home/sam/.julia/v0.3/SndLib/src/utils.jl)

---

<a id="method__fmtone.1" class="lexicon_definition"></a>
#### FMTone() [¶](#method__fmtone.1)
Generate a frequency modulated tone.

##### Parameters

* `carrierFreq`: Carrier frequency in hertz. This is the frequency of the tone at FM zero crossing.
* `MF`: Modulation frequency in Hz.
* `MI`: Modulation index (MI), also called beta. The MI is equal to deltaF/fm, where
        deltaF is the maximum deviation of the instantaneous frequency from
        the carrier frequency.
* `phase`: Starting phase in radians.
* `level`: Tone level in dB SPL. 
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone will be generated, `mono`, `right`,
  `left` or `diotic`
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of
  amplitude 1.

##### Returns

* `snd`: 2-dimensional array of floats
       
Examples

```julia
snd = FMTone(carrierFreq=1000, MF=40, MI=1, phase=0, level=55, dur=1,
     rampDur=0.01, channel="diotic", sf=48000, maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:482](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__irn.1" class="lexicon_definition"></a>
#### IRN() [¶](#method__irn.1)

Generate an iterated rippled noise

##### Parameters

* `delay`: delay in seconds
* `gain`: The gain to apply to the delayed signal
* `iterations`: The number of iterations of the delay-add cycle
* `configuration`: If `add same`, the output of iteration N-1 is added to delayed signal of the current iteration.
If `add original`, the original signal is added to the delayed signal of the current iteration.
* `spectrumLevel`: Intensity spectrum level of the noise in dB SPL.
* `dur`: Noise duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the noise will be generated (`mono`, `right`, `left`, `diotic`, `dichotic`).
* `sf`: Sampling frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns

* `snd`: array with dimensions (nSamples, nChannels).

##### Examples

```julia
irn = IRN(delay=1/440, gain=1, iterations=6, configuration="add same",
          spectrumLevel=25, dur=1, rampDur=0.01, channel="diotic",
          sf=48000, maxLevel=101)
```


*source:*
[SndLib/src/snd_generate.jl:724](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__itdshift.1" class="lexicon_definition"></a>
#### ITDShift!{T<:Real}(sig::Array{T<:Real, 2}, f1::Real, f2::Real) [¶](#method__itdshift.1)
Set the ITD of a sound within the frequency region bounded by `f1` and `f2`

##### Parameters

* `sig`: Input signal.
* `f1`: The start point in Hertz of the frequency region in which
        to apply the ITD.
* `f2`: The end point in Hertz of the frequency region in which 
        to apply the ITD.
* `ITD`: The amount of ITD shift in seconds
* `channel`: `right` or `left`. The channel in which to apply the shift.
* `sf`: The sampling frequency of the sound.
        
##### Returns

* `out` : 2-dimensional array of floats

##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, dur=1, rampDur=0.01,
     channel="diotic", sf=48000, maxLevel=100)
hp = ITDShift!(noise, 500, 600, ITD=300/1000000,
channel="left", sf=48000) #this generates a Dichotic Pitch
```


*source:*
[SndLib/src/snd_process.jl:409](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__itdtoipd.1" class="lexicon_definition"></a>
#### ITDToIPD{T<:Real}(ITD::Real, freq::Union(AbstractArray{T<:Real, 1}, T<:Real)) [¶](#method__itdtoipd.1)

Convert an interaural time difference to an equivalent interaural
phase difference for a given frequency.

##### Parameters

* `itd`: Interaural time difference in seconds.
* `freq`: Frequency in hertz.

##### Returns

`ipd`: 

##### Examples

```julia
itd = 300 #microseconds
itd = 300/1000000 #convert to seconds
ITDToIPD(itd, 1000)
```


*source:*
[SndLib/src/snd_process.jl:487](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__addsounds.1" class="lexicon_definition"></a>
#### addSounds{T<:Real, P<:Real}(snd1::Array{T<:Real, 2}, snd2::Array{P<:Real, 2}) [¶](#method__addsounds.1)
Add or concatenate two sounds.

##### Parameters

* `snd1`: First sound.
* `snd2`: Second sound.
* `delay`: Delay in seconds between the onset of `snd1` and the onset of `snd2`
* `sf`: Sampling frequency in hertz of the two sounds.

##### Returns

* `snd` : array with dimensions (nSamples, nChannels) 
       
##### Examples

```julia
snd1 = pureTone(frequency=440, phase=0, level=65, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
snd2 = pureTone(frequency=880, phase=0, level=65, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
snd = addSounds(snd1, snd2, delay=1, sf=48000)
```


*source:*
[SndLib/src/snd_process.jl:49](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__broadbandnoise.1" class="lexicon_definition"></a>
#### broadbandNoise() [¶](#method__broadbandnoise.1)
Synthetise a broadband noise.

##### Parameters:

* `spectrumLevel`: Intensity spectrum level of the noise in dB SPL. 
* `dur`: Noise duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the noise will be generated (`mono`, `right`, `left`,
             `diotic`, or `dichotic`). If `dichotic` the noise will be uncorrelated
             at the two ears.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : array. The array has dimensions (nSamples, nCh).
       
##### Examples:

```julia
noise = broadbandNoise(spectrumLevel=20, dur=180, rampDur=10,
         channel="diotic", sf=48000, maxLevel=100)
```



*source:*
[SndLib/src/snd_generate.jl:196](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__centdistance.1" class="lexicon_definition"></a>
#### centDistance(f1::Real, f2::Real) [¶](#method__centdistance.1)
Compute the distance in cents
between the frequencies `f1` and `f2`.

##### Parameters

* `f1`: frequency 1 in Hz
* `f2`: frequency 2 in Hz

##### Returns

* `deltaCents`: distance between f1 and f2 in cents.

##### Examples

```julia
centDistance(1000, 1200)
```


*source:*
[SndLib/src/utils.jl:76](file:///home/sam/.julia/v0.3/SndLib/src/utils.jl)

---

<a id="method__complextone.1" class="lexicon_definition"></a>
#### complexTone() [¶](#method__complextone.1)
Synthetise a complex tone.

##### Parameters:

* `F0`: Tone fundamental frequency in hertz.
* `harmPhase : one of `sine`, `cosine`, `alternating`, `random`, `schroeder`
        Phase relationship between the partials of the complex tone.
* `lowHarm`: Lowest harmonic component number.
* `highHarm`: Highest harmonic component number.
* `stretch`: Harmonic stretch in %F0. Increase each harmonic frequency by a fixed value
        that is equal to (F0*stretch)/100. If `stretch` is different than
        zero, an inhanmonic complex tone will be generated.
* `level`: The level of each partial in dB SPL.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: `right`, `left`, `diotic`, `odd right` or `odd left`.
        Channel in which the tone will be generated. If `channel`
        if `odd right`, odd numbered harmonics will be presented
        to the right channel and even number harmonics to the left
        channel. The opposite is true if `channel` is `odd left`.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: array of floats with dimensions (nSamples, nChannels).

##### Examples

```julia
 ct = complexTone(F0=440, harmPhase="sine", lowHarm=3, highHarm=10,
         stretch=0, level=55, dur=180, rampDur=10, channel="diotic",
         sf=48000, maxLevel=100)
```



*source:*
[SndLib/src/snd_generate.jl:308](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__delayadd.1" class="lexicon_definition"></a>
#### delayAdd!{T<:Real}(sig::Array{T<:Real, 2}) [¶](#method__delayadd.1)
Delay and add algorithm for the generation of iterated rippled noise.

##### Parameters

* `sig`: The signal to manipulate
* `delay`: delay in seconds
* `gain`: The gain to apply to the delayed signal
* `iterations`: The number of iterations of the delay-add cycle
* `configuration`: If `add same`, the output of iteration N-1 is added to delayed signal of the current iteration.
If `add original`, the original signal is added to delayed signal of the current iteration.
* `channel`: a number or a vector of numbers indicating to which columns of `sig` the delay and add process should be applied
* `sf`: Sampling frequency in Hz.

##### Returns

* `snd`:

##### References

.. [YPS1996] Yost, W. A., Patterson, R., & Sheft, S. (1996). A time domain description for the pitch strength of iterated rippled noise. J. Acoust. Soc. Am., 99(2), 1066–78. 

##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, dur=1, rampDur=0.01,
channel="diotic", sf=48000, maxLevel=100)
irn = delayAdd!(noise, delay=1/440, gain=1, iterations=6, configuration="add same", channel=[1,2], sf=48000)
```


*source:*
[SndLib/src/snd_process.jl:133](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__fir2filt.1" class="lexicon_definition"></a>
#### fir2Filt!{T<:Real}(f1::Real, f2::Real, f3::Real, f4::Real, snd::Array{T<:Real, 2}) [¶](#method__fir2filt.1)
Filter signal with a fir2 filter.

This function designs and applies a fir2 filter to a sound.
The frequency response of the ideal filter will transition
from 0 to 1 between `f1` and `f2`, and from 1 to zero
between `f3` and `f4`. The frequencies must be given in
increasing order.

##### Parameters:

* `f1`: Frequency in hertz of the point at which the transition
        for the low-frequency cutoff ends. 
* `f2`: Frequency in hertz of the point at which the transition
        for the low-frequency cutoff starts.
* `f3`: Frequency in hertz of the point at which the transition
        for the high-frequency cutoff starts.
* `f4`: Frequency in hertz of the point at which the transition
        for the high-frequency cutoff ends. 
* `snd`: The sound to be filtered.
* `nTaps`: Number of filter taps.
* `sf`: Sampling frequency of `snd`.

##### Returns:

* `snd`: 2-dimensional array of floats

##### Notes:

If `f1` and `f2` are zero the filter will be low pass.
If `f3` and `f4` are equal to or greater than the nyquist
frequency (sf/2) the filter will be high pass.
In the other cases the filter will be band pass.

This function uses internally 'scipy.signal.firwin2'.
       
##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, dur=1, rampDur=0.01,
     channel="diotic", sf=48000, maxLevel=100)
lpNoise = fir2Filt!(0, 0, 1000, 1200, 
     noise, nTaps=256, sf=48000) #lowpass filter
noise = broadbandNoise(spectrumLevel=40, dur=1, rampDur=0.01,
     channel="diotic", sf=48000, maxLevel=100)
hpNoise = fir2Filt!(0, 0, 2400, 2600, 
     noise, nTaps=256, sf=48000) #highpass filter
noise = broadbandNoise(spectrumLevel=40, dur=1, rampDur=0.01,
     channel="diotic", sf=48000, maxLevel=100)
bpNoise = fir2Filt!(400, 600, 4000, 4400, 
     noise, nTaps=256, sf=48000) #bandpass filter
```


*source:*
[SndLib/src/snd_process.jl:239](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__freqfromcentinterval.1" class="lexicon_definition"></a>
#### freqFromCentInterval{T<:Real, P<:Real}(f1::Union(AbstractArray{T<:Real, 1}, T<:Real), deltaCent::Union(AbstractArray{P<:Real, 1}, P<:Real)) [¶](#method__freqfromcentinterval.1)
Compute the frequency, in Hz, corresponding to a distance,
in equivalent cents of `deltaCents` from `f1`.

##### Parameters
 
* `f1`: frequency at one end of the interval in Hz
* `deltaCents`: distance in cents

##### Returns

* `f2`: frequency at the other end of the interval
    
##### Examples

```julia
freqFromCentInterval(100, 1.5)
freqFromCentInterval(100, -1.5)
#for several frequencies
freqFromCentInterval([100, 200, 300], 1.5)
# for several distances
freqFromCentInterval(100, [1, 1.5, 2])
```


*source:*
[SndLib/src/utils.jl:49](file:///home/sam/.julia/v0.3/SndLib/src/utils.jl)

---

<a id="method__freqfromerbinterval.1" class="lexicon_definition"></a>
#### freqFromERBInterval{T<:Real, P<:Real}(f1::Union(AbstractArray{T<:Real, 1}, T<:Real), deltaERB::Union(AbstractArray{P<:Real, 1}, P<:Real)) [¶](#method__freqfromerbinterval.1)
Compute the frequency, in Hz, corresponding to a distance,
in equivalent rectangular bandwidths (ERBs), of `deltaERB` from `f1`.

##### Parameters
 
* `f1`: frequency at one end of the interval in Hz
* `deltaERB`: distance in ERBs

##### Returns

* `f2`: frequency at the other end of the interval

##### References

.. [GM] Glasberg, B. R., & Moore, B. C. J. (1990). Derivation of auditory filter shapes from notched-noise data. Hear. Res., 47(1-2), 103–38.
    
##### Examples

```julia
freqFromERBInterval(100, 1.5)
freqFromERBInterval(100, -1.5)
#for several frequencies
freqFromERBInterval([100, 200, 300], 1.5)
# for several distances
freqFromERBInterval(100, [1, 1.5, 2])
```


*source:*
[SndLib/src/utils.jl:116](file:///home/sam/.julia/v0.3/SndLib/src/utils.jl)

---

<a id="method__gate.1" class="lexicon_definition"></a>
#### gate!{T<:Real}(sig::Array{T<:Real, 2}) [¶](#method__gate.1)
Impose onset and offset ramps to a sound.

##### Parameters:

* `sig`: The signal on which the ramps should be imposed.
* `rampDur`: The duration of the ramps.
* `sf`: The sampling frequency of `sig`

##### Returns

* `sig`: The ramped signal.

##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, dur=2, rampDur=0,
channel="diotic", sf=48000, maxLevel=100)
gate!(noise, rampDur=0.01, sf=48000)
```


*source:*
[SndLib/src/snd_process.jl:321](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__getrms.1" class="lexicon_definition"></a>
#### getRMS{T<:Real}(sig::Array{T<:Real, 2}, channel::Union(Integer, String)) [¶](#method__getrms.1)
Compute the root mean square (RMS) value of the signal.

##### Parameters

* `sig`: The signal for which the RMS needs to be computed.
* `channel`: Either an integer indicating the channel number,
or `each` for a list of the RMS values in each channel, or `all`
for the RMS across all channels.
##### Returns

* `RMS`: The RMS of `sig`

##### Examples

```julia
pt = pureTone(frequency=440, phase=0, level=65, dur=1,
     rampDur=0.01, channel="right", sf=48000, maxLevel=100)
getRMS(pt, 1)
getRMS(pt, 2)
getRMS(pt, "each")
getRMS(pt, "all")
```


*source:*
[SndLib/src/snd_process.jl:361](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__hugginspitch.1" class="lexicon_definition"></a>
#### hugginsPitch() [¶](#method__hugginspitch.1)

Synthetise a complex Huggings Pitch.

##### Parameters

* `F0`: The centre frequency of the F0 of the complex in hertz.
* `lowHarm`: Lowest harmonic component number.
* `highHarm`: Highest harmonic component number.
* `spectrumLevel`: The spectrum level of the noise from which
        the Huggins pitch is derived in dB SPL.
        If `noiseType` is `pink`, the spectrum level
        will be equal to `spectrumLevel` at 1 kHz.
* `bandwidth`: Bandwidth of the frequency regions in which the
        phase transitions occurr.
* `bandwidthUnit`: `Hz`, `Cent`, or `ERB`. Defines whether the bandwith of
  the decorrelated bands is expressed in hertz (Hz), cents (Cent), or
  equivalent rectangular bandwidths (ERB).
* `dichoticDifference`: `IPD linear`, `IPD stepped`, `IPD random`, `ITD`.
   Selects whether the decorrelation in the target regions will be achieved
   by applying an interaural phase shift that increases linearly in the
   transition regions
   (`IPD linear`), a costant interaural phase shift (`IPD stepped`),
   a costant interaural time difference (ITD), or a random IPD shift
   (`IPD random`).
* `dichoticDifferenceValue`: For `IPD linear` this is the phase difference
   between the start and the end of each transition region, in radians.
   For `IPD stepped`, this is the phase offset, in radians, between
   the correlated and the uncorrelated regions.
   For `ITD` this is the ITD in the transition region, in seconds.
   For `IPD random`, this is the range of phase shift randomization
   in the uncorrelated regions.
* `phaseRelationship`: `NoSpi` or `NpiSo`. If `NoSpi`, the phase of the regions within each frequency band will
        be shifted. If `NpiSo`, the phase of the regions between each
        frequency band will be shifted.
* `stretch`: Harmonic stretch in %`F0`. Increase each harmonic frequency by a fixed value
        that is equal to (F0*stretch)/100. If `stretch` is different than
        zero, an inhanmonic complex tone will be generated.
* `noiseType`: `white` or `pink`. The type of noise used to derive the Huggins Pitch.
* `dur`: Complex duration (excluding ramps) in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns

* `snd`: 2-dimensional array of floats.
        The array has dimensions (nSamples, nChannels).

References

.. [CH] Cramer, E. M., & Huggins, W. H. (1958). Creation of Pitch through Binaural Interaction. J. Acoust. Soc. Am., 30(5), 413. 
.. [AS] Akeroyd, M. A., & Summerfield, a Q. (2000). The lateralization of simple dichotic pitches. J. Acoust. Soc. Am., 108(1), 316–334.
.. [ZH] Zhang, P. X., & Hartmann, W. M. (2008). Lateralization of Huggins pitch. J. Acoust. Soc. Am., 124(6), 3873–87. 

##### Examples

```julia
hp = hugginsPitch(F0=300, lowHarm=1, highHarm=3, spectrumLevel=45,
bandwidth=100, bandwidthUnit="Hz", dichoticDifference="IPD stepped",
dichoticDifferenceValue=pi, phaseRelationship="NoSpi", stretch=0,
noiseType="white", dur=0.4, rampDur=0.01, sf=48000, maxLevel=101)
```
    


*source:*
[SndLib/src/snd_generate.jl:597](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__intncyclesfreq.1" class="lexicon_definition"></a>
#### intNCyclesFreq(freq::Real, dur::Real) [¶](#method__intncyclesfreq.1)

Compute the frequency closest to 'freq' that has an integer number
of cycles for the given sound duration.

##### Parameters

* `frequency`: Frequency in hertz.
* `dur` : Duration of the sound, in seconds.

##### Returns

* `adjFreq`: float
       
##### Examples

```julia
intNCyclesFreq(2.1, 1000)
intNCyclesFreq(2, 1000)
```



*source:*
[SndLib/src/utils.jl:180](file:///home/sam/.julia/v0.3/SndLib/src/utils.jl)

---

<a id="method__makepink.1" class="lexicon_definition"></a>
#### makePink!{T<:Real}(sig::Array{T<:Real, 2}) [¶](#method__makepink.1)
Convert a white noise into a pink noise.

The spectrum level of the pink noise at the frequency `ref`
will be equal to the spectrum level of the white noise input
to the function.

##### Parameters

* `sig`: The white noise to be turned into a pink noise.
* `sf`: Sampling frequency of the sound.
* `ref`: Reference frequency in Hz. The amplitude of the other
        frequencies will be scaled with respect to the amplitude
        of this frequency.

##### Returns

* `snd` : array of floats. The array has dimensions (nSamples, nChannels).

##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, dur=1, rampDur=0.01,
channel="diotic", sf=48000, maxLevel=100)
noise = makePink!(noise, sf=48000, ref=1000)
```


*source:*
[SndLib/src/snd_process.jl:613](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__phaseshift.1" class="lexicon_definition"></a>
#### phaseShift!{T<:Real}(sig::Array{T<:Real, 2}, f1::Real, f2::Real) [¶](#method__phaseshift.1)
Shift the interaural phases of a sound within a given frequency region.

##### Parameters

* `sig`: Input signal.
* `f1`: The start point of the frequency region to be
        phase-shifted in hertz.
* `f2`: The end point of the frequency region to be
        phase-shifted in hertz.
* `phaseShift`: The amount of phase shift in radians. 
* `shiftType`: If `linear` the phase changes progressively
        on a linear Hz scale from X to X+`phaseShift` from f1 to f2.
        If `step` `phaseShift` is added as a constant to the
        phases from f1 to f2.
        If `random` a random phase shift from 0 to `phaseShift`
        is added to each frequency component from `f1` to `f2`.
* `channel`: The channel(s) in which to apply the phase shift.
* `sf`: The sampling frequency of the sound.
        
##### Returns

* `out`: 2-dimensional array of floats

##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, dur=1, rampDur=0.01,
channel="diotic", sf=48000, maxLevel=100)
noise = phaseShift!(noise, 500, 600, phaseShift=pi,
channel=2, sf=48000) #this generates a Dichotic Pitch
```


*source:*
[SndLib/src/snd_process.jl:531](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__puretone.1" class="lexicon_definition"></a>
#### pureTone() [¶](#method__puretone.1)
Synthetise a pure tone.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone will be generated.  (`right`, `left` or `diotic`)
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : 2-dimensional array of floats
        The array has dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureTone(frequency=440, phase=0, level=65, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:786](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__puretoneild.1" class="lexicon_definition"></a>
#### pureToneILD() [¶](#method__puretoneild.1)
Synthetise a pure tone with an interaural level difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `level`: Tone level in dB SPL.
* `ILD`: Interaural level difference in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the ILD will be applied.  (`right` or `left`)
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneILD(frequency=440, phase=0, level=65, ILD=10, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:859](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__puretoneipd.1" class="lexicon_definition"></a>
#### pureToneIPD() [¶](#method__puretoneipd.1)
Synthetise a pure tone with an interaural phase difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `IPD`: Interaural phase difference in radians.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone IPD will be applied.  (`right` or `left`)
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneIPD(frequency=440, phase=0, IPD=pi/2, level=65, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:916](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__puretoneipdild.1" class="lexicon_definition"></a>
#### pureToneIPDILD() [¶](#method__puretoneipdild.1)
Synthetise a pure tone with an interaural phase and interaural level difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `IPD`: Interaural phase difference in radians.
* `level`: Tone level in dB SPL.
* `ILD`: Interaural level difference in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the ILD will be applied.  (`right` or `left`)
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneIPDILD(frequency=440, phase=0, IPD=pi/2, level=65, ILD=10,
dur=1, rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:1019](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__puretoneitd.1" class="lexicon_definition"></a>
#### pureToneITD() [¶](#method__puretoneitd.1)
Synthetise a pure tone with an interaural time difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `ITD`: Interaural time difference in seconds.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the ITD will be applied.  (`right` or `left`)
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneITD(frequency=440, phase=0, ITD=0.004/1000, level=65, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:970](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__puretoneitdild.1" class="lexicon_definition"></a>
#### pureToneITDILD() [¶](#method__puretoneitdild.1)
Synthetise a pure tone with an interaural time and interaural level difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `ITD`: Interaural time difference in seconds.
* `level`: Tone level in dB SPL.
* `ILD`: Interaural level difference in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the ILD will be applied.  (`right` or `left`)
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneITDILD(frequency=440, phase=0, ITD=0.004/1000, level=65,
ILD=10, dur=1, rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```


*source:*
[SndLib/src/snd_generate.jl:1075](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__scalelevel.1" class="lexicon_definition"></a>
#### scaleLevel{T<:Real}(sig::Array{T<:Real, 2}) [¶](#method__scalelevel.1)
Increase or decrease the amplitude of a sound signal.

##### Parameters:

*`level`: Desired increment or decrement in dB SPL.
* `signal`: Signal to scale.

##### Returns:

* `sig`: 2-dimensional array of floats
       
##### Examples:

```julia
noise = broadbandNoise(spectrumLevel=40, dur=1, rampDur=0.01,
channel="diotic", sf=48000, maxLevel=100)
noise = scaleLevel(noise, level=-10) #reduce level by 10 dB
```


*source:*
[SndLib/src/snd_process.jl:660](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__silence.1" class="lexicon_definition"></a>
#### silence() [¶](#method__silence.1)
Generate a silence.

This function just fills an array with zeros for the
desired duration.
    
##### Parameters:

* `dur`: Duration of the silence in seconds.
* `sf`: Samplig frequency in Hz.

##### Returns:

* `snd`: 2-dimensional array of floats
The array has dimensions (nSamples, 2).
       

##### Examples:

```julia
sil = silence(dur=2, sf=48000)
```


*source:*
[SndLib/src/snd_generate.jl:1120](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

---

<a id="method__sound.1" class="lexicon_definition"></a>
#### sound{T<:Real}(snd::Array{T<:Real, 2}) [¶](#method__sound.1)
Simple function to play sounds. Uses aplay on Linux and afplay on OSX.
Windows is not currently supported.

##### Arguments

* `snd`: The sound to be played.
* `sf`: The sampling frequency.



*source:*
[SndLib/src/snd_process.jl:684](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__sound.2" class="lexicon_definition"></a>
#### sound{T<:Real}(snd::Array{T<:Real, 2}, sf::Integer) [¶](#method__sound.2)
Simple function to play sounds. Uses aplay on Linux and afplay on OSX.
Windows is not currently supported.

##### Arguments

* `snd`: The sound to be played.
* `sf`: The sampling frequency.



*source:*
[SndLib/src/snd_process.jl:684](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__sound.3" class="lexicon_definition"></a>
#### sound{T<:Real}(snd::Array{T<:Real, 2}, sf::Integer, nbits::Integer) [¶](#method__sound.3)
Simple function to play sounds. Uses aplay on Linux and afplay on OSX.
Windows is not currently supported.

##### Arguments

* `snd`: The sound to be played.
* `sf`: The sampling frequency.



*source:*
[SndLib/src/snd_process.jl:684](file:///home/sam/.julia/v0.3/SndLib/src/snd_process.jl)

---

<a id="method__steepnoise.1" class="lexicon_definition"></a>
#### steepNoise() [¶](#method__steepnoise.1)
Synthetise band-limited noise from the addition of random-phase
sinusoids.

##### Parameters:

* `f1`: Start frequency of the noise.
* `f2`: End frequency of the noise.
* `level`: Noise spectrum level.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: `mono`, `right`, `left` or `diotic`. Channel in which the tone will be generated.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: 2-dimensional array of floats. The array has dimensions (nSamples, nChannels).
       
##### Examples:

```julia
nbNoise = steepNoise(f1=440, f2=660, level=65,
dur=1, rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```



*source:*
[SndLib/src/snd_generate.jl:1164](file:///home/sam/.julia/v0.3/SndLib/src/snd_generate.jl)

