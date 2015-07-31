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
* `phase`: Starting phase in radians.
* `AMPhase`: Starting AM phase in radians.
* `level`: Tone level in dB SPL. 
* `duration`: Tone duration (excluding ramps) in milliseconds.
* `ramp`: Duration of the onset and offset ramps in milliseconds.
        The total duration of the sound will be duration+ramp*2.
* `channel`: Channel in which the tone will be generated.
* `fs`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: 2-dimensional array of floats
       
##### Examples:

```julia
snd = AMTone(frequency=1000, AMFreq=20, AMDepth=1, phase=0, AMPhase,
level=65, duration=980, ramp=10, channel="Both", fs=48000, maxLevel=100)
```    


*source:*
[SndLib/src/SndLib.jl:57](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

---

<a id="method__broadbandnoise.1" class="lexicon_definition"></a>
#### broadbandNoise() [¶](#method__broadbandnoise.1)
Synthetise a broadband noise.

##### Parameters:

* `spectrumLevel`: Intensity spectrum level of the noise in dB SPL. 
* `duration`: Noise duration (excluding ramps) in milliseconds.
* `ramp`: Duration of the onset and offset ramps in milliseconds.
        The total duration of the sound will be duration+ramp*2.
* `channel`: Channel in which the noise will be generated.
* `fs`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : 2-dimensional array of floats
        The array has dimensions (nSamples, 2).
       
##### Examples:

```julia
noise = broadbandNoise(spectrumLevel=40, duration=180, ramp=10,
    ...     channel="Both", fs=48000, maxLevel=100)
```



*source:*
[SndLib/src/SndLib.jl:121](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

---

<a id="method__complextone.1" class="lexicon_definition"></a>
#### complexTone() [¶](#method__complextone.1)
Synthetise a complex tone.

##### Parameters:

* `F0`: Tone fundamental frequency in hertz.
* `harmPhase : one of 'Sine', 'Cosine', 'Alternating', 'Random', 'Schroeder'
        Phase relationship between the partials of the complex tone.
* `lowHarm`: Lowest harmonic component number.
* `highHarm`: Highest harmonic component number.
* `stretch`: Harmonic stretch in %F0. Increase each harmonic frequency by a fixed value
        that is equal to (F0*stretch)/100. If 'stretch' is different than
        zero, an inhanmonic complex tone will be generated.
* `level`: The level of each partial in dB SPL.
* `duration`: Tone duration (excluding ramps) in milliseconds.
* `ramp`: Duration of the onset and offset ramps in milliseconds.
        The total duration of the sound will be duration+ramp*2.
* `channel`: 'Right', 'Left', 'Both', 'Odd Right' or 'Odd Left'.
        Channel in which the tone will be generated. If `channel`
        if `Odd Right`, odd numbered harmonics will be presented
        to the right channel and even number harmonics to the left
        channel. The opposite is true if `channel` is `Odd Left`.
* `fs`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: 2-dimensional array of floats
        The array has dimensions (nSamples, 2).

##### Examples

```julia
 ct = complexTone(F0=440, harmPhase="Sine", lowHarm=3, highHarm=10,
    ...     stretch=0, level=55, duration=180, ramp=10, channel="Both",
    ...     fs=48000, maxLevel=100)
```



*source:*
[SndLib/src/SndLib.jl:212](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

---

<a id="method__fir2filt.1" class="lexicon_definition"></a>
#### fir2Filt(f1, f2, f3, f4, snd, fs) [¶](#method__fir2filt.1)
Filter signal with a fir2 filter.

This function designs and applies a fir2 filter to a sound.
The frequency response of the ideal filter will transition
from 0 to 1 between 'f1' and 'f2', and from 1 to zero
between 'f3' and 'f4'. The frequencies must be given in
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
* `snd`:The sound to be filtered.
* `fs`: Sampling frequency of 'snd'.

##### Returns:

* `snd`: 2-dimensional array of floats

##### Notes:

If 'f1' and 'f2' are zero the filter will be lowpass.
If 'f3' and 'f4' are equal to or greater than the nyquist
frequency (fs/2) the filter will be highpass.
In the other cases the filter will be bandpass.

The order of the filter (number of taps) is fixed at 256.
This function uses internally 'scipy.signal.firwin2'.
       
##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, duration=180, ramp=10,
     channel="Both", fs=48000, maxLevel=100)
lpNoise = fir2Filt(f1=0, f2=0, f3=1000, f4=1200, 
     snd=noise, fs=48000) #lowpass filter
hpNoise = fir2Filt(f1=0, f2=0, f3=24000, f4=26000, 
     snd=noise, fs=48000) #highpass filter
bpNoise = fir2Filt(f1=400, f2=600, f3=4000, f4=4400, 
.     snd=noise, fs=48000) #bandpass filter
```   



*source:*
[SndLib/src/SndLib.jl:598](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

---

<a id="method__gate.1" class="lexicon_definition"></a>
#### gate!(sig, ramps, fs) [¶](#method__gate.1)
Impose onset and offset ramps to a sound.

##### Parameters:

* `ramps`: The duration of the ramps.
* `sig`: The signal on which the ramps should be imposed.
* `fs`: The sampling frequency os 'sig'

##### Returns

* `sig`: The ramped signal.

##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, duration=200, ramp=0,
channel="Both", fs=48000, maxLevel=100)
gate!(sig=noise, ramps=10, fs=48000)
```



*source:*
[SndLib/src/SndLib.jl:357](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

---

<a id="method__makesilence.1" class="lexicon_definition"></a>
#### makeSilence() [¶](#method__makesilence.1)
Generate a silence.

This function just fills an array with zeros for the
desired duration.
    
##### Parameters:

* `duration`: Duration of the silence in milliseconds.
* `fs`: Samplig frequency in Hz.

##### Returns:

* `snd`: 2-dimensional array of floats
The array has dimensions (nSamples, 2).
       

##### Examples:

```julia
sil = makeSilence(duration=200, fs=48000)
```


*source:*
[SndLib/src/SndLib.jl:395](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

---

<a id="method__puretone.1" class="lexicon_definition"></a>
#### pureTone() [¶](#method__puretone.1)


*source:*
[SndLib/src/SndLib.jl:405](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

---

<a id="method__scalelevel.1" class="lexicon_definition"></a>
#### scaleLevel(sig, level) [¶](#method__scalelevel.1)
Increase or decrease the amplitude of a sound signal.

##### Parameters:

*`level`: Desired increment or decrement in dB SPL.
* `signal`: Signal to scale.

##### Returns:

* `sig`: 2-dimensional array of floats
       
##### Examples:

```julia
noise = broadbandNoise(spectrumLevel=40, duration=180, ramp=10,
channel="Both", fs=48000, maxLevel=100)
noise = scale(sig=noise, level=-10) #reduce level by 10 dB
```


*source:*
[SndLib/src/SndLib.jl:458](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

---

<a id="method__steepnoise.1" class="lexicon_definition"></a>
#### steepNoise() [¶](#method__steepnoise.1)
    
Synthetise band-limited noise from the addition of random-phase
sinusoids.

##### Parameters:

* `frequency1`: Start frequency of the noise.
* `frequency2`: End frequency of the noise.
* `level`: Noise spectrum level.
* `duration`: Tone duration (excluding ramps) in milliseconds.
* `ramp`: Duration of the onset and offset ramps in milliseconds.
        The total duration of the sound will be duration+ramp*2.
* `channel`: 'Right', 'Left' or 'Both'. Channel in which the tone will be generated.
* `fs`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: 2-dimensional array of floats
        The array has dimensions (nSamples, 2).
       
##### Examples:

```julia
nbNoise = steepNoise(frequency=440, frequency2=660, level=65,
duration=180, ramp=10, channel="Right", fs=48000, maxLevel=100)
```



*source:*
[SndLib/src/SndLib.jl:497](file:///home/sam/.julia/v0.3/SndLib/src/SndLib.jl)

