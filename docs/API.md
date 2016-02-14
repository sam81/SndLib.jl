Generate an amplitude modulated tone.

##### Parameters

* `carrierFreq`: Carrier frequency in hertz.
* `AMFreq`:  Amplitude modulation frequency in Hz.
* `AMDepth`: Amplitude modulation depth (a value of 1 corresponds to 100%
  modulation).
* `carrierPhase`: Starting phase in radians.
* `AMPhase`: Starting AM phase in radians.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone will be generated, one of `mono`, `right`,
  `left`, or `diotic`.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of
  amplitude 1.

##### Returns:

* `snd`: array

##### Examples:

```julia
snd = AMTone(carrierFreq=1000, AMFreq=20, AMDepth=1, carrierPhase=0,
AMPhase=0, level=65, dur=1, rampDur=0.01, channel="diotic", sf=48000,
maxLevel=100)
```

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

Generate an asynchronous chord.

This function will add a set of pure tones with a given
stimulus onset asynchrony (SOA). The temporal order of the
successive tones is random.

##### Parameters

* `freqs`: Frequencies of the chord components in hertz.
* `levels`: Level of each chord component in dB SPL.
* `phases`: Starting phase of each chord component.
* `tonesDur`: Duration of the tones composing the chord in seconds.
        All tones have the same duration.
* `tonesRampDur`: Duration of the onset and offset ramps in seconds.
        The total duration of the tones will be tonesDuration+ramp*2.
* `tonesChannel`: Channel in which the tones will be generated, one of `mono`, `right`, `left` or `diotic`.
* `SOA`: Onset asynchrony between the chord components.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: float
        Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns

* `snd` : 2-dimensional array of floats

##### Examples

```julia
freqs = [250, 500, 1000, 1500]
levels = [50, 50, 50, 50]
phases = [0, 0, 0, 0]
c1 = asynchChord(freqs=freqs, levels=levels, phases=phases,
tonesDur=0.2, tonesRampDur=0.01, tonesChannel="diotic",
SOA=0.06, sf=48000, maxLevel=100)
```

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


Generate a sinusoidally amplitude-modulated noise with an exponentially
modulated AM frequency.

##### Parameters

* `carrierFreq`: Carrier AM frequency in hertz.
* `MF`: Amplitude modulation frequency in Hz.
* `deltaCents`: AM frequency excursion in cents. The instataneous AM frequency of the noise will vary from `carrierFreq`**(-`deltaCents`/1200) to `carrierFreq`**(`deltaCents`/1200).
* `FMPhase`: Starting phase of the AM modulation in radians.
* `AMDepth`: Amplitude modulation depth.
* `spectrumLevel`: Noise spectrum level in dB SPL.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the noise will be generated, one of `right`, `left`, `diotic`, or `dichotic`.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of
  amplitude 1.

##### Returns

* `snd` : 2-dimensional array of floats

##### Examples

```julia
snd = expAMNoise(carrierFreq=150, MF=2.5, deltaCents=600, FMPhase=pi, AMDepth = 1,
     spectrumLevel=30, dur=0.4, rampDur=0.01, channel="diotic", sf=48000, maxLevel=101)
```


Generate a tone frequency modulated with an exponential sinusoid.

##### Parameters

* `carrierFreq`: Carrier frequency in hertz.
* `MF`: Modulation frequency in Hz.
* `deltaCents`: Frequency excursion in cents. The instataneous frequency of the tone
         will vary from `carrierFreq**(-deltaCents/1200)` to `carrierFreq**(deltaCents/1200)`.
* `FMPhase`: Starting FM phase in radians.
* `phase`: Starting phase in radians.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone will be generated, one of "mono", `right`, `left` or `diotic`.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of
        amplitude 1.

##### Returns

`snd` : 2-dimensional array of floats

##### Examples

```julia
tone_peak = expSinFMTone(carrierFreq=450, MF=5, deltaCents=300, FMPhase=pi, phase=0, level=60,
    dur=0.2, rampDur=0.01, channel="diotic", sf=48000, maxLevel=101)
tone_trough = expSinFMTone(carrierFreq=450, MF=5, deltaCents=300, FMPhase=0, phase=0, level=60,
    dur=0.2, rampDur=0.01, channel="diotic", sf=48000, maxLevel=101)
```

Synthetise a complex tone with an embedded frequency modulation (FM)
starting and stopping at a chosen time after the tone onset.

##### Parameters

* `midF0`: F0 at the FM zero crossing
* `harmPhase`: one of 'sine', 'cosine', 'alternating', 'random', 'schroeder'.
        Phase relationship between the partials of the complex tone.
* `lowHarm`: Lowest harmonic component number.
* `highHarm`: Highest harmonic component number.
* `level`: The level of each partial in dB SPL.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `MF`: modulation frequency in Hz.
* `FMDepth`: FM depth in %.
* `FMStartPhase`: Starting phase of FM.
* `FMStartTime`: Start of FM in ms after start of tone.
* `FMDur`: Duration of FM, in seconds.
* `levelAdj`: If `true`, scale the harmonic level so that for a complex
        tone within a bandpass filter the overall level does not
        change with F0 modulations.
* `channel`: Channel in which the tone will be generated.
        One of 'right', 'left', 'diotic', 'odd right' or 'odd left'.
         If 'channel' is 'odd right', odd numbered harmonics will be presented
        to the right channel and even number harmonics to the left
        channel. The opposite is true if 'channel' is 'odd left'.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns

* `snd`: 2-dimensional array of floats

Examples

```julia
tone_up = FMComplex2(midF0=140, harmPhase="sine",
                 lowHarm=1, highHarm=10,
                 level=60, dur=0.45, rampDur=0.01,
                 MF=1.25, FMDepth=40, FMStartPhase=1.5*pi,
                 FMStartTime=0.025, FMDur=0.4,
                 levelAdj=true, channel="diotic",
                 sf=48000, maxLevel=101)
tone_down = FMComplex2(midF0=140, harmPhase="sine",
                 lowHarm=1, highHarm=10,
                 level=60, dur=0.45, rampDur=0.01,
                 MF=1.25, FMDepth=40, FMStartPhase=0.5*pi,
                 FMStartTime=0.025, FMDur=0.4,
                 levelAdj=true, channel="diotic",
                 sf=48000, maxLevel=101)
```

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

##### References

* [Cramer, E. M., & Huggins, W. H. (1958). Creation of Pitch through Binaural Interaction. J. Acoust. Soc. Am., 30(5), 413.](http://dx.doi.org/10.1121/1.1909628)
* [Akeroyd, M. A., & Summerfield, a Q. (2000). The lateralization of simple dichotic pitches. J. Acoust. Soc. Am., 108(1), 316–334.](http://dx.doi.org/10.1121/1.429467)
* [Zhang, P. X., & Hartmann, W. M. (2008). Lateralization of Huggins pitch. J. Acoust. Soc. Am., 124(6), 3873–87.](http://dx.doi.org/10.1121/1.2977683)

##### Examples

```julia
hp = hugginsPitch(F0=300, lowHarm=1, highHarm=3, spectrumLevel=45,
bandwidth=100, bandwidthUnit="Hz", dichoticDifference="IPD stepped",
dichoticDifferenceValue=pi, phaseRelationship="NoSpi", stretch=0,
noiseType="white", dur=0.4, rampDur=0.01, sf=48000, maxLevel=101)
```



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
pt = pureToneITD(frequency=440, phase=0, ITD=0.004/10, level=65, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```

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

Compute the root mean square (RMS) value of the signal.

##### Parameters

* `sig`: The signal for which the RMS needs to be computed.
* `channel`: Either an integer indicating the channel number,
  `each` for a list of the RMS values in each channel, or `all`
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

Simple function to play sounds. Uses aplay on Linux and afplay on OSX.
Windows is not currently supported.

##### Arguments

* `snd`: The sound to be played.
* `sf`: The sampling frequency.


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


