## The MIT License (MIT)

## Copyright (c) 2013-2017 Samuele Carcagno <sam.carcagno@gmail.com>

## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:

## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.

##using DocStringExtensions
#############################
## AMTone
#############################
"""
Generate an amplitude modulated tone.

$(SIGNATURES)

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
"""

function AMTone(;carrierFreq::Real=1000, AMFreq::Real=20, AMDepth::Real=1,
                carrierPhase::Real=0, AMPhase::Real=1.5*pi, level::Real=60,
                dur::Real=1, rampDur::Real=0.01, channel::AbstractString="diotic",
                sf::Real=48000, maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["mono", "right", "left", "diotic"]) == false
        error("`channel` must be one of 'mono', 'right', 'left', or 'diotic'")
    end
    amp = 10^((level - maxLevel) / 20)

    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) / sf
    timeRamp = collect(0:nRamp-1)

    snd_mono = zeros(nTot, 1)
    snd_mono[1:nRamp, 1] = amp * (1+AMDepth*sin(2*pi*AMFreq*timeAll[1:nRamp]+AMPhase)) .* ((1-cos(pi* timeRamp/nRamp))/2) .* sin(2*pi*carrierFreq * timeAll[1:nRamp] + carrierPhase)
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp * (1 + AMDepth*sin(2*pi*AMFreq*timeAll[nRamp+1:nRamp+nSamples]+AMPhase)) .* sin(2*pi*carrierFreq * timeAll[nRamp+1:nRamp+nSamples] + carrierPhase)
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp * (1 + AMDepth*sin(2*pi*AMFreq*timeAll[nRamp+nSamples+1:nTot]+AMPhase)) .* ((1+cos(pi * timeRamp/nRamp))/2) .* sin(2*pi*carrierFreq * timeAll[nRamp+nSamples+1:nTot] + carrierPhase)

    if channel == "mono"
        snd = snd_mono
    else
        snd = zeros(nTot, 2)
    end

    if channel == "right"
        snd[:,2] = snd_mono
    elseif channel == "left"
        snd[:,1] = snd_mono
    elseif channel == "diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    end

    return snd
end

##########################
## AMToneIPD
##########################
"""
Generate an amplitude modulated tone with an IPD in the carrier and/or
modulation phase.

$(SIGNATURES)

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
"""

function AMToneIPD(;carrierFreq::Real=1000, AMFreq::Real=20, AMDepth::Real=1,
                carrierPhase::Real=0, AMPhase::Real=1.5*pi, carrierIPD::Real=0,
                AMIPD::Real=0, level::Real=60,
                dur::Real=1, rampDur::Real=0.01, channel::AbstractString="right",
                sf::Real=48000, maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["right", "left"]) == false
        error("`channel` must be either 'right', or 'left'")
    end

    fixed = AMTone(carrierFreq=carrierFreq, AMFreq=AMFreq, AMDepth=AMDepth,
                carrierPhase=carrierPhase, AMPhase=AMPhase, level=level,
                dur=dur, rampDur=rampDur, channel="mono",
                sf=sf, maxLevel=maxLevel)

    shifted = AMTone(carrierFreq=carrierFreq, AMFreq=AMFreq, AMDepth=AMDepth,
                   carrierPhase=carrierPhase+carrierIPD, AMPhase=AMPhase+AMIPD, level=level,
                   dur=dur, rampDur=rampDur, channel="mono",
                   sf=sf, maxLevel=maxLevel)
    snd = zeros(size(fixed)[1], 2)
    if channel == "right"
        snd[:,1] = fixed
        snd[:,2] = shifted
    elseif channel == "left"
        snd[:,1] = shifted
        snd[:,2] = fixed
    end
    return snd
end

#############################
## asynchChord
#############################
"""
Generate an asynchronous chord.

$(SIGNATURES)

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
"""

function asynchChord{T<:Real}(;freqs::AbstractVector{T}=[200, 400], levels::AbstractVector{T}=[60, 60], phases::AbstractVector{T}=[0, 0], tonesDur::Real=0.2, tonesRampDur::Real=0.01, tonesChannel::AbstractString="diotic", SOA::Real=0.06, sf::Real=48000, maxLevel::Real=101)

    seq = collect(1:length(freqs))
    shuffle!(seq)
    i = 1
    thisFreq = freqs[seq[i]]; thisLev = levels[seq[i]]; thisPhase = phases[seq[i]]
    snd = pureTone(frequency=thisFreq, phase=thisPhase,
                   level=thisLev, dur=tonesDur, rampDur=tonesRampDur,
                   channel=tonesChannel, sf=sf, maxLevel=maxLevel)

    for i=2:length(seq)
        thisFreq = freqs[seq[i]]; thisLev = levels[seq[i]]; thisPhase = phases[seq[i]]
        thisTone = pureTone(frequency=thisFreq, phase=thisPhase,
                            level=thisLev, dur=tonesDur,
                            rampDur=tonesRampDur,
                            channel=tonesChannel, sf=sf, maxLevel=maxLevel)

        snd = addSounds(snd, thisTone, delay=SOA*i, sf=sf)
    end
    return snd
end

#############################
## broadbandNoise
#############################
"""
Synthetise a broadband noise.

$(SIGNATURES)

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

"""
function broadbandNoise(;spectrumLevel::Real=20, dur::Real=1, rampDur::Real=0.01,
                        channel::AbstractString="diotic", sf::Real=48000, maxLevel::Real=101)

    ## Comments:
    ## The intensity spectrum level in dB is SL
    ## The peak amplitude A to achieve a desired SL is
    ## SL = 10*log10(RMS^2/NHz) that is the total RMS^2 divided by the freq band
    ## SL/10 = log10(RMS^2/NHz)
    ## 10^(SL/10) = RMS^2/NHz
    ## RMS^2 = 10^(SL/10) * NHz
    ## RMS = 10^(SL/20) * sqrt(NHz)
    ## NHz = sampRate / 2 (Nyquist)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["mono", "right", "left", "diotic", "dichotic"]) == false
        error("`channel` must be one of 'mono', 'right', 'left', 'diotic', or 'dichotic'")
    end

    amp = sqrt(sf/2)*(10^((spectrumLevel - maxLevel) / 20))
    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) ./ sf
    timeRamp = collect(0:nRamp-1)

    snd_mono = zeros(nTot, 1)
    noise = (rand(nTot) + rand(nTot)) - (rand(nTot) + rand(nTot))
    RMS = sqrt(mean(noise.*noise))
    #scale the noise so that the maxAmplitude goes from -1 to 1
    #since A = RMS*sqrt(2)
    scaled_noise = noise / (RMS * sqrt(2))

    snd_mono[1:nRamp, 1] = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* scaled_noise[1:nRamp]
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp * scaled_noise[nRamp+1:nRamp+nSamples]
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* scaled_noise[nRamp+nSamples+1:nTot]

    if channel == "mono"
        snd = snd_mono
    else
        snd = zeros(nTot, 2)
    end

    if channel == "right"
        snd[:,2] = snd_mono
    elseif channel == "left"
        snd[:,1] = snd_mono
    elseif channel == "diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    elseif channel == "dichotic"
        snd_mono2 = zeros(nTot, 1)
        noise2 = (rand(nTot) + rand(nTot)) - (rand(nTot) + rand(nTot))
        RMS2 = sqrt(mean(noise2.*noise2))
        #scale the noise so that the maxAmplitude goes from -1 to 1
        #since A = RMS*sqrt(2)
        scaled_noise2 = noise2 / (RMS2 * sqrt(2))

        snd_mono2[1:nRamp, 1] = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* scaled_noise2[1:nRamp]
        snd_mono2[nRamp+1:nRamp+nSamples, 1] = amp * scaled_noise2[nRamp+1:nRamp+nSamples]
        snd_mono2[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* scaled_noise2[nRamp+nSamples+1:nTot]
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono2

    end

    return snd
end

############################
## complexTone
############################
"""
Synthetise a complex tone.

$(SIGNATURES)

##### Parameters:

* `F0`: Tone fundamental frequency in hertz.
* `harmPhase : one of `sine`, `cosine`, `alternating`, `random`, `schroeder-`, `schroeder+`
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

"""
function complexTone(;F0::Real=220, harmPhase::AbstractString="sine", lowHarm::Integer=1, highHarm::Integer=10, stretch::Real=0,
                     level::Real=60, dur::Real=1, rampDur::Real=0.01, channel::AbstractString="diotic", sf::Real=48000,
                     maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["mono", "right", "left", "diotic", "odd left", "odd right"]) == false
        error("`channel` must be one of 'mono', 'right', 'left', 'diotic', 'odd left', 'odd right'")
    end
    if in(harmPhase, ["sine", "cosine", "alternating", "random", "schroeder-", "schroeder+"]) == false
        error("`harmPhase` must be one of 'sine', 'cosine', 'alternating', 'random', 'schroeder-', 'schroeder+'")
    end

    amp = 10^((level - maxLevel) / 20)
    stretchHz = (F0*stretch)/100

    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) / sf
    timeRamp = collect(0:nRamp-1)

    if channel == "mono"
        snd_mono = zeros(nTot, 1)
    else
        snd = zeros(nTot, 2)
    end

    if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
        tone = zeros(nTot)
    elseif channel == "odd left" || channel == "odd right"
        toneOdd = zeros(nTot)
        toneEven = zeros(nTot)
    end


    if harmPhase == "sine"
        for i=lowHarm:highHarm
            if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
                tone =  tone + sin(2 * pi * ((F0 * i) + stretchHz) * timeAll)
            elseif channel == "odd left" || channel == "odd right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                else
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                end
            end
        end
    elseif harmPhase == "cosine"
        for i=lowHarm:highHarm
            if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
                tone = tone + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
            elseif channel == "odd left" || channel == "odd right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                else
                    toneEven = toneEven + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                end
            end
        end
    elseif harmPhase == "alternating"
        for i=lowHarm:highHarm
            if i%2 > 0 #odd harmonic
                if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
                    tone = tone + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                elseif channel == "odd left" || channel == "odd right"
                    toneOdd = toneOdd + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                end
            else #even harmonic
                if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
                    tone = tone + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                elseif channel == "odd left" || channel == "odd right"
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                end
            end
        end
    elseif harmPhase == "schroeder-"
        for i=lowHarm:highHarm
            phase = -pi * i * (i - 1) / (highHarm-lowHarm+1)
            if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
                tone = tone + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
            elseif channel == "odd left" || channel == "odd right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                else
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                end
            end
        end
    elseif harmPhase == "schroeder+"
        for i=lowHarm:highHarm
            phase = pi * i * (i - 1) / (highHarm-lowHarm+1)
            if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
                tone = tone + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
            elseif channel == "odd left" || channel == "odd right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                else
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                end
            end
        end
    elseif harmPhase == "random"
        for i=lowHarm:highHarm
            phase = rand() * 2 * pi
            if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
                tone = tone + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
            elseif channel == "odd left" || channel == "odd right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                else
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                end
            end
        end
    end

    if channel == "mono" || channel == "right" || channel == "left" || channel == "diotic"
        snd_mono = zeros(nTot, 1)
        snd_mono[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* tone[1:nRamp]
        snd_mono[nRamp+1:nRamp+nSamples, 1]        = amp * tone[nRamp+1:nRamp+nSamples]
        snd_mono[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* tone[nRamp+nSamples+1:nTot]
    end
    if channel == "mono"
        snd = snd_mono
    elseif channel == "right"
        snd[:,2] = snd_mono
    elseif channel == "left"
        snd[:,1] = snd_mono
    elseif channel == "diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    elseif channel == "odd left"
        snd[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* toneOdd[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1]        = amp * toneOdd[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* toneOdd[nRamp+nSamples+1:nTot]
        snd[1:nRamp, 2]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* toneEven[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 2]        = amp * toneEven[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* toneEven[nRamp+nSamples+1:nTot]
    elseif channel == "odd right"
        snd[1:nRamp, 2]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* toneOdd[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 2]        = amp * toneOdd[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* toneOdd[nRamp+nSamples+1:nTot]
        snd[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* toneEven[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1]        = amp * toneEven[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* toneEven[nRamp+nSamples+1:nTot]
    end

    return snd
end

##############################
## expAMNoise
##############################
"""
Generate a sinusoidally amplitude-modulated noise with an exponentially
modulated AM frequency.

$(SIGNATURES)

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

"""
function expAMNoise(;carrierFreq::Real=150, MF::Real=2.5, deltaCents::Real=600, FMPhase::Real=pi, AMDepth::Real=1, spectrumLevel::Real=30, dur::Real=0.4, rampDur::Real=0.01, channel::AbstractString="diotic", sf::Real=48000, maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["mono", "right", "left", "diotic", "dichotic"]) == false
        error("`channel` must be one of 'mono', 'right', 'left', 'diotic', or 'dichotic'")
    end

    amp = sqrt(sf/2)*(10^((spectrumLevel - maxLevel) / 20))
    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) / sf
    timeRamp = collect(0:nRamp-1)
    noise = (rand(nTot) + rand(nTot)) - (rand(nTot) + rand(nTot))
    RMS = sqrt(mean(noise.*noise))
    #scale the noise so that the maxAmplitude goes from -1 to 1
    #since A = RMS*sqrt(2)
    scaled_noise = noise / (RMS * sqrt(2))

    if channel == "dichotic"
        noise2 = (rand(nTot) + rand(nTot)) - (rand(nTot) + rand(nTot))
        RMS2 = sqrt(mean(noise.*noise))
        scaled_noise2 = noise2 / (RMS2 * sqrt(2))
    end

    fArr = 2*pi*carrierFreq*2.^((deltaCents/1200)*cos(2*pi*MF*timeAll+FMPhase))
    ang = (cumsum(fArr)/sf)
    snd_mono = zeros(nTot, 1)
    snd_mono[1:nRamp, 1] = amp * (1 + AMDepth*sin(ang[1:nRamp])) .* ((1-cos(pi * timeRamp/nRamp))/2) .* scaled_noise[1:nRamp]
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp * (1 + AMDepth*sin(ang[nRamp+1:nRamp+nSamples])) .* scaled_noise[nRamp+1:nRamp+nSamples]
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp * (1 + AMDepth*sin(ang[nRamp+nSamples+1:nTot])) .* ((1+cos(pi * timeRamp/nRamp))/2) .* scaled_noise[nRamp+nSamples+1:nTot]

    if channel == "mono"
        snd = snd_mono
    else
        snd = zeros(nTot, 2)
    end
    if channel == "right"
        snd[:,2] = snd_mono
    elseif channel == "left"
        snd[:,1] = snd_mono
    elseif channel == "diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    elseif channel == "dichotic"
        snd[:,2] = snd_mono
        snd[1:nRamp, 1] = amp * (1 + AMDepth*sin(ang[1:nRamp])) .* ((1-cos(pi * timeRamp/nRamp))/2) .* scaled_noise2[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1] = amp * (1 + AMDepth*sin(ang[nRamp+1:nRamp+nSamples])) .* scaled_noise2[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * (1 + AMDepth*sin(ang[nRamp+nSamples+1:nTot])) .* ((1+cos(pi * timeRamp/nRamp))/2) .* scaled_noise2[nRamp+nSamples+1:nTot]
    end

    return snd
end

###################################
## expSinFMTone
###################################
"""
Generate a tone frequency modulated with an exponential sinusoid.

$(SIGNATURES)

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
"""

function expSinFMTone(;carrierFreq::Real=450, MF::Real=5, deltaCents::Real=600, FMPhase::Real=pi, phase::Real=0, level::Real=60, dur::Real=0.2, rampDur::Real=0.01, channel::AbstractString="diotic", sf::Real=48000, maxLevel::Real=101)
    amp = 10^((level - maxLevel) / 20)
    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) / sf
    timeRamp = collect(0:nRamp-1)
    fArr = 2*pi*carrierFreq*2.^((deltaCents/1200)*cos(2*pi*MF*timeAll+FMPhase))
    ang = (cumsum(fArr)/sf) + phase

    snd_mono = zeros(nTot, 1)
    snd_mono[1:nRamp, 1] = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* sin(ang[1:nRamp])
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp * sin(ang[nRamp+1:nRamp+nSamples])
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* sin(ang[nRamp+nSamples+1:nTot])
    if channel == "mono"
        snd = zeros(nTot, 1)
    else
        snd = zeros(nTot, 2)
    end
    if channel == "mono"
        snd[:,1] = snd_mono
    elseif channel == "right"
        snd[:,2] = snd_mono
    elseif channel == "left"
        snd[:,1] = snd_mono
    elseif channel == "diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    end
    return snd
end

####################################
## FMComplex2
####################################
"""
Synthetise a complex tone with an embedded frequency modulation (FM)
starting and stopping at a chosen time after the tone onset.

$(SIGNATURES)

##### Parameters

* `midF0`: F0 at the FM zero crossing
* `harmPhase`: one of 'sine', 'cosine', 'alternating', 'random', 'schroeder-', 'schroeder+'.
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
"""

function FMComplex2(;midF0::Real=140, harmPhase::AbstractString="sine",
                    lowHarm::Integer=1, highHarm::Integer=10,
                    level::Real=60, dur::Real=0.45, rampDur::Real=0.01,
                    MF::Real=1.25, FMDepth::Real=40,
                    FMStartPhase::Real=1.5*pi,
                    FMStartTime::Real=0.025, FMDur::Real=0.4,
                    levelAdj::Bool=true, channel::AbstractString="diotic",
                    sf::Real=48000, maxLevel::Real=101)


    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if FMDur+FMStartTime > dur
        error("`FMDur`+`FMStartTime` cannot be greater than sound duration")
    end
    if in(channel, ["mono", "right", "left", "diotic", "odd left", "odd right"]) == false
        error("`channel` must be one of 'mono', 'right', 'left', 'diotic', 'odd left', 'odd right'")
    end
    if in(harmPhase, ["sine", "cosine", "alternating", "random", "schroeder-", "schroeder+"]) == false
        error("`harmPhase` must be one of 'sine', 'cosine', 'alternating', 'random', 'schroeder-', 'schroeder+'")
    end


    amp = 10^((level - maxLevel) / 20)

    fmStartPnt = round(Int, FMStartTime*sf) #sample where FM starts
    nFMSamples = round(Int, FMDur*sf) #number of FM samples
    nSamples = round(Int, (dur-rampDur*2)* sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) / sf
    timeRamp = collect(0:nRamp-1)

    timeAllSamp = collect(0:nTot-1) #time array not scaled by sf
    time1 = timeAllSamp[1:fmStartPnt]
    time2 = timeAllSamp[fmStartPnt+1:fmStartPnt+nFMSamples]
    time3 = timeAllSamp[fmStartPnt+nFMSamples+1:nTot]
    fmTime = collect(0:nFMSamples-1)

    FMDepthHz = FMDepth*midF0/100 #convert from % to Hz
    B = FMDepthHz / MF #Beta
    fmStartDepth = FMDepthHz*sin(FMStartPhase)

    fmRadFreq = 2*pi*MF/sf

    midF0Rad = 2*pi*midF0/sf
    startF0Rad = 2*pi*(midF0 + FMDepthHz*sin(FMStartPhase))/sf
    endF0Rad = 2*pi*(midF0 + FMDepthHz*sin(FMStartPhase + nFMSamples*fmRadFreq)) / sf

    startF0 = midF0 + FMDepthHz*sin(FMStartPhase)
    endF0 = midF0 + FMDepthHz*sin(FMStartPhase + nFMSamples*fmRadFreq)

    if in(channel, ["mono", "right", "left", "diotic"])
        tone = zeros(nTot)
    elseif in(channel, ["odd left", "odd right"])
        toneOdd = zeros(nTot)
        toneEven = zeros(nTot)
    end

    #from Hartmann, WM (1997) Signals, sound, and sensation. New York: AIP Press
    #angular frequency is the time derivative of the instantaneous phase
    #if the angular frequency is given by a constant carrier `wc`, plus a
    #sinusoidal deviation `dw`:
    #eq.1: w(t) = wc + dw*cos(wm*t+phi)
    #where `wm` is the modulation frequency, then the instantaneous phase is
    #given by the integral of the above expression with respect to time
    #eq.2: PHI(t) = wc*t + (dw/wm)*sin(wm*t+phi)
    #that's the canonical form of the FM equation
    #if instead of modulating the angular freq. in eq. 1 by `cos` we modulate it by `sin`:
    #eq. 3: w(t) = wc + dw*(cos(wm*tphi)
    #and the resulting integral is:
    #eq.4: PHI(t) = wc*t - (dw/wm)*cos(wm*t+phi)
    #this is what we're actually using below

    for i=lowHarm:highHarm
        fArr = zeros(nTot)
        fArr[1:fmStartPnt] = startF0*i
        fArr[fmStartPnt+1:fmStartPnt+nFMSamples] = (midF0*i + FMDepthHz*i*sin(2*pi*MF*fmTime/sf+FMStartPhase))
        fArr[fmStartPnt+nFMSamples+1:nTot] = endF0*i
        if harmPhase == "sine"
            ang = cumsum(2*pi*fArr/sf)
            if in(channel, ["mono", "right", "left", "diotic"])
                tone = tone + sin(ang)
            elseif in(channel, ["odd left", "odd right"])
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(ang)
                else
                    toneEven = toneEven + sin(ang)
                end
            end
        elseif harmPhase == "cosine"
            ang = cumsum(2*pi*fArr/sf)
            if in(channel, ["mono", "right", "left", "diotic"])
                tone = tone + cos(ang)
            elseif in(channel, ["odd left", "odd right"])
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + cos(ang)
                else
                    toneEven = toneEven + cos(ang)
                end
            end
        elseif harmPhase == "alternating"
            ang = cumsum(2*pi*fArr/sf)
            if i%2 > 0 #odd harmonic
                if in(channel, ["mono", "right", "left", "diotic"])
                    tone = tone + cos(ang)
                elseif in(channel, ["odd left", "odd right"])
                    toneOdd = toneOdd + cos(ang)
                end
            else #even harmonic
                if in(channel, ["mono", "right", "left", "diotic"])
                    tone = tone + sin(ang)
                elseif in(channel, ["odd left", "odd right"])
                    toneEven = toneEven + sin(ang)
                end
            end
        elseif harmPhase == "schroeder-"
            phase = -pi * i * (i - 1) / (highHarm-lowHarm+1)
            ang = cumsum(2*pi*fArr/sf + phase)
            if in(channel, ["mono", "right", "left", "diotic"])
                tone = tone + sin(ang)
            elseif in(channel, ["odd left", "odd right"])
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(ang)
                else
                    toneEven = toneEven + sin(ang)
                end
            end
        elseif harmPhase == "schroeder+"
            phase = pi * i * (i - 1) / (highHarm-lowHarm+1)
            ang = cumsum(2*pi*fArr/sf + phase)
            if in(channel, ["mono", "right", "left", "diotic"])
                tone = tone + sin(ang)
            elseif in(channel, ["odd left", "odd right"])
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(ang)
                else
                    toneEven = toneEven + sin(ang)
                end
            end
        elseif harmPhase == "random"
            phase = rand() * 2 * pi
            ang = cumsum(2*pi*fArr/sf + phase)
            if in(channel, ["mono", "right", "left", "diotic"])
                tone = tone + sin(ang)
            elseif in(channel, ["odd left", "odd right"])
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(ang)
                else
                    toneEven = toneEven + sin(ang)
                end
            end
        end
    end
    #level correction --------------
    if levelAdj == true
        if in(channel, ["mono", "right", "left", "diotic"])
            tone[1:fmStartPnt] = tone[1:fmStartPnt] .* sqrt(((startF0Rad / (2*pi)) * (sf))/ midF0)
            tone[fmStartPnt+1:fmStartPnt+nFMSamples] = tone[fmStartPnt+1:fmStartPnt+nFMSamples] .* sqrt((midF0 + (FMDepthHz * sin(FMStartPhase + (fmTime * fmRadFreq)))) / midF0)
            tone[fmStartPnt+nFMSamples+1:nTot] = tone[fmStartPnt+nFMSamples+1:nTot] .* sqrt( ((endF0Rad / (2*pi)) * (sf))/ midF0)
        elseif in(channel, ["odd left", "odd right"])
            toneEven[1:fmStartPnt] = toneEven[1:fmStartPnt] * sqrt(((startF0Rad / (2*pi)) * (sf))/ midF0)
            toneEven[fmStartPnt+1:fmStartPnt+nFMSamples] = toneEven[fmStartPnt+1:fmStartPnt+nFMSamples]  .* sqrt((midF0 + (FMDepthHz * sin(FMStartPhase + (fmTime * fmRadFreq)))) / midF0)
            toneEven[fmStartPnt+nFMSamples+1:nTot] = toneEven[fmStartPnt+nFMSamples+1:nTot] .* sqrt(((endF0Rad / (2*pi)) * (sf))/ midF0)

            toneOdd[1:fmStartPnt] = toneOdd[1:fmStartPnt] * sqrt(((startF0Rad / (2*pi)) .* (sf))/ midF0)
            toneOdd[fmStartPnt+1:fmStartPnt+nFMSamples] = toneOdd[fmStartPnt+1:fmStartPnt+nFMSamples] .* sqrt((midF0 + (FMDepthHz * sin(FMStartPhase + (fmTime * fmRadFreq)))) / midF0)
            toneOdd[fmStartPnt+nFMSamples+1:nTot] = toneOdd[fmStartPnt+nFMSamples+1:nTot] .* sqrt(((endF0Rad / (2*pi)) * (sf))/ midF0)
        end
    end
#end of level correction -----------
    if channel == "mono"
        snd = zeros(nTot, 1)
    else
        snd = zeros(nTot, 2)
    end

    if in(channel, ["odd right", "odd left"]) == false
        tone[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* tone[1:nRamp]
        tone[nRamp+1:nRamp+nSamples, 1]        = amp * tone[nRamp+1:nRamp+nSamples]
        tone[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* tone[nRamp+nSamples+1:nTot]
    else
        toneOdd[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* toneOdd[1:nRamp]
        toneOdd[nRamp+1:nRamp+nSamples, 1]        = amp * toneOdd[nRamp+1:nRamp+nSamples]
        toneOdd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* toneOdd[nRamp+nSamples+1:nTot]
        toneEven[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* toneEven[1:nRamp]
        toneEven[nRamp+1:nRamp+nSamples, 1]        = amp * toneEven[nRamp+1:nRamp+nSamples]
        toneEven[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* toneEven[nRamp+nSamples+1:nTot]
    end

    if channel == "mono"
        snd[:,1] = tone
    elseif channel == "right"
        snd[:,2] = tone
    elseif channel == "left"
        snd[:,1] = tone
    elseif channel == "diotic"
        snd[:,1] = tone
        snd[:,2] = tone
    elseif channel == "odd left"
        snd[:,1] = toneOdd
        snd[:,2] = toneEven
    elseif channel == "odd right"
        snd[:,2] = toneOdd
        snd[:,1] = toneEven
    end

    return snd
end

############################
## FMTone
############################
"""
Generate a frequency modulated tone.

$(SIGNATURES)

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
"""

function FMTone(;carrierFreq::Real=1000, MF::Real=40, MI::Real=1, phase::Real=0,
                level::Real=60, dur::Real=1, rampDur::Real=0.01,
                channel::AbstractString="diotic", sf::Real=48000, maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["mono", "right", "left", "diotic"]) == false
        error("`channel` must be one of 'mono', 'right', 'left', or 'diotic'")
    end

    amp = 10^((level - maxLevel) / 20)
    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) / sf
    timeRamp = collect(0:nRamp-1)

    snd_mono = zeros(nTot, 1)
    snd_mono[1:nRamp, 1] = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* sin(2*pi*carrierFreq*timeAll[1:nRamp] + MI*sin(2*pi*MF * timeAll[1:nRamp] + phase))
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp * sin(2*pi*carrierFreq * timeAll[nRamp+1:nRamp+nSamples] +MI * sin(2*pi*MF * timeAll[nRamp+1:nRamp+nSamples] + phase))
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* sin(2*pi*carrierFreq * timeAll[nRamp+nSamples+1:nTot]+MI*sin(2*pi*MF * timeAll[nRamp+nSamples+1:nTot] + phase))

    if channel == "mono"
        snd = zeros(nTot, 1)
    else
        snd = zeros(nTot, 2)
    end

    if channel == "mono"
        snd[:,1] = snd_mono
    elseif channel == "right"
        snd[:,2] = snd_mono
    elseif channel == "left"
        snd[:,1] = snd_mono
    elseif channel == "diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    end

    return snd

end

################################
## hugginsPitch
################################

"""
Synthetise a complex Huggings Pitch.

$(SIGNATURES)

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

"""

function hugginsPitch(;F0::Real=550, lowHarm::Int=1, highHarm::Int=1,
                      spectrumLevel::Real=30, bandwidth::Real=1,
                      bandwidthUnit::AbstractString="ERB", dichoticDifference::AbstractString="IPD stepped",
                      dichoticDifferenceValue::Real=pi, phaseRelationship::AbstractString="NoSpi",
                      stretch::Real=0, noiseType::AbstractString="white", dur::Real=1, rampDur::Real=0.01,
                      sf::Real=48000, maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(bandwidthUnit, ["Hz", "Cent", "ERB"]) == false
        error("`bandwidthUnit` must be one of 'Hz', 'Cent', or 'ERB'")
    end

    if in(dichoticDifference, ["IPD linear", "IPD stepped", "IPD random", "ITD"]) == false
        error("`bandwidthUnit` must be one of 'IPD linear', 'IPD stepped', 'IPD random', or 'ITD'")
    end

    if in(phaseRelationship, ["NoSpi", "NpiSo"]) == false
        error("`phaseRelationship` must be either 'NoSpi', or 'NpiSo'")
    end

    if in(noiseType, ["white", "pink"]) == false
        error("`noiseType` must be either 'white', or 'pink'")
    end

    stretchHz = (F0*stretch)/100
    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)
    snd = zeros(nTot, 2)

    noise = broadbandNoise(spectrumLevel=spectrumLevel, dur=dur, rampDur=0,
                           channel="diotic", sf=sf, maxLevel=maxLevel)
    if noiseType == "pink"
        makePink!(noise, sf=sf, ref=1000)
    end

    cfs = collect(lowHarm:highHarm)*F0 #center frequencies
    cfs = cfs + stretchHz
    if phaseRelationship == "NoSpi"
        if bandwidthUnit == "Hz"
            shiftLo = cfs - (bandwidth/2)
            shiftHi = cfs + (bandwidth/2)
        elseif bandwidthUnit == "Cent"
            shiftLo = cfs*2.^(-(bandwidth/2)/1200)
            shiftHi = cfs*2.^((bandwidth/2)/1200)
        elseif bandwidthUnit == "ERB"
            shiftLo = freqFromERBInterval(cfs, -bandwidth/2)
            shiftHi = freqFromERBInterval(cfs, bandwidth/2)
        end
    end
    if phaseRelationship == "NpiSo"
        nHarms = length(cfs)
        shiftLo = zeros(nHarms+1)
        shiftHi = zeros(nHarms+1)

        shiftLo[1] = 10
        shiftHi[end] = sf/2
        if bandwidthUnit == "Hz"
            shiftLo[2:length(shiftLo)] = cfs + (bandwidth/2)
            shiftHi[1:length(shiftHi)-1] = cfs - (bandwidth/2)
        elseif bandwidthUnit == "Cent"
            shiftLo[2:length(shiftLo)] = cfs*2.^((bandwidth/2)/1200)
            shiftHi[1:length(shiftHi)-1] = cfs*2.^((bandwidth/2)/1200)
        elseif bandwidthUnit == "ERB"
            shiftLo[2:length(shiftLo)] = freqFromERBInterval(cfs, bandwidth/2)
            shiftHi[1:length(shiftHi)-1] = freqFromERBInterval(cfs, -bandwidth/2)
        end
    end

    for i=1:length(shiftLo)
        if dichoticDifference == "IPD linear"
            noise = phaseShift!(noise, shiftLo[i], shiftHi[i],
                                phaseShift=dichoticDifferenceValue, shiftType="linear", channel=1, sf=sf)
        elseif dichoticDifference == "IPD stepped"
            noise = phaseShift!(noise, shiftLo[i], shiftHi[i],
                               phaseShift=dichoticDifferenceValue, shiftType="step", channel=1, sf=sf)
        elseif dichoticDifference == "IPD random"
            noise = phaseShift!(noise, shiftLo[i], shiftHi[i],
                               phaseShift=dichoticDifferenceValue, shiftType="random", channel=1, sf=sf)
        elseif dichoticDifference == "ITD"
            noise = ITDShift!(noise, shiftLo[i], shiftHi[i],
                             ITD=dichoticDifferenceValue, channel="left", sf=sf)
        end
    end
    noise = gate!(noise, rampDur=rampDur, sf=sf)

    return noise

end

#################################
## IRN
#################################

"""
Generate an iterated rippled noise

$(SIGNATURES)

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
"""

function IRN(;delay::Real=0.001, gain::Real=1, iterations::Integer=6,
             configuration::AbstractString="add same", spectrumLevel::Real=25,
             dur::Real=1, rampDur::Real=0.01, channel::AbstractString="diotic",
             sf::Real=48000, maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["mono", "right", "left", "diotic", "dichotic"]) == false
        error("`channel` must be one of `mono`, `right`, `left`, `diotic`, or `dichotic`")
    end
    if in(configuration, ["add same", "add original"]) == false
        error("`configuration` must be either 'add same', or 'add original'")
    end

    snd = broadbandNoise(spectrumLevel=spectrumLevel, dur=dur, rampDur=0,
                         channel=channel, sf=sf, maxLevel=maxLevel)
    if channel == "right"
        ch = 2
    elseif channel == "left"
        ch = 1
    else
        ch = collect(1:size(snd)[2])
    end
    snd = delayAdd!(snd, delay=delay, gain=gain, iterations=iterations,
                   configuration=configuration, channel=ch, sf=sf)
    snd = gate!(snd, rampDur=rampDur, sf=sf)

    return snd
end


####################################
## pureTone
####################################
"""
Synthetise a pure tone.

$(SIGNATURES)

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone will be generated.  (`mono`, `right`, `left` or `diotic`)
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
"""
function pureTone(;frequency::Real=1000, phase::Real=0, level::Real=65,
                  dur::Real=1, rampDur::Real=0.01,
                  channel::AbstractString="diotic", sf::Real=48000,
                  maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["mono", "right", "left", "diotic"]) == false
        error("`channel` must be one of 'mono', 'right', 'left', or 'diotic'")
    end

    amp = 10^((level - maxLevel) / 20)
    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) ./ sf
    timeRamp = collect(0:nRamp-1)

    if channel == "mono"
        snd = zeros(nTot, 1)
    else
        snd = zeros(nTot, 2)
    end

    snd_mono = zeros(nTot, 1)
    snd_mono[1:nRamp, 1] = amp * ((1 - cos(pi * timeRamp/nRamp))/2) .* sin(2*pi*frequency * timeAll[1:nRamp] + phase)
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp * sin(2*pi*frequency * timeAll[nRamp+1:nRamp+nSamples] + phase)
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp * ((1 + cos(pi * timeRamp/nRamp))/2) .* sin(2*pi*frequency * timeAll[nRamp+nSamples+1:nTot] + phase)

    if channel == "mono"
        snd[:,1] = snd_mono
    elseif channel == "right"
        snd[:,2] = snd_mono
    elseif channel == "left"
        snd[:,1] = snd_mono
    elseif channel == "diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    end

    return snd
end

##################################
## pureToneILD
##################################
"""
Synthetise a pure tone with an interaural level difference.

$(SIGNATURES)

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
"""

function pureToneILD(;frequency::Real=1000, phase::Real=0,
                     level::Real=65, ILD::Real=10, dur::Real=1, rampDur::Real=0.01,
                     channel::AbstractString="right", sf::Real=48000, maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["right", "left"]) == false
        error("`channel` must be one of 'right', or 'left'")
    end


    fixed = pureTone(frequency=frequency, phase=phase, level=level, dur=dur, rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
    shifted = pureTone(frequency=frequency, phase=phase, level=level+ILD, dur=dur, rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
    snd = zeros(size(fixed)[1], 2)

    if channel == "right"
        snd[:,1] = fixed
        snd[:,2] = shifted
    elseif channel == "left"
        snd[:,1] = shifted
        snd[:,2] = fixed
    end

    return snd
end

######################################
## pureToneIPD
######################################
"""
Synthetise a pure tone with an interaural phase difference.

$(SIGNATURES)

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
"""

function pureToneIPD(;frequency::Real=1000, phase::Real=0, IPD::Real=pi,
                     level::Real=65, dur::Real=1, rampDur::Real=0.01,
                     channel::AbstractString="right", sf::Real=48000, maxLevel::Real=101)
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["right", "left"]) == false
        error("`channel` must be one of 'right', or 'left'")
    end


    fixed = pureTone(frequency=frequency, phase=phase, level=level, dur=dur,
                     rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
    shifted = pureTone(frequency=frequency, phase=phase+IPD, level=level,
                       dur=dur, rampDur=rampDur, channel="mono", sf=sf,
                       maxLevel=maxLevel)
    snd = zeros(size(fixed)[1], 2)
    if channel == "right"
        snd[:,1] = fixed
        snd[:,2] = shifted
    elseif channel == "left"
        snd[:,1] = shifted
        snd[:,2] = fixed
    end
    return snd
end

######################################
## pureToneITD
######################################
"""
Synthetise a pure tone with an interaural time difference.

$(SIGNATURES)

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
"""

function pureToneITD(;frequency::Real=1000, phase::Real=0, ITD::Real=0,
                     level::Real=65, dur::Real=1, rampDur::Real=0.01,
                     channel::AbstractString="right", sf::Real=48000, maxLevel::Real=101)
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["right", "left"]) == false
        error("`channel` must be one of 'right', or 'left'")
    end


    IPD = ITDToIPD(ITD, frequency)
    snd = pureToneIPD(frequency=frequency, phase=phase, IPD=IPD, level=level, dur=dur, rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)

    return snd
end

####################################
## pureToneIPDILD
####################################

"""
Synthetise a pure tone with an interaural phase and interaural level difference.

$(SIGNATURES)

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
"""

function pureToneIPDILD(;frequency::Real=1000, phase::Real=0, IPD::Real=0,
                     level::Real=65, ILD::Real=0, dur::Real=1, rampDur::Real=0.01,
                     channel::AbstractString="right", sf::Real=48000, maxLevel::Real=101)
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["right", "left"]) == false
        error("`channel` must be one of 'right', or 'left'")
    end

    fixed = pureTone(frequency=frequency, phase=phase, level=level, dur=dur, rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
    shifted = pureTone(frequency=frequency, phase=phase+IPD, level=level+ILD, dur=dur, rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
    snd = zeros(size(fixed)[1], 2)

    if channel == "right"
        snd[:,1] = fixed
        snd[:,2] = shifted
    elseif channel == "left"
        snd[:,1] = shifted
        snd[:,2] = fixed
    end

    return snd
end

#########################################
## pureToneITDILD
#########################################
"""
Synthetise a pure tone with an interaural time and interaural level difference.

$(SIGNATURES)

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
"""

function pureToneITDILD(;frequency::Real=1000, phase::Real=0, ITD::Real=0,
                     level::Real=65, ILD::Real=0, dur::Real=1, rampDur::Real=0.01,
                     channel::AbstractString="right", sf::Real=48000, maxLevel::Real=101)
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["right", "left"]) == false
        error("`channel` must be one of 'right', or 'left'")
    end


    IPD = ITDToIPD(ITD, frequency)
    snd = pureToneIPDILD(frequency=frequency, phase=phase, IPD=IPD, level=level,
                         ILD=ILD, dur=dur, rampDur=rampDur, channel=channel,
                         sf=sf, maxLevel=maxLevel)
    return snd
end

##################################
## silence
##################################
"""
Generate a silence.

$(SIGNATURES)

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
"""
function silence(;dur=1, channel="mono", sf=48000)
    nSamples = round(Int, dur * sf)
    if in(channel, ["mono", "diotic"]) == false
        error("`channel` must be one of 'mono', or 'diotic'")
    end

    if channel == "mono"
        snd = zeros(nSamples, 1)
    elseif channel == "diotic"
        snd = zeros(nSamples, 2)
    end
    return snd
end

#####################################
## steepNoise
#####################################
"""
Synthetise band-limited noise from the addition of random-phase
sinusoids.

$(SIGNATURES)

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

"""
function steepNoise(;f1=900, f2=1000, level=50, dur=1, rampDur=0.01,
                    channel="diotic", sf=48000, maxLevel=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end
    if in(channel, ["mono", "right", "left", "diotic", "dichotic"]) == false
        error("`channel` must be one of 'mono', 'right', 'left', 'diotic', or 'dichotic'")
    end


    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    spacing = 1 / dur
    components = 1 + floor((f2 - f1) / spacing)
    # SL = 10*log10(A^2/NHz)
    # SL/10 = log10(A^2/NHz)
    # 10^(SL/10) = A^2/NHz
    # A^2 = 10^(SL/10) * NHz
    # RMS = 10^(SL/20) * sqrt(NHz) where NHz is the spacing between harmonics
    amp =  10^((level - maxLevel) / 20) * sqrt((f2 - f1) / components)

    timeAll = collect(0:nTot-1) / sf
    timeRamp = collect(0:nRamp-1)

    if channel == "mono"
        snd = zeros(nTot, 1)
    else
        snd = zeros(nTot, 2)
    end

    noise= zeros(nTot)
    for f=f1:spacing:f2
        radFreq = 2 * pi * f
        phase = rand() * 2 * pi
        noise = noise + sin(phase + (radFreq * timeAll))
    end

    snd_mono = zeros(nTot, 1)
    snd_mono[1:nRamp, 1] = amp .* ((1-cos(pi * timeRamp/nRamp))/2) .* noise[1:nRamp]
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp * noise[nRamp+1:nRamp+nSamples]
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* noise[nRamp+nSamples+1:nTot]


    if channel == "right"
        snd[:,2] = snd_mono
    elseif channel == "left"
        snd[:,1] = snd_mono
    elseif channel == "diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    elseif channel == "dichotic"
        noise2= zeros(nTot)
        for f=f1:spacing:f2
            radFreq = 2 * pi * f
            phase = rand() * 2 * pi
            noise2 = noise2 + sin(phase + (radFreq * timeAll))
        end
        snd_mono2 = zeros(nTot, 1)
        snd_mono2[1:nRamp, 1] = amp .* ((1-cos(pi * timeRamp/nRamp))/2) .* noise2[1:nRamp]
        snd_mono2[nRamp+1:nRamp+nSamples, 1] = amp * noise2[nRamp+1:nRamp+nSamples]
        snd_mono2[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* noise2[nRamp+nSamples+1:nTot]
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono2
    end

    return snd
end

