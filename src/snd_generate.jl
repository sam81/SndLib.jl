## The MIT License (MIT)

## Copyright (c) 2013-2015 Samuele Carcagno <sam.carcagno@gmail.com>

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

#############################
## AMTone
#############################
@doc doc"""
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
""" ->

function AMTone(;carrierFreq::Real=1000, AMFreq::Real=20, AMDepth::Real=1,
                carrierPhase::Real=0, AMPhase::Real=1.5*pi, level::Real=60,
                dur::Real=1, rampDur::Real=0.01, channel::String="diotic",
                sf::Real=48000, maxLevel::Real=101)

    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
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
@doc doc"""
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
channel="diotic", sf=48000, maxLevel=100)
```
""" ->

function AMToneIPD(;carrierFreq::Real=1000, AMFreq::Real=20, AMDepth::Real=1,
                carrierPhase::Real=0, AMPhase::Real=1.5*pi, carrierIPD::Real=0,
                AMIPD::Real=0, level::Real=60,
                dur::Real=1, rampDur::Real=0.01, channel::String="right",
                sf::Real=48000, maxLevel::Real=101)

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
## broadbandNoise
#############################
@doc doc"""
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

""" ->
function broadbandNoise(;spectrumLevel::Real=20, dur::Real=1, rampDur::Real=0.01,
                        channel::String="diotic", sf::Real=48000, maxLevel::Real=101)

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
@doc doc"""
Synthetise a complex tone.

##### Parameters:

* `F0`: Tone fundamental frequency in hertz.
* `harmPhase : one of 'sine', 'cosine', 'alternating', 'random', 'schroeder'
        Phase relationship between the partials of the complex tone.
* `lowHarm`: Lowest harmonic component number.
* `highHarm`: Highest harmonic component number.
* `stretch`: Harmonic stretch in %F0. Increase each harmonic frequency by a fixed value
        that is equal to (F0*stretch)/100. If 'stretch' is different than
        zero, an inhanmonic complex tone will be generated.
* `level`: The level of each partial in dB SPL.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: 'right', 'left', 'diotic', 'odd right' or 'odd left'.
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

""" ->
function complexTone(;F0::Real=220, harmPhase::String="sine", lowHarm::Integer=1, highHarm::Integer=10, stretch::Real=0,
                     level::Real=60, dur::Real=1, rampDur::Real=0.01, channel::String="diotic", sf::Real=48000,
                     maxLevel::Real=101)
    
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
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
    elseif harmPhase == "schroeder"
        for i=lowHarm:highHarm
            phase = -pi * i * (i - 1) / highHarm
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

#################################
## IRN
#################################

@doc doc"""

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
* `channel`: Channel in which the noise will be generated (`mono`, 'right', 'left', 'diotic', 'dichotic').
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
""" ->

function IRN(;delay::Real=0.001, gain::Real=1, iterations::Integer=6,
             configuration::String="add same", spectrumLevel::Real=25,
             dur::Real=1, rampDur::Real=0.01, channel::String="diotic",
             sf::Real=48000, maxLevel::Real=101)

    if in(channel, ["mono", "right", "left", "diotic", "dichotic"]) == false
        error("`channel` must be one of `mono`, `right`, `left`, `diotic`, or `dichotic`")
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
@doc doc"""
Synthetise a pure tone.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone will be generated.  ('right', 'left' or 'diotic')
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
""" ->
function pureTone(;frequency::Real=1000, phase::Real=0, level::Real=65,
                  dur::Real=1, rampDur::Real=0.01,
                  channel::String="diotic", sf::Real=48000,
                  maxLevel::Real=101) 
    
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end

    amp = 10^((level - maxLevel) / 20)   
    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    timeAll = collect(0:nTot-1) ./ sf
    timeRamp = collect(0:nRamp-1)

    if channel == "mono"
        snd = zeros((nTot, 1))
    else
        snd = zeros((nTot, 2))
    end

    snd_mono = zeros(nTot, 1)
    snd_mono[1:nRamp, 1] = amp .* ((1.-cos(pi .* timeRamp/nRamp))./2) .* sin(2*pi*frequency .* timeAll[1:nRamp] .+ phase)
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp* sin(2*pi*frequency .* timeAll[nRamp+1:nRamp+nSamples] .+ phase)
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp .* ((1.+cos(pi .* timeRamp/nRamp))./2) .* sin(2*pi*frequency * timeAll[nRamp+nSamples+1:nTot] .+ phase)
    
    
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

##################################
## pureToneILD
##################################
@doc doc"""
Synthetise a pure tone with an interaural level difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `level`: Tone level in dB SPL.
* `ILD`: Interaural level difference in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the ILD will be applied.  ('right' or 'left')
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneILD(frequency=440, phase=0, level=65, ILD=10, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```
""" ->

function pureToneILD(;frequency::Real=1000, phase::Real=0,
                     level::Real=65, ILD::Real=0, dur::Real=1, rampDur::Real=0.01,
                     channel::String="right", sf::Real=48000, maxLevel::Real=101)

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
@doc doc"""
Synthetise a pure tone with an interaural phase difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `IPD`: Interaural phase difference in radians.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the tone IPD will be applied.  ('right' or 'left')
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneIPD(frequency=440, phase=0, IPD=pi/2, level=65, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```
""" ->

function pureToneIPD(;frequency::Real=1000, phase::Real=0, IPD::Real=0,
                     level::Real=65, dur::Real=1, rampDur::Real=0.01,
                     channel::String="right", sf::Real=48000, maxLevel::Real=101)
    fixed = pureTone(frequency=frequency, phase=phase, level=level, dur=dur, rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
    shifted = pureTone(frequency=frequency, phase=phase+IPD, level=level, dur=dur, rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
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
@doc doc"""
Synthetise a pure tone with an interaural time difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `ITD`: Interaural time difference in seconds.
* `level`: Tone level in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the ITD will be applied.  ('right' or 'left')
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneITD(frequency=440, phase=0, ITD=0.004/1000, level=65, dur=1,
rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```
""" ->

function pureToneITD(;frequency::Real=1000, phase::Real=0, ITD::Real=0,
                     level::Real=65, dur::Real=1, rampDur::Real=0.01,
                     channel::String="right", sf::Real=48000, maxLevel::Real=101)
    IPD = ITDToIPD(ITD, frequency)
    snd = pureToneIPD(frequency=frequency, phase=phase, IPD=IPD, level=level, dur=dur, rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)

    return snd
end

####################################
## pureToneIPDILD
####################################

@doc doc"""
Synthetise a pure tone with an interaural phase and interaural level difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `IPD`: Interaural phase difference in radians.
* `level`: Tone level in dB SPL.
* `ILD`: Interaural level difference in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the ILD will be applied.  ('right' or 'left')
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneIPDILD(frequency=440, phase=0, IPD=pi/2, level=65, ILD=10,
dur=1, rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```
""" ->

function pureToneIPDILD(;frequency::Real=1000, phase::Real=0, IPD::Real=0,
                     level::Real=65, ILD::Real=0, dur::Real=1, rampDur::Real=0.01,
                     channel::String="right", sf::Real=48000, maxLevel::Real=101)
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
@doc doc"""
Synthetise a pure tone with an interaural time and interaural level difference.

##### Parameters:

* `frequency`: Tone frequency in hertz.
* `phase`: Starting phase in radians.
* `ITD`: Interaural time difference in seconds.
* `level`: Tone level in dB SPL.
* `ILD`: Interaural level difference in dB SPL.
* `dur`: Tone duration in seconds.
* `ramp`: Duration of the onset and offset ramps in seconds.
* `channel`: Channel in which the ILD will be applied.  ('right' or 'left')
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : Array with dimensions (nSamples, 2).
       
##### Examples:

```julia
pt = pureToneITDILD(frequency=440, phase=0, ITD=0.004/1000, level=65,
ILD=10, dur=1, rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```
""" ->

function pureToneITDILD(;frequency::Real=1000, phase::Real=0, ITD::Real=0,
                     level::Real=65, ILD::Real=0, dur::Real=1, rampDur::Real=0.01,
                     channel::String="right", sf::Real=48000, maxLevel::Real=101)
    IPD = ITDToIPD(ITD, frequency)
    snd = pureToneIPDILD(frequency=frequency, phase=phase, IPD=IPD, level=level, ILD=ILD, dur=dur, rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
    return snd
end

##################################
## silence
##################################
@doc doc"""
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
""" ->
function silence(;dur=1000, channel="mono", sf=48000)
    nSamples = round(Int, dur * sf)
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
@doc doc"""
Synthetise band-limited noise from the addition of random-phase
sinusoids.

##### Parameters:

* `frequency1`: Start frequency of the noise.
* `frequency2`: End frequency of the noise.
* `level`: Noise spectrum level.
* `dur`: Tone duration in seconds.
* `rampDur`: Duration of the onset and offset ramps in seconds.
* `channel`: 'right', 'left' or 'diotic'. Channel in which the tone will be generated.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: 2-dimensional array of floats. The array has dimensions (nSamples, nChannels).
       
##### Examples:

```julia
nbNoise = steepNoise(frequency1=440, frequency2=660, level=65,
dur=1, rampDur=0.01, channel="right", sf=48000, maxLevel=100)
```

""" ->
function steepNoise(;frequency1=900, frequency2=1000, level=50, dur=1, rampDur=0.01, channel="diotic", sf=48000, maxLevel=101)
    
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end

    nSamples = round(Int, (dur-rampDur*2) * sf)
    nRamp = round(Int, rampDur * sf)
    nTot = nSamples + (nRamp * 2)

    spacing = 1 / dur
    components = 1 + floor((frequency2 - frequency1) / spacing)
    # SL = 10*log10(A^2/NHz) 
    # SL/10 = log10(A^2/NHz)
    # 10^(SL/10) = A^2/NHz
    # A^2 = 10^(SL/10) * NHz
    # RMS = 10^(SL/20) * sqrt(NHz) where NHz is the spacing between harmonics
    amp =  10^((level - maxLevel) / 20) * sqrt((frequency2 - frequency1) / components)
    
    timeAll = collect(0:nTot-1) / sf
    timeRamp = collect(0:nRamp-1)
    
    if channel == "mono"
        snd = zeros((nTot, 1))
    else
        snd = zeros((nTot, 2))
    end

    noise= zeros(nTot)
    for f=frequency1:spacing:frequency2
        radFreq = 2 * pi * f 
        phase = rand() * 2 * pi
        noise = noise + sin(phase + (radFreq * timeAll))
    end

    if channel == "mono"
        snd = zeros((nTot, 1))
    else
        snd = zeros((nTot, 2))
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
        for f=frequency1:spacing:frequency2
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
 
