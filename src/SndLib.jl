#   Copyright (C) 2013-2015 Samuele Carcagno <sam.carcagno@gmail.com>
#   This file is part of SndLib.jl

#    SndLib.jl is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    SndLib.jl is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with SndLib.jl.  If not, see <http://www.gnu.org/licenses/>.
#

#SndLib.jls is a module to generate sounds in julia
module SndLib

export addSounds, AMTone, AMToneIPD,
broadbandNoise,
complexTone,
fir2Filt!,
gate!,
ITDToIPD,
makePink!, makeSilence,
pureTone, pureToneILD, pureToneIPD, pureToneIPDILD, pureToneITD, pureToneITDILD,
scaleLevel, sound, steepNoise

VERSION < v"0.4-" && using Docile
using PyCall, WAV
#pyinitialize("python3")
@pyimport scipy.signal as scisig

@doc doc"""
Add or concatenate two sounds.

##### Parameters

* `snd1`: First sound.
* `snd2`: Second sound.
* `delay`: Delay in seconds between the onset of 'snd1' and the onset of 'snd2'
* `sf`: Sampling frequency in hertz of the two sounds.

##### Returns

* `snd` : array with dimensions (nSamples, nChannels) 
       
##### Examples

```julia
snd1 = pureTone(frequency=440, phase=0, level=65, duration=180,
ramp=10, channel="right", sf=48000, maxLevel=100)
snd2 = pureTone(frequency=880, phase=0, level=65, duration=180,
ramp=10, channel="right", sf=48000, maxLevel=100)
snd = addSounds(snd1=snd1, snd2=snd2, delay=100, sf=48000)
```
""" ->

function addSounds{T<:Real}(snd1::Array{T, 2}, snd2::Array{T, 2}; delay::Real=0, sf::Real=48000)
    nChans1 = size(snd1)[1]
    nChans2 = size(snd2)[1]
    if nChans1 != nChans2
        error("snd1 and snd2 must have the same number of channels")
    end
    nChans = nChans1
    nSnd1 = length(snd1[:,1])
    nSnd2 = length(snd2[:,1])
    snd1Duration = nSnd1/sf
    snd2Duration = nSnd2/sf

    #------------------------
    #            ...............................

    # Seg1           Seg2              Seg3
    
    nSampSeg1 = round(delay * sf)
    if nSampSeg1 < nSnd1
        nSampSeg2 = nSnd1 - nSampSeg1
        nSampSeg3 = nSnd2 - nSampSeg2
        seg1 = snd1[1:nSampSeg1,:]
        seg2a = snd1[nSampSeg1+1:nSnd1,:]
        if nSampSeg2 > nSnd2 # snd2 shorter than seg2, fill with zeros
            ldiff = nSampSeg2 - nSnd2
            diffSeg = zeros((ldiff, nChans))
            seg2b = vcat(snd2, diffSeg)
        else
            seg2b = snd2[1:nSampSeg2,:]
            seg3 = snd2[nSampSeg2+1:nSnd2,:]
        end

        seg2 = seg2a+seg2b
        snd = vcat(seg1, seg2)

        if nSampSeg2 < nSnd2
            snd = vcat(snd, seg3)
        end
            
    else
        seg1 = snd1
        seg2 = makeSilence(dur=delay-snd1Duration, channel=ifelse(nChans==1, "mono", "diotic"), sf=sf)
        seg3 = snd2
        snd = vcat(seg1, seg2)
        snd = vcat(snd, seg3)
    end
        
    return snd
end

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
snd = AMTone(carrierFreq=1000, AMFreq=20, AMDepth=1, carrierPhase=0, AMPhase,
level=65, dur=1, rampDur=0.01, channel="diotic", sf=48000, maxLevel=100)
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

    nSamples = int(round((dur-rampDur*2) * sf))
    nRamp = int(round(rampDur * sf))
    nTot = nSamples + (nRamp * 2)

    timeAll = [0:nTot-1] / sf
    timeRamp = [0:nRamp-1]

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
snd = AMToneIPD(carrierFreq=1000, AMFreq=20, AMDepth=1, carrierPhase=0, AMPhase,
carrierIPD=0, AMIPD=pi/2,
level=65, dur=1, rampDur=0.01, channel="diotic", sf=48000, maxLevel=100)
```    
""" ->

function AMToneIPD(;carrierFreq::Real=1000, AMFreq::Real=20, AMDepth::Real=1,
                carrierPhase::Real=0, AMPhase::Real=1.5*pi, carrierIPD::Real=0,
                AMIPD::Real=0, level::Real=60,
                dur::Real=1, rampDur::Real=0.01, channel::String="right",
                sf::Real=48000, maxLevel::Real=101)

    fixed = AMTone(carrierFreq=carrierFreq, AMFreq=AMfreq, AMDepth=AMDepth,
                carrierPhase=carrierPhase, AMPhase=AMPhase, level=level,
                dur=dur, rampDur=rampDur, channel="mono",
                sf=sf, maxLevel=maxLevel)

    shifted = AMTone(carrierFreq=carrierFreq, AMFreq=AMfreq, AMDepth=AMDepth,
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
    nSamples = int(round((dur-rampDur*2) * sf))
    nRamp = int(round(rampDur * sf))
    nTot = nSamples + (nRamp * 2)

    timeAll = [0:nTot-1] ./ sf
    timeRamp = [0:nRamp-1]

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
function complexTone(;F0::Real=220, harmPhase::String="sine", lowHarm::Int=1, highHarm::Int=10, stretch::Real=0,
                     level::Real=60, dur::Real=1, rampDur::Real=0.01, channel::String="diotic", sf::Real=48000,
                     maxLevel::Real=101)
    
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end

    amp = 10^((level - maxLevel) / 20)
    stretchHz = (F0*stretch)/100
    
    nSamples = int(round((duration-rampDur*2) * sf))
    nRamp = int(round(rampDur * sf))
    nTot = nSamples + (nRamp * 2)

    timeAll = [0:nTot-1] / sf
    timeRamp = [0:nRamp-1]

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
        snd_mono[1:nRamp, 2]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* tone[1:nRamp]
        snd_mono[nRamp+1:nRamp+nSamples, 2]        = amp * tone[nRamp+1:nRamp+nSamples]
        snd_mono[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* tone[nRamp+nSamples+1:nTot]
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

@doc doc"""
Impose onset and offset ramps to a sound.

##### Parameters:

* `rampDur`: The duration of the ramps.
* `sig`: The signal on which the ramps should be imposed.
* `sf`: The sampling frequency of 'sig'

##### Returns

* `sig`: The ramped signal.

##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, dur=2, rampDur=0,
channel="diotic", sf=48000, maxLevel=100)
gate!(sig=noise, rampDur=0.01, sf=48000)
```

""" ->
function gate!{T<:Real}(sig::Array{T, 2}; rampDur::Real=0.01, sf::Real=48000)

    nRamp = int(round(rampDur * sf))
    timeRamp = [0:nRamp-1] 
    nTot = length(sig[:,1])
    nStartSecondRamp = nTot - nRamp

    nChans = size(sig)[2]
    for i = 1:nChans
        sig[1:nRamp, i] = sig[1:nRamp, i] .*  ((1-cos(pi * timeRamp/nRamp))/2)
        sig[nStartSecondRamp+1:nTot, i] = sig[nStartSecondRamp+1:nTot, i] .* ((1+cos(pi * timeRamp/nRamp))/2)
    end

    return sig
end

@doc doc"""

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

""" ->


function ITDToIPD(ITD::Real, freq::Real)
    
    IPD = (ITD / (1/freq)) * 2 * pi
    
    return IPD
end

@doc doc"""
Convert a white noise into a pink noise.

The spectrum level of the pink noise at the frequency "ref"
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
noise = makePink(noise, sf=48000, ref=1000)
```
""" ->


function makePink!{T<:Real}(sig::Array{T, 2}; sf::Real=48000, ref::Real=1000)
    
    nSamples = size(sig)[1]
    nChans = size(sig)[2]
    ref = 1 + (ref * nSamples/sf)

    for i=1:nChans
        x = rfft(sig[:,i])
        n = length(x)
        idx = [1:n]#arange(1, len(x))
        mag = zeros(n)
        mag[1:n] = abs(x[1:n]) .* sqrt(ref./idx)
        mag[1] = abs(x[1])
        ph = angle(x)
        x = mag .* (cos(ph) + 1im * sin(ph))
        sig[:,i] = irfft(x, nSamples)
    end

    return sig
end


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
sil = makeSilence(dur=2, sf=48000)
```
""" ->
function makeSilence(;dur=1000, channel="mono", sf=48000)
    nSamples = int(round(dur * sf))
    if channel == "mono"
        snd = zeros(nSamples, 1)
    elseif channel == "diotic"
        snd = zeros(nSamples, 2)
    end
    return snd
end

@doc doc"""
Simple function to play sounds. Uses aplay on Linux and afplay on OSX.
Windows is not currently supported.

##### Arguments

* `snd`: The sound to be played.
* `sf`: The sampling frequency.

""" ->
function sound(snd, sf::Int=48000, nbits::Int=32)
    tmp = tempname()
    wavwrite(snd, tmp, Fs=sf, nbits=nbits)
    if OS_NAME == :Linux
        run(`aplay $tmp`)
    elseif OS_NAME == :OSX
        run(`afplay $tmp`)
    else
        println("Sorry sound function does not currently support your operating system")
    end
end


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
pt = pureTone(frequency=440, phase=0, level=65, duration=180,
ramp=10, channel="right", sf=48000, maxLevel=100)
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
    nSamples = int(round((dur-rampDur*2) * sf))
    nRamp = int(round(rampDur * sf))
    nTot = nSamples + (nRamp * 2)

    timeAll = [0:nTot-1] ./ sf
    timeRamp = [0:nRamp-1]

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
pt = pureToneILD(frequency=440, phase=0, level=65, ILD=10, duration=180,
ramp=10, channel="right", sf=48000, maxLevel=100)
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
pt = pureToneIPD(frequency=440, phase=0, IPD=pi/2, level=65, duration=180,
ramp=10, channel="right", sf=48000, maxLevel=100)
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
pt = pureToneITD(frequency=440, phase=0, ITD=0.004/1000, level=65, duration=180,
ramp=10, channel="right", sf=48000, maxLevel=100)
```
""" ->

function pureToneITD(;frequency::Real=1000, phase::Real=0, ITD::Real=0,
                     level::Real=65, dur::Real=1, rampDur::Real=0.01,
                     channel::String="right", sf::Real=48000, maxLevel::Real=101)
    IPD = ITDToIPD(ITD, frequency)
    snd = pureToneIPD(frequency=frequency, phase=phase, IPD=IPD, level=level, dur=dur, rampDur=rampDur, channel=channel, sf=sf, maxLevel=maxLevel)

    return snd
end


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
pt = pureToneIPDILD(frequency=440, phase=0, IPD=pi/2, level=65, ILD=10, duration=180,
ramp=10, channel="right", sf=48000, maxLevel=100)
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
pt = pureToneITDILD(frequency=440, phase=0, ITD=0.004/1000, level=65, ILD=10, duration=180,
ramp=10, channel="right", sf=48000, maxLevel=100)
```
""" ->

function pureToneITDILD(;frequency::Real=1000, phase::Real=0, ITD::Real=0,
                     level::Real=65, ILD::Real=0, dur::Real=1, rampDur::Real=0.01,
                     channel::String="right", sf::Real=48000, maxLevel::Real=101)
    IPD = ITDToIPD(ITD, frequency)
    snd = pureToneIPDILD(frequency=frequency, phase=phase, IPD=IPD, level=level, ILD=ILD, dur=dur, rampDur=rampDur, channel="mono", sf=sf, maxLevel=maxLevel)
    return snd
end


@doc doc"""
Increase or decrease the amplitude of a sound signal.

##### Parameters:

*`level`: Desired increment or decrement in dB SPL.
* `signal`: Signal to scale.

##### Returns:

* `sig`: 2-dimensional array of floats
       
##### Examples:

```julia
noise = broadbandNoise(spectrumLevel=40, duration=180, ramp=10,
channel="diotic", sf=48000, maxLevel=100)
noise = scale(sig=noise, level=-10) #reduce level by 10 dB
```
""" ->

function scaleLevel{T<:Real}(sig::Array{T, 2}; level::Real=10)
    
    #10**(level/20) is the amplitude corresponding to level
    #by multiplying the amplitudes we're adding the decibels
    # 20*log10(a1*a2) = 20*log10(a1) + 20*log10(a2)
    sig = sig * 10^(level/20)
    return sig
end

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
nbNoise = steepNoise(frequency=440, frequency2=660, level=65,
dur=180, rampDur=10, channel="right", sf=48000, maxLevel=100)
```

""" ->
function steepNoise(;frequency1=900, frequency2=1000, level=50, dur=1, rampDur=0.01, channel="diotic", sf=48000, maxLevel=101)
    
    if dur < rampDur*2
        error("Sound duration cannot be less than total duration of ramps")
    end

    nSamples = int(round((dur-rampDur*2) * sf))
    nRamp = int(round(rampDur * sf))
    nTot = nSamples + (nRamp * 2)

    spacing = 1 / dur
    components = 1 + floor((frequency2 - frequency1) / spacing)
    # SL = 10*log10(A^2/NHz) 
    # SL/10 = log10(A^2/NHz)
    # 10^(SL/10) = A^2/NHz
    # A^2 = 10^(SL/10) * NHz
    # RMS = 10^(SL/20) * sqrt(NHz) where NHz is the spacing between harmonics
    amp =  10^((level - maxLevel) / 20) * sqrt((frequency2 - frequency1) / components)
    
    timeAll = [0:nTot-1] / sf
    timeRamp = [0:nRamp-1]
    
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


@doc doc"""
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
* `snd`: The sound to be filtered.
* `nTaps`: Number of filter taps.
* `sf`: Sampling frequency of 'snd'.

##### Returns:

* `snd`: 2-dimensional array of floats

##### Notes:

If 'f1' and 'f2' are zero the filter will be low pass.
If 'f3' and 'f4' are equal to or greater than the nyquist
frequency (sf/2) the filter will be high pass.
In the other cases the filter will be band pass.

The order of the filter (number of taps) is fixed at 256.
This function uses internally 'scipy.signal.firwin2'.
       
##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, duration=180, ramp=10,
     channel="diotic", sf=48000, maxLevel=100)
lpNoise = fir2Filt(f1=0, f2=0, f3=1000, f4=1200, 
     snd=noise, sf=48000) #lowpass filter
hpNoise = fir2Filt(f1=0, f2=0, f3=24000, f4=26000, 
     snd=noise, sf=48000) #highpass filter
bpNoise = fir2Filt(f1=400, f2=600, f3=4000, f4=4400, 
.     snd=noise, sf=48000) #bandpass filter
```   

""" ->
function fir2Filt!{T<:Real}(f1::Real, f2::Real, f3::Real, f4::Real, snd::Array{T, 2}; nTaps::Int=256, sf::Real=48000)

    f1 = (f1 * 2) / sf
    f2 = (f2 * 2) / sf
    f3 = (f3 * 2) / sf
    f4 = (f4 * 2) / sf

    if f2 == 0 #low pass
        f = [0, f3, f4, 1]
        m = [1, 1, 0.00003, 0]
        
    elseif f3 < 1 #bandpass
        f = [0, f1, f2, ((f2+f3)/2), f3, f4, 1]
        m = [0, 0.00003, 1, 1, 1, 0.00003, 0]
        
    else #high pass
        f = [0, f1, f2, 0.999999, 1] #scipy wants that gain at the Nyquist is 0
        m = [0, 0.00003, 1, 1, 0]
    end
        
    b = convert(Array{eltype(snd),1}, scisig.firwin2(nTaps,f,m))  
    nChans = size(snd)[2]
    for i=1:nChans
        snd[:, i] = fftconvolve(snd[:,i], b, "same")
    end
    
    return snd
end

function _centered(arr, newsize)
    # Return the center newsize portion of the array.
    currsize = size(arr)[1]
    startind = div((currsize - newsize), 2) + 1
    endind = startind + newsize -1 #check indexing is the same pyth julia?
    return arr[startind:endind]
end

function fftconvolve(x, y, mode)
    s1 = size(x)[1]#check if array has two dim?
    s2 = size(y)[1]
    #println(typeof(x), typeof(y))
    convArray = conv(x,y)
    if mode == "full"
        return convArray
    elseif mode == "same"
        return _centered(convArray, s1)
    elseif mode == "valid"
        return _centered(convArray, abs(s1 - s2) + 1)
    end
end


end #module
