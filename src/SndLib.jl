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

export AMTone, broadbandNoise, complexTone, gate!, fir2Filt, makeSilence, pureTone, scaleLevel, sound, steepNoise

VERSION < v"0.4-" && using Docile
using PyCall, WAV
#using Devectorize
#pyinitialize("python3")
@pyimport scipy.signal as scisig

@doc doc"""
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
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: 2-dimensional array of floats
       
##### Examples:

```julia
snd = AMTone(frequency=1000, AMFreq=20, AMDepth=1, phase=0, AMPhase,
level=65, duration=980, ramp=10, channel="Diotic", sf=48000, maxLevel=100)
```    
""" ->

function AMTone(;frequency::Real=1000, AMFreq::Real=20, AMDepth::Real=1,
                phase::Real=0, AMPhase::Real=1.5*pi, level::Real=60,
                duration::Real=980, ramp::Real=10, channel::String="Diotic",
                sf::Real=48000, maxLevel::Real=101)

    amp = 10^((level - maxLevel) / 20)
    duration = duration / 1000 #convert from ms to sec
    ramp = ramp / 1000

    nSamples = int(round(duration * sf))
    nRamp = int(round(ramp * sf))
    nTot = nSamples + (nRamp * 2)

    timeAll = [0:nTot-1] / sf
    timeRamp = [0:nRamp-1]

    snd_mono = zeros(nTot, 1)
    snd_mono[1:nRamp, 1] = amp * (1+AMDepth*sin(2*pi*AMFreq*timeAll[1:nRamp]+AMPhase)) .* ((1-cos(pi* timeRamp/nRamp))/2) .* sin(2*pi*frequency * timeAll[1:nRamp] + phase)
    snd_mono[nRamp+1:nRamp+nSamples, 1] = amp * (1 + AMDepth*sin(2*pi*AMFreq*timeAll[nRamp+1:nRamp+nSamples]+AMPhase)) .* sin(2*pi*frequency * timeAll[nRamp+1:nRamp+nSamples] + phase)
    snd_mono[nRamp+nSamples+1:nTot, 1] = amp * (1 + AMDepth*sin(2*pi*AMFreq*timeAll[nRamp+nSamples+1:nTot]+AMPhase)) .* ((1+cos(pi * timeRamp/nRamp))/2) .* sin(2*pi*frequency * timeAll[nRamp+nSamples+1:nTot] + phase)
    
    if channel == "Mono"
        snd = snd_mono
    else
        snd = zeros(nTot, 2)
    end

    if channel == "Right"
        snd[:,2] = snd_mono
    elseif channel == "Left"
        snd[:,1] = snd_mono
    elseif channel == "Diotic"
        snd[:,1] = snd_mono
        snd[:,2] = snd_mono
    end
       
    return snd
end

@doc doc"""
Synthetise a broadband noise.

##### Parameters:

* `spectrumLevel`: Intensity spectrum level of the noise in dB SPL. 
* `duration`: Noise duration (excluding ramps) in milliseconds.
* `ramp`: Duration of the onset and offset ramps in milliseconds.
        The total duration of the sound will be duration+ramp*2.
* `channel`: Channel in which the noise will be generated.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd` : 2-dimensional array of floats
        The array has dimensions (nSamples, 2).
       
##### Examples:

```julia
noise = broadbandNoise(spectrumLevel=40, duration=180, ramp=10,
         channel="Both", sf=48000, maxLevel=100)
```

""" ->
function broadbandNoise(;spectrumLevel::Real=40, duration::Real=980, ramp::Real=10, channel::String="Both", sf::Real=48000, maxLevel::Real=101)

    ## Comments:
    ## The intensity spectrum level in dB is SL
    ## The peak amplitude A to achieve a desired SL is
    ## SL = 10*log10(RMS^2/NHz) that is the total RMS^2 divided by the freq band
    ## SL/10 = log10(RMS^2/NHz)
    ## 10^(SL/10) = RMS^2/NHz
    ## RMS^2 = 10^(SL/10) * NHz
    ## RMS = 10^(SL/20) * sqrt(NHz)
    ## NHz = sampRate / 2 (Nyquist)

    amp = sqrt(sf/2)*(10^((spectrumLevel - maxLevel) / 20))
    duration = duration / 1000 #convert from ms to sec
    ramp = ramp / 1000

    nSamples = int(round(duration * sf))
    nRamp = int(round(ramp * sf))
    nTot = nSamples + (nRamp * 2)

    timeAll = [0:nTot-1] ./ sf
    timeRamp = [0:nRamp-1]

    snd = zeros(nTot, 2)
    #random is a numpy module
    noise = (rand(nTot) + rand(nTot)) - (rand(nTot) + rand(nTot))
    RMS = sqrt(mean(noise.*noise))
    #scale the noise so that the maxAmplitude goes from -1 to 1
    #since A = RMS*sqrt(2)
    scaled_noise = noise / (RMS * sqrt(2))


    if channel == "Right"
        snd[1:nRamp, 2] = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* scaled_noise[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 2] = amp * scaled_noise[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* scaled_noise[nRamp+nSamples+1:nTot]
    elseif channel == "Left"
        snd[1:nRamp, 1] = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* scaled_noise[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1] = amp * scaled_noise[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* scaled_noise[nRamp+nSamples+1:nTot]
    elseif channel == "Both"
        snd[1:nRamp, 2] = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* scaled_noise[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 2] = amp * scaled_noise[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* scaled_noise[nRamp+nSamples+1:nTot]
        snd[1:nRamp, 1] = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* scaled_noise[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1] = amp * scaled_noise[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* scaled_noise[nRamp+nSamples+1:nTot]
    end

    return snd
end

@doc doc"""
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
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: 2-dimensional array of floats
        The array has dimensions (nSamples, 2).

##### Examples

```julia
 ct = complexTone(F0=440, harmPhase="Sine", lowHarm=3, highHarm=10,
         stretch=0, level=55, duration=180, ramp=10, channel="Both",
         sf=48000, maxLevel=100)
```

""" ->
function complexTone(;F0=220, harmPhase="Sine", lowHarm=1, highHarm=10, stretch=0, level=60, duration=980, ramp=10, channel="Both", sf=48000, maxLevel=101)
    
    amp = 10^((level - maxLevel) / 20)
    duration = duration / 1000 #convert from ms to sec
    ramp = ramp / 1000
    stretchHz = (F0*stretch)/100
    
    nSamples = int(round(duration * sf))
    nRamp = int(round(ramp * sf))
    nTot = nSamples + (nRamp * 2)

    timeAll = [0:nTot-1] / sf
    timeRamp = [0:nRamp-1]

    snd = zeros(nTot, 2)
    if channel == "Right" || channel == "Left" || channel == "Both"
        tone = zeros(nTot)
    elseif channel == "Odd Left" || channel == "Odd Right"
        toneOdd = zeros(nTot)
        toneEven = zeros(nTot)
    end

    if harmPhase == "Sine"
        for i=lowHarm:highHarm
            if channel == "Right" || channel == "Left" || channel == "Both"
                tone =  tone + sin(2 * pi * ((F0 * i) + stretchHz) * timeAll)
            elseif channel == "Odd Left" || channel == "Odd Right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                else
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                end
            end
        end
    elseif harmPhase == "Cosine"
        for i=lowHarm:highHarm
            if channel == "Right" || channel == "Left" || channel == "Both"
                tone = tone + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
            elseif channel == "Odd Left" || channel == "Odd Right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                else
                    toneEven = toneEven + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                end
            end
        end
    elseif harmPhase == "Alternating"
        for i=lowHarm:highHarm
            if i%2 > 0 #odd harmonic
                if channel == "Right" || channel == "Left" || channel == "Both"
                    tone = tone + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                elseif channel == "Odd Left" || channel == "Odd Right"
                    toneOdd = toneOdd + cos(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                end
            else #even harmonic
                if channel == "Right" || channel == "Left" || channel == "Both"
                    tone = tone + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                elseif channel == "Odd Left" || channel == "Odd Right"
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll)
                end
            end
        end
    elseif harmPhase == "Schroeder"
        for i=lowHarm:highHarm
            phase = -pi * i * (i - 1) / highHarm
            if channel == "Right" || channel == "Left" || channel == "Both"
                tone = tone + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
            elseif channel == "Odd Left" || channel == "Odd Right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                else
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                end
            end
        end
    elseif harmPhase == "Random"
        for i=lowHarm:highHarm
            phase = rand() * 2 * pi
            if channel == "Right" || channel == "Left" || channel == "Both"
                tone = tone + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
            elseif channel == "Odd Left" || channel == "Odd Right"
                if i%2 > 0 #odd harmonic
                    toneOdd = toneOdd + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                else
                    toneEven = toneEven + sin(2 * pi * ((F0 * i)+stretchHz) * timeAll + phase)
                end
            end
        end
    end
  

    if channel == "Right"
        snd[1:nRamp, 2]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* tone[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 2]        = amp * tone[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* tone[nRamp+nSamples+1:nTot]
    elseif channel == "Left"
        snd[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* tone[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1]        = amp * tone[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* tone[nRamp+nSamples+1:nTot]
    elseif channel == "Both"
        snd[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* tone[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1]        = amp * tone[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* tone[nRamp+nSamples+1:nTot]
        snd[:, 2] = snd[:, 1]
    elseif channel == "Odd Left"
        snd[1:nRamp, 1]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* toneOdd[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1]        = amp * toneOdd[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* toneOdd[nRamp+nSamples+1:nTot]
        snd[1:nRamp, 2]                     = amp * ((1-cos(pi * timeRamp/nRamp))/2) .* toneEven[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 2]        = amp * toneEven[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* toneEven[nRamp+nSamples+1:nTot]
    elseif channel == "Odd Right"
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

* `ramps`: The duration of the ramps.
* `sig`: The signal on which the ramps should be imposed.
* `sf`: The sampling frequency os 'sig'

##### Returns

* `sig`: The ramped signal.

##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, duration=200, ramp=0,
channel="Both", sf=48000, maxLevel=100)
gate!(sig=noise, ramps=10, sf=48000)
```

""" ->
function gate!(sig, ramps, sf)
    ramps = ramps / 1000
    nRamp = int(round(ramps * sf))
    timeRamp = [0:nRamp-1] 
    nTot = length(sig[:,1])
    nStartSecondRamp = nTot - nRamp
    
    sig[1:nRamp, 1] = sig[1:nRamp, 1] .*  ((1-cos(pi * timeRamp/nRamp))/2)
    sig[1:nRamp, 2] = sig[1:nRamp, 2] .*  ((1-cos(pi * timeRamp/nRamp))/2)
    sig[nStartSecondRamp+1:nTot, 1] = sig[nStartSecondRamp+1:nTot, 1] .* ((1+cos(pi * timeRamp/nRamp))/2)
    sig[nStartSecondRamp+1:nTot, 2] = sig[nStartSecondRamp+1:nTot, 2] .* ((1+cos(pi * timeRamp/nRamp))/2)

    return sig
end

@doc doc"""
Generate a silence.

This function just fills an array with zeros for the
desired duration.
    
##### Parameters:

* `duration`: Duration of the silence in milliseconds.
* `sf`: Samplig frequency in Hz.

##### Returns:

* `snd`: 2-dimensional array of floats
The array has dimensions (nSamples, 2).
       

##### Examples:

```julia
sil = makeSilence(duration=200, sf=48000)
```
""" ->
function makeSilence(;duration=1000, sf=48000)
    duration = duration / 1000 #convert from ms to sec
    nSamples = int(round(duration * sf))
    snd = zeros(nSamples, 2)
    
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
""" ->
function pureTone(;frequency=1000, phase=0, level=65, duration=980, ramp=10, channel="Both", sf=48000, maxLevel=101) 

    amp = 10^((level - maxLevel) / 20)
    duration = duration / 1000 #convert from ms to sec
    ramp = ramp / 1000
    
    nSamples = int(round(duration * sf))
    nRamp = int(round(ramp * sf))
    nTot = nSamples + (nRamp * 2)

    timeAll = [0:nTot-1] ./ sf
    timeRamp = [0:nRamp-1]
    
    snd = zeros((nTot, 2))
    if channel == "Right"
        snd[1:nRamp, 2] = amp .* ((1.-cos(pi .* timeRamp/nRamp))./2) .* sin(2*pi*frequency .* timeAll[1:nRamp] .+ phase)
        snd[nRamp+1:nRamp+nSamples, 2] = amp* sin(2*pi*frequency .* timeAll[nRamp+1:nRamp+nSamples] .+ phase)
        snd[nRamp+nSamples+1:nTot, 2] = amp .* ((1.+cos(pi .* timeRamp/nRamp))./2) .* sin(2*pi*frequency * timeAll[nRamp+nSamples+1:nTot] .+ phase)
    elseif channel == "Left"
        snd[1:nRamp, 1] = amp .* ((1.-cos(pi .* timeRamp/nRamp))./2) .* sin(2*pi*frequency .* timeAll[1:nRamp] .+ phase)
        snd[nRamp+1:nRamp+nSamples, 1] = amp* sin(2*pi*frequency .* timeAll[nRamp+1:nRamp+nSamples] .+ phase)
        snd[nRamp+nSamples+1:nTot, 1] = amp .* ((1.+cos(pi .* timeRamp/nRamp))./2) .* sin(2*pi*frequency * timeAll[nRamp+nSamples+1:nTot] .+ phase)
    elseif channel == "Both"
        snd[1:nRamp, 1] = amp .* ((1 .-cos(pi .* timeRamp/nRamp))./2) .* sin(2*pi*frequency .* timeAll[1:nRamp] .+ phase)
        snd[nRamp+1:nRamp+nSamples, 1] = amp* sin(2*pi*frequency .* timeAll[nRamp+1:nRamp+nSamples] .+ phase)
        snd[nRamp+nSamples+1:nTot, 1] = amp .* ((1.+cos(pi .* timeRamp/nRamp))./2) .* sin(2*pi*frequency * timeAll[nRamp+nSamples+1:nTot] .+ phase)
        snd[:, 2] = snd[:, 1]
    end
    
        
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
channel="Both", sf=48000, maxLevel=100)
noise = scale(sig=noise, level=-10) #reduce level by 10 dB
```
""" ->
function scaleLevel(sig, level)
    
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
* `duration`: Tone duration (excluding ramps) in milliseconds.
* `ramp`: Duration of the onset and offset ramps in milliseconds.
        The total duration of the sound will be duration+ramp*2.
* `channel`: 'Right', 'Left' or 'Both'. Channel in which the tone will be generated.
* `sf`: Samplig frequency in Hz.
* `maxLevel`: Level in dB SPL output by the soundcard for a sinusoid of amplitude 1.

##### Returns:

* `snd`: 2-dimensional array of floats
        The array has dimensions (nSamples, 2).
       
##### Examples:

```julia
nbNoise = steepNoise(frequency=440, frequency2=660, level=65,
duration=180, ramp=10, channel="Right", sf=48000, maxLevel=100)
```

""" ->
function steepNoise(;frequency1=900, frequency2=1000, level=50, duration=980, ramp=10, channel="Both", sf=48000, maxLevel=101)

    duration = duration/1000 #convert from ms to sec
    ramp = ramp/1000

    totDur = duration + (2 * ramp)
    nSamples = int(round(duration * sf))
    nRamp = int(round(ramp * sf))
    nTot = nSamples + (nRamp * 2)

    spacing = 1 / totDur
    components = 1 + floor((frequency2 - frequency1) / spacing)
    # SL = 10*log10(A^2/NHz) 
    # SL/10 = log10(A^2/NHz)
    # 10^(SL/10) = A^2/NHz
    # A^2 = 10^(SL/10) * NHz
    # RMS = 10^(SL/20) * sqrt(NHz) where NHz is the spacing between harmonics
    amp =  10^((level - maxLevel) / 20) * sqrt((frequency2 - frequency1) / components)
    
    timeAll = [0:nTot-1] / sf
    timeRamp = [0:nRamp-1]
    snd = zeros(nTot, 2)

    noise= zeros(nTot)
    for f=frequency1:spacing:frequency2
        radFreq = 2 * pi * f 
        phase = rand() * 2 * pi
        noise = noise + sin(phase + (radFreq * timeAll))
    end
    
    if channel == "Right"
        snd[1:nRamp, 2] = amp .* ((1-cos(pi * timeRamp/nRamp))/2) .* noise[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 2] = amp * noise[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* noise[nRamp+nSamples+1:nTot]
    elseif channel == "Left"
        snd[1:nRamp, 1] = amp .* ((1-cos(pi * timeRamp/nRamp))/2) .* noise[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1] = amp * noise[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* noise[nRamp+nSamples+1:nTot]
    elseif channel == "Both"
        snd[1:nRamp, 2] = amp .* ((1-cos(pi * timeRamp/nRamp))/2) .* noise[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 2] = amp * noise[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 2] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* noise[nRamp+nSamples+1:nTot]
        snd[1:nRamp, 1] = amp .* ((1-cos(pi * timeRamp/nRamp))/2) .* noise[1:nRamp]
        snd[nRamp+1:nRamp+nSamples, 1] = amp * noise[nRamp+1:nRamp+nSamples]
        snd[nRamp+nSamples+1:nTot, 1] = amp * ((1+cos(pi * timeRamp/nRamp))/2) .* noise[nRamp+nSamples+1:nTot]
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
* `snd`:The sound to be filtered.
* `sf`: Sampling frequency of 'snd'.

##### Returns:

* `snd`: 2-dimensional array of floats

##### Notes:

If 'f1' and 'f2' are zero the filter will be lowpass.
If 'f3' and 'f4' are equal to or greater than the nyquist
frequency (sf/2) the filter will be highpass.
In the other cases the filter will be bandpass.

The order of the filter (number of taps) is fixed at 256.
This function uses internally 'scipy.signal.firwin2'.
       
##### Examples

```julia
noise = broadbandNoise(spectrumLevel=40, duration=180, ramp=10,
     channel="Both", sf=48000, maxLevel=100)
lpNoise = fir2Filt(f1=0, f2=0, f3=1000, f4=1200, 
     snd=noise, sf=48000) #lowpass filter
hpNoise = fir2Filt(f1=0, f2=0, f3=24000, f4=26000, 
     snd=noise, sf=48000) #highpass filter
bpNoise = fir2Filt(f1=400, f2=600, f3=4000, f4=4400, 
.     snd=noise, sf=48000) #bandpass filter
```   

""" ->
function fir2Filt(f1, f2, f3, f4, snd, sf)

 
       
    ## if filterType == "lowpass"
    ##     f3 = cutoffs[1]
    ##     f4 = cutoffs[1] * (1+transitionWidth)
    ##     f3 = (f3*2) / sampRate
    ##     f4 = (f4*2) / sampRate
    ##     f = [0, f3, f4, 1]
    ##     m = [1, 1, 0.00003, 0]
    ## elseif filterType == "highpass"
    ##     f1 = cutoffs[1] * (1-transitionWidth)
    ##     f2 = cutoffs[1]
    ##     f1 = (f1*2) / sampRate
    ##     f2 = (f2*2) / sampRate
    ##     f = [0, f1, f2, 0.999999, 1] #high pass
    ##     m = [0, 0.00003, 1, 1, 0]
    ## elseif filterType == "bandpass"
    ##     f1 = cutoffs[1] * (1-transitionWidth)
    ##     f2 = cutoffs[1]
    ##     f3 = cutoffs[2]
    ##     f4 = cutoffs[2] * (1+transitionWidth)
    ##     f1 = (f1*2) / sampRate
    ##     f2 = (f2*2) / sampRate
    ##     f3 = (f3*2) / sampRate
    ##     f4 = (f4*2) / sampRate
    ##     f = [0, f1, f2, ((f2+f3)/2), f3, f4, 1]
    ##     m = [0, 0.00003, 1, 1, 1, 0.00003, 0]
    ## end

    ## b = convert(Array{eltype(rec),1}, scisig.firwin2(nTaps,f,m))
    ## #b = ones(Float32,nTaps)
    ## nChannels = size(rec)[1]
    ## if channels == nothing
    ##     channels = [1:nChannels]
    ## end
   
    ## for i=1:nChannels
    ##     if in(i, channels) == true
    ##         rec[i,:] = fftconvolve(reshape(rec[i,:], size(rec[i,:])[2]), b, "same")
    ##     end
    ## end
    ## return rec




    f1 = (f1 * 2) / sf
    f2 = (f2 * 2) / sf
    f3 = (f3 * 2) / sf
    f4 = (f4 * 2) / sf

    n = 256

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
        
    b = convert(Array{eltype(snd),1}, scisig.firwin2(n,f,m))  
    #b = firwin2 (n,f,m);
    x = copy(snd)
    x[:, 2] = fftconvolve(x[:,2], b, "same")
    x[:, 1] = fftconvolve(x[:,1], b, "same")
    #x[:, 0] = convolve(snd[:,0], b, 1)
    #x[:, 1] = convolve(snd[:,1], b, 1)
    
    return x





    
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
