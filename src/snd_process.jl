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

################################
## addSounds
################################
@doc doc"""
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
""" ->

function addSounds{T<:Real, P<:Real}(snd1::Array{T, 2}, snd2::Array{P, 2}; delay::Real=0, sf::Real=48000)
    nChans1 = size(snd1)[2]
    nChans2 = size(snd2)[2]
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
    
    nSampSeg1 = round(Int, delay * sf)
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
        seg2 = silence(dur=delay-snd1Duration, channel=ifelse(nChans==1, "mono", "diotic"), sf=sf)
        seg3 = snd2
        snd = vcat(seg1, seg2)
        snd = vcat(snd, seg3)
    end
        
    return snd
end

############################
## delayAdd
############################
@doc doc"""
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
""" ->

function delayAdd!{T<:Real, P<:Integer}(sig::Array{T,2}; delay::Real=0.01,
                                        gain::Real=1, iterations::Integer=1,
                                        configuration::String="add same",
                                        channel::Union(P, AbstractVector{P})=[1:size(sig)[2]],
                                        sf::Real=48000)

    if in(configuration, ["add same", "add original"]) == false
        error("`configuration` must be either 'add same', or 'add original'")
    end

    delayPnt = round(Int, delay * sf)
    nChans = size(sig)[2]
    nSamples = length(sig[:,1])
    if channel == "all"
        chans = [1:nChans]
    else
        chans = channel
    end
    #snd = zeros(nSamples, nChans)

    if configuration == "add original"
        original_sig = copy(sig)
    end

    for ch in chans
        delayed_sig = zeros(nSamples, 1)
        en_input = sqrt(sum(sig[:,ch].^2))
        if configuration == "add same"
            for i=1:iterations
                delayed_sig = vcat(sig[delayPnt+1:nSamples,ch], sig[1:delayPnt,ch])
                delayed_sig = delayed_sig * gain
                sig[:,ch] = sig[:,ch] + delayed_sig
            end
        elseif configuration == "add original"
            for i=1:iterations
                delayed_sig = vcat(sig[delayPnt+1:nSamples,ch], sig[1:delayPnt,ch])
                delayed_sig = delayed_sig * gain
                sig[:,ch] = original_sig[:,ch] + delayed_sig
            end
        end
        en_output = sqrt(sum(sig[:,ch].^2))
        scale_factor = en_input / en_output
        sig[:,ch] = sig[:,ch] * scale_factor
    end

    return sig
end


#####################################
## fir2Filt!
#####################################
@doc doc"""
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
""" ->
function fir2Filt!{T<:Real}(f1::Real, f2::Real, f3::Real, f4::Real, snd::Array{T, 2}; nTaps::Integer=256, sf::Real=48000)

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

###################
## _centered
###################
function _centered(arr, newsize)
    # Return the center newsize portion of the array.
    currsize = size(arr)[1]
    startind = div((currsize - newsize), 2) + 1
    endind = startind + newsize -1 #check indexing is the same pyth julia?
    return arr[startind:endind]
end

#######################
## fftconvolve
#######################
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


####################################
## gate!
####################################
@doc doc"""
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
""" ->
function gate!{T<:Real}(sig::Array{T, 2}; rampDur::Real=0.01, sf::Real=48000)

    nRamp = round(Int, rampDur * sf)
    timeRamp = collect(0:nRamp-1)
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
""" ->
function getRMS{T<:Real}(sig::Array{T,2}, channel::Union(String, Integer))
    RMS = (FloatingPoint)[]
    if channel == "all"
        push!(RMS, sqrt(mean(sig.*sig)))
        #RMS = sqrt(mean(sig.*sig))
    elseif channel == "each"
        nChans = size(sig)[2]
        for i=1:nChans
            push!(RMS, sqrt(mean(sig[:,i].*sig[:,i])))
        end
    else
        push!(RMS, sqrt(mean(sig[:,channel].*sig[:,channel])))
        #RMS = sqrt(mean(sig[:,channel].*sig[:,channel]))
    end
    return RMS
    end

########################
## ITDShift!
########################

@doc doc"""
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
"""->
function ITDShift!{T<:Real}(sig::Array{T,2}, f1::Real, f2::Real; ITD::Real=300/1000000, channel::String="left", sf::Real=48000)

    if in(channel, ["right", "left"]) == false
        error("Channel must be one of 'right', or 'left'")
    end

    nSamples = size(sig)[1]
    fftPoints = nSamples#nextpow2(nSamples)
    nUniquePnts = ceil((fftPoints+1)/2)
    #compute the frequencies of the first half of the FFT
    freqArray1 = collect(0:nUniquePnts) * (sf / fftPoints)
    #remove DC offset and nyquist for the second half of the FFT
    freqArray2 = -flipdim(collect(1:(nUniquePnts-1)),1) * (sf / fftPoints) 
    #find the indexes of the frequencies for which to set the ITD for the first half of the FFT
    sh1 = find((freqArray1 .>= f1) & (freqArray1 .<= f2))
    #same as above for the second half of the FFT
    sh2 = find((freqArray2 .<= -f1) & (freqArray2 .>= -f2))
    #compute IPSs for the first half of the FFT
    phaseShiftArray1 = ITDToIPD(ITD/1000000, freqArray1[sh1])
    #same as above for the second half of the FFT
    phaseShiftArray2 = ITDToIPD(ITD/1000000, freqArray2[sh2])
    #get the indexes of the first half of the FFT
    p1Start = 1; p1End = length(freqArray1)
    #get the indexes of the second half of the FFT
    p2Start = length(freqArray1)+1; p2End = fftPoints 
        
    if channel == "left"
        x = fft(sig[:,1])#, fftPoints)
    elseif channel == "right"
        x = fft(sig[:,2])#, fftPoints)
    end
    
    x1 = x[p1Start:p1End] #first half of the FFT
    x2 = x[p2Start:p2End] #second half of the FFT
    x1mag = abs(x1); x2mag = abs(x2) 
    x1Phase =  angle(x1); x2Phase =  angle(x2);
    x1Phase[sh1] = x1Phase[sh1] + phaseShiftArray1 #change phases
    x2Phase[sh2] = x2Phase[sh2] + phaseShiftArray2
    x1 = x1mag .* (cos(x1Phase) + (1im * sin(x1Phase))) #rebuild FFTs
    x2 = x2mag .* (cos(x2Phase) + (1im * sin(x2Phase)))
    x = vcat(x1, x2)
    x = real(ifft(x)) #inverse transform to get the sound back
    
    if channel == "left"
        sig[:,1] = x[1:nSamples]
    elseif channel == "right"
        sig[:,2] = x[1:nSamples]
    end

    return sig
end


#################################
## ITDToIPD
#################################
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

function ITDToIPD{T<:Real}(ITD::Real, freq::Union(T, AbstractVector{T}))
    
    IPD = (ITD ./ (1./freq)) * 2 * pi
    
    return IPD
end

#########################
## phaseShift!
#########################
@doc doc"""
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
"""->

function phaseShift!{T<:Real, P<:Integer}(sig::Array{T, 2}, f1::Real, f2::Real; phaseShift::Real=pi, shiftType::String="step", channel::Union(P, AbstractVector{P})=1, sf::Real=48000)

    if in(shiftType, ["linear", "random", "step"]) == false
        error("`shiftType`must be one of 'linear', 'random', 'step'")
    end
    
    nSamples = size(sig)[1]
    nChans = size(sig)[2]
    fftPoints = nSamples#nextpow2(nSamples)
    nUniquePnts = ceil((fftPoints+1)/2)
    freqArray1 = collect(0:nUniquePnts) * (sf / fftPoints)
    freqArray2 = -flipdim(collect(1:(nUniquePnts-1)), 1) * (sf / fftPoints) #remove DC offset and nyquist
    ## sh1 = where((freqArray1>f1) & (freqArray1<f2))
    ## sh2 = where((freqArray2<-f1) & (freqArray2>-f2))
    sh1 = find((freqArray1 .>= f1) & (freqArray1 .<= f2))
    sh2 = find((freqArray2 .<= -f1) & (freqArray2 .>= -f2))

    p1Start = 1; p1End = length(freqArray1)
    p2Start = length(freqArray1)+1; p2End = fftPoints
    #println(string(p1Start, " ", p1End, " ", p2Start, " ", p2End))

    if shiftType == "linear"
        phaseShiftArray1 = linspace(0, phaseShift, length(sh1))
        phaseShiftArray2 = - linspace(phaseShift, 0, length(sh2))
    elseif shiftType == "step"
        phaseShiftArray1 = [float(phaseShift) for ll=1:length(sh1)]
        phaseShiftArray2 = [float(-phaseShift) for ll=1:length(sh1)]
    elseif shiftType == "random"
        phaseShiftArray1 = rand(length(sh1))*phaseShift
        phaseShiftArray2 = -flipdim(phaseShiftArray1, 1)
    end

    for c=1:length(channel)
        ch = channel[c]
        x = fft(sig[:, ch])#, fftPoints)
        x1 = x[p1Start:p1End]
        x2 = x[p2Start:p2End]
        x1mag = abs(x1); x2mag = abs(x2)
        x1Phase =  angle(x1); x2Phase =  angle(x2);
        x1Phase[sh1] = x1Phase[sh1] + phaseShiftArray1
        x2Phase[sh2] =  x2Phase[sh2] + phaseShiftArray2
        x1 = x1mag .* (cos(x1Phase) + (1im * sin(x1Phase)))
        x2 = x2mag .* (cos(x2Phase) + (1im * sin(x2Phase)))
        x = vcat(x1, x2)
        x = real(ifft(x))
        sig[:, ch] = x[1:nSamples]
    end

    return sig
end

###############################
## makePink!
###############################
@doc doc"""
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
""" ->

function makePink!{T<:Real}(sig::Array{T, 2}; sf::Real=48000, ref::Real=1000)
    
    nSamples = size(sig)[1]
    nChans = size(sig)[2]
    ref = 1 + (ref * nSamples/sf)

    for i=1:nChans
        x = rfft(sig[:,i])
        n = length(x)
        idx = collect(1:n)#arange(1, len(x))
        mag = zeros(n)
        mag[1:n] = abs(x[1:n]) .* sqrt(ref./idx)
        mag[1] = abs(x[1])
        ph = angle(x)
        x = mag .* (cos(ph) + 1im * sin(ph))
        sig[:,i] = irfft(x, nSamples)
    end

    return sig
end



###################################
## scaleLevel
###################################
@doc doc"""
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
""" ->

function scaleLevel{T<:Real}(sig::Array{T, 2}; level::Real=10)
    
    #10**(level/20) is the amplitude corresponding to level
    #by multiplying the amplitudes we're adding the decibels
    # 20*log10(a1*a2) = 20*log10(a1) + 20*log10(a2)
    sig = sig * 10^(level/20)
    return sig
end


#############################
## sound
#############################
@doc doc"""
Simple function to play sounds. Uses aplay on Linux and afplay on OSX.
Windows is not currently supported.

##### Arguments

* `snd`: The sound to be played.
* `sf`: The sampling frequency.

""" ->
function sound{T<:Real}(snd::Array{T,2}, sf::Integer=48000, nbits::Integer=32)
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

