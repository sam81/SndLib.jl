################################
## addSounds
################################
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
snd = addSounds(snd1=snd1, snd2=snd2, delay=1, sf=48000)
```
""" ->

function addSounds{T<:Real, P<:Real}(snd1::Array{T, 2}, snd2::Array{P, 2}; delay::Real=0, sf::Real=48000)
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

#####################################
## fir2Filt!
#####################################
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

function ITDToIPD(ITD::Real, freq::Real)
    
    IPD = (ITD / (1/freq)) * 2 * pi
    
    return IPD
end



###############################
## makePink!
###############################
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
function sound{T<:Real}(snd::Array{T,2}, sf::Int=48000, nbits::Int=32)
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

