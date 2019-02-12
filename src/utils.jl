## The MIT License (MIT)

## Copyright (c) 2013-2019 Samuele Carcagno <sam.carcagno@gmail.com>

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
## freqFromCentInterval
################################
"""
Compute the frequency, in Hz, corresponding to a distance,
in equivalent cents of `deltaCents` from `f1`.

$(SIGNATURES)

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
"""
function freqFromCentInterval(f1::Union{T, AbstractVector{T}}, deltaCent::Union{P, AbstractVector{P}}) where {T<:Real, P<:Real}
    f2 = f1*2 .^ (deltaCent/1200)
    return f2
end

#######################################
## centDistance
#######################################
"""
Compute the distance in cents
between the frequencies `f1` and `f2`.

$(SIGNATURES)

##### Parameters

* `f1`: frequency 1 in Hz
* `f2`: frequency 2 in Hz

##### Returns

* `deltaCents`: distance between f1 and f2 in cents.

##### Examples

```julia
centDistance(1000, 1200)
```
"""
function centDistance(f1::Real, f2::Real)

    deltaCents = 1200*log2(f2/f1)

    return deltaCents
end


################################
## freqFromERBInterval
################################
"""
Compute the frequency, in Hz, corresponding to a distance,
in equivalent rectangular bandwidths (ERBs), of `deltaERB` from `f1`.

$(SIGNATURES)

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
"""
function freqFromERBInterval(f1::Union{T, AbstractVector{T}}, deltaERB::Union{P, AbstractVector{P}}) where {T<:Real, P<:Real}
    f2 = (10 .^ ((deltaERB .+ 21.4*log10.(0.00437*f1 .+ 1))/21.4) .-1) / 0.00437
    return f2
end

#######################################
## ERBDistance
#######################################
"""
Compute the distance in equivalent rectangular bandwiths (ERBs)
between the frequencies `f1` and `f2`.

$(SIGNATURES)

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
"""
function ERBDistance(f1::Real, f2::Real)

    deltaERB = 21.4*log10(0.00437*f2+1) - 21.4*log10(0.00437*f1+1)

    return deltaERB
end

###################
## intNCyclesFreq
###################
"""
Compute the frequency closest to 'freq' that has an integer number
of cycles for the given sound duration.

$(SIGNATURES)

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

"""
function intNCyclesFreq(freq::Real, dur::Real)
    adjFreq = round(freq*dur)/dur

    return adjFreq
end
