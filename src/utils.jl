
################################
## freqFromERBInterval
################################
@doc doc"""
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
""" ->
function freqFromERBInterval{T<:Real, P<:Real}(f1::Union(T, AbstractVector{T}), deltaERB::Union(P, AbstractVector{P}))
    f2 = (10.^((deltaERB + 21.4*log10(0.00437*f1 +1))/21.4)-1) / 0.00437 
    return f2
end

#######################################
## ERBDistance
#######################################
@doc doc"""
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
""" ->

function ERBDistance(f1::Real, f2::Real)

    deltaERB = 21.4*log10(0.00437*f2+1) - 21.4*log10(0.00437*f1+1)

    return deltaERB
end
 
