using Test, SndLib, WAV


wavDir = "wavDir/"
sf = 48000
nbits = 32
maxLevel = 101

F0 = 300
lowHarm = 1
highHarm = 2
spectrumLevel = 30
dur = 1
rampDur = 0.01
stretch = 0

channelOpts = ["right", "left"]
bandwidthUnitOpts = ["Hz", "Cent", "ERB"]
dichoticDifferenceOpts = ["IPD linear", "IPD stepped", "IPD random", "ITD"]
phaseRelationshipOpts = ["NoSpi", "NpiSo"]
noiseTypeOpts = ["white", "pink"]

for channel in channelOpts
    for bandwidthUnit in bandwidthUnitOpts
        for dichoticDifference in dichoticDifferenceOpts
            for phaseRelationship in phaseRelationshipOpts
                for noiseType in noiseTypeOpts
                    if bandwidthUnit == "Hz"
                        bandwidth=100
                    elseif bandwidthUnit == "Cent"
                        bandwidth = 300
                    elseif bandwidthUnit == "ERB"
                        bandwidth = 1
                    end
                    if dichoticDifference == "ITD"
                        dichoticDifferenceValue = 300/1000000
                    else
                        dichoticDifferenceValue = pi
                    end
                    hp = hugginsPitch(F0=F0, lowHarm=lowHarm,
                                      highHarm=highHarm,
                                      spectrumLevel=spectrumLevel,
                                      bandwidth=bandwidth,
                                      bandwidthUnit=bandwidthUnit,
                                      dichoticDifference=dichoticDifference,
                                      dichoticDifferenceValue=dichoticDifferenceValue,
                                      phaseRelationship=phaseRelationship,
                                      stretch=stretch,
                                      noiseType=noiseType,
                                      dur=dur, rampDur=rampDur,
                                      sf=sf, maxLevel=maxLevel)
                    wavwrite(hp, wavDir*"huggins_"*channel*"_"*bandwidthUnit*"_"*dichoticDifference*"_"*phaseRelationship*"_"*noiseType*".wav", Fs=sf, nbits=nbits)
                end
            end
        end
    end
end


## Test error condition
@test_throws(ErrorException, hugginsPitch(dur=1, rampDur=0.6))

@test_throws(ErrorException, hugginsPitch(bandwidthUnit="foo"))

@test_throws(ErrorException, hugginsPitch(phaseRelationship="foo"))

@test_throws(ErrorException, hugginsPitch(noiseType="foo"))
