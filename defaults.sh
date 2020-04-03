#=====================================
# Flags and options
#=====================================
# Options for codes
export runPythia=true
#export useParallel=true
export useArbitraryParticle=false
export centralitySelectionInPythia=false
export runHBTEG=true
export runFitCF=true
export runSV=true
#export runBF=false

# system specifications
export projectile="p"
export target="p"
export beamEnergy="13000.0"                     # GeV
export Nevents="10"
export chosenHBTparticle="211"                  # pion(+)
export centralityClass="0-100%"
export bMin="0.0"
export bMax="20.0"

# BE and related specifications
export QRefValue="0.2"                          # GeV
export BEeffects='off'
export BEEnhancementMode='1'                    # 0 - use fixed QRef
                                         # 1 - use ST interval with spherical Bessel BE form
export SetFragmentationVertices='on'
export SetPartonVertices='off'
export ThermalOnly='false'

export useInvariantSourceSize='off'             # Lorentz-invariant size vs. spatial size only
export useDistribution='off'                    # Estimate QRef vs. take as input parameter
export useRelativeDistance='on'                 # Use relative distances or absolute sizes
export useRestFrame='on'                        # Use rest frame vs. lab frame
export includePhaseSpace='on'                   # Include phase-space factor
export linearInterpolateCDF='on'                # Estimate pair density via linear interpolation
export usePositiveShiftsForCompensation='on'    # Pairs shifted apart used to compensate pairs shifted together
export computeBEEnhancementExactly='on'         # Whether to evaluate BE enhancement approximately or exactly


# other flags
export UseColorReconnection='off'
export UseRopeHadronization='off'
export IncludeStringShoving='off'
export IncludeFlavourRopesMechanism='off'
export storeBjorkenCoordinates='false'
#=====================================
