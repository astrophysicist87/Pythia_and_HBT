#! /usr/bin/env bash

#=====================================
# Flags and options
#=====================================
# Options for codes
export runPythia=true
export eventClassSelectionInPythia=false
export runHBTEG=true
export runFitCF=true
export runSV=false
#export runBF=false

# system specifications
export projectile="p"
export target="p"
export beamEnergy="13000.0"                     # GeV
export Nevents="10"
export chosenHBTparticle="211"                  # pion(+)
export eventClassSelectionMode="centrality"
export eventClass="0-100%"
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
export computeBEEnhancementExactly='on'         # Whether to evaluate BE enhancement approximately or exactly

export shiftingSet='1'							# Which pairs to shift
export compensationSet='0'						# Which pairs to use in compensation
export compensationMode='1'						# How to compute compensation
export compensationVersion='0'					# Version for implementing compensation

# other flags
export UseColorReconnection='on'
export UseRopeHadronization='off'
export IncludeStringShoving='off'
export IncludeFlavourRopesMechanism='off'
export storeBjorkenCoordinates='false'
#=====================================
