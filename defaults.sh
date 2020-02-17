#=====================================
# Flags and options
#=====================================
# Options for codes
runPythia=true
#useParallel=true
useArbitraryParticle=false
centralitySelectionInPythia=false
runHBTEG=true
runFitCF=true
runSV=true
#runBF=false

# system specifications
projectile="p"
target="p"
beamEnergy="13000.0"	# GeV
Nevents="10"
chosenHBTparticle="211"	# pion(+)
centralityClass="0-100%"
bMin="0.0"
bMax="20.0"

# BE and related specifications
QRefValue="0.2"			# GeV
BEeffects='off'
BEEnhancementMode='1'	# 0 - use fixed QRef
						# 1 - use ST interval with spherical Bessel BE form
SetFragmentationVertices='on'
SetPartonVertices='off'
ThermalOnly='false'

useInvariantSize='false'				# Lorentz-invariant size vs. spatial size only
useDistribution='false'					# Estimate QRef vs. take as input parameter
useRelativeDistance='true'				# Use relative distances or absolute sizes
useRestFrame='true'						# Use rest frame vs. lab frame
includePhaseSpace='true'				# Include phase-space factor
linearInterpolateCDF='true'				# Estimate pair density via linear interpolation
usePositiveShiftsForCompensation='true'	# Pairs shifted apart used to compensate pairs shifted together
computeBEEnhancementExactly='true'		# Whether to evaluate BE enhancement approximately or exactly


# other flags
UseColorReconnection='off'
UseRopeHadronization='off'
IncludeStringShoving='off'
IncludeFlavourRopesMechanism='off'
storeBjorkenCoordinates='false'
#=====================================
