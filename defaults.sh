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

# BE and related specifications
QRefValue="0.2"			# GeV
BEeffects='off'
BEEnhancementMode='2'	# 0 - use fixed QRef
						# 1 - use ST interval with Gaussian BE form
						# 2 - use ST interval with spherical Bessel BE form
SetFragmentationVertices='on'
SetPartonVertices='off'
ThermalOnly='true'

# other flags
UseColorReconnection='off'
UseRopeHadronization='off'
IncludeStringShoving='off'
IncludeFlavourRopesMechanism='off'
#=====================================
