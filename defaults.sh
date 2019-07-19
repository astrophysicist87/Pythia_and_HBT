#=====================================
# Flags and options
#=====================================
# Options for codes
runPythia=true
useParallel=true
centralitySelectionInPythia=false

# system specifications
projectile="p"
target="p"
beamEnergy="13000.0"	#GeV
Nevents="10"

# BE and related specifications
QRefValue="0.2"			#GeV
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
