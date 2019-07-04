#=====================================
# Flags and options
#=====================================
runPythia=true

# system specifications
projectile="p"
target="p"
beamEnergy="13000.0"	#GeV
Nevents="1"

# BE and related specifications
QRefValue="0.2"			#GeV
BEeffects='on'
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
