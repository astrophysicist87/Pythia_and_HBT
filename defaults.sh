#=====================================
# Flags and options
#=====================================
# Options for codes
runPythia=true
useParallel=true
centralitySelectionInPythia=false
runHBTEG=true
runFitCF=true
runSV=true
#runBF=false

# system specifications
projectile="p"
target="p"
beamEnergy="13000.0"	#GeV
Nevents="10"
centralityClass="0-100%"

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

output_settings () {

echo 'runPythia='$runPythia
echo 'useParallel='$useParallel
echo 'centralitySelectionInPythia='$centralitySelectionInPythia
echo
echo 'projectile='$projectile
echo 'target='$target
echo 'beamEnergy='$beamEnergy
echo 'Nevents='$Nevents
echo 'centralityClass='$centralityClass
echo
echo 'QRefValue='$QRefValue
echo 'BEeffects='$BEeffects
echo 'BEEnhancementMode='$BEEnhancementMode
echo 'SetFragmentationVertices='$SetFragmentationVertices
echo 'SetPartonVertices='$SetPartonVertices
echo 'ThermalOnly='$ThermalOnly
echo
echo 'UseColorReconnection='$UseColorReconnection
echo 'UseRopeHadronization='$UseRopeHadronization
echo 'IncludeStringShoving='$IncludeStringShoving
echo 'IncludeFlavourRopesMechanism='$IncludeFlavourRopesMechanism

}
