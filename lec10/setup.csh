#!/bin/csh

if(! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH

setenv BASEDIR ${HOME}

#change to your root setup!
setenv ROOTSYS /usr/local/root_v5.32_binary_m64/
setenv PATH ${ROOTSYS}/bin:${PATH}
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
setenv DYLD_LIBRARY_PATH ${ROOTSYS}/lib:${DYLD_LIBRARY_PATH}

setenv GSL /opt/local
setenv LD_LIBRARY_PATH ${GSL}/lib:${LD_LIBRARY_PATH}
setenv DYLD_LIBRARY_PATH ${GSL}/lib:${DYLD_LIBRARY_PATH}

if ($?TERM == 0 || $?prompt == 0) exit 0

echo ''
echo 'Setup GSL and ROOT'
echo '====================================='
echo ''
echo "<I>---------------Info--------------------<I>"
echo "Setting up the following environments: "
echo "ROOT: " $ROOTSYS
echo "GSL : " $GSL
#echo "PYTHIA8: " $PYTHIA8DIR
#echo "ktJet: " $ktJetDIR
#echo "FastJet: " $FASTJETDIR
echo "<I>---------------Info--------------------<I>"
echo ""
