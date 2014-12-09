#! /bin/csh -f

setenv JOB
setenv DECAY
setenv FILE1 $JOB-$DECAY.root
setenv STORAGE1
setenv MCATNLO 150
setenv HEPMCBASE "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/ursl/strebeli/ggHyy/store/mcanlo_$MCATNLO"


setenv SRMCP

echo "============================="
echo "--> Running PSIT3 job wrapper"
echo "============================="

# BATCH START

# ----------------------------------------------------------------------
# -- Setup production
# ----------------------------------------------------------------------
echo "--> Setup production"

echo $VO_CMS_SW_DIR
ls -l $VO_CMS_SW_DIR
source $VO_CMS_SW_DIR/cmsset_default.csh
echo "-> which lcg-ls"
which lcg-ls
echo "-> which globus-url-copy"
which globus-url-copy
echo "-> which srmcp"
which srmcp

pwd
date
setenv ROOTSYS /shome/ursl/root
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
setenv PATH ${ROOTSYS}/bin:${PATH}

echo "--> environment"
df -kl 
printenv

echo "--> Extract tar file"
date
tar zxf $JOB.tar.gz

echo "--> running DELPHES"
date
cd ./delphes
lcg-cp -b -D srmv2 "$HEPMCBASE"/$JOB-mcanlo_$MCATNLO/$JOB-MC.hepmc file:///`pwd`/$JOB-MC.hepmc

./DelphesHepMC ./delphes_card_CMS.tcl $JOB-delphes.root ./$JOB-MC.hepmc
ls -l 

# ----------------------------------------------------------------------
# -- Save Output to NFS, not the SE
# ----------------------------------------------------------------------

# several files possible if tree grows to big. copy them all

echo "--> Save output to SE: $STORAGE1/$FILE1"
echo $SRMCP 

set FILES=`ls $JOB*.root`
echo "Found the following output root files: $FILES"
foreach f ($FILES)
  echo lcg-del  "$STORAGE1/$f"
  lcg-del -b -D srmv2 -l "$STORAGE1/$f"
  echo lcg-cp    file:///`pwd`/$f "$STORAGE1/$f"
  lcg-cp -b -D srmv2  file:///`pwd`/$f "$STORAGE1/$f"
  echo lcg-ls     "$STORAGE1/$f"
  lcg-ls -b -D srmv2 -l "$STORAGE1/$f"
end

date

# BATCH END

echo "run: This is the end, my friend"
