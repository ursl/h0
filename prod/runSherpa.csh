#! /bin/csh -f

setenv JOB
setenv DECAY
setenv FILE1 $JOB-$DECAY.root
setenv STORAGE1

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
/bin/cp /shome/ursl/h/prod/jobs/sherpa/dat/$JOB.dat .

echo "--> running SHERPA"
date
/shome/ursl/h/sherpa/SHERPA-MC-1.4.0/bin/Sherpa -f $JOB.dat HEPMC2_SHORT_OUTPUT=$JOB
ls -l 

echo "--> running DELPHES"
date
cd ./delphes
mv ../$JOB.hepmc .
./DelphesHepMC ./delphes_card_CMS.tcl $JOB.root ./$JOB.hepmc
ls -l 

# ----------------------------------------------------------------------
# -- Save Output to NFS, not the SE
# ----------------------------------------------------------------------

# several files possible if tree grows to big. copy them all

echo "--> Save hepmc to SE: $STORAGE1/hepmc/$FILE1"
echo "--> Save root  to SE: $STORAGE1/delphes/$FILE1"
echo $SRMCP 

set FILES=`ls $JOB*.root`
echo "Found the following output root files: $FILES"
foreach f ($FILES)
  echo lcg-del  "$STORAGE1/delphes/$f"
  lcg-del -b -D srmv2 -l "$STORAGE1/delphes/$f"
  echo lcg-cp    file:///`pwd`/$f "$STORAGE1/delphes/$f"
  lcg-cp -b -D srmv2  file:///`pwd`/$f "$STORAGE1/delphes/$f"
  echo lcg-ls     "$STORAGE1/delphes/$f"
  lcg-ls -b -D srmv2 -l "$STORAGE1/delphes/$f"
end

set FILES=`ls $JOB*.hepmc`
echo "Found the following output hepmc files: $FILES"
foreach f ($FILES)
  echo lcg-del  "$STORAGE1/hepmc/$f"
  lcg-del -b -D srmv2 -l "$STORAGE1/hepmc/$f"
  echo lcg-cp    file:///`pwd`/$f "$STORAGE1/hepmc/$f"
  lcg-cp -b -D srmv2  file:///`pwd`/$f "$STORAGE1/hepmc/$f"
  echo lcg-ls     "$STORAGE1/hepmc/$f"
  lcg-ls -b -D srmv2 -l "$STORAGE1/hepmc/$f"
end

date

# BATCH END

echo "run: This is the end, my friend"
