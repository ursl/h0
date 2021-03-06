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
echo "--> Setup production for job $JOB"
date

echo $VO_CMS_SW_DIR
ls -l $VO_CMS_SW_DIR
source $VO_CMS_SW_DIR/cmsset_default.csh
echo "-> which lcg-ls"
which lcg-ls
echo "-> which globus-url-copy"
which globus-url-copy
echo "-> which srmcp"
which srmcp
echo "-> which lcg-cp"
setenv BLA  `which lcg-cp`
echo "-> ldd $BLA"
ldd $BLA

pwd
ls -l 

setenv HEPMCLOCATION /shome/ursl/h/HepMC-2.06.0
setenv LD_LIBRARY_PATH /swshare/cms/slc5_amd64_gcc462/external/gcc/4.6.2/lib64:/usr/local/lib:/usr/lib:/usr/X11/lib:
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/shome/ursl/h/yaml-cpp-0.5.1/lib

setenv ROOTSYS /shome/ursl/root
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
setenv PATH ${ROOTSYS}/bin:${PATH}


echo "--> environment"
df -kl 
printenv

echo "--> running SHERPA"
/bin/cp /shome/ursl/h/prod/jobs/sherpa/dat/$JOB.dat .
ls -l 
date
/shome/ursl/h/sherpa/SHERPA-MC-1.4.0/bin/Sherpa -f $JOB.dat HEPMC2_SHORT_OUTPUT=$JOB
date
ls -l 

echo "--> Extract tar file"
date
ls -l 
tar zxvf $JOB.tar.gz
ls -l 
if (! -d "delphes") then
  echo "ERROR could not extract directory delphes from tar file"
  exit
else 
  echo "extracted directory delphes from tar file, contents:"
  ls -l delphes
endif

echo "--> running DELPHES"
date
cd ./delphes
cat ./delphes_card_CMS.tcl
mv ../$JOB.hepmc .
./DelphesHepMC ./delphes_card_CMS.tcl $JOB.root ./$JOB.hepmc
ls -l 

# ----------------------------------------------------------------------
# -- Save Output to NFS, not the SE
# ----------------------------------------------------------------------
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/swshare/emi/emi-wn/usr/lib64
echo "--> LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

echo "--> Save hepmc to SE: $STORAGE1/hepmc"
echo "--> Save root  to SE: $STORAGE1/delphes"
echo "--> srmcp command: $SRMCP"

# several files possible if tree grows to big. copy them all
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
