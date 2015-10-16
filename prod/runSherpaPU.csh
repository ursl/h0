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
echo "%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%"
echo "cat $JOB.dat"
cat $JOB.dat
echo "%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%"
ls -l 
date
/shome/ursl/h/sherpa/SHERPA-MC-1.4.0/bin/Sherpa -f $JOB.dat HEPMC2_SHORT_OUTPUT=$JOB
#test /bin/cp /shome/ursl/h/prod/jobs/testpu/gg-132-10000.hepmc .
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
echo "%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%"
echo "cat ./delphes_card_CMS_PileUp.tcl"
cat ./delphes_card_CMS_PileUp.tcl
echo "%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%"
echo "cat ./delphes_card_CMS.tcl"
cat ./delphes_card_CMS.tcl
echo "%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%"

mv ../$JOB.hepmc .
lcg-cp -b -D srmv2  "$STORAGE1/../MinBias.pileup" file:///`pwd`/MinBias.pileup
ls -l 
./DelphesHepMC ./delphes_card_CMS_PileUp.tcl $JOB-wpu.root ./$JOB.hepmc
ls -l 
./DelphesHepMC ./delphes_card_CMS.tcl $JOB-npu.root ./$JOB.hepmc
ls -l 

# ----------------------------------------------------------------------
# -- Save Output to SE
# ----------------------------------------------------------------------
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/swshare/emi/emi-wn/usr/lib64
echo "--> LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

echo "--> Save root  to SE: $STORAGE1/sherpa/wpu/"
echo "--> Save root  to SE: $STORAGE1/sherpa/npu/"
echo "--> srmcp command: $SRMCP"

set FILES=`ls $JOB*-wpu.root`
echo "Found the following output root files: $FILES"
foreach f ($FILES)
  echo lcg-del  "$STORAGE1/sherpa/wpu/$f"
  lcg-del -b -D srmv2 -l "$STORAGE1/sherpa/wpu/$f"
  echo lcg-cp    file:///`pwd`/$f "$STORAGE1/sherpa/wpu/$f"
  lcg-cp -b -D srmv2  file:///`pwd`/$f "$STORAGE1/sherpa/wpu/$f"
  echo lcg-ls     "$STORAGE1/sherpa/wpu/$f"
  lcg-ls -b -D srmv2 -l "$STORAGE1/sherpa/wpu/$f"
end

set FILES=`ls $JOB*-npu.root`
echo "Found the following output root files: $FILES"
foreach f ($FILES)
  echo lcg-del  "$STORAGE1/sherpa/npu/$f"
  lcg-del -b -D srmv2 -l "$STORAGE1/sherpa/npu/$f"
  echo lcg-cp    file:///`pwd`/$f "$STORAGE1/sherpa/npu/$f"
  lcg-cp -b -D srmv2  file:///`pwd`/$f "$STORAGE1/sherpa/npu/$f"
  echo lcg-ls     "$STORAGE1/sherpa/npu/$f"
  lcg-ls -b -D srmv2 -l "$STORAGE1/sherpa/npu/$f"
end

gzip $JOB*.hepmc
set FILES=`ls $JOB*.hepmc.gz`
echo "Found the following output hepmc files: $FILES"
foreach f ($FILES)
  echo lcg-del  "$STORAGE1/sherpa/hepmc/$f"
  lcg-del -b -D srmv2 -l "$STORAGE1/sherpa/hepmc/$f"
  echo lcg-cp    file:///`pwd`/$f "$STORAGE1/sherpa/hepmc/$f"
  lcg-cp -b -D srmv2  file:///`pwd`/$f "$STORAGE1/sherpa/hepmc/$f"
  echo lcg-ls     "$STORAGE1/sherpa/hepmc/$f"
  lcg-ls -b -D srmv2 -l "$STORAGE1/sherpa/hepmc/$f"
end

date

# BATCH END

echo "run: This is the end, my friend"
