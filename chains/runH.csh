#! /bin/csh -f

setenv JOB
setenv EXECUTABLE
setenv FILE1 $JOB.root
setenv STORAGE1

setenv SRMCP


echo "========================"
echo "====> PSIT3 wrapper <===="
echo "========================"

echo "--> Running PSIT3 job wrapper"

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
# -- 
setenv LD_LIBRARY_PATH /swshare/cms/slc5_amd64_gcc462/external/gcc/4.6.2/lib64:/shome/ursl/h/HepMC-2.06.0/lib
setenv LD_LIBRARY_PATH /swshare/emi-wn-2.5.1-1_v1/emi-wn/usr/lib64:${LD_LIBRARY_PATH}
source $VO_CMS_SW_DIR/cmsset_default.csh

setenv ROOTSYS /shome/ursl/root
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:/shome/ursl/macros/lib:${LD_LIBRARY_PATH}
setenv PATH ${ROOTSYS}/bin:${PATH}

echo "--> environment"
df -kl 
printenv

echo "--> Extract tar file"
date
tar zxf $JOB.tar.gz

echo "--> running analysis"
pwd
ls -l 
date
echo "$EXECUTABLE -c $JOB  -o $FILE1 |& tee $JOB.log"
$EXECUTABLE -c $JOB -o $FILE1 |& tee $JOB.log
pwd 
ls -l 
date

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

# copy the log file.
# echo lcg-del  "$STORAGE1/$JOB.log"
# lcg-del -b -D srmv2 -l "$STORAGE1/$JOB.log"
# echo lcg-cp file:///`pwd`/$JOB.log "$STORAGE1/$JOB.log"
# lcg-cp -b -D srmv2     file:///`pwd`/$JOB.log "$STORAGE1/$JOB.log"
# echo lcg-ls     "$STORAGE1/$JOB.log"
# lcg-ls -b -D srmv2 -l  "$STORAGE1/$JOB.log"

date

# BATCH END


# -- cleanup
#/bin/rm -rf /scratch/ursl/dana*

echo "run: This is the end, my friend"
