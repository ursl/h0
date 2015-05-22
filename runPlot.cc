#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>

#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TRandom.h"
#include "TUnixSystem.h"

#include "plotHpt.hh"
#include "hstat.hh"

using namespace std;

// run plotHpt:          bin/runPlot -s 1 -c "PTLO=400;PTHI=10000;G0PT=100;g1PT=40"
// create all plots:     bin/runPlot -s m -m 0
// run run2D vanilla:    bin/runPlot -s 10 -a hstat -m 0 -n 1000
// run run1D vanilla:    bin/runPlot -s 0  -a hstat -m 0 -n 1000
// Run run1D sys no top: bin/runPlot -s 0  -a hstat -m 1 -n 1000
// run run1D sys scale:  bin/runPlot -s 0  -a hstat -m 3 -n 1000
// run run1D sys bg:     bin/runPlot -s 0  -a hstat -m 10 -n 1000
// run sigStudies nbg:   bin/runPlot -s 2  -a hstat -m 0 -n 1000
// run sigStudies lumi:  bin/runPlot -s 2  -a hstat -m 1 -n 1000
// run sigStudies repr:  bin/runPlot -s 2  -a hstat -m 2 -n 1000
// run ALL systematics:  bin/runPlot -s 1  -a hstat -n 2000
/*
 run default, all systematics and LUMI sig study:
 ------------------------------------------------
 bin/runPlot -s 1 -a hstat -m 10 -n 500 >& s1m10.log & 
 bin/runPlot -s 1 -a hstat -m 11 -n 500 >& s1m11.log & 
 bin/runPlot -s 1 -a hstat -m 12 -n 500 >& s1m12.log & 
 bin/runPlot -s 1 -a hstat -m 13 -n 500 >& s1m13.log & 
 bin/runPlot -s 2 -a hstat -m 11 -n 500 >& s2m11.log & 
 bin/runPlot -s 2 -a hstat -m 12 -n 500 >& s2m12.log & 
 bin/runPlot -s 2 -a hstat -m 13 -n 500 >& s2m13.log & 

*/
// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  string progName  = argv[0]; 
  string dir("hpt0"), 
    ifiles("plotHpt.files"), 
    setup(""), 
    cuts("nada"),
    ana("plotHpt");
  
  int mode(-1), ntoy(2000), rndms(111); 
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-a"))  {ana    = string(argv[++i]);}
    if (!strcmp(argv[i], "-c"))  {cuts   = string(argv[++i]);}
    if (!strcmp(argv[i], "-d"))  {dir    = string(argv[++i]);}
    if (!strcmp(argv[i], "-f"))  {ifiles = string(argv[++i]);}
    if (!strcmp(argv[i], "-m"))  {mode   = atoi(argv[++i]);}
    if (!strcmp(argv[i], "-n"))  {ntoy   = atoi(argv[++i]);}
    if (!strcmp(argv[i], "-r"))  {rndms  = atoi(argv[++i]);}
    if (!strcmp(argv[i], "-s"))  {setup  = string(argv[++i]);}
  }


  gROOT->Clear();  gROOT->DeleteAll();

  if (ana == "plotHpt") {
    plotHpt a(dir, ifiles, setup);
    if (cuts != "nada") {
      cout << "setting cuts: " << cuts << endl;
      a.setCuts(cuts);
    } 
    a.setNtoy(ntoy); 
    a.setRndmSeed(rndms); 
    if (mode > -1) {
      a.makeAll(mode);
    } else {
      a.makeAll(1); 
    } 
  } else if (ana == "hstat") {
    hstat a; 
    if ("0" == setup) {
      a.setRndmSeed(rndms); 
      a.run1D(ntoy, mode);
    } else if ("1" == setup) {
      a.setRndmSeed(rndms); 
      a.systematics(mode, ntoy);
    } else if ("2" == setup) {
      a.setRndmSeed(rndms); 
      a.sigStudies(mode, ntoy);
    } else if ("10" == setup) {
      a.setRndmSeed(rndms); 
      a.run2D(ntoy, mode);
    }
  }
}
