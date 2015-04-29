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

// bin/runPlot -s 1 -c "PTLO=400;PTHI=10000;G0PT=100;g1PT=40"
// bin/runPlot -s 0 -a hstat -m 0 -n 200

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
      a.run1D(ntoy, mode);
    }
  }
}
