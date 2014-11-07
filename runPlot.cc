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

using namespace std;

// bin/runPlot -s 1 -c "PTLO=400;PTHI=10000;G0PT=100;g1PT=40"

// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  string progName  = argv[0]; 
  string dir("hpt0"), 
    ifiles("plotHpt.files"), 
    setup(""), 
    cuts("nada");
  
  int mode(-1), ntoy(200); 
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-c"))  {cuts   = string(argv[++i]);}
    if (!strcmp(argv[i], "-d"))  {dir    = string(argv[++i]);}
    if (!strcmp(argv[i], "-f"))  {ifiles = string(argv[++i]);}
    if (!strcmp(argv[i], "-m"))  {mode   = atoi(argv[++i]);}
    if (!strcmp(argv[i], "-n"))  {ntoy   = atoi(argv[++i]);}
    if (!strcmp(argv[i], "-s"))  {setup  = string(argv[++i]);}
  }


  gROOT->Clear();  gROOT->DeleteAll();
  plotHpt a(dir, ifiles, setup);
  if (cuts != "nada") {
    cout << "setting cuts: " << cuts << endl;
    a.setCuts(cuts);
  } 
  a.setNtoy(ntoy); 
  if (mode > -1) {
    a.makeAll(mode);
  } else {
    a.makeAll(1); 
  } 

}
