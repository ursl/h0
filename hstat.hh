#ifndef STAT_h
#define STAT_h

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"

#include "plotClass.hh"
#include "dataset.hh"
#include "redTreeData.hh"


// ----------------------------------------------------------------------
class hstat: public TObject {

public :
  hstat(std::string dir = "hpt0", std::string files = "plotHpt.files", std::string setup = "m");
  virtual        ~hstat();

  RooAbsPdf* genModel(int i, int ini);
  void sigNoBackground(); 
  void genData(); 

private: 

  // -- essential analysis numbers
  double fSg0, fSg1, fBg; 
  double fHiggsMpeak, fHiggsMres;
  double fBgMp1, fBgMp1E; // pol1
  double fBgMc0, fBgMc0E; // chebychev
  double fBgTau, fBgTauE;
  double fSgTau, fSg0Tau, fSg0TauE;
  double fSg1Tau, fSg1TauE;

  // -- and the corresponding RooVars
  RooRealVar *fRm, *fRpt; // the variables
  RooRealVar *fRsgP, *fRsgS; // signal mass peak and sigma

  // -- variables for various models
  RooRealVar *fRsg0N, *fRbg0N, *fRbg0Slope;;
  RooRealVar *fRsg1N, *fRsg1Tau, *fRbg1N, *fRbg1Tau, *fRbg1Slope;
  RooRealVar *fRsg2N, *fRsg2Tau, *fRbg2N, *fRbg2Tau, *fRbg2Slope;
  RooRealVar *fRsg3N, *fRsg3Tau, *fRbg3N, *fRbg3Tau, *fRbg3Slope;
  RooRealVar *fRsg4N, *fRsg4Tau, *fRbg4N, *fRbg4Tau, *fRbg4Slope;

  RooDataSet *bgData, *sg0Data, *sg1Data, *data1, *data0;


  int fNtoy; 

  std::string fSetup, fHistFileName, fTexFileName; 

  TFile *fHistFile; 

  // ----------------------------------------------------------------------
  ClassDef(hstat,1) 

};


#endif
