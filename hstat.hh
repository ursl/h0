#ifndef STAT_h
#define STAT_h

#include "TH1.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"

#include "plotClass.hh"
#include "dataset.hh"
#include "redTreeData.hh"

struct model {
  RooRealVar *m, *pt;   // variables: mass and pt
  // -- fit (fixed) parameters:
  RooRealVar *massP, *massS, *bgSlope;
  RooRealVar *sgTau, *bgTau;
  RooRealVar *sgN, *bgN;

  RooAbsPdf  *sgPdf, *bgPdf; 
  RooAbsPdf  *modelPdf; 

  RooArgSet poi;

};


// ----------------------------------------------------------------------
class hstat: public TObject {

public :
  hstat(std::string dir = "hpt0", std::string files = "plotHpt.files", std::string setup = "m");
  virtual        ~hstat();

  model* genModel1(int mode, double nsg, double tau);
  void   delModel(model*); 

  // -- various studies
  void run1(); 

  // -- 1D 
  void run1D(int ntoys = 1000, int mode = 0); 
  void toy1D(); 


  void findMidPoint(TH1D* hq0, TH1D* hq1, double &midpoint, double &tailprob);
  double oneSidedGaussianSigma(double prob);
  
private: 
  int NBINS;
  double MGGLO, MGGHI; 

  // -- essential analysis numbers
  double fSg0, fSg1, fBg; 
  double fHiggsMpeak, fHiggsMres;
  double fBgMp1, fBgMp1E; // pol1
  double fBgMc0, fBgMc0E; // chebychev
  double fBgTau, fBgTauE;
  double fSgTau, fSg0Tau, fSg0TauE;
  double fSg1Tau, fSg1TauE;

  int fNtoy; 

  double fD0M0, fD0M1; 
  double fD1M0, fD1M1; 

  std::string fSetup, fHistFileName, fTexFileName; 

  TFile *fHistFile; 

  // ----------------------------------------------------------------------
  ClassDef(hstat,1) 

};


#endif
