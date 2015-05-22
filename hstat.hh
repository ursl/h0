#ifndef STAT_h
#define STAT_h

#include "TH1.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooPolynomial.h"

#include "plotClass.hh"
#include "dataset.hh"
#include "redTreeData.hh"

#define EXPO
#undef EXPO

struct model {
  RooRealVar *m, *pt;   // variables: mass and pt
  // -- fit (fixed) parameters:
  RooRealVar *massP, *massS, *bgSlope;
  RooRealVar *sgTau, *bgTau, *sgMu, *bgMu;
  RooRealVar *sgN, *bgN;

#ifdef EXPO
  RooExponential *sgPt, *bgPt;
#else 
  RooLognormal *sgPt, *bgPt;
#endif
  RooPolynomial  *bgM;
  RooGaussian    *sgM; 

  RooAbsPdf  *sgPdf, *bgPdf; 
  RooAbsPdf  *modelPdf; 

  RooArgSet poi;

};


class TGraph;
// ----------------------------------------------------------------------
class hstat: public TObject {

public :
  hstat(std::string dir = "hpt0", std::string files = "plotHpt.files", std::string setup = "m");
  virtual        ~hstat();

  model* genModel1(int mode, double nsg, double nbg);
  model* genModel2(int mode, double nsg, double nbg, double tsg, double tbg, double msg = 0, double mbg = 0);
  void   delModel(model*); 

  // -- various studies
  void run1(); 

  // -- 1D 
  void systematics(int mode = 0, int n = 1000); 
  void sigStudies(int mode = 0, int n = 1000);
  void run1D(int ntoys = 1000, int mode = 0); 
  void toy1D(); 

  // -- 2D 
  void run2D(int ntoys = 1000, int mode = 0); 
  void toy2D(); 
  void plotResults(); 
  void setGraph(TGraph *, int color, int fillStyle = 1000); 

  void setRndmSeed(int rndms) {fRndmSeed = rndms;}

  void testTools(int nbins = 400); 
  void findMidPoint(TH1D* hq0, TH1D* hq1, double &midpoint, double &tailprob, double &lo, double &hi);
  double oneSidedGaussianSigma(double prob);
  void expModification(int evt = 10000);
  
private: 
  int NBINS;
  double MGGLO, MGGHI; 

  int fRndmSeed; 

  // -- essential analysis numbers
  double fSG0, fSG1, fBG; //unmutable const default parameters!
  double fSG0TAU, fSG1TAU, fBGTAU; //unmutable const default parameters!
  double fSG0MU, fSG1MU, fBGMU; //unmutable const default parameters!
  
  double fSg0, fSg1, fBg; 
  double fHiggsMpeak, fHiggsMres;
  double fBgMp1, fBgMp1E; // pol1
  double fBgMc0, fBgMc0E; // chebychev
  double fBgTau, fBgTauE;
  double fBgMu, fBgMuE;
  double fSgTau, fSg0Tau, fSg0TauE;
  double fSg1Tau, fSg1TauE;
  double fSg0Mu, fSg0MuE;
  double fSg1Mu, fSg1MuE;

  int fNtoy; 

  double fD0M0, fD0M1; 
  double fD1M0, fD1M1; 
  double fSeparation, fSeparationE; 

  std::string fSetup, fHistFileName, fTexFileName; 

  bool fDoPlot; 
  std::string fDoPlotName;
  std::string fDirectory; 
  TFile *fHistFile; 
  ofstream fTEX; 

  // ----------------------------------------------------------------------
  ClassDef(hstat,1) 

};


#endif
