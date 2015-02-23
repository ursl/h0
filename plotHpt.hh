#ifndef PLOTHPT_h
#define PLOTHPT_h

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"

#include "plotClass.hh"
#include "dataset.hh"
#include "redTreeData.hh"

struct ptbins {
  double peak, sigma;
  double pt, ptlo, pthi;
};

// ----------------------------------------------------------------------
class plotHpt: public plotClass {

public :
  plotHpt(std::string dir = "hpt0", std::string files = "plotHpt.files", std::string setup = "m");
  virtual        ~plotHpt();

  void   loadFiles(std::string afiles);
  void   setCuts(std::string cuts);


  // -- Main analysis methods 
  void   makeAll(int bitmask = 0);
  virtual void treeAnalysis(int nevts = -1); 
  void   massResolution(std::string hname = "mpt_mcatnlo5_nopt");
  // -- these two are historic 
  void   optimizeCuts(std::string fname = "opt.root", double lumi = 1000., int nevts = -1);
  void   optAnalysis(int mode = 1, std::string filename = "opt.root", std::string treename = "opt");
  // -- this parses the output of running the complete analysis
  void   anaOptimization(); 


  // toy1 is an UML implementation: 2d UML fit of Higgs signal against a gg background
  //      - fixed m resolution
  //      - pT is modeled with RooHistPdf
  //      - no reference to normalization anywhere
  void   toy1(int nsg = 50, int nbg = 500); 

  // toy2 is another UML implementation: 2d UML fit of Higgs signal against a gg background
  //     - fixed m resolution
  //     - pT is modeled with RooHistPdf
  //     - extended UML
  void   toy2(int nsg = 50, int nbg = 500); 

  // toy3 is another UML implementation: 2d UML fit of contact Higgs signal against a gg+SM Higgs bg
  //     - signal is     Higgs(mtop -> infty)
  //     - background is diphoton + Higgs(mtop = 173.5GeV)
  void   toy3(int nsg = 100, double nbgH = 50, double nbg = 450); 

  // toy4 adds a sensitivity calculation
  //     - signal is     Higgs(mtop -> infty)
  //     - background is diphoton + Higgs(mtop = 173.5GeV)
  void   toy4(double nsg0 = 50, double nsg1 = 100, double nbg = 450, int ntoy = -1); 

  // toy5 introduces the delta LogLikelihood as a significance estimate using RooFit tools
  //     - extended UML for mass ONLY
  //     - signal is     Higgs(mtop -> infty)
  //     - background is diphoton + Higgs(mtop = 173.5GeV)
  void   toy5(double nsg0 = 50, double nsg1 = 100, double nbg = 450, int ntoy = 200); 

  // -- this has nothing to do with Higgs!
  void   toy6();   
  void   toy7(); 
  void   toy8(); 

  void   toy9(int nsig = 50, int nbkg = 1000); 
  void   toy10(int nsig0 = 50, int nbkg0 = 1000);

  void   toy11(int nb0 = 1000, double start = 0.2);
  void   toy12(int nb0 = 1000, double start = 0.2);
  
  // validation produces the simple plots
  //     - pT for various components
  //     - mass for various selections
  //     > all necessary plots are produced. This emulates overlayAll, in a sense...
  void   validation(); 

  // study the background mass shape in different pT bins
  void   bgShape(int nevts = -1);
  void   bgSyst(std::string ds0, std::string ds1, std::string ds2 = "nada");

  // analysis: Everything after treeAnalysis and before toyX()
  //     - event yields in lumi
  //     - hi-pt diphoton resolution
  //     - all required histograms
  void allNumbers(int ntoy = 100); 

  // -- 1D eUML mass fit with DLL
  void allNumbers1(int ntoy = -1); 

  // -- 2D eUML mass-pT fit wtih DLL
  void allNumbers2(int ntoy = -1); 

  // -- now start to factorize the preprocessing
  void fitMass(TH1D *hs, double &peak, double &reso); 
  void ptDistributions();
  void displayCuts();
  void massHistograms(TH1D *hb, TH1D *hs0, TH1D *hs1);

  // -- 2D eUML mass-pT fit with LLR
  void allNumbers3(int ntoy = -1); 
  // -- 2D eUML mass-pT fit with PROFILED LLR, try the setup of 1310.1397
  void allNumbers4(int ntoy = -1); 

  // -- 2D eUML mass-pT with the official RooStats tools (this follows now after toy10()
  void allNumbers5(int ntoy = -1); 


  virtual void   bookHist(std::string name, std::string cuts); 
  void   readHistograms();
  void   setupTree(TTree *t); 
  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  void   candAnalysis(); 
  void   loopFunction1(); // treeAnalysis
  void   loopFunction2(); // optimizeCuts
  void   loopFunction3(); // bgShape
  // hq0 = SM (model0), hq1 = model1
  void   findMidPoint(TH1D* hq0, TH1D* hq1, double &midpoint, double &tailprob);
  double oneSidedGaussianSigma(double prob); 

  TH1*   combMCAtNLOHist(TH1 *h0, TH1 *h1);

  void setNtoy(int ntoy) {fNtoy = ntoy;} 
  void shutup();

private: 

  // -- essential analysis numbers
  double fMu, fSg, fSg0, fSg1, fBg; 
  double fHiggsMpeak, fHiggsMres;
  double fBgMp1, fBgMp1E; // pol1
  double fBgMc0, fBgMc0E; // chebychev
  double fBgTau, fBgTauE;
  double fSgTau, fSg0Tau, fSg0TauE;
  double fSg1Tau, fSg1TauE;

  // -- and the corresponding RooVars
  void iniRooVars(); 
  RooRealVar *fRm, *fRpt; // the variables
  RooRealVar *fRsgP, *fRsgS; // signal mass peak and sigma
  RooRealVar *fRsg0Tau, *fRsg1Tau; // signal pT slope
  RooRealVar *fRC0, *fRbg0Tau, *fRbg1Tau; // bg mass and pt shapes
  RooRealVar *fRbg0Slope, *fRbg1Slope, *fRbg2Slope; // bg mass and pt shapes
  RooRealVar *fRmu, *fRsg0N, *fRsg1N, *fRsg2N, *fRbg0N, *fRbg1N, *fRbg2N; // the event numbers for the extended LH
  RooRealVar *fRbgN, *fRbgTau, *fRbgSlope;
  RooFormulaVar *fRsgN, *fRsgTau;

  RooRealVar *fRsg3N, *fRsg3Tau, *fRbg3N, *fRbg3Tau, *fRbg3Slope;
  RooRealVar *fRsg4N, *fRsg4Tau, *fRbg4N, *fRbg4Tau, *fRbg4Slope;

  
  double fNormSg0, fNormSg0E; 
  double fNormSg1, fNormSg1E; 
  double fNormBg;
  double fNormHiggsMpeak, fNormHiggsMres;
  
  // -- setup
  int    NBINS, NPTBINS; 

  // -- cuts
  double GETA;
  double G0ISO, G1ISO; 
  double G0ISOR, G1ISOR; 
  double G0PT, G1PT; 
  double G0PTLO, G1PTLO; 
  double PTNL, PTNH; 
  double PTLO, PTHI; 
  double MGGLO, MGGHI; 

  // -- variables
  double fptnl, fptnh, fptlo, fpthi, fg0ptlo, fg1ptlo, fg0pt, fg1pt; 

  bool fGoodCand, fGoodCandNoPtCuts; 

  TTree* fTree;

  struct redTreeData fb; 

  RooWorkspace fWorkspace; 

  // -- optimization
  int                   fOptMode; 
  std::vector<selpoint> fSelPoints;

  // -- pt bins
  std::vector<ptbins>   fPtBins;


  int fNtoy; 

  // ----------------------------------------------------------------------
  ClassDef(plotHpt,1) 

};


#endif
