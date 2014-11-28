#ifndef PLOTHPT_h
#define PLOTHPT_h

#include "RooWorkspace.h"

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
  void   optimizeCuts(std::string fname = "opt.root", double lumi = 1000., int nevts = -1);
  void   optAnalysis(int mode = 1, std::string filename = "opt.root", std::string treename = "opt");
  void   massResolution(std::string hname = "mpt_mcatnlo5_nopt");

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
  void   toy4(double nsg0 = 50, double nsg1 = 100, double nbg = 450, int ntoy = 200); 

  // toy5 introduces the LLR
  //     - extended UML
  //     - signal is     Higgs(mtop -> infty)
  //     - background is diphoton + Higgs(mtop = 173.5GeV)
  void   toy5(double nsg0 = 50, double nsg1 = 100, double nbg = 450, int ntoy = 200); 

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
  void allNumbers(); 

  virtual void   bookHist(std::string name, std::string cuts); 
  void   readHistograms();
  void   setupTree(TTree *t); 
  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  void   candAnalysis(); 
  void   loopFunction1(); // treeAnalysis
  void   loopFunction2(); // optimizeCuts
  void   loopFunction3(); // bgShape

  TH1*   combMCAtNLOHist(TH1 *h0, TH1 *h1);

  void setNtoy(int ntoy) {fNtoy = ntoy;} 


private: 

  // -- essential analysis numbers
  double fSg0, fSg1, fBg; 
  double fNormSg0, fNormSg0E; 
  double fNormSg1, fNormSg1E; 
  double fNormBg;
  double fHiggsMpeak, fHiggsMres;
  double fBgMp1, fBgMp1E;
  double fBgTau, fBgTauE;
  double fSg0Tau, fSg0TauE;
  double fSg1Tau, fSg1TauE;
  
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
