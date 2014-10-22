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


  // -- Main analysis methods 
  void   makeAll(int bitmask = 0);
  void   treeAnalysis(int nevts = -1); 
  void   optimizeCuts(string fname = "opt.root", double lumi = 1000., int nevts = -1);
  void   optAnalysis(int mode = 1, std::string filename = "opt.root", std::string treename = "opt");
  void   massResolution(std::string hname = "mpt_mcatnlo5_nopt");

  // toy1 is an UML implementation: 2d UML fit of Higgs signal against a gg background
  //      - fixed m resolution
  //      - pT is modeled with RooHistPdf
  //      - no reference to normalization anywhere
  void   toy1(int nsg = 60, int nbg = 150); 

  // toy2 is another UML implementation: 2d UML fit of Higgs signal against a gg background
  //     - fixed m resolution
  //     - pT is modeled with RooHistPdf
  //     - extended UML
  void   toy2(int nsg = 60, int nbg = 150); 

  // toy3 is another UML implementation: 2d UML fit of Higgs signal against a gg background
  //     - fixed m resolution
  //     - pT is modeled with RooHistPdf
  //     - extended UML
  //  void   toy3(int nsg = 60, int nbg = 150); 

  // validation produces the simple plots
  //     - pT for various components
  //     - mass for various selections
  //     > all necessary plots are produced. This emulates overlayAll, in a sense...
  void   validation(); 


  void   bookHist(std::string name, std::string cuts); 
  void   readHistograms();
  void   setupTree(TTree *t); 
  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  void   candAnalysis(); 
  void   loopFunction1(); 
  void   loopFunction2(); 

  TH1*   combMCAtNLOHist(TH1 *h0, TH1 *h1);

private: 

  int    NBINS; 
  double GETA;
  double G0ISO, G1ISO; 
  double G0ISOR, G1ISOR; 
  double G0PT, G1PT; 
  double PTNL, PTNH; 
  double PTLO, PTHI; 
  double MGGLO, MGGHI; 

  double fSg0Xs, fSg1Xs, fBgXs, fBgRf, fSgBF;

  bool fGoodCand, fGoodCandNoPtCuts; 

  TH1D *fHistSg0, *fHistSg1, *fHistBg; 
  TH2D *fHistBg2;

  TTree* fTree;

  struct redTreeData fb; 

  RooWorkspace fWorkspace; 

  // -- optimization
  int                   fOptMode; 
  std::vector<selpoint> fSelPoints;

  // -- pt bins
  std::vector<ptbins>   fPtBins;

  // ----------------------------------------------------------------------
  ClassDef(plotHpt,1) 

};


#endif
