#ifndef PLOTHPT_h
#define PLOTHPT_h

#include "RooWorkspace.h"

#include "plotClass.hh"
#include "dataset.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotHpt: public plotClass {

public :
  plotHpt(std::string dir = "hpt0", std::string files = "plotHpt.files", std::string setup = "m");
  virtual        ~plotHpt();

  virtual void   loadFiles(std::string afiles);


  // -- Main analysis methods 
  virtual void   makeAll(int bitmask = 0);
  virtual void   treeAnalysis(); 

  void           toy1(int nsg = 60, int nbg = 150); 

  virtual void   overlayAll();

  virtual void   bookHist(std::string name); 
  virtual void   setupTree(TTree *t); 
  virtual void   loopOverTree(TTree *t, int nevts = -1, int nstart = 0); 
  void           candAnalysis(); 
  void           loopFunction(); 

private: 

  int    NBINS; 
  double GETA;
  double G0ISO, G1ISO; 
  double G0PT, G1PT; 
  double PTLO, PTHI; 
  double MGGLO, MGGHI; 

  double fSg0Xs, fSg1Xs, fBgXs, fBgRf, fSgBF;

  bool fGoodCand; 

  TH1D *fHistSg0, *fHistSg1, *fHistBg; 
  TH2D *fHistBg2;

  TTree* fTree;

  struct redTreeData fb; 

  RooWorkspace fWorkspace; 

  int             fOptMode; 

  // ----------------------------------------------------------------------
  ClassDef(plotHpt,1) 

};


#endif
