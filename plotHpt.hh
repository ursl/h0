#ifndef PLOTHPT_h
#define PLOTHPT_h


#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include "RooWorkspace.h"

#include "dataset.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotHpt: public TObject {

public :
  plotHpt(std::string dir = "hpt0", std::string files = "plotHpt.files", std::string setup = "m");
  ~plotHpt();

  void   loadFiles(std::string afiles);
  TFile* loadFile(std::string afiles);


  // -- Main analysis methods 
  void makeAll(int bitmask = 0);
  void treeAnalysis(); 
  void toy1(int nsg = 60, int nbg = 150); 

  void overlayAll();
  void overlay(TH1D* h1, std::string f1, TH1D *h2, std::string f2, bool legend = true);
  void overlay(std::string f1, std::string h1name, std::string f2, std::string h2name, bool legend = true);

  void bookHist(std::string name); 
  TTree* getTree(std::string ds); 
  void setupTree(TTree *t); 
  void loopOverTree(TTree *t, int nevts = -1, int nstart = 0); 
  void candAnalysis(); 
  void loopFunction(); 

  void cd(std::string dataset) {fDS[dataset]->cd("");}
  void replaceAll(std::string &sInput, const std::string &oldString, const std::string &newString);
  void newLegend(double x1, double y1, double x2, double y2, std::string title = "");
  void makeCanvas(int i = 3);
  void normHist(TH1D *, double integral = -1., std::string type=""); 

private: 
  int    fVerbose; 
  double fEpsilon; 

  std::string fDirectory, fSetup, fSuffix;   

  int    NBINS; 
  double GETA;
  double G0ISO, G1ISO; 
  double G0PT, G1PT; 
  double PTLO, PTHI; 
  double MGGLO, MGGHI; 

  double fLumi, fSg0Xs, fSg1Xs, fBgXs, fBgRf, fSgBF;

  bool fGoodCand; 

  TH1D *fHistSg0, *fHistSg1, *fHistBg; 
  TH2D *fHistBg2;

  TTree* fTree;

  struct redTreeData fb; 

  RooWorkspace fWorkspace; 

  int             fOptMode; 
  std::map<std::string, TH1*> fHists;

  // -- datasets (files and associated information)
  std::map<std::string, dataset*> fDS; 
  // -- current dataset for analysis
  std::string fCds; 

  // -- Display utilities
  int fFont; 
  double fSize; 
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5;
  TLatex *tl; 
  TBox *box;
  TArrow *pa; 
  TLine *pl; 
  TLegend *legg;
  TLegendEntry *legge;


  // ----------------------------------------------------------------------
  ClassDef(plotHpt,1) 

};


#endif
