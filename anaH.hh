#ifndef ANAH_H
#define ANAH_H

#include <iostream>
#include <string>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "delphes/ExRootTreeReader.hh"
#include "delphes/DelphesClasses.h"

#include "redTreeData.hh"


// ----------------------------------------------------------------------
// -- class for analysis of Delphes trees
//    output is displayed by plotH
// ----------------------------------------------------------------------


class anaH {
public:
  anaH(TChain *tree);
  virtual      ~anaH();
  virtual void init();
  
  virtual void openHistFile(std::string filename);
  virtual void closeHistFile();
  virtual void bookHist();
  virtual void setupReducedTree();
  
  virtual void startAnalysis();
  virtual void endAnalysis();
  virtual int  loop(int nevents = 1, int start = -1);
  virtual void eventProcessing();
  virtual void initVariables(); 
  virtual void genLevelAnalysis();
  virtual void ggAnalysis();
  virtual void llAnalysis();
  virtual void ggReco();
  virtual void llReco();
  virtual void fillHist(int type = 0);
  virtual void fillRedTreeData(int type = 0); 
  virtual void dumpGenBlock(bool withGluons=false); 

  virtual void muonEfficiency(); 

  virtual bool isAncestor(GenParticle *mo, GenParticle *dau);  
  virtual void printParticle(GenParticle *); 
  virtual void dumpDaughters(GenParticle *); 

  double       iso(Photon *, double radius, double ptmin); 

  GenParticle* getParticle(int i) {return (GenParticle*)fbParticles->At(i);}
  Photon*      getPhoton(int i) {return (Photon*)fbPhotons->At(i);}
  Electron*    getElectron(int i) {return (Electron*)fbElectrons->At(i);}
  Muon*        getMuon(int i) {return (Muon*)fbMuons->At(i);}
  HepMCEvent*  getEvent(int i) {return (HepMCEvent*)fbEvent->At(i);}

  
  ExRootTreeReader *fTR;
  TClonesArray     
  *fbParticles, 
    *fbEvent,
    *fbJets,
    *fbPhotons,  
    *fbTracks,
    *fbPFtowers, 
    *fbPFtracks, 
    *fbElectrons,
    *fbMuons;


  TFile       *fpHistFile; 
  TTree       *fTree;
  redTreeData fRtd;

  TChain      *fpChain; 
  int         fNentries;

  GenParticle *fH0; // NLO Higgs
  GenParticle *fH1; // `final' SMC Higgs
  GenParticle *fG0, *fG1; // photons

  TLorentzVector fp4genH0, fp4genH1, fp4genG0, fp4genG1,
    fp4H, fp4G0, fp4G1; 
  
  int         fClass; // 0 gamma gamma, 1 (mu mu)  (mu mu), 2 (e e) (e e), 3 (e e) (mu mu), 4 l l 
  double      fW8;
  int         fNRecoPhotons;

  double      fG0Iso, fG1Iso, fgenG0Iso, fgenG1Iso;

  static const int PTN   = 40; 
  static const int PTMAX = 800; 
};

// ----------------------------------------------------------------------
inline void mk4Vector(TLorentzVector &p4, const Double_t p, const Double_t t, const Double_t f, const Double_t m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

#endif
