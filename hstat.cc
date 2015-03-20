#include "hstat.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TFitResult.h"

#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooDataSet.h"

#include "RooAbsPdf.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooMCStudy.h"
#include "RooMsgService.h"
#include "RooConstVar.h"
#include "RooMinuit.h"
#include "RooProfileLL.h"


#include "RooDLLSignificanceMCSModule.h"

#include "RooGlobalFunc.h"


#include "RooStats/ModelConfig.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestPlot.h"

#include "dataset.hh"
#include "util.hh"

#include "RooOneSidedProfileLL.hh"

ClassImp(hstat)

using namespace std; 
using namespace RooFit; 
using namespace RooStats; 

// ----------------------------------------------------------------------
hstat::hstat(string dir,  string files, string setup)  {

  fNtoy = 1000; 
  fSetup = setup;

  if (setup == "") {
    fHistFileName = Form("%s/hstat.root", dir.c_str()); 
  } else {
    fHistFileName = Form("%s/hstat-%s.root", dir.c_str(), setup.c_str()); 
  }

  fTexFileName = fHistFileName; 
  replaceAll(fTexFileName, ".root", ".tex"); 
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));


  double MGGLO = 70.;
  double MGGHI = 180.;

  fSg0 =   60; 
  fSg1 =  115;
  fBg  = 1360; 
  fHiggsMpeak = 125.0; 
  fHiggsMres  =   9.4;
  // pol1
  fBgMp1   = 0.071276;
  fBgMp1E  = 0.00529807;
  // chebychev
  fBgMc0   = 0.159472;
  fBgMc0E  = 0.0119849;
  fBgTau   = -0.0135662;
  fBgTauE  =  0.000110543;
  fSg0Tau  = -0.0117118;
  fSg0TauE =  0.000176586;
  fSg1Tau  = -0.00791874;
  fSg1TauE = -9.58177e-05;

  bgData = sg0Data = sg1Data = data1 = data0 = 0; 

  fRm = new RooRealVar("m", "m", MGGLO, MGGHI, "GeV"); 
  fRsgP = new RooRealVar("sgP", "signal peak mass", 125., MGGLO, MGGHI);  
  fRsgP->setConstant(kTRUE);
  fRsgS = new RooRealVar("sgS", "signal sigma mass", fHiggsMres, 0., 15.);  
  fRsgS->setConstant(kTRUE);

  fRpt = new RooRealVar("pt", "pt", 200., 1000., "GeV"); 
    
}


// ----------------------------------------------------------------------
hstat::~hstat() {
  if (fHistFile) fHistFile->Close();
}


// ----------------------------------------------------------------------
void hstat::genData() {


  
}


// ----------------------------------------------------------------------
void hstat::sigNoBackground() {

}




// ----------------------------------------------------------------------
RooAbsPdf* hstat::genModel(int i, int ini) {
  /*
  double nsg(0.), ntau(0.); 
  nsg = fSg0; 
  ntau = fSg0Tau;
  if (1 == ini) {
    nsg = fSg1; 
    ntau = fSg1Tau;
  }     
  RooRealVar  *sgN = new RooRealVar(Form("sgN%d", i), Form("Number of Higgs signal events %d", i), nsg, 0., 1.e4);
  RooGaussian *sgM = new RooRealVar(Form("sgM%d", i), Form("signal mass %d", i), *fRm, *fRsgP, *fRsgS); 
  RooRealVar  *bgS = new RooRealVar(Form("bgS%d", i), Form("bg slope %d", i), fBgMc0, -10., 10.); 
  RooChebychev *bgM = new RooChebychev(Form("bgM%d", i), Form("background gamma gamma mass %d", i), 
				       *fRm, RooArgList(*bgS)); 
  RooRealVar  *bgN  = new RooRealVar(Form("bgN%d", i), Form("Number of background events %d", i), fBg, 0., 1.e5);
  RooAddPdf   *model = new RooAddPdf(Form("model%d"), "model2", RooArgList(sgM, bgM), RooArgList(*fRsgN, *fRbgN));
  */

  return 0; 
}
