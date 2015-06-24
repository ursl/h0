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
#include "TGraphAsymmErrors.h"

#include "RooPlot.h"
#include "RooDataSet.h"

#include "RooAbsPdf.h"
#include "RooRandom.h"
#include "RooFitResult.h"
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
// par[0] -> const
// par[1] -> mean
// par[2] -> sigmaL
// par[3] -> sigmaR
double gaussBifurcated(double *x, double *par) {
  
  double arg = x[0] - par[1];
  double sigmaL = par[2];
  double sigmaR = par[3];

  double coef(0.0);

  if (arg < 0.0){
    if (TMath::Abs(sigmaL) > 1e-30) {
      coef = -0.5/(sigmaL*sigmaL);
    }
  } else {
    if (TMath::Abs(sigmaR) > 1e-30) {
      coef = -0.5/(sigmaR*sigmaR);
    }
  }
    
  return par[0]*exp(coef*arg*arg);
}



// ----------------------------------------------------------------------
hstat::hstat(double lumi, string dir,  string files, string setup)  {

  fLumi = lumi; 
  fDoPlotName = Form("2dfit-lumi-%2.1f-setup-%s", fLumi, setup.c_str()); 
  fNtoy = 1000; 
  fSetup = setup;

  fRndmSeed = 111; 

  fHistFile = 0; 

  if (setup == "") {
    fHistFileName = Form("%s/hstat.root", dir.c_str()); 
  } else {
    fHistFileName = Form("%s/hstat-%s.root", dir.c_str(), setup.c_str()); 
  }

  fTexFileName = fHistFileName; 
  replaceAll(fTexFileName, ".root", ".tex"); 
  //  system(Form("/bin/rm -f %s", fTexFileName.c_str()));

  fDirectory = "hpt0";

  NBINS = 55; 
  MGGLO = 70.;
  MGGHI = 180.;

  // -- const unmutable setup for 1/ab
  fCONSTSG0 =   86; 
  fCONSTSG1 =  150;
  fCONSTBG  = 1363; 

#ifdef EXPO
  fCONSTSG0TAU = -0.0117118; 
  fCONSTSG1TAU = -0.00791874;
  fCONSTBGTAU  = -0.0135662;
#else
  // this is "k"
  fCONSTSG0TAU = 2.20; 
  fCONSTSG1TAU = 2.11;
  fCONSTBGTAU  = 2.04;

  // this is "mu"
  fCONSTSG0MU = 40.;
  fCONSTSG1MU = 95.;
  fCONSTBGMU  = 39.;
#endif

  // -- hypothesis setup
  fSG0 =  fLumi*fCONSTSG0; 
  fSG1 =  fLumi*fCONSTSG1;
  fBG  =  fLumi*fCONSTBG; 

#ifdef EXPO
  fSG0TAU = -0.0117118; 
  fSG1TAU = -0.00791874;
  fBGTAU  = -0.0135662;
#else
  // this is "k"
  fSG0TAU = 2.20; 
  fSG1TAU = 2.11;
  fBGTAU  = 2.04;

  // this is "mu"
  fSG0MU = 40.;
  fSG1MU = 95.;
  fBGMU  = 39.;
#endif

  fSg0 = fSG0; 
  fSg1 = fSG1;
  fBg  = fBG; 

  fHiggsMpeak = 125.0; 
  fHiggsMres  =   9.4;
  // mass pol1
  fBgMp1   = 0.071276;
  fBgMp1E  = 0.00529807;
  // mass chebychev
  fBgMc0   = 0.159472;
  fBgMc0E  = 0.0119849;

#ifdef EXPO
  fBgTau   = -0.0135662;
  fBgTauE  =  0.000110543;
  fSg0Tau  = -0.0117118;
  fSg0TauE =  0.000176586;
  fSg1Tau  = -0.00791874;
  fSg1TauE = -9.58177e-05;
#else
  fBgTau   = 2.04;
  fBgTauE  = 0.20;
  fBgMu    = 39.0;
  fBgMuE   =  4.0;

  fSg0Tau  = 2.20; 
  fSg0TauE = 0.02;
  fSg0Mu   = 40.0;
  fSg0MuE  =  4.0;

  fSg1Tau  = 2.11; 
  fSg1TauE = 0.02;
  fSg1Mu   = 95.0;
  fSg1MuE  =  9.0;
#endif  
}


// ----------------------------------------------------------------------
hstat::~hstat() {
  if (fHistFile) fHistFile->Close();
}


// ----------------------------------------------------------------------
void hstat::run1() {

  model *m0 = genModel1(0, fSg0, fBg); 
  model *m1 = genModel1(1, 90, fBg); 
  model *m2 = genModel1(1, fSg0, fBg); 

  RooDataSet *d0(0), *bgData(0), *sg0Data(0);

  bgData  = m0->bgPdf->generate(RooArgSet(*m0->m), fBg); 
  sg0Data = m0->sgPdf->generate(RooArgSet(*m0->m), fSg0); 

  // -- model 0 plus background
  d0 = new RooDataSet(*bgData);
  d0->append(*sg0Data);

  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "Input: sg0N = " << fSg0 << endl;
  cout << "Input: sg0N'= " << 90 << endl;
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  
  cout << "XXXXXXXXX create Nll0" << endl;
  RooAbsReal *m0NLL = m0->modelPdf->createNLL(*d0, /*ExternalConstraints(bgSlopeC0),*/ Extended(kTRUE), Verbose(kFALSE));
  RooAbsReal *m2NLL = m2->modelPdf->createNLL(*d0, /*ExternalConstraints(bgSlopeC0),*/ Extended(kTRUE), Verbose(kFALSE));
  cout << endl << "========= sg = " << m0->sgN->getVal() << " +/- " << m0->sgN->getError() << endl;
  cout << endl << "XXXXXXXXX calling m0NLL->getVal() " << endl;
  double m0nll = m0NLL->getVal(); 
  cout << endl << "XXXXXXXXXXXXXXXXXXXX m0nll = " << m0nll << endl;
  cout << endl << "========= sg = " << m0->sgN->getVal() << " +/- " << m0->sgN->getError() << endl;
  cout << "XXXXXXXXX create profile from m0NLL" << endl;
  RooProfileLL *pm0NLL = dynamic_cast<RooProfileLL*>(m0NLL->createProfile(m0->poi));
  cout << endl << "========= sg = " << m0->sgN->getVal() << " +/- " << m0->sgN->getError() << endl;
  cout << endl << "XXXXXXXXX calling pm0NLL0->getVal() = " << endl;
  double pm0nll      = pm0NLL->getVal(); 
  cout << endl << "XXXXXXXXX pm0nll = " << pm0nll << endl;
  cout << endl << "========= sg = " << m0->sgN->getVal() << " +/- " << m0->sgN->getError() << endl;
  cout << endl << "XXXXXXXXXXXXXXXXXXXX Nll0->getVal() = " << m0NLL->getVal() << " after profile setup" << endl;

  cout << "XXXXXXXXX RooMinuit(t0)" << endl;
  RooMinuit t0(*m0NLL);
  t0.setPrintLevel(-1);
  t0.setNoWarn();
  cout << "XXXXXXXXX t0.migrad" << endl;
  t0.migrad();
  cout << endl << "XXXXXXXXXXXXXXXXXXXX Nll0->getVal() = " << m0NLL->getVal() << " after migrad" << endl;
  cout << endl << "========= sg = " << m0->sgN->getVal() << " +/- " << m0->sgN->getError() << endl;

  
  cout << "XXXXXXXXX t0.minos" << endl;
  t0.minos();
  cout << endl << "XXXXXXXXXXXXXXXXXXXX Nll0->getVal() = " << m0NLL->getVal() << " after minos" << endl;
  cout << endl << "========= sg = " << m0->sgN->getVal() << " +/- " << m0->sgN->getError() << endl;

  cout << "XXXXXXXXX create another profile from m0NLL that minimizes that" << endl;
  RooProfileLL *pm2NLL = dynamic_cast<RooProfileLL*>(m2NLL->createProfile(m2->poi));
  cout << endl << "========= sg = " << m2->sgN->getVal() << " +/- " << m2->sgN->getError() << endl;
  cout << endl << "XXXXXXXXX calling pm2NLL1->getVal() = " << endl;
  double pm2nll      = pm2NLL->getVal(); 
  cout << endl << "XXXXXXXXX pm2nll = " << pm2nll << endl;
  cout << endl << "========= sg = " << m2->sgN->getVal() << " +/- " << m2->sgN->getError() << endl;
  cout << endl << "XXXXXXXXXXXXXXXXXXXX m2Nll->getVal() = " << m2NLL->getVal() << " after another profile setup" << endl;

  cout << "XXXXXXXXX create another profile from M1" << endl;
  cout << endl << "========= sg = " << m1->sgN->getVal() << " +/- " << m1->sgN->getError() << endl;
  RooAbsReal *m1NLL = m1->modelPdf->createNLL(*d0,  Extended(kTRUE), Verbose(kFALSE));
  cout << endl << "XXXXXXXXX m1nll = " << m1NLL->getVal() << endl;
  cout << endl << "========= sg = " << m1->sgN->getVal() << " +/- " << m1->sgN->getError() << endl;
  RooProfileLL *pm1NLL = dynamic_cast<RooProfileLL*>(m1NLL->createProfile(m1->poi));
  cout << endl << "========= sg = " << m1->sgN->getVal() << " +/- " << m1->sgN->getError() << endl;
  cout << endl << "XXXXXXXXX calling pm1NLL->getVal() = " << endl;
  double pm1nll    = pm1NLL->getVal(); 
  cout << endl << "XXXXXXXXX pm1nll = " << pm1nll << endl;
  cout << endl << "========= sg = " << m1->sgN->getVal() << " +/- " << m1->sgN->getError() << endl;

  RooMinuit t1(*m1NLL);
  t1.setPrintLevel(-1);
  t1.setNoWarn();
  t1.migrad();



  RooPlot *frame2 = m0->sgN->frame(Name("plot2"), Range(86., 110.));
  cout << endl << "========= sg = " << m0->sgN->getVal() << " +/- " << m0->sgN->getError() << endl;
  m0NLL->plotOn(frame2, LineColor(kBlack), ShiftToZero());
  cout << endl << "========= sg = " << m0->sgN->getVal() << " +/- " << m0->sgN->getError() << endl;
  pm0NLL->plotOn(frame2, LineColor(kBlue), LineStyle(kSolid));
  frame2->SetMinimum(0) ;
  frame2->SetMaximum(0.2) ;
  frame2->Draw();


  RooPlot *frame3 = m1->sgN->frame(Name("plot3"), Range(86., 110.));
  m1NLL->plotOn(frame3, LineColor(kRed), LineStyle(kDashed), ShiftToZero());
  pm1NLL->plotOn(frame3, LineColor(kMagenta), LineStyle(kDashed));
  frame3->SetMinimum(0) ;
  frame3->SetMaximum(0.2) ;
  frame3->Draw("same");


  return;


//   RooPlot *frame2 = m0->sgN->frame(Name("plot2"), Range(30., 200.));
//   m0NLL->plotOn(frame2, LineColor(kBlack));
//   pm0NLL->plotOn(frame2, LineColor(kBlue), LineStyle(kDashed));
//   m1NLL->plotOn(frame2, LineColor(kRed), LineStyle(kDashDotted));
//   pm1NLL->plotOn(frame2, LineColor(kMagenta), LineStyle(kDotted));
//   frame2->Draw();


}


// ----------------------------------------------------------------------
void hstat::systematics(int mode, int n) {

  fTEX.open(fTexFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;

  int ntoys(n); 

  fDoPlotName = Form("2dfit-syst-lumi-%2.1f-setup-%s-mode-%d-n-%d", fLumi, fSetup.c_str(), mode, n); 

  // -- no systematics 2D
  if (10 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, 0); 
    double sRaw = fSeparation; 
    fTEX << "2d sNoSyst: " << Form("%4.3f +/- %4.3f", sRaw, fSeparationE) << endl;
    cout << "2d sNoSyst: " << Form("%4.3f +/- %4.3f", sRaw, fSeparationE) << endl;
  }

  if (11 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, 101); 
    double sTopLo = fSeparation; 
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, -101); 
    double sTopHi = fSeparation; 
    fTEX << "2d sTop: " << Form("%4.3f .. %4.3f", sTopLo, sTopHi) << endl;
    cout << "2d sTop: " << Form("%4.3f .. %4.3f", sTopLo, sTopHi) << endl;
  }

  if (12 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, 102); 
    double sBgLo = fSeparation; 
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, -102); 
    double sBgHi = fSeparation; 
    fTEX << "2d sBg: " << Form("%4.3f .. %4.3f", sBgLo, sBgHi) << endl;
    cout << "2d sBg: " << Form("%4.3f .. %4.3f", sBgLo, sBgHi) << endl;
  }

  // -- scale
  if (13 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, 103); 
    double sScaleLo = fSeparation; 
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, -103); 
    double sScaleHi = fSeparation; 
    fTEX << "2d sScale: " << Form("%4.3f .. %4.3f", sScaleLo, sScaleHi) << endl;
    cout << "2d sScale: " << Form("%4.3f .. %4.3f", sScaleLo, sScaleHi) << endl;
  }

  // -- signal efficiency
  if (14 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, 104); 
    double sEffLo = fSeparation; 
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run2D(ntoys, -104); 
    double sEffHi = fSeparation; 
    fTEX << "2d sEff: " << Form("%4.3f .. %4.3f", sEffLo, sEffHi) << endl;
    cout << "2d sEff: " << Form("%4.3f .. %4.3f", sEffLo, sEffHi) << endl;
  }



  // -- no systematics
  if (0 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run1D(ntoys, 0); 
    double sRaw = fSeparation; 
    fTEX << "1d sNoSyst: " << Form("%4.3f +/- %4.3f", sRaw, fSeparationE) << endl;
    cout << "1d sNoSyst: " << Form("%4.3f +/- %4.3f", sRaw, fSeparationE) << endl;
  }

  // -- missing top systematics
  if (1 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run1D(ntoys, 101); 
    double sTopLo = fSeparation; 
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run1D(ntoys, -101); 
    double sTopHi = fSeparation; 
    fTEX << "1d sTop: " << Form("%4.3f .. %4.3f", sTopLo, sTopHi) << endl;
    cout << "1d sTop: " << Form("%4.3f .. %4.3f", sTopLo, sTopHi) << endl;
  }

  // -- background systematics
  if (2 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run1D(ntoys, 102); 
    double sBgLo = fSeparation; 
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run1D(ntoys, -102); 
    double sBgHi = fSeparation; 
    fTEX << "1d sBg: " << Form("%4.3f .. %4.3f", sBgLo, sBgHi) << endl;
    cout << "1d sBg: " << Form("%4.3f .. %4.3f", sBgLo, sBgHi) << endl;
  }

  // -- scale
  if (3 == mode) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run1D(ntoys, 103); 
    double sScaleLo = fSeparation; 
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    run1D(ntoys, -103); 
    double sScaleHi = fSeparation; 
    fTEX << "1d sScale: " << Form("%4.3f .. %4.3f", sScaleLo, sScaleHi) << endl;
    cout << "1d sScale: " << Form("%4.3f .. %4.3f", sScaleLo, sScaleHi) << endl;
  }


  // -- pdf eigenvectors
  if (4 == mode) {
    run1D(ntoys, 104); 
    double sPdfLo = fSeparation; 
    run1D(ntoys, -104); 
    double sPdfHi = fSeparation; 
    fTEX << "1d sPdf: " << Form("%4.3f .. %4.3f", sPdfLo, sPdfHi) << endl;
    cout << "1d sPdf: " << Form("%4.3f .. %4.3f", sPdfLo, sPdfHi) << endl;
  }

  // -- all 
  if (5 == mode) {
    run1D(ntoys, 20); 
    double sSyst = fSeparation; 
    fTEX << "1d sSyst: " << Form("%4.3f +/- %4.3f", sSyst, fSeparationE) << endl;
    cout << "1d sSyst: " << Form("%4.3f +/- %4.3f", sSyst, fSeparationE) << endl;
  }

  fTEX.close(); 
}


// ----------------------------------------------------------------------
void hstat::sigStudies(int mode, int n) {

  fDoPlotName = Form("2dfit-sigStudies-lumi-%2.1f-setup-%s-mode-%d-n-%d", fLumi, fSetup.c_str(), mode, n); 

  fTEX.open(fTexFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  
  TH1D *hsig(0); 
  if (10 == mode) {
    // -- significance vs background numbers
    int nBINS(9); 
    hsig = new TH1D("hsig", "nBG = 0 .. 2000", 40, 0., 2000.); 
    double bg[] = {fBG, 200., 500., 750., 1000., 1250., 1500., 1750., 1950.};
    hsig->GetXaxis()->SetTitle("N_{BG}");
      for (int i = 0; i < nBINS; ++i) {
	fBg = bg[i];
	fBG = bg[i];
	fSg0 = fSG0;
	fSg1 = fSG1; 
	cout << "XXXXXXXXXX run " << i << "  " << fBg << " .. " << fSg0 << " .. " << fSg1 << endl;
	run2D(n, -1); 
	if (fSeparation < 100.) {
	  hsig->SetBinContent(hsig->FindBin(fBg), fSeparation); 
	  hsig->SetBinError(hsig->FindBin(fBg), fSeparationE); 
	}
	fTEX << "sigStudies::nbg : " << fBg << " "  << fSg0 << " " << fSg1 << " " 
	     << fSeparation << " +/- " << fSeparationE << endl;
      }      
  } else if (11 == mode) {
    // -- significance vs lumi
    int nBINS(6); 
    hsig = new TH1D("hsig", "lumi = 500 .. 3000", 30, 500., 3500.); 
    hsig->GetXaxis()->SetTitle("luminosity [/fb]");
    double scale[] = {1.0, 0.5, 1.5, 2.0, 2.5, 3.0};
    //    double scale[] = {1.0, 2.0, 3.0};
    for (int i = 0; i < nBINS; ++i) {
      fBg =  scale[i]*fCONSTBG; 
      fBG =  scale[i]*fCONSTBG; 
      fSg0 = scale[i]*fCONSTSG0;
      fSG0 = scale[i]*fCONSTSG0;
      fSg1 = scale[i]*fCONSTSG1; 
      fSG1 = scale[i]*fCONSTSG1; 
      cout << "XXXXXXXXXX run " << i << "  " << fBg << " .. " << fSg0 << " .. " << fSg1 << endl;
      run2D(n, -1); 
      if (fSeparation < 100.) {
	hsig->SetBinContent(hsig->FindBin(scale[i]*1000.), fSeparation); 
	hsig->SetBinError(hsig->FindBin(scale[i]*1000.), fSeparationE); 
      }
      fTEX << "sigStudies::lumi : " << scale[i] << " " << fBg << " "  << fSg0 << " " << fSg1 
	   << " " << fSeparation << " +/- " << fSeparationE << endl;
    }      
  } else  if (12 == mode) {
    // -- significance error from repeated trials
    int nBINS(20); 
    hsig = new TH1D("hsig", "run = 0 .. 9 ", 400, 1., 5.); 
    hsig->GetXaxis()->SetTitle("significance");
    for (int i = 0; i < nBINS; ++i) {
      fBg =  fCONSTBG; 
      fBG =  fCONSTBG; 
      fSg0 = fCONSTSG0;
      fSG0 = fCONSTSG0;
      fSg1 = fCONSTSG1; 
      fSG1 = fCONSTSG1; 
      cout << "XXXXXXXXXX run " << i << "  " << fBg << " .. " << fSg0 << " .. " << fSg1 << endl;

      run2D(n, -1); 
      if (fSeparation < 100.) {
	hsig->Fill(fSeparation); 
      }
      fTEX << "sigStudies::repeat : " << i << " " << fBg << " "  << fSg0 << " " << fSg1 
	   << " " << fSeparation << " +/- " << fSeparationE << endl;
    }      
  } else  if (13 == mode) {
    // -- significance vs signal error
    int nBINS(4); 
    TH1D *hLo = new TH1D("hLo", "error", 25, 1., 25.); 
    TH1D *hHi = new TH1D("hHi", "error", 25, 1., 25.); 
    double error[] = {3., 5., 10., 20.};
    for (int i = 0; i < nBINS; ++i) {
      fBg  = fCONSTBG; 
      fBG  = fCONSTBG; 
      fSg0 = fCONSTSG0 * (1. + 0.01*error[i]);
      fSG0 = fCONSTSG0 * (1. + 0.01*error[i]);
      fSg1 = fCONSTSG1 * (1. + 0.01*error[i]);
      fSG1 = fCONSTSG1 * (1. + 0.01*error[i]);
      cout << "XXXXXXXXXX run " << i << "  " << fBg << " .. " << fSg0 << " .. " << fSg1 << endl;
      run2D(n, -1); 
      if (fSeparation < 100.) {
	hHi->SetBinContent(hHi->FindBin(error[i]), fSeparation); 
      }

      fTEX << "sigStudies::sigError : +" << error[i] << " " << fBg << " "  << fSg0 << " " << fSg1 
	   << " " << fSeparation << " +/- " << fSeparationE << endl;

      fBg  = fBG; 
      fBG  = fBG; 
      fSg0 = fSG0 * (1. - 0.01*error[i]);
      fSG0 = fSG0 * (1. - 0.01*error[i]);
      fSg1 = fSG1 * (1. - 0.01*error[i]);
      fSG1 = fSG1 * (1. - 0.01*error[i]);
      cout << "XXXXXXXXXX run " << i << "  " << fBg << " .. " << fSg0 << " .. " << fSg1 << endl;
      run2D(n, -1); 
      if (fSeparation < 100.) {
	hLo->SetBinContent(hLo->FindBin(error[i]), fSeparation); 
      }
      fTEX << "sigStudies::sigError : -" << error[i] << " " << fBg << " "  << fSg0 << " " << fSg1 
	   << " " << fSeparation << " +/- " << fSeparationE << endl;
    }      
    
    TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0"); 
    if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);
    c0->Clear();
    gPad->SetLogy(0); 
    hHi->SetMarkerStyle(20); 
    hHi->Draw("p][");
    hLo->SetMarkerStyle(21); 
    hLo->Draw("p][same");
    
    c0->SaveAs(Form("hstat-sigStudies-%d-%d.pdf", mode, n)); 
    fTEX.close();
    
    fHistFile = TFile::Open(Form("hstat-sigStudies-%d-%d.root", mode, n), "RECREATE"); 
    hHi->SetDirectory(fHistFile); 
    hHi->Write();
    hLo->SetDirectory(fHistFile); 
    hLo->Write();
    fHistFile->Close();
    fHistFile = 0; 
    
    return;
    
  } else if (0 == mode) {
    // -- significance vs background numbers
    int nBINS(9); 
    hsig = new TH1D("hsig", "nBG = 0 .. 2000", 40, 0., 2000.); 
    double bg[] = {fBG, 200., 500., 750., 1000., 1250., 1500., 1750., 1950.};
    hsig->GetXaxis()->SetTitle("N_{BG}");
    for (int i = 0; i < nBINS; ++i) {
      fBg = bg[i];
      fSg0 = fSG0;
      fSg1 = fSG1; 
      cout << "XXXXXXXXXX run " << i << "  " << fBg << " .. " << fSg0 << " .. " << fSg1 << endl;
      run1D(n, -1); 
      if (fSeparation < 100.) {
	hsig->SetBinContent(hsig->FindBin(fBg), fSeparation); 
	hsig->SetBinError(hsig->FindBin(fBg), fSeparationE); 
      }
      fTEX << "sigStudies::nbg : " << fBg << " "  << fSg0 << " " << fSg1 << " " 
	   << fSeparation << " +/- " << fSeparationE << endl;
    }
  } else if (1 == mode) {
    // -- significance vs lumi
    int nBINS(6); 
    hsig = new TH1D("hsig", "lumi = 500 .. 3000", 30, 500., 3500.); 
    hsig->GetXaxis()->SetTitle("luminosity [/fb]");
    double scale[] = {1.0, 0.5, 1.5, 2.0, 2.5, 3.0};
    for (int i = 0; i < nBINS; ++i) {
      fBg =  scale[i]*fBG; 
      fSg0 = scale[i]*fSG0;
      fSg1 = scale[i]*fSG1; 
      cout << "XXXXXXXXXX run " << i << "  " << fBg << " .. " << fSg0 << " .. " << fSg1 << endl;
      run1D(n, -1); 
      if (fSeparation < 100.) {
	hsig->SetBinContent(hsig->FindBin(scale[i]*1000.), fSeparation); 
	hsig->SetBinError(hsig->FindBin(scale[i]*1000.), fSeparationE); 
      }
      fTEX << "sigStudies::lumi : " << scale[i] << " " << fBg << " "  << fSg0 << " " << fSg1 
	   << " " << fSeparation << " +/- " << fSeparationE << endl;
    }      
  } else if (2 == mode) {
    // -- significance error from repeated trials
    int nBINS(20); 
    hsig = new TH1D("hsig", "run = 0 .. 9 ", 400, 1., 5.); 
    hsig->GetXaxis()->SetTitle("significance");
    for (int i = 0; i < nBINS; ++i) {
      fBg =  fBG; 
      fSg0 = fSG0;
      fSg1 = fSG1; 
      cout << "XXXXXXXXXX run " << i << "  " << fBg << " .. " << fSg0 << " .. " << fSg1 << endl;

      run1D(n, -1); 
      if (fSeparation < 100.) {
	hsig->Fill(fSeparation); 
      }
      fTEX << "sigStudies::repeat : " << i << " " << fBg << " "  << fSg0 << " " << fSg1 
	   << " " << fSeparation << " +/- " << fSeparationE << endl;
    }      
  }
    
  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);
  c0->Clear();
  hsig->SetMarkerStyle(20); 
  hsig->Draw("p][");
  if (12 == mode || 2 == mode)   hsig->Draw("hist");
  c0->SaveAs(Form("hstat-sigStudies-%d-%d.pdf", mode, n)); 
  fTEX.close();

  fHistFile = TFile::Open(Form("hstat-sigStudies-%d-%d.root", mode, n), "RECREATE"); 
  hsig->SetDirectory(fHistFile); 
  hsig->Write();
  fHistFile->Close();
  fHistFile = 0; 

}

// ----------------------------------------------------------------------
void hstat::run1D(int ntoys, int mode) {

  TH1D *ht0 = new TH1D("ht0", "prof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 1000, -40., 40.); 
  ht0->SetLineColor(kBlue); 
  TH1D *ht1 = new TH1D("ht1", "prof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 1000, -40., 40.); 
  ht1->SetLineColor(kRed); 

  TH1D *hsg0 = (TH1D*)gROOT->Get("hsg0"); 
  if (0 == hsg0) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    hsg0 = new TH1D("hsg0", "hsg0", 50, 0., 200.); 
  } else { 
    hsg0->Reset();
  }
  TH1D *hsg1 = (TH1D*)gROOT->Get("hsg1"); 
  if (0 == hsg1) {
    hsg1 = new TH1D("hsg1", "hsg1", 50, 0., 400.); 
  } else {
    hsg1->Reset();
  }
  TH1D *hrat = (TH1D*)gROOT->Get("hrat");
  if (0 == hrat) {
    hrat = new TH1D("hrat", "hrat", 50, 0.5, 3.0); 
  } else {
    hrat->Reset(); 
  }

  if ((0 == mode || -1 == mode)) {
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "No systematics: " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      cout << "XXXXXX run " << i << endl;
      toy1D();
    
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (101 == TMath::Abs(mode)) {
    double syst(0.17); 
    if (mode > 0) {
      syst = 1. + syst;
    } else {
      syst = 1. - syst;
    }
    fSg0 = fSG0*syst; 
    fSg1 = fSG1;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from missing top calculations, syst = " << syst 
	 << ", " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (102 == TMath::Abs(mode)) {
    double syst(0.40); 
    if (mode > 0) {
      syst = 1. + syst;
    } else {
      syst = 1. - syst;
    }
    fSg0 = fSG0; 
    fSg1 = fSG1;
    fBg  = fBG*syst;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from background uncertainties, syst = " << syst 
	 << ", " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (103 == TMath::Abs(mode)) {
    double syst(0.20); 
    if (mode > 0) {
      syst = 1. + syst;
    } else {
      syst = 1. - syst;
    }
    fSg0 = fSG0*syst; 
    fSg1 = fSG1*syst;
    fBg  = fBG;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from scale uncertainties, syst = " << syst 
	 << ", " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (104 == TMath::Abs(mode)) {
    double syst(0.03); 
    if (mode > 0) {
      syst = 1. + syst;
    } else {
      syst = 1. - syst;
    }
    fSg0 = fSG0*syst; 
    fSg1 = fSG1/syst;
    fBg  = fBG;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from pdf/eigenvector uncertainties, syst = " << syst << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (1 == TMath::Abs(mode)) {
    double syst(0.17); 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from missing top calculations" << endl;
    cout << "Input: sg0N = " << fSg0 << ", varied by " << 100*syst << "%" << endl;
    cout << "Input: sg1N = " << fSg1 << ", varied through delta(fSg1/fSg0) = " << 100*syst << "%" <<  endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    double randFact(1.0); 
    for (int i = 0; i < ntoys; ++i) {
      randFact = 1. + gRandom->Gaus(0., syst);
      while (randFact < 0.) {
	randFact = 1. + gRandom->Gaus(0., syst);
      }
      fSg0 = fSG0*randFact; 
      fSg1 = fSG1;
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (2 == TMath::Abs(mode)) {
    double systHi(0.22); 
    double systLo(0.06); 
    TF1 *f1 = new TF1("f1BifGauss", gaussBifurcated, -2., 2., 4); 
    f1->SetParameters(1., 0., systLo, systHi); 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from pdf eigenvectors" << endl;
    cout << "Input: sg0N = " << fSg0 << ", varied by +" << systHi << " -" << systLo<< endl;
    cout << "Input: sg1N = " << fSg1 << ", varied through delta(fSg1/fSg0) = +" << 1.-systHi << " -" << 1.-systLo << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    double randFact(1.0);
    for (int i = 0; i < ntoys; ++i) {
      randFact = 1. + f1->GetRandom(); 
      while (randFact < 0.) {
	randFact = 1. + f1->GetRandom(); 
      }
      fSg0 = fSG0*randFact; 
      fSg1 = fSG1;
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (3 == TMath::Abs(mode)) {
    double syst(0.20); 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from scale uncertainties" << endl;
    cout << "Input: sg0N = " << fSg0 << ", varied by " << 100*syst << "%" << endl;
    cout << "Input: sg1N = " << fSg1 << ", varied through delta(fSg1/fSg0) = " << 100*syst << "%" <<  endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    double randFact(1.0); 
    for (int i = 0; i < ntoys; ++i) {
      randFact = 1. + gRandom->Gaus(0., syst);
      while (randFact < 0.) {
	randFact = 1. + gRandom->Gaus(0., syst);
      }
      fSg0 = fSG0*randFact; 
      fSg1 = fSG1*randFact;
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (10 == TMath::Abs(mode)) {
    double syst(0.4); 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from background uncertainty" << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 <<  endl;
    cout << "Input: bgN  = " << fBg << ", varied by " << 100*syst << "%, truncated to be > 0" << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    double randFact(1.0); 
    for (int i = 0; i < ntoys; ++i) {
      fSg0 = fSG0; 
      fSg1 = fSG1;
      randFact = 1. + gRandom->Gaus(0., syst);
      while (randFact < 0.) {
	randFact = 1. + gRandom->Gaus(0., syst);
      }
      fBg = fBG*randFact; 
      if (fBg < 0) fBg = 10.; 
      cout << "XXXXXX run " << i << ", fBg = " << fBg << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (20 == TMath::Abs(mode)) {
    double syst(0.25); 
    double systBg(0.40); 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from complete uncertainty" << endl;
    cout << "Input: sg0N = " << fSg0 << ", varied by " << 100*syst << "%" << endl;
    cout << "Input: sg1N = " << fSg1 << ", varied by " << 100*syst << "%" << endl;
    cout << "Input: bgN  = " << fBg << ", varied by " << 100*systBg << "%, truncated to be > 0" << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    double randFact(1.0); 
    for (int i = 0; i < ntoys; ++i) {
      randFact = 1. + gRandom->Gaus(0., systBg);
      while (randFact < 0.) {
	randFact = 1. + gRandom->Gaus(0., systBg);
      }
      fBg = fBG*randFact; 

      randFact = 1. + gRandom->Gaus(0., syst);
      while (randFact < 0.) {
	randFact = 1. + gRandom->Gaus(0., syst);
      }
      fSg0 = fSG0*randFact; 
      randFact = 1. + gRandom->Gaus(0., syst);
      while (randFact < 0.) {
	randFact = 1. + gRandom->Gaus(0., syst);
      }
      fSg1 = fSG1*randFact; 
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << endl;

      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy1D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  }

  TLatex *tl = new TLatex(); 
  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.04);

  TArrow *ta = new TArrow(); 

  gStyle->SetOptStat(111111); 
  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);
  c0->Clear(); 
  c0->Divide(2,2);
  c0->cd(1);
  hsg0->Draw();
  tl->DrawLatex(0.15, 0.92, Form("rel rms: %4.3f", hsg0->GetRMS()/hsg0->GetMean())); 
  ta->DrawArrow(fSG0, 0.5*hsg0->GetMaximum(), fSG0, 0.); 
  c0->cd(2); 
  hsg1->Draw();
  tl->DrawLatex(0.15, 0.92, Form("rel rms: %4.3f", hsg1->GetRMS()/hsg1->GetMean())); 
  ta->DrawArrow(fSG1, 0.5*hsg1->GetMaximum(), fSG1, 0.); 
  c0->cd(3); 
  hrat->Draw();
  tl->DrawLatex(0.15, 0.92, Form("rel rms: %4.3f", hrat->GetRMS()/hrat->GetMean())); 
  ta->DrawArrow(fSG1/fSG0, 0.5*hrat->GetMaximum(), fSG1/fSG0, 0.); 

  c0->cd(4);

  double midpoint, tailprob, sigma, lo, hi;
  findMidPoint(ht0, ht1, midpoint, tailprob, lo, hi); 
  //  sigma = 2.*oneSidedGaussianSigma(tailprob);
  sigma = 2.*RooStats::PValueToSignificance(tailprob);
  fSeparation = sigma; 
  fSeparationE = 2.*(RooStats::PValueToSignificance(hi) - RooStats::PValueToSignificance(lo)); 
  if (fSeparationE < 1.e-3) fSeparationE = 1.e-3; 
  // -- new error on separation according to 
  // evernote:///view/40567474/s295/local/x-coredata://B740E98C-51E1-432D-B0F3-002D230CF640/ENNote/p7145/
  fSeparationE = 0.02;
  if (ntoys < 2200) fSeparationE = 0.02*fSeparation; 
  if (ntoys < 1200) fSeparationE = 0.02*fSeparation; 
  if (ntoys < 700) fSeparationE = 0.03*fSeparation; 
  if (ntoys < 200) fSeparationE = 0.08*fSeparation; 
  

  ht0->Rebin(10)->Draw();
  ht1->Rebin(10)->Draw("same");
  
  tl->DrawLatex(0.17, 0.85, Form("ntoys/mode: %d/%d", ntoys, mode)); 
  tl->DrawLatex(0.17, 0.80, Form("midpoint: %4.3f", midpoint)); 
  tl->DrawLatex(0.17, 0.75, Form("tailprob: %4.3f", tailprob)); 
  tl->DrawLatex(0.17, 0.70, Form("Sep: %4.3f +/- %4.3f", sigma, fSeparationE)); 
  
  cout << "----------------------------------------------------------------------" << endl;
  cout << "midpoint: " << midpoint << endl;
  cout << "sigma:    " << sigma << " +/- " << fSeparationE << endl;
  cout << "sigma(lo): " << 2.*RooStats::PValueToSignificance(lo) << endl;
  cout << "sigma(hi): " << 2.*RooStats::PValueToSignificance(hi) << endl;
  cout << "----------------------------------------------------------------------" << endl;
  
  if (mode != -1) c0->SaveAs(Form("hstat-%d-%d.pdf", mode, ntoys)); 
  
}


// ----------------------------------------------------------------------
void hstat::toy1D() {

  // -- the model is fixed to fSG0/fSG1/fBG, while the data below may be affected by systematic effects
  model *m0 = genModel1(0, fSG0, fBG); 
  model *m1 = genModel1(1, fSG1, fBG); 

  RooDataSet *d0(0), *d1(0), *bgData(0), *sg0Data(0), *sg1Data(0);

  bgData  = m0->bgPdf->generate(RooArgSet(*m0->m), fBg); 
  sg0Data = m0->sgPdf->generate(RooArgSet(*m0->m), fSg0); 
  sg1Data = m1->sgPdf->generate(RooArgSet(*m1->m), fSg1); 

  // -- signal model plus background
  d0 = new RooDataSet(*bgData);
  d0->append(*sg0Data);
  d1 = new RooDataSet(*bgData);
  d1->append(*sg1Data);

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooAbsReal *d0m0NLL = m0->modelPdf->createNLL(*d0, Extended(kTRUE), Verbose(kFALSE));
  RooAbsReal *d0m1NLL = m1->modelPdf->createNLL(*d0, Extended(kTRUE), Verbose(kFALSE));
  RooAbsReal *d1m0NLL = m0->modelPdf->createNLL(*d1,  Extended(kTRUE), Verbose(kFALSE));
  RooAbsReal *d1m1NLL = m1->modelPdf->createNLL(*d1,  Extended(kTRUE), Verbose(kFALSE));
  
  RooProfileLL *pd0m0NLL = dynamic_cast<RooProfileLL*>(d0m0NLL->createProfile(m0->poi));
  fD0M0 = pd0m0NLL->getVal(); 
  RooProfileLL *pd0m1NLL = dynamic_cast<RooProfileLL*>(d0m1NLL->createProfile(m1->poi));
  fD0M1 = pd0m1NLL->getVal(); 
  RooProfileLL *pd1m0NLL = dynamic_cast<RooProfileLL*>(d1m0NLL->createProfile(m0->poi));
  fD1M0 = pd1m0NLL->getVal(); 
  RooProfileLL *pd1m1NLL = dynamic_cast<RooProfileLL*>(d1m1NLL->createProfile(m1->poi));
  fD1M1 = pd1m1NLL->getVal();   

  delModel(m0); 
  delModel(m1); 

  delete d0;
  delete d1;
  delete bgData;
  delete sg0Data; 
  delete sg1Data;

  delete d0m0NLL;
  delete d0m1NLL;
  delete d1m0NLL;
  delete d1m1NLL;

  delete pd0m0NLL;
  delete pd0m1NLL;
  delete pd1m0NLL;
  delete pd1m1NLL;

  return;

}


// ----------------------------------------------------------------------
void hstat::run2D(int ntoys, int mode) {

  TH1D *ht0 = new TH1D("ht0", "prof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 2000, -80., 80.); 
  ht0->SetLineColor(kBlue); 
  TH1D *ht1 = new TH1D("ht1", "prof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 2000, -80., 80.); 
  ht1->SetLineColor(kRed); 

  TH1D *hsg0 = (TH1D*)gROOT->Get("hsg0"); 
  if (0 == hsg0) {
    RooRandom::randomGenerator()->SetSeed(fRndmSeed);
    hsg0 = new TH1D("hsg0", "hsg0", 50, 0., 200.); 
  } else { 
    hsg0->Reset();
  }
  TH1D *hsg1 = (TH1D*)gROOT->Get("hsg1"); 
  if (0 == hsg1) {
    hsg1 = new TH1D("hsg1", "hsg1", 50, 0., 400.); 
  } else {
    hsg1->Reset();
  }
  TH1D *hrat = (TH1D*)gROOT->Get("hrat");
  if (0 == hrat) {
    hrat = new TH1D("hrat", "hrat", 50, 0.5, 3.0); 
  } else {
    hrat->Reset(); 
  }

  if ((0 == mode || -1 == mode)) {
    //     fSg1Tau = fSg0Tau; 
    //     fBgTau = fSg0Tau; 
    fSg0Tau = fSG0TAU;
    fSg1Tau = fSG1TAU;
    fBgTau  = fBGTAU;
    fSg0Mu  = fSG0MU;
    fSg1Mu  = fSG1MU;
    fBgMu   = fBGMU;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "No systematics: " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << " SG0N = " << fSG0 << endl;
    cout << "Input: sg1N = " << fSg1 << " SG1N = " << fSG1 << endl;
    cout << "Input: bgN  = " << fBg << " BGN = " << fBG << endl;
    cout << "Input: sg0T = " << fSg0Tau << " SG0T = " << fSG0TAU << endl;
    cout << "Input: sg0M = " << fSg0Mu << " SG0MU = " << fSG0MU << endl;
    cout << "Input: sg1T = " << fSg1Tau << " SG1T = " << fSG1TAU << endl;
    cout << "Input: sg1M = " << fSg1Mu << " SG1MU = " << fSG1MU << endl;
    cout << "Input: bgT  = " << fBgTau << " BGT = " << fBGTAU << endl;
    cout << "Input: bgM  = " << fBgMu << " BGMU = " << fBGMU << endl; 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << ", fSg0/Sg1/BgTau = " << fSg0Tau << "/" << fSg0Tau << "/" << fBgTau
	   << ", fSg0/Sg1/BgMu = " << fSg0Mu << "/" << fSg0Mu << "/" << fBgMu
	   << endl;
      fDoPlot = false; 
      toy2D();
    
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
    // -- make plot
    int bacSeed = RooRandom::randomGenerator()->GetSeed();
    cout << "create plot" << endl;
    RooRandom::randomGenerator()->SetSeed(1234);
    fDoPlot = true; 
    toy2D();
    fDoPlot = false; 
    RooRandom::randomGenerator()->SetSeed(bacSeed);
  } else if (101 == TMath::Abs(mode)) {
    double syst(0.17); 
    if (mode > 0) {
      syst = 1. + syst;
      // evernote:///view/40567474/s295/local/x-coredata://B740E98C-51E1-432D-B0F3-002D230CF640/ENNote/p7213/
      fSg0Tau = 2.20;
      fSg0Mu  = 40.1;
    } else {
      syst = 1. - syst;
      fSg0Tau = 2.22;
      fSg0Mu  = 40.5;
    }
    fSg0 = fSG0*syst; 
    fSg1 = fSG1;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from missing top calculations, syst = " << syst 
	 << ", " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "Input: sg0T = " << fSg0Tau << endl;
    cout << "Input: sg0M = " << fSg0Mu << endl;
    cout << "Input: sg1T = " << fSg1Tau << endl;
    cout << "Input: sg1M = " << fSg1Mu << endl;
    cout << "Input: bgT  = " << fBgTau << endl;
    cout << "Input: bgM  = " << fBgMu << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << ", fSg0/Sg1/BgTau = " << fSg0Tau << "/" << fSg0Tau << "/" << fBgTau
	   << ", fSg0/Sg1/BgMu = " << fSg0Mu << "/" << fSg0Mu << "/" << fBgMu
	   << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      fDoPlot = false; 
      toy2D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
    // -- make plot
    cout << "create plot" << endl;
    int bacSeed = RooRandom::randomGenerator()->GetSeed();
    RooRandom::randomGenerator()->SetSeed(1234);
    fDoPlot = true; 
    toy2D();
    fDoPlot = false; 
    RooRandom::randomGenerator()->SetSeed(bacSeed);
  } else if (102 == TMath::Abs(mode)) {
    double syst(0.4); 
    double systTau(0.0), systMu(0); 
    if (mode > 0) {
      syst = 1. + syst;
      fBgTau = fBGTAU + systTau;
      fBgMu  = fBGMU + systMu;
    } else {
      syst = 1. - syst;
      fBgTau = fBGTAU - systTau;
      fBgMu  = fBGMU - systMu;
    }
    fSg0    = fSG0;
    fSg1    = fSG1;
    fSg0Tau = fSG0TAU;
    fSg1Tau = fSG1TAU;
    fBg     = fBG*syst;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from background uncertainty, syst = " << syst 
	 << ", " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "Input: sg0T = " << fSg0Tau << endl;
    cout << "Input: sg0M = " << fSg0Mu << endl;
    cout << "Input: sg1T = " << fSg1Tau << endl;
    cout << "Input: sg1M = " << fSg1Mu << endl;
    cout << "Input: bgT  = " << fBgTau << endl;
    cout << "Input: bgM  = " << fBgMu << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << ", fSg0/Sg1/BgTau = " << fSg0Tau << "/" << fSg0Tau << "/" << fBgTau
	   << ", fSg0/Sg1/BgMu = " << fSg0Mu << "/" << fSg0Mu << "/" << fBgMu
	   << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      fDoPlot = false; 
      toy2D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
    // -- make plot
    cout << "create plot" << endl;
    int bacSeed = RooRandom::randomGenerator()->GetSeed();
    RooRandom::randomGenerator()->SetSeed(1234);
    fDoPlot = true; 
    toy2D();
    fDoPlot = false; 
    RooRandom::randomGenerator()->SetSeed(bacSeed);
  } else if (103 == TMath::Abs(mode)) {
    double syst(0.2); 
    double systTau(0.0), systMu(0.0); 
    if (mode > 0) {
      syst = 1. + syst;
      fSg0Tau = fSG0TAU + systTau;
      fSg0Mu  = fSG0MU + systMu;
      fSg1Tau = fSG1TAU + systTau;
      fSg1Mu  = fSG1MU + systMu;
    } else {
      syst = 1. - syst;
      fSg0Tau = fSG0TAU - systTau;
      fSg0Mu  = fSG0MU - systMu;
      fSg1Tau = fSG1TAU - systTau;
      fSg1Mu  = fSG1MU - systMu;
    }
    fSg0    = fSG0*syst;
    fSg1    = fSG1*syst;
    fBg     = fBG;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from scale uncertainty, syst = " << syst 
	 << ", " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "Input: sg0T = " << fSg0Tau << endl;
    cout << "Input: sg0M = " << fSg0Mu << endl;
    cout << "Input: sg1T = " << fSg1Tau << endl;
    cout << "Input: sg1M = " << fSg1Mu << endl;
    cout << "Input: bgT  = " << fBgTau << endl;
    cout << "Input: bgM  = " << fBgMu << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << ", fSg0/Sg1/BgTau = " << fSg0Tau << "/" << fSg0Tau << "/" << fBgTau
	   << ", fSg0/Sg1/BgMu = " << fSg0Mu << "/" << fSg0Mu << "/" << fBgMu
	   << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      fDoPlot = false; 
      toy2D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
    int bacSeed = RooRandom::randomGenerator()->GetSeed();
    RooRandom::randomGenerator()->SetSeed(1234);
    fDoPlot = true; 
    toy2D();
    RooRandom::randomGenerator()->SetSeed(bacSeed);
  } else if (104 == TMath::Abs(mode)) {
    double syst(0.20); 
    fSg0Tau = fSG0TAU;
    fSg1Tau = fSG1TAU;
    fBgTau  = fBGTAU;
    fSg0Mu  = fSG0MU;
    fSg1Mu  = fSG1MU;
    fBgMu   = fBGMU;
    if (mode > 0) {
      syst = 1. + syst;
    } else {
      syst = 1. - syst;
    }
    fSg0 = fSG0*syst; 
    fSg1 = fSG1*syst;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from efficiency systematics, syst = " << syst 
	 << ", " << RooRandom::randomGenerator()->GetSeed() << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "Input: sg0T = " << fSg0Tau << endl;
    cout << "Input: sg0M = " << fSg0Mu << endl;
    cout << "Input: sg1T = " << fSg1Tau << endl;
    cout << "Input: sg1M = " << fSg1Mu << endl;
    cout << "Input: bgT  = " << fBgTau << endl;
    cout << "Input: bgM  = " << fBgMu << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 
	   << ", fBg = " << fBg 
	   << ", fSg0/Sg1/BgTau = " << fSg0Tau << "/" << fSg0Tau << "/" << fBgTau
	   << ", fSg0/Sg1/BgMu = " << fSg0Mu << "/" << fSg0Mu << "/" << fBgMu
	   << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      fDoPlot = false; 
      toy2D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
    // -- make plot
    cout << "create plot" << endl;
    int bacSeed = RooRandom::randomGenerator()->GetSeed();
    RooRandom::randomGenerator()->SetSeed(1234);
    fDoPlot = true; 
    toy2D();
    fDoPlot = false; 
    RooRandom::randomGenerator()->SetSeed(bacSeed);
  } else if (1 == mode) {
    double syst(0.17); 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from missing top calculations" << endl;
    cout << "Input: sg0N = " << fSg0 << ", varied by " << 100*syst << "%" << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "Input: sg0T = " << fSg0Tau << endl;
    cout << "Input: sg1T = " << fSg1Tau << endl;
    cout << "Input: bgT  = " << fBgTau << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    double randFact(1.0); 
    for (int i = 0; i < ntoys; ++i) {
      randFact = 1. + gRandom->Gaus(0., syst);
      while (randFact < 0.) {
	randFact = 1. + gRandom->Gaus(0., syst);
      }
      fSg0 = fSG0*randFact; 
      fSg1 = fSG1;
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy2D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } else if (3 == mode) {
    double syst(0.20); 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "systematics from scale uncertainties" << endl;
    cout << "Input: sg0N = " << fSg0 << ", varied by " << 100*syst << "%" << endl;
    cout << "Input: sg1N = " << fSg1 << ", varied through delta(fSg1/fSg0) = " << 100*syst << "%" <<  endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    double randFact(1.0); 
    for (int i = 0; i < ntoys; ++i) {
      randFact = 1. + gRandom->Gaus(0., syst);
      while (randFact < 0.) {
	randFact = 1. + gRandom->Gaus(0., syst);
      }
      fSg0 = fSG0*randFact; 
      fSg1 = fSG1*randFact;
      cout << "XXXXXX run " << i << ", fSg0 = " << fSg0 << ", fSg1 = " << fSg1 << ", fSg1/fSg0 = " << fSg1/fSg0 << endl;
      hsg0->Fill(fSg0); 
      hsg1->Fill(fSg1); 
      hrat->Fill(fSg1/fSg0); 
      toy2D();
      
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  } 


  TLatex *tl = new TLatex(); 
  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.04);

  TArrow *ta = new TArrow(); 

  gStyle->SetOptStat(111111); 
  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);
  c0->Clear(); 
  c0->Divide(2,2);
  c0->cd(1);
  hsg0->Draw();
  tl->DrawLatex(0.15, 0.92, Form("rel rms: %4.3f", hsg0->GetRMS()/hsg0->GetMean())); 
  ta->DrawArrow(fSG0, 0.5*hsg0->GetMaximum(), fSG0, 0.); 
  c0->cd(2); 
  hsg1->Draw();
  tl->DrawLatex(0.15, 0.92, Form("rel rms: %4.3f", hsg1->GetRMS()/hsg1->GetMean())); 
  ta->DrawArrow(fSG1, 0.5*hsg1->GetMaximum(), fSG1, 0.); 
  c0->cd(3); 
  hrat->Draw();
  tl->DrawLatex(0.15, 0.92, Form("rel rms: %4.3f", hrat->GetRMS()/hrat->GetMean())); 
  ta->DrawArrow(fSG1/fSG0, 0.5*hrat->GetMaximum(), fSG1/fSG0, 0.); 

  c0->cd(4);
  //  gStyle->SetOptStat(0); 

  double midpoint, tailprob, sigma, lo, hi;
  findMidPoint(ht0, ht1, midpoint, tailprob, lo, hi); 
  //  sigma = 2.*oneSidedGaussianSigma(tailprob);
  sigma = 2.*RooStats::PValueToSignificance(tailprob);
  fSeparation = sigma; 
  fSeparationE = 2.*(RooStats::PValueToSignificance(hi) - RooStats::PValueToSignificance(lo)); 
  if (fSeparationE < 1.e-3) fSeparationE = 1.e-3; 

  // -- new error on separation according to 
  // evernote:///view/40567474/s295/local/x-coredata://B740E98C-51E1-432D-B0F3-002D230CF640/ENNote/p7145/
  fSeparationE = 0.02;
  if (ntoys < 2200) fSeparationE = 0.02*fSeparation; 
  if (ntoys < 1200) fSeparationE = 0.02*fSeparation; 
  if (ntoys < 700) fSeparationE = 0.03*fSeparation; 
  if (ntoys < 200) fSeparationE = 0.08*fSeparation; 

  ht0->Rebin(10)->Draw();
  ht1->Rebin(10)->Draw("same");
  
  tl->DrawLatex(0.17, 0.85, Form("ntoys/mode: %d/%d", ntoys, mode)); 
  tl->DrawLatex(0.17, 0.80, Form("midpoint: %4.3f", midpoint)); 
  tl->DrawLatex(0.17, 0.75, Form("tailprob: %4.3f", tailprob)); 
  tl->DrawLatex(0.17, 0.70, Form("Sep: %4.3f +/- %4.3f", sigma, fSeparationE)); 
  
  cout << "----------------------------------------------------------------------" << endl;
  cout << "midpoint:  " << midpoint << endl;
  cout << "sigma:     " << sigma << " +/- " << fSeparationE << endl;
  cout << "sigma(lo): " << 2.*RooStats::PValueToSignificance(lo) << endl;
  cout << "sigma(hi): " << 2.*RooStats::PValueToSignificance(hi) << endl;
  cout << "----------------------------------------------------------------------" << endl;
  
  if (mode != -1) c0->SaveAs(Form("hstat2d-%s-%d-%d.pdf", fSetup.c_str(), mode, ntoys)); 
  
}


// ----------------------------------------------------------------------
void hstat::toy2D() {

#ifdef EXPO
  model *m0 = genModel2(0, fSG0, fBG, fSG0TAU, fBGTAU); 
  model *m1 = genModel2(1, fSG1, fBG, fSG1TAU, fBGTAU); 
#else
  model *m0 = genModel2(0, fSG0, fBG, fSG0TAU, fBGTAU, fSG0MU, fBGMU); 
  model *m1 = genModel2(1, fSG1, fBG, fSG1TAU, fBGTAU, fSG1MU, fBGMU); 
#endif

  RooDataSet *d0(0), *d1(0), *bgData(0), *sg0Data(0), *sg1Data(0);

  // -- set pT slopes for the toy data generation
  m0->bgTau->setVal(fBgTau);
  m0->bgMu->setVal(fBgMu);

  m0->sgTau->setVal(fSg0Tau);
  m0->sgMu->setVal(fSg0Mu);

  m1->sgTau->setVal(fSg1Tau);
  m1->sgMu->setVal(fSg1Mu);

  bgData  = m0->bgPdf->generate(RooArgSet(*m0->m, *m0->pt), fBg); 
  sg0Data = m0->sgPdf->generate(RooArgSet(*m0->m, *m0->pt), fSg0); 
  sg1Data = m1->sgPdf->generate(RooArgSet(*m1->m, *m1->pt), fSg1); 

  // -- reset sg0tau and sg1tau to the constant default values!
  m0->sgTau->setVal(fSG0TAU);
  m0->sgMu->setVal(fSG0MU);

  m1->sgTau->setVal(fSG1TAU);
  m1->sgMu->setVal(fSG1MU);

  // -- model 0 plus background
  d0 = new RooDataSet(*bgData);
  d0->append(*sg0Data);
  d1 = new RooDataSet(*bgData);
  d1->append(*sg1Data);

  RooAbsReal *d0m0NLL = m0->modelPdf->createNLL(*d0, Extended(kTRUE), Verbose(kFALSE), NumCPU(2));
  RooAbsReal *d0m1NLL = m1->modelPdf->createNLL(*d0, Extended(kTRUE), Verbose(kFALSE), NumCPU(2));
  RooAbsReal *d1m0NLL = m0->modelPdf->createNLL(*d1,  Extended(kTRUE), Verbose(kFALSE), NumCPU(2));
  RooAbsReal *d1m1NLL = m1->modelPdf->createNLL(*d1,  Extended(kTRUE), Verbose(kFALSE), NumCPU(2));
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooProfileLL *pd0m0NLL = dynamic_cast<RooProfileLL*>(d0m0NLL->createProfile(m0->poi));
  fD0M0 = pd0m0NLL->getVal(); 
  RooProfileLL *pd0m1NLL = dynamic_cast<RooProfileLL*>(d0m1NLL->createProfile(m1->poi));
  fD0M1 = pd0m1NLL->getVal(); 
  RooProfileLL *pd1m0NLL = dynamic_cast<RooProfileLL*>(d1m0NLL->createProfile(m0->poi));
  fD1M0 = pd1m0NLL->getVal(); 
  RooProfileLL *pd1m1NLL = dynamic_cast<RooProfileLL*>(d1m1NLL->createProfile(m1->poi));
  fD1M1 = pd1m1NLL->getVal(); 


  if (fDoPlot) {

    TLatex *tl = new TLatex(); 
    tl->SetNDC(kTRUE);
    tl->SetTextSize(0.04);
    
    TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0"); 
    if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);

    // -- data 0 
    c0->Clear();
    c0->Divide(2,2);

    c0->cd(1);
    gPad->SetLogy(0);
    RooPlot *fd0m0 = m0->m->frame(Title("d0"), Name("mass"), Range(MGGLO, MGGHI));
    d0->plotOn(fd0m0);
    m0->modelPdf->plotOn(fd0m0);
    m0->modelPdf->plotOn(fd0m0, Components("bgPdf"), LineStyle(kDashed)) ;
    m0->modelPdf->paramOn(fd0m0, Layout(0.5, 0.8, 0.45));
    fd0m0->Draw();
    tl->DrawLatex(0.2, 0.8, Form("fD0M0 = %5.4f", fD0M0)); 

    c0->cd(2);
    gPad->SetLogy(1);
    RooPlot *fd0pt0 = m0->pt->frame(Title("d0"), Name("pt"), Range(300., 1000.));
    d0->plotOn(fd0pt0);
    m0->modelPdf->plotOn(fd0pt0);
    m0->modelPdf->plotOn(fd0pt0, Components("bgPdf"), LineStyle(kDashed)) ;
    m0->modelPdf->paramOn(fd0pt0, Layout(0.5, 0.8, 0.45));
    fd0pt0->Draw();
    tl->DrawLatex(0.2, 0.8, Form("fD0M0 = %5.4f", fD0M0)); 

    c0->cd(3);
    gPad->SetLogy(0);
    RooPlot *fd0m1 = m1->m->frame(Title("d0"), Name("mass"), Range(MGGLO, MGGHI));
    d0->plotOn(fd0m1);
    m1->modelPdf->plotOn(fd0m1);
    m1->modelPdf->plotOn(fd0m1, Components("bgPdf"), LineStyle(kDashed)) ;
    m1->modelPdf->paramOn(fd0m1, Layout(0.5, 0.8, 0.45));
    fd0m1->Draw();
    tl->DrawLatex(0.2, 0.8, Form("fD0M1 = %5.4f", fD0M1)); 

    c0->cd(4);
    gPad->SetLogy(1);
    RooPlot *fd0pt1 = m1->pt->frame(Title("d0"), Name("pt"), Range(300., 1000.));
    d0->plotOn(fd0pt1);
    m1->modelPdf->plotOn(fd0pt1);
    m1->modelPdf->plotOn(fd0pt1, Components("bgPdf"), LineStyle(kDashed)) ;
    m1->modelPdf->paramOn(fd0pt1, Layout(0.5, 0.8, 0.45));
    fd0pt1->Draw();
    tl->DrawLatex(0.2, 0.8, Form("fD0M1 = %5.4f", fD0M1)); 
    c0->SaveAs(Form("%s/d0-%s.pdf", fDirectory.c_str(), fDoPlotName.c_str())); 


    // -- data 1 
    c0->Clear();
    c0->Divide(2,2);

    c0->cd(1);
    gPad->SetLogy(0);
    RooPlot *fd1m0 = m0->m->frame(Title("d1"), Name("mass"), Range(MGGLO, MGGHI));
    d1->plotOn(fd1m0);
    m0->modelPdf->plotOn(fd1m0);
    m0->modelPdf->plotOn(fd1m0, Components("bgPdf"), LineStyle(kDashed)) ;
    m0->modelPdf->paramOn(fd1m0, Layout(0.5, 0.8, 0.45));
    fd1m0->Draw();
    tl->DrawLatex(0.2, 0.8, Form("fD1M0 = %5.4f", fD1M0)); 

    c0->cd(2);
    gPad->SetLogy(1);
    RooPlot *fd1pt0 = m0->pt->frame(Title("d1"), Name("pt"), Range(300., 1000.));
    d1->plotOn(fd1pt0);
    m0->modelPdf->plotOn(fd1pt0);
    m0->modelPdf->plotOn(fd1pt0, Components("bgPdf"), LineStyle(kDashed)) ;
    m0->modelPdf->paramOn(fd1pt0, Layout(0.5, 0.8, 0.45));
    fd1pt0->Draw();
    tl->DrawLatex(0.2, 0.8, Form("fD1M0 = %5.4f", fD1M0)); 

    c0->cd(3);
    gPad->SetLogy(0);
    RooPlot *fd1m1 = m1->m->frame(Title("d1"), Name("mass"), Range(MGGLO, MGGHI));
    d1->plotOn(fd1m1);
    m1->modelPdf->plotOn(fd1m1);
    m1->modelPdf->plotOn(fd1m1, Components("bgPdf"), LineStyle(kDashed)) ;
    m1->modelPdf->paramOn(fd1m1, Layout(0.5, 0.8, 0.45));
    fd1m1->Draw();
    tl->DrawLatex(0.2, 0.8, Form("fD1M1 = %5.4f", fD1M1)); 

    c0->cd(4);
    gPad->SetLogy(1);
    RooPlot *fpt1 = m1->pt->frame(Title("d1"), Name("pt"), Range(300., 1000.));
    d1->plotOn(fpt1);
    m1->modelPdf->plotOn(fpt1);
    m1->modelPdf->plotOn(fpt1, Components("bgPdf"), LineStyle(kDashed)) ;
    m1->modelPdf->paramOn(fpt1, Layout(0.5, 0.8, 0.45));
    fpt1->Draw();
    tl->DrawLatex(0.2, 0.8, Form("fD1M1 = %5.4f", fD1M1)); 
    c0->SaveAs(Form("%s/d1-%s.pdf", fDirectory.c_str(), fDoPlotName.c_str())); 

    
    // -- frame for FITTED models
    int nbinsM(55), nbinsPt(70);
    c0->Clear();
    gPad->SetLogy(0);
    RooFitResult *r0 = m0->modelPdf->fitTo(*d0, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
    shrinkPad(0.15, 0.15);
    RooPlot *fM0 = m0->m->frame(Title(" "), Range(MGGLO, MGGHI), Bins(nbinsM));
    d0->plotOn(fM0, MarkerColor(kBlue), MarkerStyle(24), LineColor(kBlue), MarkerSize(2));
    m0->modelPdf->plotOn(fM0, LineColor(kBlack));
    m0->modelPdf->plotOn(fM0, Components("bgPdf"), LineColor(kRed), LineStyle(kDashed)) ;
    //    m0->modelPdf->paramOn(fM0, Layout(0.5, 0.8, 0.45));
    setTitles(fM0, "m [GeV]", Form("Events / %d GeV", static_cast<int>((MGGHI-MGGLO)/nbinsM)));
    fM0->Draw();
    tl->SetTextFont(42);
    tl->DrawLatex(0.15, 0.93, "Background plus Higgs m_{top} = 173.5 GeV");
    tl->DrawLatex(0.2, 0.82, "#sqrt{s} = 14 TeV");
    tl->DrawLatex(0.2, 0.77, "L = 1000 fb^{-1}");
    c0->SaveAs(Form("%s/fit0-m-%s.pdf", fDirectory.c_str(), fDoPlotName.c_str())); 

    gPad->SetLogy(1);
    shrinkPad(0.15, 0.15);
    RooPlot *fPt0 = m0->pt->frame(Title(" "), Range(300., 1000.), Bins(nbinsPt));
    // -- display the background data as well
    //    bgData->plotOn(fPt0, MarkerColor(kRed), MarkerStyle(25), LineColor(kRed), MarkerSize(1));
    d0->plotOn(fPt0, MarkerColor(kBlue), MarkerStyle(24), LineColor(kBlue), MarkerSize(2));
    // fPt0->Print();
    RooHist *rh = fPt0->getHist("h_bgPdfData");
    removeEmptyBins(rh); 

    m0->modelPdf->plotOn(fPt0, LineColor(kBlack));
    m0->modelPdf->plotOn(fPt0, Components("bgPdf"), LineColor(kRed), LineStyle(kDashed)) ;
    setTitles(fPt0, "p_{T} [GeV]", Form("Events / %d GeV", static_cast<int>((1000.-300.)/nbinsPt)), 0.05, 1.1, 1.3);
    fPt0->SetNdivisions(405);
    fPt0->Draw();
    tl->SetTextFont(42);
    tl->DrawLatex(0.15, 0.93, "Background plus Higgs m_{top} = 173.5 GeV");
    tl->DrawLatex(0.5, 0.82, "#sqrt{s} = 14 TeV");
    tl->DrawLatex(0.5, 0.77, "L = 1000 fb^{-1}");
    c0->SaveAs(Form("%s/fit0-pt-%s.pdf", fDirectory.c_str(), fDoPlotName.c_str())); 

    gPad->SetLogy(0);
    RooFitResult *r1 = m1->modelPdf->fitTo(*d1, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
    shrinkPad(0.15, 0.15);
    RooPlot *fM1 = m1->m->frame(Title(" "), Range(MGGLO, MGGHI), Bins(nbinsM));
    d1->plotOn(fM1, MarkerColor(kGreen+2), MarkerStyle(26), LineColor(kGreen+2), MarkerSize(2));
    removeEmptyBins(rh); 
    m1->modelPdf->plotOn(fM1, LineColor(kBlack));
    m1->modelPdf->plotOn(fM1, Components("bgPdf"), LineColor(kRed), LineStyle(kDashed)) ;
    //    m1->modelPdf->paramOn(fM1, Layout(0.5, 0.8, 0.45));
    setTitles(fM1, "m [GeV]", Form("Events / %d GeV", static_cast<int>((MGGHI-MGGLO)/nbinsM)));
    fM1->Draw();
    tl->SetTextFont(42);
    tl->DrawLatex(0.15, 0.93, "Background plus Higgs m_{top} #rightarrow #infty");
    tl->DrawLatex(0.2, 0.82, "#sqrt{s} = 14 TeV");
    tl->DrawLatex(0.2, 0.77, "L = 1000 fb^{-1}");
    c0->SaveAs(Form("%s/fit1-m-%s.pdf", fDirectory.c_str(), fDoPlotName.c_str())); 

    gPad->SetLogy(1);
    shrinkPad(0.15, 0.15);
    RooPlot *fPt1 = m1->pt->frame(Title(" "), Range(300., 1000.), Bins(nbinsPt));
    d1->plotOn(fPt1, MarkerColor(kGreen+2), MarkerStyle(26), LineColor(kGreen+2), MarkerSize(2));
    // -- display the background data as well
    //    bgData->plotOn(fPt1, MarkerColor(kRed), MarkerStyle(25), LineColor(kRed), MarkerSize(1));
    rh = fPt1->getHist("h_bgPdfData");
    removeEmptyBins(rh); 
    m1->modelPdf->plotOn(fPt1, LineColor(kBlack));
    m1->modelPdf->plotOn(fPt1, Components("bgPdf"), LineColor(kRed), LineStyle(kDashed)) ;
    //    m1->modelPdf->paramOn(fPt1, Layout(0.5, 0.8, 0.45));
    setTitles(fPt1, "p_{T} [GeV]", Form("Events / %d GeV", static_cast<int>((1000.-300.)/nbinsPt)), 0.05, 1.1, 1.3);
    fPt1->SetNdivisions(405);
    fPt1->Draw();
    tl->SetTextFont(42);
    tl->DrawLatex(0.15, 0.93, "Background plus Higgs m_{top} #rightarrow #infty");
    tl->DrawLatex(0.5, 0.82, "#sqrt{s} = 14 TeV");
    tl->DrawLatex(0.5, 0.77, "L = 1000 fb^{-1}");
    c0->SaveAs(Form("%s/fit1-pt-%s.pdf", fDirectory.c_str(), fDoPlotName.c_str())); 

  }



  delModel(m0); 
  delModel(m1); 

  delete d0;
  delete d1;
  delete bgData;
  delete sg0Data; 
  delete sg1Data;

  delete d0m0NLL;
  delete d0m1NLL;
  delete d1m0NLL;
  delete d1m1NLL;

  delete pd0m0NLL;
  delete pd0m1NLL;
  delete pd1m0NLL;
  delete pd1m1NLL;

  return;

}


// ----------------------------------------------------------------------
void hstat::delModel(model *m) {
  delete m->m;
  delete m->pt; 
  // -- fit (fixed) parameters:
  delete m->massP;
  delete m->massS;
  delete m->bgSlope;
  delete m->sgTau;
  delete m->bgTau;
  delete m->sgN;
  delete m->bgN;
  delete m->sgPdf;
  delete m->bgPdf; 
  delete m->modelPdf; 

  if (m->sgPt) delete m->sgPt;
  if (m->bgPt) delete m->bgPt;
  if (m->sgM) delete m->sgM;
  if (m->bgM) delete m->bgM;

  if (m->sgMu) delete m->sgMu;
  if (m->bgMu) delete m->bgMu;

  delete m;
}


// ----------------------------------------------------------------------
model* hstat::genModel1(int mode, double nsg, double nbg) {

  model *aModel = new model(); 

  // -- variables
  aModel->m  = new RooRealVar("m", "m", MGGLO, MGGHI, "GeV"); 

  // -- parameters
  aModel->massP = new RooRealVar(Form("m%d_mP", mode), "signal peak mass", fHiggsMpeak, MGGLO, MGGHI);  
  aModel->massP->setConstant(kTRUE);
  aModel->massS = new RooRealVar(Form("m%d_mS", mode), "signal sigma mass", fHiggsMres, 0., 15.);  
  aModel->massS->setConstant(kTRUE); 
  aModel->bgSlope = new RooRealVar(Form("m%d_bgSlope", mode), "coefficient #0 for bg", fBgMp1, 0., 1.); 

  aModel->sgN = new RooRealVar(Form("m%d_sgN", mode), "signal events", nsg, 0., 1.e5);
  aModel->bgN = new RooRealVar(Form("m%d_bgN", mode), "background events", nbg, 0., 1.e5);

  // -- create PDF
  aModel->bgPdf = new RooPolynomial(Form("m%d_bgM", mode), "background gamma gamma mass", 
				    *aModel->m, RooArgList(*aModel->bgSlope)); 
  aModel->sgPdf = new RooGaussian(Form("m%d_sgM", mode), "signal mass", 
				  *aModel->m, *aModel->massP, *aModel->massS); 
  aModel->modelPdf = new RooAddPdf(Form("m%d_model", mode), "model", 
				   RooArgList(*aModel->sgPdf, *aModel->bgPdf), 
				   RooArgList(*aModel->sgN, *aModel->bgN));

  aModel->sgM  = 0;   	     
  aModel->bgM  = 0; 
  aModel->sgPt = 0;   	     
  aModel->bgPt = 0; 
  
  // -- define POI
  aModel->poi.add(*aModel->sgN);

  return aModel; 
}



// ----------------------------------------------------------------------
model* hstat::genModel2(int mode, double nsg, double nbg, double tsg, double tbg, double msg, double mbg) {
  
  model *aModel = new model(); 
  
  // -- variables
  aModel->m  = new RooRealVar("m", "m", MGGLO, MGGHI, "GeV"); 
  aModel->pt = new RooRealVar("pt", "pt", 300., 1000., "GeV"); 

  // -- parameters
  aModel->massP = new RooRealVar(Form("m%d_mP", mode), "signal peak mass", fHiggsMpeak, MGGLO, MGGHI);  
  aModel->massP->setConstant(kTRUE);
  aModel->massS = new RooRealVar(Form("m%d_mS", mode), "signal sigma mass", fHiggsMres, 0., 15.);  
  aModel->massS->setConstant(kTRUE); 
  aModel->bgSlope = new RooRealVar(Form("m%d_bgSlope", mode), "coefficient #0 for bg", fBgMp1, 0., 1.); 

  aModel->sgN = new RooRealVar(Form("m%d_sgN", mode), "signal events", nsg, 0., 1.e5);
  aModel->bgN = new RooRealVar(Form("m%d_bgN", mode), "background events", nbg, 0., 1.e5);

  aModel->sgTau = new RooRealVar(Form("m%d_sgTau", mode), "signal tau", tsg, -10., 10.);
  aModel->bgTau = new RooRealVar(Form("m%d_bgTau", mode), "background tau", tbg, -10., 10.);

  aModel->sgMu = new RooRealVar(Form("m%d_sgMu", mode), "signal mu", msg, 0., 1000.);  
  aModel->bgMu = new RooRealVar(Form("m%d_bgMu", mode), "background mu", mbg, 0., 1000.);  

  // -- create PDF
#ifdef EXPO
  //  cout << "creating pT model based on EXPO, mu should be zero: " << msg << " " << mbg << endl;
  aModel->sgPt = new RooExponential(Form("m%d_sgPt", mode), "signal pT", *aModel->pt, *aModel->sgTau);   
  aModel->bgPt = new RooExponential(Form("m%d_bgPt", mode), "background pT", *aModel->pt, *aModel->bgTau);   
#else
  //  cout << "creating pT model based on LOGNORMAL, mu should be non-zero: " << msg << " " << mbg << endl;
  aModel->sgPt = new RooLognormal(Form("m%d_sgPt", mode), "signal pT", *aModel->pt, *aModel->sgMu, *aModel->sgTau);   
  aModel->bgPt = new RooLognormal(Form("m%d_bgPt", mode), "background pT", *aModel->pt, *aModel->bgMu, *aModel->bgTau);
#endif

  aModel->sgM = new RooGaussian(Form("m%d_sgM", mode), "signal gamma gamma mass", 
				*aModel->m, *aModel->massP, *aModel->massS);
  aModel->bgM = new RooPolynomial(Form("m%d_bgM", mode), "background mass", 
				  *aModel->m, RooArgList(*aModel->bgSlope)); 
  
  aModel->sgPdf = new RooProdPdf("sgPdf", "sgPdf", RooArgSet(*aModel->sgM, *aModel->sgPt));
  aModel->bgPdf = new RooProdPdf("bgPdf", "bgPdf", RooArgSet(*aModel->bgM, *aModel->bgPt));

  aModel->modelPdf = new RooAddPdf(Form("m%d_model", mode), "model", 
				   RooArgList(*aModel->sgPdf, *aModel->bgPdf), 
				   RooArgList(*aModel->sgN, *aModel->bgN));
  
  // -- define POI
  aModel->poi.add(*aModel->sgN);
  aModel->poi.add(*aModel->sgTau);
#ifdef EXPO
#else
  aModel->poi.add(*aModel->sgMu);
#endif

  return aModel; 
}


// ----------------------------------------------------------------------
void  hstat::findMidPoint(TH1D* hq0, TH1D* hq1, double &midpoint, double &tailprob, double &lo, double &hi) {
  cout << endl;
  double iq0(0), iq1(0), delta(0.); 
  int nbinsx = hq0->GetNbinsX()+1;
  TH1D *nq0 = (TH1D*)hq0->Clone("nq0");
  //  nq0->Scale(1./nq0->GetSumOfWeights());
  nq0->Scale(1./nq0->Integral());
  TH1D *nq1 = (TH1D*)hq1->Clone("nq1");
  //  nq1->Scale(1./nq1->GetSumOfWeights());
  nq1->Scale(1./nq1->Integral());
  for (int i = 0; i < nbinsx; ++i) {
    iq0 = nq0->Integral(0, i); 
    iq1 = nq1->Integral(i, nbinsx+1); 
    if (iq0 > iq1) {
      midpoint = nq0->GetBinCenter(i); 
      if ((iq0 - iq1) > delta) {
	cout << "taking previous, delta = " << delta << " while this one has iq0 = " << iq0 << " and iq1 = " << iq1 << endl;
	midpoint = nq0->GetBinCenter(i-1); 
	iq0 = nq0->Integral(0, i-1); 
	iq1 = nq1->Integral(i-1, nbinsx+1); 
      }
      break;
    }
    delta = iq1 - iq0; 
  }
    
  // -- take the average to account a bit for binning effects
  tailprob = 0.5*(iq0+iq1); 
  lo = iq0; 
  hi = iq1; 

  cout << "found "
       << Form(" iq0: %5.4f, iq1: %5.4f, diff: %5.4f", iq0, iq1, iq1-iq0)
       << " at midpoint = " << Form("%5.4f", midpoint) 
       << Form(" tail prob: %5.4f, lo: %5.4f, hi: %5.4f", tailprob, lo, hi)
       << endl;



  // -- printout for confirmation
  iq0 = nq0->Integral(0, nq0->FindBin(midpoint)+1); 
  iq1 = nq1->Integral(nq0->FindBin(midpoint)+1, nbinsx+1); 
  cout << " next "
       << Form(" iq0: %5.4f, iq1: %5.4f, diff: %5.4f", iq0, iq1, iq1-iq0)
       << " at " << nq0->GetBinCenter(nq0->FindBin(midpoint)+1)
       << endl;
  
  iq0 = nq0->Integral(0, nq0->FindBin(midpoint)-1); 
  iq1 = nq1->Integral(nq0->FindBin(midpoint)-1, nbinsx); 
  cout << " prev "
       << Form(" iq0: %5.4f, iq1: %5.4f, diff: %5.4f", iq0, iq1, iq1-iq0)
       << " at " << nq0->GetBinCenter(nq0->FindBin(midpoint)-1)
       << endl;

}


// ----------------------------------------------------------------------
double hstat::oneSidedGaussianSigma(double prob) {

  RooRealVar x("x", "x", -10, 10) ;
  RooGaussian gx("gx", "gx", x, RooConst(0), RooConst(1));
  
  int i(0); 
  double qold(0.),  eps(1.e-9); 
  double delta(0.1);
  double sigma(0.); 
  double q(0.); 
  while (1) {
    
    q = 0.5*TMath::Erfc(sigma/TMath::Sqrt(2));
    if (q > prob) {
      sigma += delta;
    } else {
      sigma -= delta;
      delta /= 10.; 
    }

    //    cout << "tail prob " << Form("%4d", i) << ": " << q << " from sigma = " << sigma << endl;

    if (TMath::Abs(qold - q) < eps) {
      //      cout << "qold = " << qold << " q  = " << q << " ->  TMath::Abs(qold - q) = " << TMath::Abs(qold - q) << endl;
      break;
    }
    if (i > 100) {
      //      cout << "i = " << i << endl;
      break;
    }
    ++i;
    qold = q;
  }
  
  cout << "tail probability " << prob << " with local prob = " << q
       << " yields significance " << sigma << " " <<  0.5*TMath::Erfc(sigma/TMath::Sqrt(2)) << endl;

  //   x.setRange("signal", sigma, 10);
  //   RooAbsReal* igx(0);   
  //   igx = gx.createIntegral(x, NormSet(x), Range("signal")) ;
  //   cout << gx.createIntegral(x, NormSet(x), Range("signal"))->getVal() << endl;
  //  cout << "with formula: " <<  << endl;

  return sigma;

}


// ----------------------------------------------------------------------
void hstat::testTools(int nbins) {
  int NBINS(nbins); 
  int ntoys(10000); 
  double rms(8.); 

  TH1D *h0 = new TH1D("h0", "h0", NBINS, -40., 40.); h0->SetLineColor(kBlue); 
  TH1D *h1 = new TH1D("h1", "h1", NBINS, -40., 40.); h1->SetLineColor(kRed); 

  TF1 *f0 = new TF1("f0", "[0]*exp(-0.5*((x-[1])/[2])**2)", -40., 40.); 
  f0->SetParameters(1, 10., rms); 

  TF1 *f1 = new TF1("f1", "[0]*exp(-0.5*((x-[1])/[2])**2)", -40., 40.); 
  f1->SetParameters(1, -11., rms); 

  h0->FillRandom("f0", ntoys); 
  h1->FillRandom("f1", ntoys); 

  h0->Draw();
  h1->Draw("same");

  double midpoint, tailprob, sigma0, sigma1, sigma2, lo, hi;
  findMidPoint(h0, h1, midpoint, tailprob, lo, hi); 
  sigma0 = 2.*oneSidedGaussianSigma(tailprob);
  // the following is also in http://agenda.nikhef.nl/getFile.py/access?sessionId=36&resId=0&materialId=0&confId=1488
  sigma1 = 2.*ROOT::Math::gaussian_quantile_c(tailprob, 1.);
  sigma2 = 2.*RooStats::PValueToSignificance(tailprob); 

  cout << "midpoint: " << Form("%4.3f", midpoint) << endl;
  cout << "tailprob: " << Form("%4.3f", tailprob) << endl;
  cout << "integral(h0; -40, midpoint)* = " << h0->Integral(0, h0->FindBin(midpoint)) << endl;
  cout << "integral(h0; midpoint, 40)   = " << h0->Integral(h0->FindBin(midpoint), h0->GetNbinsX()+1) << endl;
  cout << "integral(h1; -40, midpoint)  = " << h1->Integral(0, h1->FindBin(midpoint)) << endl;
  cout << "integral(h1; midpoint, 40)*  = " << h1->Integral(h1->FindBin(midpoint), h1->GetNbinsX()+1) << endl;
  cout << "sigma0: " << sigma0 << endl;
  cout << "sigma1: " << sigma1 << endl;
  cout << "sigma2: " << sigma2 << endl;


}


// ----------------------------------------------------------------------
// determine modification of fSg0Tau by MS' missing top parametrization
void hstat::expModification(int evt) {

  TF1 *fExp = new TF1("fExp",  "expo", 0., 1000.);
  fExp->SetLineColor(kBlack);
  fExp->SetParameters(1., -0.0117118);
  TF1 *f1 = new TF1("f1",  "1.0 + 0.01*(x-40)*0.01*(x-40)*0.015", 0., 1000.);
  TF1 *f2 = new TF1("f2",  "1.0 - 0.01*(x-40)*0.01*(x-40)*0.015", 0., 1000.);

  TH1D *h0 = new TH1D("h0", "", 80, 200., 1000.); 
  h0->FillRandom("fExp", evt); 

  TH1D *h1 = (TH1D*)h0->Clone("h1");
  h1->Multiply(f1); 
  h1->Scale(evt/h1->Integral());

  TH1D *h2 = (TH1D*)h0->Clone("h2");
  h2->Multiply(f2); 
  h2->Scale(evt/h2->Integral());

  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);

  zone(2,2);
  c0->cd(1);
  gPad->SetLogy(1);
  h0->Fit("expo");

  c0->cd(2);
  gPad->SetLogy(1);
  h1->Fit("expo");

  c0->cd(3);
  gPad->SetLogy(1);
  h2->Fit("expo");

} 

// ----------------------------------------------------------------------
void hstat::setGraph(TGraph *g, int color, int fillStyle) {
  g->SetTitle("");
  g->SetLineWidth(3); 
  g->SetMinimum(0.); 
  g->SetFillStyle(fillStyle); 

  g->GetXaxis()->SetTitle("Luminosity #kern[-0.3]{[}fb^{-1}#kern[0.3]{]}");
  g->GetYaxis()->SetTitle("expected separation [#sigma]");

  g->SetLineColor(color); 
  g->SetMarkerColor(color); 
  g->SetFillColor(color); 
}

// ----------------------------------------------------------------------
// This is the LOGNORMAL version
// show experimental and theoretical errors separately
void hstat::plotResults() {

  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);
  c0->Clear();

  double sig0 = 2.99744;

  // -- various signficance ranges for index
  //    0: sampling error
  //    1: efficiency uncertainty (14 == mode)
  //    2: background uncertainty (12 == mode)
  //    3: scale uncertainty (13 == mode)
  //    4: missing top calculations (11 == mode)

  const int nexp(3), nthy(2); 
  const int nerr(nexp+nthy);
  const double sig0Lo[nerr] = {2.90817, 2.70   , 2.56711, 2.47062, 2.51942};
  const double sig0Hi[nerr] = {3.06369, 3.28   , 3.11513, 3.45626, 3.16759};
  
  double errLo[nerr]; 
  double errHi[nerr];

  for (int i = 0; i < nerr; ++i) {
    errLo[i] = (sig0 - sig0Lo[i])/sig0;
    errHi[i] = (sig0Hi[i] - sig0)/sig0;
  }

  // -- Significances vs lumi:  grep sigma: s10m10-*.log
  const int nlumi(6);
  double x[nlumi] =   {500.,    1000.,   1500.,   2000.,   2500.,   3000.};
  double sig[nlumi] = {1.99009, 2.99744, 3.34701, 3.63724, 3.73092, 4.02801};

  double eBare[nlumi] = {0., 0., 0., 0., 0., 0.};

  double eSysLo[nerr][nlumi], eSysHi[nerr][nlumi];
  
  // -- running quadratic sum for experimental systematic errors
  for (int i = 0; i < nlumi; ++i) {
    for (int j = 0; j < nexp; ++j) {
      if (0 == j) {
	eSysLo[j][i] = errLo[j]*errLo[j];
	eSysHi[j][i] = errHi[j]*errHi[j];
      } else {
	eSysLo[j][i] = eSysLo[j-1][i] + errLo[j]*errLo[j];
	eSysHi[j][i] = eSysHi[j-1][i] + errHi[j]*errHi[j];
      }
      if (0) cout << "a lumi " << i << "/" << j << " eSysLo[" << j << "][" << i << "] = " << eSysLo[j][i] 
	   << " sqrt = (" << TMath::Sqrt(eSysLo[j][i]) << ") " 
	   << " eSysHi[" << j << "][" << i << "] = " << eSysHi[j][i] 
	   << " sqrt = (" << TMath::Sqrt(eSysHi[j][i]) << ") " 
	   << endl;
    }
  }    

  // -- running linear sum for theory errors
  for (int i = 0; i < nlumi; ++i) {
    for (int j = nexp; j < nerr; ++j) {
      if (j == nexp) {
	eSysLo[j][i] = errLo[j];
	eSysHi[j][i] = errHi[j];
      } else {
	eSysLo[j][i] = eSysLo[j-1][i] + errLo[j];
	eSysHi[j][i] = eSysHi[j-1][i] + errHi[j];
      }
      if (0) cout << "b lumi " << i << "/" << j << " eSysLo[" << j << "][" << i << "] = " << eSysLo[j][i] 
	   << " sqrt = (" << TMath::Sqrt(eSysLo[j][i]) << ") " 
	   << " eSysHi[" << j << "][" << i << "] = " << eSysHi[j][i] 
	   << " sqrt = (" << TMath::Sqrt(eSysHi[j][i]) << ") " 
	   << endl;
    }
  }

  // -- sqrt for systematic errors
  for (int i = 0; i < nlumi; ++i) {
    for (int j = 0; j < nexp; ++j) {
      eSysLo[j][i] = sig[i]*TMath::Sqrt(eSysLo[j][i]); 
      eSysHi[j][i] = sig[i]*TMath::Sqrt(eSysHi[j][i]); 
      if (0) cout << "c lumi " << i << " eSysLo[" << j << "][" << i << "] = " << eSysLo[j][i] 
	   << " (" << TMath::Sqrt(eSysLo[j][i]) << ") "
	   << " eSysHi[" << j << "][" << i << "] = " << eSysHi[j][i] 
	   << " (" << TMath::Sqrt(eSysHi[j][i]) << ") "
	   << endl;
    }
  }

  // -- no sqrt for theory errors
  for (int i = 0; i < nlumi; ++i) {
    for (int j = nexp; j < nerr; ++j) {
      eSysLo[j][i] = sig[i]*eSysLo[j][i]; 
      eSysHi[j][i] = sig[i]*eSysHi[j][i]; 
      if (0) cout << "c lumi " << i << " eSysLo[" << j << "][" << i << "] = " << eSysLo[j][i] 
	   << " (" << TMath::Sqrt(eSysLo[j][i]) << ") "
	   << " eSysHi[" << j << "][" << i << "] = " << eSysHi[j][i] 
	   << " (" << TMath::Sqrt(eSysHi[j][i]) << ") "
	   << endl;
    }
  }


  // -- display systematic errors
  TGraphAsymmErrors *gBare = new TGraphAsymmErrors(nlumi, x, sig, eBare, eBare, eBare, eBare); 
  setGraph(gBare, kBlack, 1000); 

  double errGrLo[nlumi], errGrHi[nlumi];
  for (int i = 0; i < nlumi; ++i) {
    errGrLo[i] = eSysLo[0][i];
    errGrHi[i] = eSysHi[0][i];
  }
  TGraphAsymmErrors *gSys0 = new TGraphAsymmErrors(nlumi, x, sig, eBare, eBare, errGrLo, errGrHi); 
  setGraph(gSys0, kCyan+3, 1000); 

  for (int i = 0; i < nlumi; ++i) {
    errGrLo[i] = eSysLo[1][i];
    errGrHi[i] = eSysHi[1][i];
  }
  TGraphAsymmErrors *gSys1 = new TGraphAsymmErrors(nlumi, x, sig, eBare, eBare, errGrLo, errGrHi); 
  setGraph(gSys1, kCyan+2, 1000); 

  for (int i = 0; i < nlumi; ++i) {
    errGrLo[i] = eSysLo[2][i];
    errGrHi[i] = eSysHi[2][i];
  }
  TGraphAsymmErrors *gSys2 = new TGraphAsymmErrors(nlumi, x, sig, eBare, eBare, errGrLo, errGrHi); 
  setGraph(gSys2, kCyan+1, 1000); 

  double YMAX(5.2);
  gSys1->SetMaximum(YMAX);

  gSys2->Draw("a3");
  gSys1->Draw("3");
  gSys0->Draw("3");
  gBare->Draw("l");

  TLatex *tl = new TLatex(); 
  tl->SetNDC(kTRUE);
  tl->SetTextFont(42);
  tl->DrawLatex(0.2, 0.8, "#sqrt{s} = 14 TeV");

  TLegend *legg = new TLegend(0.3, 0.2, 0.8, 0.4); 
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 

  legg->AddEntry(gBare, "central expectation", "l"); 
  legg->AddEntry(gSys0, "sampling uncertainty", "f"); 
  legg->AddEntry(gSys1, "efficiency uncertainty", "f"); 
  legg->AddEntry(gSys2, "background uncertainty", "f"); 
  
  legg->Draw();

  c0->SaveAs(Form("%s/result-sys-lumi.pdf", fDirectory.c_str())); 



  // -- display theory errors
  for (int i = 0; i < nlumi; ++i) {
    errGrLo[i] = eSysLo[nexp][i];
    errGrHi[i] = eSysHi[nexp][i];
  }
  TGraphAsymmErrors *gThy0 = new TGraphAsymmErrors(nlumi, x, sig, eBare, eBare, errGrLo, errGrHi); 
  setGraph(gThy0, kGreen-1, 1000); 

  for (int i = 0; i < nlumi; ++i) {
    errGrLo[i] = eSysLo[nexp+1][i];
    errGrHi[i] = eSysHi[nexp+1][i];
  }
  TGraphAsymmErrors *gThy1 = new TGraphAsymmErrors(nlumi, x, sig, eBare, eBare, errGrLo, errGrHi); 
  setGraph(gThy1, kGreen-2, 1000); 


  gThy1->SetMaximum(YMAX);
  gThy1->Draw("a3");
  gThy0->Draw("3");
  gBare->Draw("l");

  tl->SetNDC(kTRUE);
  tl->SetTextFont(42);
  tl->DrawLatex(0.2, 0.8, "#sqrt{s} = 14 TeV");

  legg = new TLegend(0.3, 0.2, 0.8, 0.4); 
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 

  legg->AddEntry(gBare, "central expectation", "l"); 
  legg->AddEntry(gThy0, "scale uncertainties", "f"); 
  legg->AddEntry(gThy1, "missing top mass effects", "f"); 
  
  legg->Draw();

  c0->SaveAs(Form("%s/result-thy-lumi.pdf", fDirectory.c_str())); 



}




// ----------------------------------------------------------------------
// this is with the EXPO paramtrization
void hstat::plotResults1() {

  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);
  c0->Clear();

  double sig0 = 3.92765;

  // -- various signficance ranges for index
  //    0: background uncertainty (12 == mode)
  //    1: scale uncertainty (13 == mode)
  //    2: missing top calculations (11 == mode)

  const double sig0Lo[] = {3.60458, 3.41723, 3.66894};
  const double sig0Hi[] = {4.81397, 4.62856, 4.5463};

  double errLo[] = {sig0 - sig0Lo[0], sig0 - sig0Lo[1], sig0 - sig0Lo[2]};
  double errHi[] = {sig0Hi[0] - sig0, sig0Hi[1] - sig0, sig0Hi[2] - sig0};

  double erelLo[] = {errLo[0]/sig0, errLo[1]/sig0,  errLo[2]/sig0};
  double erelHi[] = {errHi[0]/sig0, errHi[1]/sig0,  errHi[2]/sig0};

  // -- Significances: mixed 3 and 4!
  double x[] =   {500.,    1000.,   1500.,   2000.,   2500.,   3000.};
  double sig[] = {2.98989, 3.92765, 4.69809, 5.11973, 5.49504, 5.7};
  double eBare[] = {0., 0., 0., 0., 0., 0.};
  double e0Lo[6], e0Hi[6];
  double e1Lo[6], e1Hi[6];
  double e2Lo[6], e2Hi[6];
  double eALo[6], eAHi[6];

  for (int i = 0; i < 6; ++i) {
    e0Lo[i] = sig[i]*erelLo[0];
    e0Hi[i] = sig[i]*erelHi[0];

    e1Lo[i] = sig[i]*TMath::Sqrt(erelLo[0]*erelLo[0] + erelLo[1]*erelLo[1]);
    e1Hi[i] = sig[i]*TMath::Sqrt(erelHi[0]*erelHi[0] + erelHi[1]*erelHi[1]);

    eALo[i] = sig[i]*TMath::Sqrt(erelLo[0]*erelLo[0] + erelLo[1]*erelLo[1] + erelLo[2]*erelLo[2]);
    eAHi[i] = sig[i]*TMath::Sqrt(erelHi[0]*erelHi[0] + erelHi[1]*erelHi[1] + erelHi[2]*erelHi[2]);
  }
  
  TGraphAsymmErrors *gBare = new TGraphAsymmErrors(6, x, sig, eBare, eBare, eBare, eBare); 
  setGraph(gBare, kBlack, 1000); 

  TGraphAsymmErrors *g0 = new TGraphAsymmErrors(6, x, sig, eBare, eBare, e0Lo, e0Hi); 
  setGraph(g0, kCyan+2, 1000); 

  TGraphAsymmErrors *g1 = new TGraphAsymmErrors(6, x, sig, eBare, eBare, e1Lo, e1Hi); 
  setGraph(g1, kCyan+1, 1000); 

  TGraphAsymmErrors *gA = new TGraphAsymmErrors(6, x, sig, eBare, eBare, eALo, eAHi); 
  setGraph(gA, kYellow-3, 1000); 

  gA->SetMinimum(0.);
  
  //  zone(2,2);
  gA->Draw("a3");
  g1->Draw("3");
  g0->Draw("3");
  gBare->Draw("c");

  TLatex *tl = new TLatex(); 
  tl->SetNDC(kTRUE);
  tl->SetTextFont(42);
  tl->DrawLatex(0.2, 0.8, "#sqrt{s} = 14 TeV");

  TLegend *legg = new TLegend(0.3, 0.2, 0.8, 0.4); 
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 

  legg->AddEntry(gBare, "central expectation", "l"); 
  legg->AddEntry(g0, "background uncertainty", "f"); 
  legg->AddEntry(g1, "scale uncertainty", "f"); 
  legg->AddEntry(gA, "missing top calculations", "f"); 
  
  legg->Draw();

  c0->SaveAs(Form("%s/result-lumi.pdf", fDirectory.c_str())); 

}


