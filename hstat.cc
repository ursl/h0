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

  fHistFile = 0; 

  if (setup == "") {
    fHistFileName = Form("%s/hstat.root", dir.c_str()); 
  } else {
    fHistFileName = Form("%s/hstat-%s.root", dir.c_str(), setup.c_str()); 
  }

  fTexFileName = fHistFileName; 
  replaceAll(fTexFileName, ".root", ".tex"); 
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));


  NBINS = 55; 
  MGGLO = 70.;
  MGGHI = 180.;

  fSg0 =   86; 
  fSg1 =  150;
  fBg  = 1363; 
  fHiggsMpeak = 125.0; 
  fHiggsMres  =   9.4;
  // mass pol1
  fBgMp1   = 0.071276;
  fBgMp1E  = 0.00529807;
  // mass chebychev
  fBgMc0   = 0.159472;
  fBgMc0E  = 0.0119849;
  // pT
  fBgTau   = -0.0135662;
  fBgTauE  =  0.000110543;
  fSg0Tau  = -0.0117118;
  fSg0TauE =  0.000176586;
  fSg1Tau  = -0.00791874;
  fSg1TauE = -9.58177e-05;
   
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
void hstat::run1D(int ntoys, int mode) {

  TH1D *ht0 = new TH1D("ht0", "prof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 400, -40., 40.); ht0->SetLineColor(kBlue); 
  TH1D *ht1 = new TH1D("ht1", "prof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 400, -40., 40.); ht1->SetLineColor(kRed); 

  if (0 == mode) {
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "Input: sg0N = " << fSg0 << endl;
    cout << "Input: sg1N = " << fSg1 << endl;
    cout << "Input: bgN  = " << fBg << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    for (int i = 0; i < ntoys; ++i) {
      cout << "XXXXXX run " << i << endl;
      toy1D();
   
      ht0->Fill(2.*(fD0M1 - fD0M0));    
      ht1->Fill(2.*(fD1M1 - fD1M0));
    }
  }

  zone(); 
  gStyle->SetOptStat(0); 
  ht0->Draw();
  ht1->Draw("same");
  
  double midpoint, tailprob, sigma;
  findMidPoint(ht0, ht1, midpoint, tailprob); 
  sigma = 2.*oneSidedGaussianSigma(tailprob);
  
  TLatex *tl = new TLatex(); 
  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.04);
  tl->DrawLatex(0.6, 0.85, Form("ntoys: %d", ntoys)); 
  tl->DrawLatex(0.6, 0.80, Form("midpoint: %4.3f", midpoint)); 
  tl->DrawLatex(0.6, 0.75, Form("tailprob: %4.3f", tailprob)); 
  tl->DrawLatex(0.6, 0.70, Form("Sep: %4.3f", sigma)); 

  cout << "----------------------------------------------------------------------" << endl;
  cout << "midpoint: " << midpoint << endl;
  cout << "sigma:    " << sigma << endl;
  cout << "----------------------------------------------------------------------" << endl;

  gPad->SaveAs(Form("hstat-%d-%d.pdf", mode, ntoys)); 

}

// ----------------------------------------------------------------------
void hstat::toy1D() {

  model *m0 = genModel1(0, fSg0, fBg); 
  model *m1 = genModel1(1, fSg1, fBg); 

  RooDataSet *d0(0), *d1(0), *bgData(0), *sg0Data(0), *sg1Data(0);

  bgData  = m0->bgPdf->generate(RooArgSet(*m0->m), fBg); 
  sg0Data = m0->sgPdf->generate(RooArgSet(*m0->m), fSg0); 
  sg1Data = m1->sgPdf->generate(RooArgSet(*m1->m), fSg1); 

  // -- model 0 plus background
  d0 = new RooDataSet(*bgData);
  d0->append(*sg0Data);
  d1 = new RooDataSet(*bgData);
  d1->append(*sg1Data);

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

  delete m0; 
  delete m1;
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
  
  // -- define POI
  aModel->poi.add(*aModel->sgN);

  return aModel; 
}


// ----------------------------------------------------------------------
void  hstat::findMidPoint(TH1D* hq0, TH1D* hq1, double &midpoint, double &tailprob) {
  double iq0(0), iq1(0), delta(0.); 
  int nbinsx = hq0->GetNbinsX();
  TH1D *nq0 = (TH1D*)hq0->Clone("nq0");
  nq0->Scale(1./nq0->GetSumOfWeights());
  TH1D *nq1 = (TH1D*)hq1->Clone("nq1");
  nq1->Scale(1./nq1->GetSumOfWeights());
  for (int i = 0; i < nbinsx; ++i) {
    iq0 = nq0->Integral(0, i); 
    iq1 = nq1->Integral(i, nbinsx); 
    if (iq0 > iq1) {
      midpoint = nq0->GetBinCenter(i); 
      if ((iq0 - iq1) > delta) {
	cout << "taking previous, delta = " << delta << " while this one has iq0 = " << iq0 << " and iq1 = " << iq1 << endl;
	midpoint = nq0->GetBinCenter(i-1); 
	iq0 = nq0->Integral(0, i-1); 
	iq1 = nq1->Integral(i-1, nbinsx); 
      }
      break;
    }
    delta = iq1 - iq0; 
  }
    
  cout << "found iq0 = " << iq0 << " iq1 = " << iq1 << " at midpoint = " << midpoint << endl;
  iq0 = nq0->Integral(0, nq0->FindBin(midpoint)+1); 
  iq1 = nq1->Integral(nq0->FindBin(midpoint)+1, nbinsx); 
  cout << " next iq0 = " << iq0 << " iq1 = " << iq1 << endl;
  iq0 = nq0->Integral(0, nq0->FindBin(midpoint)-1); 
  iq1 = nq1->Integral(nq0->FindBin(midpoint)-1, nbinsx); 
  cout << " prev iq0 = " << iq0 << " iq1 = " << iq1 << endl;

  tailprob = 0.5*(iq0+iq1); 
 
}


// ----------------------------------------------------------------------
double hstat::oneSidedGaussianSigma(double prob) {

  RooRealVar x("x", "x", -10, 10) ;
  RooGaussian gx("gx", "gx", x, RooConst(0), RooConst(1));
  
  RooAbsReal* igx(0);   
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
       << " yields significance " << sigma << 0.5*TMath::Erfc(sigma/TMath::Sqrt(2)) << endl;

  //   x.setRange("signal", sigma, 10);
  //   igx = gx.createIntegral(x, NormSet(x), Range("signal")) ;
  //   cout << gx.createIntegral(x, NormSet(x), Range("signal"))->getVal() << endl;
  //  cout << "with formula: " <<  << endl;

  return sigma;

}
