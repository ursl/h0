#include "plotHpt.hh"

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

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooDataSet.h"

#include "RooAbsPdf.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"

#include "RooStats/NumberCountingUtils.h"
#include "RooGlobalFunc.h"
#include "RooStats/RooStatsUtils.h"

#include "dataset.hh"
#include "util.hh"

ClassImp(plotHpt)

using namespace std; 
using namespace RooFit; 
using namespace RooStats; 

// ----------------------------------------------------------------------
plotHpt::plotHpt(string dir,  string files, string setup): plotClass(dir, files, setup), fWorkspace("w") {
  loadFiles(files);

  fNtoy = 200; 

  if (setup == "") {
    fHistFileName = Form("%s/plotHpt.root", dir.c_str()); 
  } else {
    fHistFileName = Form("%s/plotHpt-%s.root", dir.c_str(), setup.c_str()); 
  }

  fTexFileName = fHistFileName; 
  replaceAll(fTexFileName, ".root", ".tex"); 
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));

  fSg0 = fSg1 =  fBg = 
    fNormSg0 = fNormSg0E = 
    fHiggsMpeak =  fHiggsMres = 
    fBgMp1 = fBgMp1E = 
    fBgTau = fBgTauE = 
    fSg0Tau = fSg0TauE = 
    fSg1Tau = fSg1TauE = 
    0.;
  
  NPTBINS = 10; 
  GETA = 2.5;
  G0ISO = G1ISO = 0.2; 
  G0ISOR = 0.005; 
  G1ISOR = 0.010; 

  PTNL = 150.;
  PTNH = 250.;
  PTLO = 400.;
  PTHI = 10000.;

  NBINS = 35;  
  MGGLO = 90.;
  MGGHI = 160.;

  fLumi = 1000.;
  
  G0PT = 150.;
  G1PT =  70.;
  
  G0PTLO = 60.; 
  G1PTLO = 40.;
}


// ----------------------------------------------------------------------
plotHpt::~plotHpt() {

}


// ----------------------------------------------------------------------
void plotHpt::readHistograms() {
  fHists.clear();
  cout << "readHistograms from " << fHistFileName << endl;
  fHistFile = TFile::Open(fHistFileName.c_str()); 

  vector<string> ds;
  ds.push_back("sherpa");
  ds.push_back("mcatnlo5");
  ds.push_back("mcatnlo");
  string dataset("");

  vector<string> hists;
  hists.push_back("pt");
  hists.push_back("m");
  hists.push_back("g0pt");
  hists.push_back("g1pt");

  vector<string> cuts;
  cuts.push_back("lopt");
  cuts.push_back("hipt");
  cuts.push_back("nopt");
  cuts.push_back("goodcand");


  TH1* h1(0);
  char hist[200];
  for (unsigned int ih = 0; ih < hists.size(); ++ih) {
    for (unsigned int ic = 0; ic < cuts.size(); ++ic) {
      for (unsigned int id = 0; id < ds.size(); ++id) {
	sprintf(hist, "%s_%s_%s", hists[ih].c_str(), ds[id].c_str(), cuts[ic].c_str());
	h1 = (TH1*)fHistFile->Get(hist);
	fHists.insert(make_pair(hist, h1));
	dataset = ds[id];
	if (!dataset.compare("mcatnlo")) dataset = "mcatnlo0";
	setHistTitles(fHists[hist], fDS[dataset], h1->GetXaxis()->GetTitle(), "Entries/bin");
      }
    }
  }
  
//   fHistFile->Close();

}

// ----------------------------------------------------------------------
void plotHpt::bookHist(string name, string cuts) {
  char hist[200], thist[200], ahist[200];

  // -- gen pT for efficiency
  sprintf(hist, "%s_%s_%s", "genpt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "genpt");
  sprintf(ahist, "%s", "p_{T}(#gamma#gamma)^{gen} [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  // -- gen pT vs g1pt
  sprintf(hist, "genpt_g1pt_%s_%s", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "pt");
  sprintf(ahist, "%s", "p_T^{gen}(#gamma#gamma) [GeV]"); 
  fHists.insert(make_pair(hist, new TH2D(hist, thist, 100, 0, 1000., 100, 0., 300.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "pt_{#gamma1}^{gen} [GeV]");

  // -- reco overlays
  sprintf(hist, "%s_%s_%s", "mpt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "mpt");
  sprintf(ahist, "%s", "p_T(#gamma#gamma) [GeV]"); 
  fHists.insert(make_pair(hist, new TH2D(hist, thist, 100, 0., 1000., NBINS, MGGLO, MGGHI))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "m_{#gamma #gamma} [GeV]");

  sprintf(hist, "%s_%s_%s", "m", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "m");
  sprintf(ahist, "%s", "m(#gamma#gamma) [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, NBINS, MGGLO, MGGHI))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "pt");
  sprintf(ahist, "%s", "pT(#gamma#gamma) [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  
  sprintf(hist, "%s_%s_%s", "eta", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "eta");
  sprintf(ahist, "%s", "#eta(#gamma#gamma)"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, -5., 5.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  
  sprintf(hist, "%s_%s_%s", "g0pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0pt");
  sprintf(ahist, "%s", "p_{T}(#gamma_{0}) [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 800, 0., 800.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g1pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1pt");
  sprintf(ahist, "%s", "p_T(#gamma_{1}) [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 300, 0., 300.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g0iso", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0iso");
  sprintf(ahist, "%s", "I(#gamma_{0})"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 2.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g1iso", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1iso");
  sprintf(ahist, "%s", "I(#gamma_{1})"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 2.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g0isor", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0isor");
  sprintf(ahist, "%s", "I^{rel}(#gamma_{0})"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 0.006))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g1isor", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1isor");
  sprintf(ahist, "%s", "I^{rel}(#gamma_{1})"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 0.012))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

}


// ----------------------------------------------------------------------
void plotHpt::makeAll(int bitmask) {
  if (bitmask & 0x1) {
    treeAnalysis();
    allNumbers(); 
  }
  if (bitmask & 0x2) {
    optimizeCuts(); 
  }
  if (bitmask & 0x4) {
    optAnalysis(1); 
    optAnalysis(2); 
    optAnalysis(3); 
  }
  if (bitmask & 0x10) {
    validation(); 
  }

}


// ----------------------------------------------------------------------
void plotHpt::bgShape(int nevts) {
  string ds("sherpa");
  fCds = ds; 

  char hist[200], thist[200];
  sprintf(hist, "bgshape_default_m");
  sprintf(thist, "bgshape default m");
  fHists.insert(make_pair(hist, new TH1D(hist, thist, NBINS, MGGLO, MGGHI))); 
  fHists[hist]->Sumw2();
  sprintf(hist, "bgshape_default_pt");
  sprintf(thist, "bgshape default pt");
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
  fHists[hist]->Sumw2();

  for (int i = 0; i < NPTBINS; ++i) {
    sprintf(hist, "%s_ptbin%d_m", "bgshape", i);
    sprintf(thist, "ptbin %d m", i);
    fHists.insert(make_pair(hist, new TH1D(hist, thist, NBINS, MGGLO, MGGHI))); 
    fHists[hist]->Sumw2();
    sprintf(hist, "%s_ptbin%d_pt", "bgshape", i);
    sprintf(thist, "ptbin %d pt", i);
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
    fHists[hist]->Sumw2();
  }
  
  
  TTree *t = getTree(ds); 
  setupTree(t); 
  loopOverTree(t, 3, nevts); 

  TH1D *h1 = new TH1D("h1", "", NPTBINS+1, 0., NPTBINS+1); 
  
  zone(1,2); 
  fHists[Form("bgshape_default_m")]->Fit("pol1", "wl");
  TF1 *f1 = fHists[Form("bgshape_default_m")]->GetFunction("pol1");
  double defx = f1->GetParameter(1)/fHists[Form("bgshape_default_m")]->Integral();
  double defxe = f1->GetParError(1)/fHists[Form("bgshape_default_m")]->Integral();
  h1->SetBinContent(1, f1->GetParameter(1)/fHists[Form("bgshape_default_m")]->Integral()); 
  h1->SetBinError(1, f1->GetParError(1)/fHists[Form("bgshape_default_m")]->Integral()); 
  h1->GetXaxis()->SetBinLabel(1, "default");
  c0->cd(2); 
  fHists[Form("bgshape_default_pt")]->Draw(); 
  c0->SaveAs(Form("%s/bgshape-default.pdf", fDirectory.c_str())); 
  
  
  tl->SetNDC(kTRUE); 
  TH1D *h0(0); 
  for (int i = 0; i < 10; ++i) {
    zone(1,3);
    h0 = (TH1D*)fHists[Form("bgshape_ptbin%d_m", i)]->Clone("h0"); 
    h0->Scale(1./h0->Integral()); 
    fHists[Form("bgshape_ptbin%d_m", i)]->Fit("pol1", "l"); 
    f1 = fHists[Form("bgshape_ptbin%d_m", i)]->GetFunction("pol1");
    h1->SetBinContent(i+2, f1->GetParameter(1)/fHists[Form("bgshape_ptbin%d_m", i)]->Integral()); 
    h1->SetBinError(i+2, f1->GetParError(1)/fHists[Form("bgshape_ptbin%d_m", i)]->Integral()); 
    h1->GetXaxis()->SetBinLabel(i+2, Form("pt bin %d", i));
    tl->DrawLatex(0.25, 0.8, Form("p1 = %8.7f +/- %8.7f", 
				  f1->GetParameter(1)/fHists[Form("bgshape_ptbin%d_m", i)]->Integral(),  
				  f1->GetParError(1)/fHists[Form("bgshape_ptbin%d_m", i)]->Integral())
		  ); 
    
    c0->cd(2); 
    h0->Fit("pol1"); 
    tl->DrawLatex(0.25, 0.8, Form("p1 = %8.7f +/- %8.7f", 
				  h0->GetFunction("pol1")->GetParameter(1), 
				  h0->GetFunction("pol1")->GetParError(1))
		  ); 

    c0->cd(3);
    fHists[Form("bgshape_ptbin%d_pt", i)]->Draw(); 
    c0->SaveAs(Form("%s/bgshape-ptbin%d.pdf", fDirectory.c_str(), i)); 
  }

  zone(); 
  h1->Fit("pol1", "r", "", h1->GetBinLowEdge(4), h1->GetBinLowEdge(h1->GetNbinsX()+1)); 
  pl->SetLineColor(kBlue);
  pl->SetLineStyle(kSolid);
  pl->DrawLine(0., defx, NPTBINS+1, defx); 
  pl->SetLineStyle(kDashed);
  pl->DrawLine(0., defx+defxe, NPTBINS+1, defx+defxe); 
  pl->SetLineStyle(kDashed);
  pl->DrawLine(0., defx-defxe, NPTBINS+1, defx-defxe); 
  
  c0->SaveAs(Form("%s/bgshape-summary.pdf", fDirectory.c_str())); 

}


// ----------------------------------------------------------------------
void plotHpt::treeAnalysis(int nevts) {

  cout << "treeAnalysis: open " << fHistFileName << endl;
  cout << " cuts: " << endl;
  cout << "  PT > " << PTLO << endl;
  cout << "  G0PT > " << G0PT << endl;
  cout << "  G1PT > " << G1PT << endl;
  
  TFile *f = TFile::Open(fHistFileName.c_str(), "RECREATE"); 
  
  string ds("sherpa");
  fCds = ds; 
  bookHist(ds, "nopt"); 
  bookHist(ds, "goodcand"); 
  bookHist(ds, "lopt"); 
  bookHist(ds, "hipt"); 
  TTree *t = getTree(ds); 
  setupTree(t); 
  loopOverTree(t, 1, nevts); 


  ds = "mcatnlo5"; 
  fCds = ds; 
  bookHist(ds, "nopt"); 
  bookHist(ds, "goodcand"); 
  bookHist(ds, "lopt"); 
  bookHist(ds, "hipt"); 
  t = getTree(ds); 
  setupTree(t);
  loopOverTree(t, 1, nevts); 


  ds = "mcatnlo0"; 
  fCds = ds; 
  bookHist(ds, "nopt"); 
  bookHist(ds, "goodcand"); 
  bookHist(ds, "lopt"); 
  bookHist(ds, "hipt"); 
  t = getTree(ds); 
  setupTree(t);
  loopOverTree(t, 1, nevts); 
  TH1D *hmcatnlo0 = new TH1D("hmcatnlo0", "hmcatnlo0", 100, 0., 1000.); 
  t->Draw("gpt>>hmcatnlo0"); 

  ds = "mcatnlo1"; 
  fCds = ds; 
  bookHist(ds, "nopt"); 
  bookHist(ds, "goodcand"); 
  bookHist(ds, "lopt"); 
  bookHist(ds, "hipt"); 
  t = getTree(ds); 
  setupTree(t);
  loopOverTree(t, 1, nevts); 
  TH1D *hmcatnlo1 = new TH1D("hmcatnlo1", "hmcatnlo1", 100, 0., 1000.); 
  t->Draw("gpt>>hmcatnlo1"); 

  map<string, TH1*>::iterator hit = fHists.begin();
  map<string, TH1*>::iterator hite = fHists.end();
  string h0name(""), h1name("");
  TH1 *h1(0), *h2(0); 
  map<string, TH1*> tlist;
  for (; hit != hite; ++hit) {
    h0name = hit->second->GetName(); 
    //    cout << h0name << endl;
    if (string::npos == h0name.find("mcatnlo0")) continue;
    h1name = h0name;
    replaceAll(h1name, "mcatnlo0", "mcatnlo1"); 
    h1 = fHists[h1name];
    h2 = combMCAtNLOHist(hit->second, h1); 
    tlist.insert(make_pair(h2->GetName(), h2)); 
  }
  fHists.insert(tlist.begin(), tlist.end()); 


  f->Write();
  f->Close();


}


// ----------------------------------------------------------------------
void plotHpt::massResolution(string hname) {
  fPtBins.clear();

  int NBINS(10);
  double aptbins[] = {0., 100., 150., 200., 250., 300., 350., 400., 600., 800., 1000.};

  fHistFile = TFile::Open(fHistFileName.c_str()); 
  TH2D *hmpt = (TH2D*)fHistFile->Get(hname.c_str());
  double xmax(hmpt->GetXaxis()->GetXmax());
  TH1D *hmpeak = new TH1D("hmpeak", "", 10*NBINS, 0., xmax);   hmpeak->Sumw2();
  TH1D *hmres = new TH1D("hmres", "", 10*NBINS, 0., xmax);     hmres->Sumw2();
  setTitles(hmpeak, "pT_{#gamma#gamma} [GeV]", "peak value and resolution [GeV]",  0.05, 1.1, 1.3); 
  int totbins(hmpt->GetNbinsX());
  TH1D *hy(0);
  
  int ibin0(0), ibin1(0);
  c0->Clear();
  zone(3,4);
  int ipad(1);
  for (int i = 0; i < NBINS; ++i) {
    c0->cd(ipad);
    //     hy = hmpt->ProjectionY(Form("ptbin%d", i), 
    // 			   i*totbins/NBINS+1, 
    // 			   (i+1)*totbins/NBINS+1
    // 			   ); 
    ibin0 = hmpt->GetXaxis()->FindBin(aptbins[i]);
    ibin1 = hmpt->GetXaxis()->FindBin(aptbins[i+1]);
    hy = hmpt->ProjectionY(Form("ptbin%d", i), ibin0, ibin1);
    cout << "x bins: " << ibin0 << " .. " << ibin1 << endl;
    hy->Fit("gaus");    
    ptbins a; 
    a.peak  = hy->GetFunction("gaus")->GetParameter(1);
    a.sigma = hy->GetFunction("gaus")->GetParameter(2);
    a.ptlo  = hmpt->GetXaxis()->GetBinLowEdge(ibin0);
    a.pthi  = hmpt->GetXaxis()->GetBinLowEdge(ibin1);
    a.pt    = 0.5*(a.ptlo + a.pthi);
    ibin0 = hmpeak->FindBin(a.pt);
    hmpeak->SetBinContent(ibin0, a.peak); 
    hmpeak->SetBinError(ibin0, a.sigma); 
    hmres->SetBinContent(ibin0, a.sigma); 
    hmres->SetBinError(ibin0, hy->GetFunction("gaus")->GetParError(2)); 
    cout << "========> " << a.peak << " " << a.sigma 
	 << " " << a.pt << ", " << a.ptlo << " .. " << a.pthi << endl;
    fPtBins.push_back(a);
    ++ipad;
  }
  c0->SaveAs(Form("%s/massResolution-binfits.pdf", fDirectory.c_str())); 
  
  c0->Clear();
  gStyle->SetOptStat(0);
  hmpeak->SetMinimum(100.);
  hmpeak->Draw();
  pl->SetLineColor(kRed);
  pl->SetLineWidth(2);
  pl->DrawLine(0., 125., xmax, 125.);
  c0->SaveAs(Form("%s/massResolution.pdf", fDirectory.c_str()));


  c0->Clear();
  hmres->Fit("pol1");
  c0->SaveAs(Form("%s/massResolution-fit.pdf", fDirectory.c_str()));

  fHistFile->Close();

}

// ----------------------------------------------------------------------
void plotHpt::optimizeCuts(string fname, double lumi, int nevts) {

  // -- output file
  TFile *f(0); 
  cout << "open file " << Form("%s/%s", 
			       fDirectory.c_str(), fname.c_str()) << " RECREATE" << endl;
  f = TFile::Open(Form("%s/%s", 
		       fDirectory.c_str(), fname.c_str()), "RECREATE");     

  // -- setup selection points
  static const double ptloArr[] = {150., 160., 170., 180., 190., 200, 220., 240., 250., 260., 280., 300., 350., 400.
  };
  
  vector<double> ptloCuts(ptloArr, ptloArr + sizeof(ptloArr)/sizeof(ptloArr[0]));
  
  static const double ptg0Arr[] = {90., 100., 110., 120., 125., 130., 140., 150., 160., 170., 180., 200.
  };
  vector<double> ptg0Cuts(ptg0Arr, ptg0Arr + sizeof(ptg0Arr)/sizeof(ptg0Arr[0]));
  
  static const double ptg1Arr[] = {40., 45., 50., 55., 60., 65., 70., 80.
  };
  vector<double> ptg1Cuts(ptg1Arr, ptg1Arr + sizeof(ptg1Arr)/sizeof(ptg1Arr[0]));
  
  for (unsigned int i = 0; i < ptloCuts.size(); ++i) {
    for (unsigned int j = 0; j < ptg0Cuts.size(); ++j) {
      for (unsigned int k = 0; k < ptg1Cuts.size(); ++k) {
	selpoint s; 
	s.fLargerThan.push_back(make_pair(&fb.pt, ptloCuts[i])); 
	s.fLargerThan.push_back(make_pair(&fb.g0pt, ptg0Cuts[j])); 
	s.fLargerThan.push_back(make_pair(&fb.g1pt, ptg1Cuts[k])); 
	fSelPoints.push_back(s); 
      }
    }
  }
  
  cout << "optimizing " << fSelPoints.size() << " selection points " 
       << endl;
  
  // -- loop over signal
  fCds = "mcatnlo5"; 
  fOptMode = 0; 
  TTree *ts = getTree(fCds); 
  double sgScale = lumi/fDS[fCds]->fLumi;
  setupTree(ts); 
  loopOverTree(ts, 2, nevts); 

  // -- loop over background
  fCds = "sherpa"; 
  fOptMode = 1; 
  double bg1Scale = lumi/fDS[fCds]->fLumi;
  TTree *tb1 = getTree(fCds); 
  setupTree(tb1); 
  loopOverTree(tb1, 2, nevts); 

  // -- loop over background from SM Higgs production
  fCds = "mcatnlo0"; 
  fOptMode = 2; 
  double bg2Scale = lumi/fDS[fCds]->fLumi;
  TTree *tb2 = getTree(fCds); 
  setupTree(tb2); 
  loopOverTree(tb2, 2, nevts); 


  f->cd();
  TTree *t = new TTree("opt", "opt");
  double s, b, b0, b1, rs, rb0, rb1; 
  double ptlo, ptg0, ptg1;
  double ssb, sb;
  t->Branch("s",    &s,     "s/D");
  t->Branch("b",    &b,     "b/D");
  t->Branch("b0",   &b0,    "b0/D");
  t->Branch("b1",   &b1,    "b1/D");

  t->Branch("rs",   &rs,    "rs/D");
  t->Branch("rb0",  &rb0,   "rb0/D");
  t->Branch("rb1",  &rb1,   "rb1/D");

  t->Branch("ssb",  &ssb,   "ssb/D");
  t->Branch("sb",   &sb,    "sb/D");

  t->Branch("ptlo",  &ptlo,    "ptlo/D");
  t->Branch("ptg0",  &ptg0,    "ptg0/D");
  t->Branch("ptg1",  &ptg1,    "ptg1/D");
  
  for (unsigned int i = 0; i < fSelPoints.size(); ++i) {
    rs  = fSelPoints[i].fCnt[0];
    rb0 = fSelPoints[i].fCnt[1];
    rb1 = fSelPoints[i].fCnt[2];
    s   = sgScale * fSelPoints[i].fCnt[0];
    b0  = bg1Scale*fSelPoints[i].fCnt[1]; 
    b1  = bg2Scale*fSelPoints[i].fCnt[2]; 
    b   = b0 + b1; 
    ssb = (s+b>0? s/TMath::Sqrt(s+b):0.);
    sb  = (b>0?   s/TMath::Sqrt(b)  :0.);
    ptlo = fSelPoints[i].fLargerThan[0].second;
    ptg0 = fSelPoints[i].fLargerThan[1].second;
    ptg1 = fSelPoints[i].fLargerThan[2].second;
    t->Fill();
  }

  t->Write();
  f->Close();
}


// ----------------------------------------------------------------------
// if (1 == mode) fom = ssb;
// if (2 == mode) fom = (s-b1)/TMath::Sqrt(s+0.3*b);
// if (3 == mode) fom = (s-b1)/TMath::Sqrt(s+b);
void plotHpt::optAnalysis(int mode, string filename, string treename) {
  TFile *f(0); 
  if (filename.compare("")) {
    f = TFile::Open(Form("%s/%s", fDirectory.c_str(), filename.c_str())); 
  }
  TTree *t = (TTree*)f->Get(treename.c_str()); 
  if (0 == t) {
    cout << "no tree with name " << treename << " found in file " << filename << endl;
    return;
  }

  // -- setup tree
  double ssb, sb, s, b, b1; 
  double ptlo, ptg0, ptg1; 
  t->SetBranchAddress("ssb",   &ssb);
  t->SetBranchAddress("s",     &s);
  t->SetBranchAddress("b",     &b);
  t->SetBranchAddress("b1",    &b1);
  t->SetBranchAddress("ptlo",  &ptlo);
  t->SetBranchAddress("ptg0",  &ptg0);
  t->SetBranchAddress("ptg1",  &ptg1);

  TH1D *hptlo = new TH1D("hptlo", "", 400, 100., 500);
  TH1D *hptg0 = new TH1D("hptg0", "", 200, 50., 250);
  TH1D *hptg1 = new TH1D("hptg1", "", 100, 0., 100);
  map<string, TH2D*> hg0g1; 
  TH2D *h2(0); 
  t->Draw("ptlo>>hptlo"); 
  float ptl(0.); 
  for (int i = 1; i < hptlo->GetNbinsX(); ++i) {
    if (hptlo->GetBinContent(i) > 0) {
      ptl = hptlo->GetBinLowEdge(i);
      h2 = new TH2D(Form("ptlo_%3.0f", ptl), Form("ptlo_%3.0f", ptl), 
		    40, 50., 250., 20, 0., 100.); 
      setTitles(h2, "ptg0", "ptg1");
      hg0g1.insert(make_pair(Form("ptlo_%3.0f", ptl), h2)); 
    }
  }
  hptlo->Reset();
  
  // -- loop over tree
  int ibin(0);
  int nentries = Int_t(t->GetEntries());
  double fom(0.); 
  for (int jentry = 0; jentry < nentries; jentry++) {
    t->GetEntry(jentry);
    if (1 == mode) fom = ssb;
    if (2 == mode) fom = (s-b1)/TMath::Sqrt(s+0.3*b);
    if (3 == mode) fom = (s-b1)/TMath::Sqrt(s+b);
    if (fom > hptlo->GetBinContent(hptlo->FindBin(ptlo))) {
      hptlo->SetBinContent(hptlo->FindBin(ptlo), fom); 
    }

    if (fom > hptg0->GetBinContent(hptg0->FindBin(ptg0))) {
      hptg0->SetBinContent(hptg0->FindBin(ptg0), fom); 
    }

    if (fom > hptg1->GetBinContent(hptg1->FindBin(ptg1))) {
      hptg1->SetBinContent(hptg1->FindBin(ptg1), fom); 
    }

    // -- fill 2D histogram specifically for each ptlo
    h2 = hg0g1[Form("ptlo_%3.0f", ptlo)]; 
    if (0 == h2) {
      cout << "histogram " << Form("ptlo_%3.0f", ptlo)  << " not found" << endl;
    }
    ibin = h2->FindBin(ptg0, ptg1); 
    if (fom > h2->GetBinContent(ibin)) {
      h2->SetBinContent(ibin, fom); 
    }

  }

  zone(2,2);
  hptlo->Draw();
  c0->cd(2);
  hptg0->Draw();
  c0->cd(3);
  hptg1->Draw();
  c0->SaveAs(Form("%s/opt-%d.pdf", fDirectory.c_str(), mode)); 
    
  map<string, TH2D*>::iterator hit  = hg0g1.begin();   
  map<string, TH2D*>::iterator hite = hg0g1.end();   

  zone(4,4);
  gStyle->SetOptStat(0); 
  int ix, iy, iz; 
  double ax, ay;
  int ipad(1); 
  for (; hit != hite; ++hit) {
    c0->cd(ipad);
    shrinkPad(0.2, 0.1, 0.2);
    cout << hit->second->GetName() << endl;
    ibin = hit->second->GetMaximumBin(); 
    hit->second->GetBinXYZ(ibin, ix, iy, iz); 
    ax = hit->second->GetXaxis()->GetBinLowEdge(ix); 
    ay = hit->second->GetYaxis()->GetBinLowEdge(iy); 
    hit->second->SetTitle(Form("%s, max = %4.3f at ptg0 > %4.3f, ptg1 > %4.3f", 
			       hit->second->GetTitle(), 
			       hit->second->GetMaximum(), 
			       ax, 
			       ay));
    hit->second->SetMinimum(hit->second->GetMinimum(1.0));
    hit->second->Draw("colz");
    ++ipad;
  }
  c0->SaveAs(Form("%s/opt-ptranges-%d.pdf", fDirectory.c_str(), mode));

}




// ----------------------------------------------------------------------
void plotHpt::toy1(int nsg, int nbg) {

  fHistFile = TFile::Open(fHistFileName.c_str()); 

  // =======
  // -- mass
  // =======

  // -- Signal
  RooRealVar m("m", "m", 100., 150.); 
  RooRealVar meanM("meanM", "meanM", 125., 100, 150);  
  RooRealVar sigmaM("sigmaM", "sigmaM", 5., 3., 10.);  
  RooGaussian sgM("sgM", "signal mass", m, meanM, sigmaM); 

  // -- Background
  RooRealVar C0("C0", "coefficient #0", 1.0, -1., 1.); 
  RooRealVar C1("C1", "coefficient #1", 0.1, -1., 1.); 
  RooPolynomial bgM("bgM", "background mass", m, RooArgList(C0, C1)); 


  // =====
  // -- pT
  // =====

  RooRealVar pt("pt", "pt", 0., 1000); 
  RooRealVar tau("tau","tau", -0.01412, -10, -1.e-4); 
  RooExponential sgPt("sgPt", "signal pT", pt, tau);   

  RooRealVar bgTau("bgTau", "background tau", -0.019, -10, -1.e-4); 
  RooExponential bgPt("bgPt", "background pT", pt, bgTau);   

  TH1D *h1 = (TH1D*)fHistFile->Get("pt_mcatnlo_hipt")->Clone("h1");
  //  h1->Smooth();
  RooDataHist hSgPt("hSgPt", "signal histogram pT", pt, h1); 
  RooHistPdf sgHistPt("sgHistPt","signal histogam pT", pt, hSgPt, 2) ;

  h1 = (TH1D*)fHistFile->Get("pt_sherpa_hipt")->Clone("h1");
  //  h1->Smooth();
  RooDataHist hBgPt("hBgPt", "background histogram pT", pt, h1); 
  RooHistPdf bgHistPt("bgHistPt","background histogam pT", pt, hBgPt, 2) ;

  // -- determine Chebychev coefficients
  if (0) {

    // -- OLD Background
    RooRealVar C0("C0", "coefficient #0", 0.086, -1., 1.); 
    RooRealVar C1("C1", "coefficient #1", 0.006, -1., 1.); 
    RooRealVar C2("C2", "coefficient #2", 0.008, -1., 1.); 
    RooChebychev bgM("bgM", "background mass", m, RooArgList(C0, C1, C2)); 

    TH1D *h1 = (TH1D*)fHistFile->Get("m_sherpa_hipt");
    RooDataHist dataBg("data", "background mass distribution", m, h1);
    bgM.fitTo(dataBg); 
    
    RooPlot *frameBg = m.frame(); 
    dataBg.plotOn(frameBg);  
    bgM.plotOn(frameBg); 
    
    bgM.paramOn(frameBg);
    frameBg->Draw();
  }

  RooRealVar fsig("fsig", "signal fraction", 0.25, 0., 1.) ;

  RooAddPdf modelM("modelM", "model for mass", RooArgList(sgM, bgM), fsig); 
  RooAddPdf modelPt("modelPt", "model for pT", RooArgList(sgHistPt, bgHistPt), fsig) ;
  RooProdPdf model("model", "complete model", RooArgSet(modelM, modelPt)); 
  

  RooDataSet *data = model.generate(RooArgSet(m, pt), nsg+nbg);  

  sigmaM.setConstant(kTRUE);
  meanM.setConstant(kTRUE);

  C0.setRange(0.9, 1.1);
   C1.setRange(0.0, 0.2);

  model.fitTo(*data); 
  RooPlot *frameM = m.frame(); 
  data->plotOn(frameM);  
  model.plotOn(frameM); 
  model.plotOn(frameM, Components("bgM"), LineStyle(kDashed)) ;
  
  RooPlot *framePt = pt.frame();  
  data->plotOn(framePt);  
  model.plotOn(framePt); 
  model.plotOn(framePt, Components("bgHistPt"), LineStyle(kDashed)) ;


  zone(1, 2);
  frameM->Draw();
  c0->cd(2);
  gPad->SetLogy(1);
  framePt->Draw();

  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotHpt::toy2(int nsg, int nbg) {

  fHistFile = TFile::Open(fHistFileName.c_str()); 

  double ggSigma(5.0); 

  // =======
  // -- mass
  // =======

  RooRealVar nSg("nSg", "signal fraction", nsg, 0., 100*nsg);
  RooRealVar nBg("nBg", "background fraction", nbg, 0., 100.*nbg);

  // -- Signal
  RooRealVar m("m", "m", 100., 150.); 
  RooRealVar meanM("meanM", "meanM", 125., 100, 150);  
  RooRealVar sigmaM("sigmaM", "sigmaM", ggSigma, 3., 10.);  
  RooGaussian sgM("sgM", "signal mass", m, meanM, sigmaM); 
  RooExtendPdf esgM("esgM", "extended signal mass", sgM, nSg);

  // -- Background
  RooRealVar C0("C0", "coefficient #0", 1.0, -1., 2.); 
  RooRealVar C1("C1", "coefficient #1", 0.1, -1., 2.); 
  RooPolynomial bgM("bgM", "background mass", m, RooArgList(C0, C1)); 
  RooExtendPdf ebgM("ebgM", "extended background mass", bgM, nBg);
  

  // =====
  // -- pT
  // =====

  RooRealVar pt("pt", "pt", 0., 1000); 
  RooRealVar tau("tau","tau", -0.01412, -10, -1.e-4); 
  RooExponential sgPt("sgPt", "signal pT", pt, tau);   
  RooExtendPdf esgPt("esgPt", "extended signal pT", sgPt, nSg); 

  RooRealVar bgTau("bgTau", "background tau", -0.019, -10, 0.); 
  RooExponential bgPt("bgPt", "background pT", pt, bgTau);   

  //  TH1D *h1 = (TH1D*)fHistFile->Get("pt_mcatnlo_hipt")->Clone("h1");
  TH1D *h1 = (TH1D*)fHistFile->Get("pt_mcatnlo5_hipt")->Clone("h1");
  RooDataHist hSgPt("hSgPt", "signal histogram pT", pt, h1); 
  RooHistPdf sgHistPt("sgHistPt","signal histogam pT", pt, hSgPt, 2) ;
  //  RooExtendPdf esgHistPt("esgHistPt", "extended signal histogram pT", sgHistPt, nSg); 

  TH1D *h2 = (TH1D*)fHistFile->Get("pt_sherpa_hipt")->Clone("h2");
  RooDataHist hBgPt("hBgPt", "background histogram pT", pt, h2); 
  RooHistPdf bgHistPt("bgHistPt","background histogam pT", pt, hBgPt, 2) ;
  //  RooExtendPdf ebgHistPt("ebgHistPt", "extended background histogram pT", bgHistPt, nBg); 

  RooAddPdf modelM("modelM", "model for mass", RooArgList(sgM, bgM), RooArgList(nSg, nBg)); 
  //  RooAddPdf modelPt("modelPt", "model for pT", RooArgList(sgHistPt, bgHistPt), RooArgList(nSg, nBg));
  RooAddPdf modelPt("modelPt", "model for pT", RooArgList(sgPt, bgPt), RooArgList(nSg, nBg));
  RooProdPdf model("model", "complete model", RooArgSet(modelM, modelPt)); 
  
  //  RooDataSet *data = model.generate(RooArgSet(m, pt), nsg+nbg);  
  RooDataSet *data = model.generate(RooArgSet(m, pt));  

  sigmaM.setConstant(kTRUE);
  meanM.setConstant(kTRUE);

  C0.setRange(0.9, 1.1);
  C1.setRange(0.0, 0.2);

  model.fitTo(*data, Extended(kTRUE)); 
  RooPlot *frameM = m.frame(); 
  data->plotOn(frameM);  
  model.plotOn(frameM); 
  model.plotOn(frameM, Components("bgM"), LineStyle(kDashed)) ;
  
  RooPlot *framePt = pt.frame();  
  data->plotOn(framePt);  
  model.plotOn(framePt); 
  model.plotOn(framePt, Components("bgHistPt"), LineStyle(kDashed)) ;


  zone(2, 2);
  frameM->Draw();
  c0->cd(2);
  gPad->SetLogy(1);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);
  h1->Draw(); 
  h2->Draw("same"); 
  c0->cd(3);
  gPad->SetLogy(0);
  framePt->Draw();
  c0->cd(4);
  gPad->SetLogy(1);
  framePt->Draw();

  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotHpt::toy3(int nsg0, double nsg1, double  nbg) {

  fHistFile = TFile::Open(fHistFileName.c_str()); 

  // -- determine function parameters
  if (0) { 
    double xmin(300.), xmax(1000.); 
    gStyle->SetOptStat(0); 
    gStyle->SetOptFit(111); 
    zone(2,2);
    TH1D *h1 = (TH1D*)fHistFile->Get("m_sherpa_hipt")->Clone("h1");
    //    h1->Scale(1./h1->GetSumOfWeights()); 
    h1->SetMinimum(0.);
    h1->Fit("pol1"); 
    h1->GetFunction("pol1")->SetLineWidth(2); 

    c0->cd(2);
    gPad->SetLogy(1);
    TH1D *h2 = (TH1D*)fHistFile->Get("pt_sherpa_hipt")->Clone("h2");
    //    h2->Scale(1./h2->GetSumOfWeights()); 
    h2->SetTitle("Sherpa gg"); 
    h2->Fit("expo", "rl", "", xmin, xmax); 
    h2->GetFunction("expo")->SetLineWidth(2); 

    c0->cd(3);
    gPad->SetLogy(1);
    TH1D *h3 = (TH1D*)fHistFile->Get("pt_mcatnlo5_hipt")->Clone("h3");
    //    h3->Scale(1./h3->GetSumOfWeights()); 
    h3->SetTitle("m(top) #rightarrow #infty"); 
    h3->Fit("expo", "rl", "", xmin, xmax); 
    h3->GetFunction("expo")->SetLineWidth(2); 

    c0->cd(4);
    gPad->SetLogy(1);
    TH1D *h4 = (TH1D*)fHistFile->Get("pt_mcatnlo_hipt")->Clone("h4");
    //    h4->Scale(1./h4->GetSumOfWeights()); 
    h4->SetTitle("m(top) = 173.5 GeV"); 
    h4->Fit("expo", "rl", "", xmin, xmax); 
    h4->GetFunction("expo")->SetLineWidth(2); 

    //  ****************************************
    //  Minimizer is Linear
    //  Chi2                      =      14.1314
    //  NDf                       =           23
    //  p0                        =      166.971   +/-   28.4753     
    //  p1                        =     0.837631   +/-   0.227452    

    //  NORMALIZED TO 1:
    //  ****************************************
    //  Minimizer is Linear
    //  Chi2                      =   0.00207631
    //  NDf                       =           23
    //  p0                        =    0.0245329   +/-   0.345162    
    //  p1                        =  0.000123072   +/-   0.00275705  
   
    //  SHERPA
    //  FCN=80.7532 FROM MIGRAD    STATUS=CONVERGED      83 CALLS          84 TOTAL
    //                      EDM=1.57571e-10    STRATEGY= 1      ERROR MATRIX ACCURATE
    //   EXT PARAMETER                                   STEP         FIRST
    //   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
    //    1  Constant     1.05700e+01   6.00323e-02   7.54946e-05   8.85296e-04
    //    2  Slope       -1.26957e-02   1.55210e-04   1.95185e-07   3.77248e-01
    //                                ERR DEF= 0.5
    //  m(top) -> infty
    //  FCN=44.2806 FROM MIGRAD    STATUS=CONVERGED      42 CALLS          43 TOTAL
    //                      EDM=8.07575e-08    STRATEGY= 1      ERROR MATRIX ACCURATE
    //   EXT PARAMETER                                   STEP         FIRST
    //   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
    //    1  Constant     8.84227e+00   4.01906e-02   5.02440e-05   1.74232e-02
    //    2  Slope       -7.81368e-03   9.10171e-05   1.13778e-07   9.73403e+00
    //                                ERR DEF= 0.5
    // m(top) = 173.5 GeV
    //  FCN=35.0345 FROM MIGRAD    STATUS=CONVERGED      70 CALLS          71 TOTAL
    //                      EDM=2.92507e-12    STRATEGY= 1      ERROR MATRIX ACCURATE
    //   EXT PARAMETER                                   STEP         FIRST
    //   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
    //    1  Constant     9.45657e+00   6.72258e-02   6.20582e-05   1.01431e-04
    //    2  Slope       -1.09325e-02   1.67445e-04   1.54573e-07   4.42596e-02
    //                                ERR DEF= 0.5
    return;
  }




  // mass
  // ====
  
  // -- Signal
  RooRealVar m("m", "m", 100., 150.); 
  RooRealVar sgP("sgP", "signal peak mass", 125., 100, 150);  
  sgP.setConstant(kTRUE);
  RooRealVar sgS("sgS", "signal sigma mass", 10., 5., 15.);  
  sgS.setConstant(kTRUE);
  RooGaussian sg0M("sg0M", "signal mass", m, sgP, sgS); 
  RooGaussian sg1M("sg1M", "signal mass", m, sgP, sgS); 

  // -- Background
  double bgc0  = 0.0245329;
  double bgc0e = 0.00418385;  // derived from the normalized unscaled fit! 
  RooRealVar C0("C0", "coefficient #0", bgc0, -1., 1.); 
  //  C0.setRange(bgc0 - bgc0e, bgc0 + bgc0e);
  //   double bgc1  = 0.0001231;
  //   double bgc1e = 3.341922e-5; // derived from the normalized unscaled fit! 
  //   RooRealVar C1("C1", "coefficient #1", bgc1, -1., 1.); 
  //  C1.setRange(bgc1 - bgc1e, bgc1 + bgc1e);
  RooPolynomial bg0M("bg0M", "background 0 gamma gamma mass", m, RooArgList(C0)); 
  RooPolynomial bg1M("bg1M", "background 1 gamma gamma mass", m, RooArgList(C0)); 


  // pT
  // ==
  RooRealVar pt("pt", "pt", 300., 1000.); 
  // -- SM Higgs
  double sg0tau  = -1.09325e-02; 
  double sg0taue = 1.67445e-04;
  RooRealVar sg0Tau("sg0Tau", "signal 0 tau", sg0tau, -10., 10.); 
  sg0Tau.setRange(sg0tau - sg0taue, sg0tau + sg0taue); 
  RooExponential sg0Pt("sg0Pt", "signal 0 pT", pt, sg0Tau);   

  // -- "contact" Higgs
  double sg1tau  = -7.81368e-03; 
  double sg1taue = 9.10171e-05;
  RooRealVar sg1Tau("sg1Tau", "signal 1 tau", sg1tau, -10., 10.); 
  sg1Tau.setRange(sg1tau - sg1taue, sg1tau + sg1taue); 
  RooExponential sg1Pt("sg1Pt", "signal 1 pT", pt, sg1Tau);   


  double bgtau  = -1.26957e-02; 
  double bgtaue = 1.55210e-04; 
  RooRealVar bgTau("bgTau", "background gamma gamma tau", bgtau, -10., 10.); 
  //  bg0Tau.setRange(bg0tau - bg0taue, bg0tau + bg0taue); 
  RooExponential bg0Pt("bg0Pt", "background 0 gamma gamma pT", pt, bgTau);   
  RooExponential bg1Pt("bg1Pt", "background 1 gamma gamma pT", pt, bgTau);   


  RooRealVar sg0frac("sg0frac","fraction of Higgs 0 signal", nsg0/nbg, 0., 1.);
  RooRealVar sg1frac("sg1frac","fraction of Higgs 1 signal", nsg1/nbg, 0., 1.);
  RooAddPdf model0M("model0M", "model 0 for mass", RooArgList(sg0M, bg0M), sg0frac);
  RooAddPdf model0Pt("model0Pt", "model 0 for pT", RooArgList(sg0Pt, bg0Pt), sg0frac);

  RooAddPdf model1M("model1M", "model 1 for mass", RooArgList(sg1M, bg1M), sg1frac);
  RooAddPdf model1Pt("model1Pt", "model 1 for pT", RooArgList(sg1Pt, bg1Pt), sg1frac);
  
  RooProdPdf model0("model0", "complete model 0", RooArgSet(model0M, model0Pt)); 
  RooProdPdf model1("model1", "complete model 1", RooArgSet(model1M, model1Pt)); 
  RooDataSet *data0 = model0.generate(RooArgSet(m, pt), nsg0 + nbg);  
  RooDataSet *data1 = model1.generate(RooArgSet(m, pt), nsg1 + nbg);  

  cout << "XXXXX initial signal fractions, 0 = " << sg0frac.getVal() 
       << ", 1 = " << sg1frac.getVal() 
       << endl;


  RooFitResult *resData0Model0 = model0.fitTo(*data0, Save()); 
  RooPlot *f0M = m.frame(); 
  data0->plotOn(f0M);  
  model0.plotOn(f0M); 
  model0.plotOn(f0M, Components("bg0M"), LineStyle(kDashed)) ;
  
  RooPlot *f0Pt = pt.frame();  
  data0->plotOn(f0Pt);  
  model0.plotOn(f0Pt); 
  model0.plotOn(f0Pt, Components("bg0Pt"), LineStyle(kDashed)) ;

  RooFitResult *resData0Model1 = model1.fitTo(*data0, Save()); 
  RooPlot *f1M = m.frame(); 
  data0->plotOn(f1M);  
  model1.plotOn(f1M); 
  model1.plotOn(f1M, Components("bg1M"), LineStyle(kDashed)) ;
  
  RooPlot *f1Pt = pt.frame();  
  data0->plotOn(f1Pt);  
  model1.plotOn(f1Pt); 
  model1.plotOn(f1Pt, Components("bg1Pt"), LineStyle(kDashed)) ;


  zone(2, 2);
  f0M->Draw();
  c0->cd(2);
  gPad->SetLogy(1);
  f0Pt->Draw();

  c0->cd(3);
  f1M->Draw();
  c0->cd(4);
  gPad->SetLogy(1);
  f1Pt->Draw();

  c0->SaveAs("toy3-data0.pdf");


  RooFitResult *resData1Model1 =  model1.fitTo(*data1, Save()); 
  RooPlot *g1M = m.frame(); 
  data1->plotOn(g1M);  
  model1.plotOn(g1M); 
  model1.plotOn(g1M, Components("bg1M"), LineStyle(kDashed)) ;
  
  RooPlot *g1Pt = pt.frame();  
  data1->plotOn(g1Pt);  
  model1.plotOn(g1Pt); 
  model1.plotOn(g1Pt, Components("bg1Pt"), LineStyle(kDashed)) ;

  RooFitResult *resData1Model0 = model0.fitTo(*data1, Save()); 
  RooPlot *g0M = m.frame(); 
  data1->plotOn(g0M);  
  model0.plotOn(g0M); 
  model0.plotOn(g0M, Components("bg0M"), LineStyle(kDashed)) ;
  
  RooPlot *g0Pt = pt.frame();  
  data1->plotOn(g0Pt);  
  model0.plotOn(g0Pt); 
  model0.plotOn(g0Pt, Components("bg0Pt"), LineStyle(kDashed)) ;


  zone(2, 2);
  g1M->Draw();
  c0->cd(2);
  gPad->SetLogy(1);
  g1Pt->Draw();

  c0->cd(3);
  g0M->Draw();
  c0->cd(4);
  gPad->SetLogy(1);
  g0Pt->Draw();

  cout << "Fitting model 0 on data 0" << endl;
  resData0Model0->Print("v");
  cout << "Fitting model 1 on data 0" << endl;
  resData0Model1->Print("v");

  cout << "Fitting model 0 on data 1" << endl;
  resData1Model0->Print("v");
  cout << "Fitting model 1 on data 1" << endl;
  resData1Model1->Print("v");

  c0->SaveAs("toy3-data1.pdf");

  fHistFile->Close(); 

}


// ----------------------------------------------------------------------
void plotHpt::toy4(double nsg0, double nsg1, double nbg, int ntoy) {

  fHistFile = TFile::Open(fHistFileName.c_str()); 

  int NBINS(10); 

  // mass
  // ====
  
  // -- Signal
  RooRealVar m("m", "m", 100., 150., "GeV"); 
  RooRealVar sgP("sgP", "signal peak mass", 125., 100, 150);  
  sgP.setConstant(kTRUE);
  RooRealVar sgS("sgS", "signal sigma mass", fHiggsMres, 0., 15.);  
  sgS.setConstant(kTRUE);
  RooGaussian sg0M("sg0M", "signal mass", m, sgP, sgS); 
  RooGaussian sg1M("sg1M", "signal mass", m, sgP, sgS); 

  // -- Background
  double bgc0  = 0.0245329; // FIXME!!!
  //  double bgc0e = 0.00418385;  // derived from the normalized unscaled fit! 
  RooRealVar C0("C0", "coefficient #0", bgc0, -1., 1.); 
  //  C0.setRange(bgc0 - bgc0e, bgc0 + bgc0e);
  //   double bgc1  = 0.0001231;
  //   double bgc1e = 3.341922e-5; // derived from the normalized unscaled fit! 
  //   RooRealVar C1("C1", "coefficient #1", bgc1, -1., 1.); 
  //  C1.setRange(bgc1 - bgc1e, bgc1 + bgc1e);
  RooPolynomial bg0M("bg0M", "background 0 gamma gamma mass", m, RooArgList(C0)); 
  RooPolynomial bg1M("bg1M", "background 1 gamma gamma mass", m, RooArgList(C0)); 


  // pT
  // ==
  RooRealVar pt("pt", "pt", 300., 1000., "GeV"); 
  //  SM Higgs count:           17.2289 fitted exponential slope = -0.0108046+/-0.000305826
  // -- SM Higgs
  //  double sg0tau  = -1.09325e-02; // FIXME!!!
  double sg0tau = fSg0Tau; 
  //  double sg0taue = 1.67445e-04;
  double sg0taue = fSg0TauE;
  RooRealVar sg0Tau("sg0Tau", "signal 0 tau", sg0tau, -10., 10.); 
  sg0Tau.setRange(sg0tau - sg0taue, sg0tau + sg0taue); 
  RooExponential sg0Pt("sg0Pt", "signal 0 pT", pt, sg0Tau);   

  //  Contact Higgs count:      45.3107 fitted exponential slope = -0.00707933+/-0.000137285
  // -- "contact" Higgs
  //  double sg1tau  = -7.81368e-03; // FIXME!!!
  double sg1tau = fSg1Tau; 
  //  double sg1taue = 9.10171e-05;
  double sg1taue = fSg1TauE;
  RooRealVar sg1Tau("sg1Tau", "signal 1 tau", sg1tau, -10., 10.); 
  sg1Tau.setRange(sg1tau - sg1taue, sg1tau + sg1taue); 
  RooExponential sg1Pt("sg1Pt", "signal 1 pT", pt, sg1Tau);   

  //  SHERPA background:        168.859 fitted exponential slope = -0.0113285+/-0.000262761
  //  double bgtau  = -1.26957e-02; // FIXME!!!
  double bgtau = fBgTau; 
  //  double bgtaue = 1.55210e-04; 
  RooRealVar bgTau("bgTau", "background gamma gamma tau", bgtau, -10., 10.); 
  RooExponential bg0Pt("bg0Pt", "background 0 gamma gamma pT", pt, bgTau);   
  RooExponential bg1Pt("bg1Pt", "background 1 gamma gamma pT", pt, bgTau);   


  RooRealVar sg0frac("sg0frac","fraction of Higgs 0 signal", nsg0/nbg, 0., 1.);
  RooRealVar sg1frac("sg1frac","fraction of Higgs 1 signal", nsg1/nbg, 0., 1.);
  RooAddPdf model0M("model0M", "model 0 for mass", RooArgList(sg0M, bg0M), sg0frac);
  RooAddPdf model0Pt("model0Pt", "model 0 for pT", RooArgList(sg0Pt, bg0Pt), sg0frac);

  RooAddPdf model1M("model1M", "model 1 for mass", RooArgList(sg1M, bg1M), sg1frac);
  RooAddPdf model1Pt("model1Pt", "model 1 for pT", RooArgList(sg1Pt, bg1Pt), sg1frac);
  
  RooProdPdf model0("model0", "complete model 0", RooArgSet(model0M, model0Pt)); 
  RooProdPdf model1("model1", "complete model 1", RooArgSet(model1M, model1Pt)); 

  // -- alternative way to assemble data sets
  RooProdPdf modelBg("modelBg", "complete background model", RooArgSet(bg0M, bg0Pt)); 
  RooProdPdf model0Sg("model0Sg", "complete signal model 0", RooArgSet(sg0M, sg0Pt)); 
  RooProdPdf model1Sg("model1Sg", "complete signal model 1", RooArgSet(sg1M, sg1Pt)); 

//   RooDataSet *data0 = model0.generate(RooArgSet(m, pt), nsg0 + nbg);  
//   RooDataSet *data1 = model1.generate(RooArgSet(m, pt), nsg1 + nbg);  

  RooDataSet *bg0 = modelBg.generate(RooArgSet(m, pt), nbg);  
  RooDataSet *data0 = new RooDataSet(*bg0); 
  RooDataSet *data1 = new RooDataSet(*data0); 

  RooDataSet *s0 = model0Sg.generate(RooArgSet(m, pt), nsg0);  
  RooDataSet *s1 = model1Sg.generate(RooArgSet(m, pt), nsg1);  
  data0->append(*s0); 
  data1->append(*s1); 


  TVirtualPad *pad = gPad; 
  makeCanvas(1); 
  zone(2,2, c1);
  
  RooPlot *f0M = m.frame(); 
  f0M->SetTitle("m(top) = 173.5 GeV");
  data0->plotOn(f0M, Binning(NBINS));  
  model0.plotOn(f0M); 
  model0.plotOn(f0M, Components("bg0M"), LineStyle(kDashed)) ;
  model0.paramOn(f0M, Layout(0.15, 0.85)) ;
  f0M->Draw();

  c1->cd(2);
  RooPlot *a = m.frame(); 
  a->SetTitle("signal (blue) and background (red) for m(top) = 173.5 GeV");
  s0->plotOn(a, Binning(NBINS), LineColor(kBlue), MarkerColor(kBlue));  
  bg0->plotOn(a, Binning(NBINS), LineColor(kRed), MarkerColor(kRed)); 
  a->Draw();

  c1->cd(4);
  RooPlot *a1 = m.frame(); 
  a1->SetTitle("signal (blue) and background (red) for m(top) #rightarrow #infty");
  s1->plotOn(a1, Binning(NBINS), LineColor(kBlue), MarkerColor(kBlue));  
  bg0->plotOn(a1, Binning(NBINS), LineColor(kRed), MarkerColor(kRed)); 
  a1->Draw();
  
  double bExpected     = sg0frac.getVal()*data0->numEntries();
  double bExpectedRelE = sg0frac.getError()/sg0frac.getVal();

  cout << "--> bExpected = " << bExpected << " +/- " << bExpected*bExpectedRelE 
       << " (rel: " << bExpectedRelE << "), sg0frac.getError()= " << sg0frac.getError() 
       << " sg0frac.getVal() = " << sg0frac.getVal() << ", data0 entries = " << data0->numEntries()
       << endl;

 
  RooPlot *g1M = m.frame(); 
  g1M->SetTitle("m(top) #rightarrow #infty");
  data1->plotOn(g1M, Binning(NBINS));  
  model1.plotOn(g1M); 
  model1.plotOn(g1M, Components("bg1M"), LineStyle(kDashed)) ;
  c1->cd(3);
  model1.paramOn(g1M, Layout(0.15, 0.85)) ;
  g1M->Draw();
  
  
  double sExpected     = sg1frac.getVal()*data1->numEntries();
  double sExpectedRelE = sg1frac.getError()/sg1frac.getVal();
  
  cout << "--> bExpected = " << bExpected << " +/- " << bExpected*bExpectedRelE 
       << " (rel: " << bExpectedRelE << "), sg0frac.getError()= " << sg0frac.getError() 
       << " sg0frac.getVal() = " << sg0frac.getVal() << ", data0 entries = " << data0->numEntries()
       << endl;

  cout << "--> sExpected = " << sExpected << " +/- " << sExpected*sExpectedRelE 
       << " (rel: " << sExpectedRelE << "), sg1frac.getError() = " << sg1frac.getError() 
       << " sg1frac.getVal() = " << sg1frac.getVal() << ", data1 entries = " << data1->numEntries()
       << endl;
  
  double zExp = NumberCountingUtils::BinomialExpZ(sExpected, bExpected, bExpectedRelE);
  cout << "--> zExp: " << zExp << endl;

  
  double tau = 6.; 
  double zExpWithTau = NumberCountingUtils::BinomialWithTauExpZ(sExpected, bExpected, tau);
  cout << "--> zExp(tau): " << zExpWithTau << endl;

  // -- now setup toy loop
  RooCmdArg fitargs; 
  if (1) {
  fitargs.addArg(Minos(kTRUE)); 
  fitargs.addArg(PrintLevel(-1));
  fitargs.addArg(PrintEvalErrors(-1));
  fitargs.addArg(Verbose(kFALSE));
  fitargs.addArg(Save()); 
  fitargs.addArg(Minimizer("Minuit2")); 
  }
  TH1D *h1z = new TH1D("h1z", " ", 100, 0., 10.); 
  TH1D *h1zTh = new TH1D("h1zTh", "significance wrt theory expectation", 100, 0., 10.); 
  TH1D *h1zZero = new TH1D("h1zZero", "significance wrt zero expectation ", 100, 0., 10.); 
  RooFitResult *fr0(0), *fr1(0); 
  zone(2, 1, c1);
  double zExpTh(0.); 
  double zExpZero(0.); 
  for (int imc = 0; imc < fNtoy; ++imc) {
    data0->reset();
    data1->reset();
    cout << "==> starting run " << imc << endl;
    C0.setVal(bgc0);
    bgTau.setVal(bgtau);

    sg0Tau.setVal(sg0tau);
    sg0frac.setVal(nsg0/nbg);

    sg1Tau.setVal(sg1tau);
    sg1frac.setVal(nsg1/nbg);

    if (0) {
      modelBg.Print("t");
      model0Sg.Print("t");
      model1Sg.Print("t");
    }

    bg0 = modelBg.generate(RooArgSet(m, pt), nbg);  
    s0 = model0Sg.generate(RooArgSet(m, pt), nsg0);  
    s1 = model1Sg.generate(RooArgSet(m, pt), nsg1);  
    
    data0->append(*bg0); 
    data0->append(*s0); 

    data1->append(*bg0); 
    data1->append(*s1); 

    cout << "----------------------------------------------------------------------" << endl;
    cout << "Fit model 0 to data 0" << endl;
    cout << "----------------------------------------------------------------------" << endl;
    fr0 = model0.fitTo(*data0, fitargs); 
    cout << "----------------------------------------------------------------------" << endl;
    cout << "Fit model 1 to data 1" << endl;
    cout << "----------------------------------------------------------------------" << endl;
    fr1 = model1.fitTo(*data1, fitargs); 

    bExpected     = sg0frac.getVal()*data0->numEntries();
    bExpectedRelE = sg0frac.getError()/sg0frac.getVal();

    sExpected     = sg1frac.getVal()*data1->numEntries();
    sExpectedRelE = sg1frac.getError()/sg1frac.getVal();
 
    int reset(0); 
    if (0 && bExpected < 1.) {
      reset = 1; 
      bExpected = 2.3; 
      bExpectedRelE = 0.99;
    }

    if (sg0frac.getVal() < sg0frac.getError()) {
      bExpected = nsg0; 
      bExpectedRelE = TMath::Sqrt(nsg0)/nsg0; 
    }
    
    zExp = NumberCountingUtils::BinomialExpZ(sExpected, bExpected, bExpectedRelE);
    zExpTh = NumberCountingUtils::BinomialExpZ(sExpected, nsg0, TMath::Sqrt(nsg0)/nsg0);
    zExpZero = NumberCountingUtils::BinomialExpZ(sExpected, 2.3, 0.99);
    cout << Form("%3d zExp = %5.4f from bExpected = %4.1f+/-%4.1f and sExpected = %4.1f (reset = %d, zExpTh = %4.1f, zExpZero = %4.1f)", 
		 imc, zExp, bExpected, bExpected*bExpectedRelE, sExpected, reset, zExpTh, zExpZero) 
	 << endl;
    if (zExp < 0) {
      --imc;
      continue;
    }

    h1z->Fill(zExp); 
    h1zTh->Fill(zExpTh); 
    h1zZero->Fill(zExpZero); 

    if (zExp < 0.1) {
      RooPlot *g0M = m.frame(); 
      g0M->SetTitle("m(top) = 173.5GeV");
      data0->plotOn(g0M, Binning(NBINS));  
      model0.plotOn(g0M); 
      model0.plotOn(g0M, Components("bg0M"), LineStyle(kDashed)) ;
      c1->cd(1);
      model0.paramOn(g0M, Layout(0.15, 0.85)) ;
      g0M->Draw();

      RooPlot *g1M = m.frame(); 
      g1M->SetTitle("m(top) #rightarrow #infty");
      data1->plotOn(g1M, Binning(NBINS));  
      model1.plotOn(g1M); 
      model1.plotOn(g1M, Components("bg1M"), LineStyle(kDashed)) ;
      c1->cd(2);
      model1.paramOn(g1M, Layout(0.15, 0.85)) ;
      g1M->Draw();
      tl->SetNDC(kTRUE);
      tl->DrawLatex(0.2, 0.2, Form("zExp = %4.3f", zExp)); 
      tl->DrawLatex(0.2, 0.15, Form("S = %4.1f, B = %4.1f +/- %4.1f", sExpected, bExpected, bExpectedRelE*bExpected)); 
      c1->SaveAs(Form("%s/toy4-%d.pdf", fDirectory.c_str(), imc)); 
    }
  }
  pad->cd();
  gStyle->SetOptStat(0); 
  h1z->Draw();
  tl->SetTextColor(kBlack); 
  tl->DrawLatex(0.6, 0.80, Form("Mean = %4.3f", h1z->GetMean())); 
  tl->DrawLatex(0.6, 0.75, Form("RMS = %4.3f", h1z->GetRMS())); 
  tl->DrawLatex(0.6, 0.70, Form("Mean(TH) =  %4.3f", h1zTh->GetMean())); 
  tl->DrawLatex(0.6, 0.65, Form("Mean(zero) =  %4.3f", h1zZero->GetMean())); 
  tl->DrawLatex(0.6, 0.92, Form("Setup: %s", fSetup.c_str())); 
  fHistFile->Close();

  fTEX.open(fTexFileName.c_str(), ios::app);
  fTEX << formatTex(h1z->GetMean(), Form("%s:zexp:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(h1z->GetRMS(), Form("%s:zexp:rms", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(h1zTh->GetMean(), Form("%s:zexpTh:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(h1zTh->GetRMS(), Form("%s:zexpTh:rms", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(h1zZero->GetMean(), Form("%s:zexpZero:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(h1zZero->GetRMS(), Form("%s:zexpZero:rms", fSetup.c_str()), 2) << endl;
  fTEX.close();
}






// ----------------------------------------------------------------------
void plotHpt::loadFiles(string afiles) {
  
  string files = fDirectory + "/" + afiles;
  cout << "==> Loading files listed in " << files << endl;

  // -- mH = 125.0, from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3
  //  fBF["H2GammaGamma"] = 2.28E-03;

  char buffer[1000];
  ifstream is(files.c_str());
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    
    string sbuffer = string(buffer); 

    string::size_type m1 = sbuffer.find("xsec="); 
    string stype = sbuffer.substr(5, m1-6); 

    string::size_type m2 = sbuffer.find("file="); 
    string sxsec = sbuffer.substr(m1+5, m2-m1-6); 
    string sfile = sbuffer.substr(m2+5); 
    string sname, sdecay; 

    TFile *pF(0); 
    if (string::npos != stype.find("data")) {
      // -- DATA
      cout << "XXXX do not know what to do with data?!" << endl;
    } else {
      // -- MC
      pF = loadFile(sfile); 
      TTree *t = (TTree*)pF->Get("events"); 
      int nevt = t->GetEntries();
      if (string::npos != stype.find("mcatnlo")) {
	dataset *ds = new dataset(); 
	if (string::npos != stype.find("0")) {
	  sname = "mcatnlo0"; 
	  sdecay = "m_{t}=173.5GeV";
	  ds->fColor = kBlue; 
	  ds->fLcolor = kBlue; 
	  ds->fFcolor = kBlue; 
	  ds->fSymbol = 24; 
	  ds->fFillStyle = 3305; 
	} else if (string::npos != stype.find("1")) {
	  sname = "mcatnlo1";
	  sdecay = "blurp";
	  ds->fColor = kBlue+2; 
	  ds->fLcolor = kBlue+2; 
	  ds->fFcolor = kBlue+2; 
	  ds->fSymbol = 25; 
	  ds->fFillStyle = 3395; 
	} else if (string::npos != stype.find("5")) {
	  sname = "mcatnlo5";
	  sdecay = "m_{t}#rightarrow#infty";
	  ds->fColor = kBlack; 
	  ds->fLcolor = kBlack; 
	  ds->fFcolor = kBlack; 
	  ds->fSymbol = 26; 
	  ds->fFillStyle = 3356; 
	} 
	ds->fF = pF; 
	ds->fXsec = atof(sxsec.c_str());
	ds->fBf   = 2.28E-03;
	ds->fLumi = nevt/ds->fXsec/ds->fBf/1000.;
	ds->fName = "MC@NLO " + sdecay; 
	ds->fSize = 1; 
	ds->fWidth = 2; 
	fDS.insert(make_pair(sname, ds)); 
      }


      if (string::npos != stype.find("sherpa")) {
	dataset *ds = new dataset(); 
	sname = "sherpa";
	sdecay = "#gamma #gamma j {j}";
	ds->fColor = kRed; 
	ds->fLcolor = kRed; 
	ds->fFcolor = kRed; 
	ds->fSymbol = 27; 
	ds->fF = pF; 
	ds->fXsec = atof(sxsec.c_str());
	ds->fBf   = 1.;
	ds->fLumi = nevt/ds->fXsec/ds->fBf/1000.;
	ds->fName = "SHERPA " + sdecay; 
	ds->fFillStyle = 3365; 
	ds->fSize = 1; 
	ds->fWidth = 2; 
	fDS.insert(make_pair(sname, ds)); 
      } 

      // mb ub nb pb fb 
      cout << "MC "  << sfile  << " as " << sname << " (" << stype << ") with xsec = " << sxsec
	   << " eq. lumi = " << fDS[sname]->fLumi << "/fb"
	   << " (Nevts = " << nevt << ")"
	   << endl;

    }
  }
}

// ----------------------------------------------------------------------
void plotHpt::candAnalysis() {
  fGoodCand = true; 

  if (fb.m < MGGLO) fGoodCand = false; 
  if (fb.m > MGGHI) fGoodCand = false; 
  if (TMath::Abs(fb.g0eta) > GETA) fGoodCand = false; 
  if (TMath::Abs(fb.g1eta) > GETA) fGoodCand = false; 

  if (fb.g0iso > G0ISO) fGoodCand = false; 
  if (fb.g1iso > G1ISO) fGoodCand = false; 

  fGoodCandNoPtCuts = fGoodCand;

  if (fb.g0pt < G0PTLO) fGoodCand = false; 
  if (fb.g1pt < G1PTLO) fGoodCand = false; 
//   if (fb.g0pt < G0PT) fGoodCand = false; 
//   if (fb.g1pt < G1PT) fGoodCand = false; 
}

// ----------------------------------------------------------------------
void plotHpt::loopFunction2() {

  // -- cand selection cuts
  if (fb.m < MGGLO) return; 
  if (fb.m > MGGHI) return;
  if (TMath::Abs(fb.g0eta) > GETA) return;
  if (TMath::Abs(fb.g1eta) > GETA) return;
  if (fb.g0iso > G0ISO) return;
  if (fb.g1iso > G1ISO) return;

  for (unsigned int i = 0; i < fSelPoints.size(); ++i) {
    fSelPoints[i].eval(fOptMode, 1.);
  }

}


// ----------------------------------------------------------------------
void plotHpt::loopFunction3() {

  if (!fGoodCand) return;

  if (fb.pt > PTLO /*&& fb.pt < PTHI*/ && fb.g0pt > G0PT && fb.g1pt > G1PT) {
    fHists["bgshape_default_pt"]->Fill(fb.pt); 
    fHists["bgshape_default_m"]->Fill(fb.m); 
  }

  int ptbin = fb.pt/100; 
  if (ptbin > NPTBINS - 1) ptbin = NPTBINS - 1; 
  //  cout << "pt = " << fb.pt << " -> ptbin = " << ptbin << endl;
  fHists[Form("bgshape_ptbin%d_pt", ptbin)]->Fill(fb.pt); 
  fHists[Form("bgshape_ptbin%d_m", ptbin)]->Fill(fb.m); 
  
}

// ----------------------------------------------------------------------
void plotHpt::loopFunction1() {
  char cds[100], cut[100];

  sprintf(cds, "%s", fCds.c_str());
  
  fHists[Form("genpt_%s_goodcand", cds)]->Fill(fb.gpt); 
  fHists[Form("genpt_g1pt_%s_goodcand", cds)]->Fill(fb.gpt, fb.gg1pt); 

  if (fGoodCandNoPtCuts) { 
    fHists[Form("genpt_%s_nopt", cds)]->Fill(fb.gpt); 
    fHists[Form("genpt_g1pt_%s_nopt", cds)]->Fill(fb.gpt, fb.gg1pt); 
    fHists[Form("pt_%s_nopt", cds)]->Fill(fb.pt); 
    fHists[Form("m_%s_nopt", cds)]->Fill(fb.m); 
    fHists[Form("mpt_%s_nopt", cds)]->Fill(fb.pt, fb.m); 
  }

  if (fGoodCand) { 
    fHists[Form("pt_%s_goodcand", cds)]->Fill(fb.pt); 
    fHists[Form("m_%s_goodcand", cds)]->Fill(fb.m); 
    fHists[Form("mpt_%s_goodcand", cds)]->Fill(fb.pt, fb.m); 
    
    if (fb.pt < PTNL) return;
    if (fb.pt > PTNL && fb.pt < PTNH && fb.g0pt > G0PTLO && fb.g1pt > G1PTLO)  { 
      sprintf(cut, "lopt"); 
      fHists[Form("m_%s_%s", cds, cut)]->Fill(fb.m); 
      fHists[Form("pt_%s_%s", cds, cut)]->Fill(fb.pt); 
      fHists[Form("mpt_%s_%s", cds, cut)]->Fill(fb.pt, fb.m); 
      fHists[Form("eta_%s_%s", cds, cut)]->Fill(fb.eta); 
      fHists[Form("g0pt_%s_%s", cds, cut)]->Fill(fb.g0pt); 
      fHists[Form("g1pt_%s_%s", cds, cut)]->Fill(fb.g1pt); 
      fHists[Form("g0iso_%s_%s", cds, cut)]->Fill(fb.g0iso); 
      fHists[Form("g1iso_%s_%s", cds, cut)]->Fill(fb.g1iso); 
      fHists[Form("g0isor_%s_%s", cds, cut)]->Fill(fb.g0iso/fb.g0pt); 
      fHists[Form("g1isor_%s_%s", cds, cut)]->Fill(fb.g1iso/fb.g1pt); 
    } 
    if (fb.pt > PTLO && fb.pt < PTHI && fb.g0pt > G0PT && fb.g1pt > G1PT) {
      sprintf(cut, "hipt"); 
      fHists[Form("m_%s_%s", cds, cut)]->Fill(fb.m); 
      fHists[Form("pt_%s_%s", cds, cut)]->Fill(fb.pt); 
      fHists[Form("mpt_%s_%s", cds, cut)]->Fill(fb.pt, fb.m); 
      fHists[Form("eta_%s_%s", cds, cut)]->Fill(fb.eta); 
      fHists[Form("g0pt_%s_%s", cds, cut)]->Fill(fb.g0pt); 
      fHists[Form("g1pt_%s_%s", cds, cut)]->Fill(fb.g1pt); 
      fHists[Form("g0iso_%s_%s", cds, cut)]->Fill(fb.g0iso); 
      fHists[Form("g1iso_%s_%s", cds, cut)]->Fill(fb.g1iso); 
      fHists[Form("g0isor_%s_%s", cds, cut)]->Fill(fb.g0iso/fb.g0pt); 
      fHists[Form("g1isor_%s_%s", cds, cut)]->Fill(fb.g1iso/fb.g1pt); 
    }
  }

}


// ----------------------------------------------------------------------
void plotHpt::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
  int nentries = Int_t(t->GetEntries());
  int nbegin(0), nend(nentries); 
  if (nevts > 0 && nentries > nevts) {
    nentries = nevts;
    nbegin = 0; 
    nend = nevts;
  }
  if (nevts > 0 && nstart > 0) {
    nentries = nstart + nevts;
    nbegin = nstart; 
    if (nstart + nevts < t->GetEntries()) {
      nend = nstart + nevts; 
    } else {
      nend = t->GetEntries();
    }
  }
  
  nentries = nend - nstart; 
  
  int step(1000000); 
  if (nentries < 5000000)  step = 500000; 
  if (nentries < 1000000)  step = 100000; 
  if (nentries < 100000)   step = 10000; 
  if (nentries < 10000)    step = 1000; 
  if (nentries < 1000)     step = 100; 
  step = 500000; 
  cout << "==> plotHpt::loopOverTree> loop over dataset " << fCds << " in file " 
       << t->GetDirectory()->GetName() 
       << " with " << nentries << " entries" 
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotHpt::*pF)(void);
  if (ifunc == 1) pF = &plotHpt::loopFunction1;
  if (ifunc == 2) pF = &plotHpt::loopFunction2;
  if (ifunc == 3) pF = &plotHpt::loopFunction3;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
   
    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotHpt::setupTree(TTree *t) {
  t->SetBranchAddress("type", &fb.type);
  t->SetBranchAddress("m", &fb.m); 
  t->SetBranchAddress("w8", &fb.w8); 
  t->SetBranchAddress("pt", &fb.pt);
  t->SetBranchAddress("eta", &fb.eta);
  t->SetBranchAddress("phi", &fb.phi);
  t->SetBranchAddress("gm", &fb.gm);
  t->SetBranchAddress("gpt", &fb.gpt);
  t->SetBranchAddress("geta", &fb.geta);
  t->SetBranchAddress("gphi", &fb.gphi);
  t->SetBranchAddress("g0pt", &fb.g0pt);
  t->SetBranchAddress("g0eta", &fb.g0eta);
  t->SetBranchAddress("g0phi", &fb.g0phi);
  t->SetBranchAddress("g0iso", &fb.g0iso);
  t->SetBranchAddress("g1pt", &fb.g1pt);
  t->SetBranchAddress("g1eta", &fb.g1eta);
  t->SetBranchAddress("g1phi", &fb.g1phi);
  t->SetBranchAddress("g1iso", &fb.g1iso);
  t->SetBranchAddress("gg0pt", &fb.gg0pt);
  t->SetBranchAddress("gg0eta", &fb.gg0eta);
  t->SetBranchAddress("gg0phi", &fb.gg0phi);
  t->SetBranchAddress("gg0iso", &fb.gg0iso);
  t->SetBranchAddress("gg1pt", &fb.gg1pt);
  t->SetBranchAddress("gg1eta", &fb.gg1eta);
  t->SetBranchAddress("gg1phi", &fb.gg1phi);
  t->SetBranchAddress("gg1iso", &fb.gg1iso);
}



// ----------------------------------------------------------------------
TH1* plotHpt::combMCAtNLOHist(TH1 *h0, TH1 *h1) {
  double xs0(36.44), xs1(-2.37); 

  string h0name = h0->GetName(); 
  replaceAll(h0name, "mcatnlo0", "mcatnlo"); 

  TH1 *h(0);
  if (h0->InheritsFrom(TH2D::Class())) {
    h = (TH2D*)h0->Clone(h0name.c_str());
    h->Reset();
  } else {
    h = (TH1D*)h0->Clone(h0name.c_str());
    h->Reset();
    //     h = new TH1D(h0name.c_str(), h0->GetTitle(), 
    // 		 h0->GetNbinsX(), 
    // 		 h0->GetBinLowEdge(1), h0->GetBinLowEdge(h0->GetNbinsX()+1));
  }

  h->Add(h0, h1, 1., xs1/xs0); 

  return h;
}


// ----------------------------------------------------------------------
void plotHpt::validation() {

  readHistograms();

  gStyle->SetOptTitle(0);
  
  string what("pt");
  string sel("nopt");
  
  int OTYPE = LUMI;
  zone();
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
 	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
 	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
 	  OTYPE, true, true); 

  c0->SaveAs(Form("%s/plot-%s-%s.pdf", fDirectory.c_str(), what.c_str(), sel.c_str()));

  zone();
  what = "m";
  sel = "hipt"; 
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
 	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
 	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
 	  OTYPE, false, false); 
  c0->SaveAs(Form("%s/plot-%s-%s.pdf", fDirectory.c_str(), what.c_str(), sel.c_str()));

  what = "m";
  sel = "lopt"; 
  zone();
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
 	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
 	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
 	  OTYPE, false, false); 
  c0->SaveAs(Form("%s/plot-%s-%s.pdf", fDirectory.c_str(), what.c_str(), sel.c_str()));

  what = "pt";
  sel = "hipt"; 
  zone();
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
 	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
 	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
 	  OTYPE, true, false); 
  c0->SaveAs(Form("%s/plot-%s-%s.pdf", fDirectory.c_str(), what.c_str(), sel.c_str()));

  zone();
  c0->Clear(); 
  what = "pt";
  sel = "hipt"; 
  TH1D *h1 = (TH1D*)fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())]->Clone("h1");
  h1->Scale(1./h1->Integral()); 
  TH1D *h2 = (TH1D*)fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())]->Clone("h2");
  h2->Scale(1./h2->Integral()); 
  h1->Divide(h2); 
  gStyle->SetOptFit(1);
  TFitResultPtr frp1 = h1->Fit("pol1", "R", "", 300., 1000.);
  TFitResultPtr frp2 = h1->Fit("pol1", "R+", "same", 400., 850.);

  c0->SaveAs(Form("%s/plot-ratio-%s-%s.pdf", fDirectory.c_str(), what.c_str(), sel.c_str()));

}


// ----------------------------------------------------------------------
void plotHpt::allNumbers() {

  tl->SetNDC(kTRUE);
  tl->SetTextColor(kBlack);

  TH1D *h1(0); 

  readHistograms();

  zone(2,3);

  double ptlo(0.), pthi(0.), ptnl(0.), ptnh(0.), g0pt(0.), g1pt(0.); 

  cout << "fullAnalysis: " << endl;
  h1 = (TH1D*)fHists["pt_sherpa_hipt"]; 
  cout << h1 << endl;
  ptlo = h1->GetBinLowEdge(h1->FindFirstBinAbove(1.)); 
  pthi = h1->GetBinLowEdge(h1->GetNbinsX()+1);
  cout << " signal region: " << ptlo << " < pT < " << pthi << endl;
  h1 = (TH1D*)fHists["pt_sherpa_lopt"]; 
  ptnl = h1->GetBinLowEdge(h1->FindFirstBinAbove(1.));
  ptnh = h1->GetBinLowEdge(h1->FindLastBinAbove(1.)+1);
  cout << " norm   region: " << ptnl << " < pT < " << ptnh << endl;

  // -- mass fits in hipt and lopt region
  h1 = (TH1D*)fHists["m_mcatnlo_hipt"]; 
  h1->SetTitle("HIPT"); 
  h1->Fit("gaus", "", "e"); 
  cout << " SM Higgs high-mass fitted resolution: " << h1->GetFunction("gaus")->GetParameter(2) 
       << " +/- " << h1->GetFunction("gaus")->GetParError(2) 
       << " RMS: " << h1->GetRMS() << " +/- " << h1->GetRMSError()
       << endl;
  cout << " SM Higgs mass peak: "       << h1->GetFunction("gaus")->GetParameter(1) 
       << " +/- " << h1->GetFunction("gaus")->GetParError(1) << endl;
  tl->DrawLatex(0.16, 0.8, Form("resolution: %4.1f GeV", h1->GetFunction("gaus")->GetParameter(2))); 
  fHiggsMres  = h1->GetFunction("gaus")->GetParameter(2);
  fHiggsMpeak = h1->GetFunction("gaus")->GetParameter(1);

  // -- printout single photon ET cuts
  h1 = (TH1D*)fHists["g0pt_mcatnlo_hipt"]->Clone("h1");
  g0pt = h1->GetBinLowEdge(h1->FindFirstBinAbove(1.));
  cout << "  Gamma0 pT > " << g0pt << endl;
  tl->SetTextColor(kBlack); tl->DrawLatex(0.16, 0.60, Form("pT(G0) > %4.0f", h1->GetBinLowEdge(h1->FindFirstBinAbove(0.5)))); 
  h1 = (TH1D*)fHists["g1pt_mcatnlo_hipt"]->Clone("h1");
  g1pt = h1->GetBinLowEdge(h1->FindFirstBinAbove(1.)); 
  cout << "  Gamma1 pT > " << g1pt << endl;
  tl->SetTextColor(kBlack); tl->DrawLatex(0.16, 0.55, Form("pT(G1) > %4.0f", h1->GetBinLowEdge(h1->FindFirstBinAbove(0.5)))); 
  

  
  c0->cd(2);
  h1 = (TH1D*)fHists["m_mcatnlo_lopt"]; 
  h1->SetTitle("LOPT"); 
  h1->Fit("gaus", "", "e"); 
  cout << " SM Higgs low-mass fitted resolution:  " << h1->GetFunction("gaus")->GetParameter(2) 
       << " +/- " << h1->GetFunction("gaus")->GetParError(2) 
       << " RMS: " << h1->GetRMS() << " +/- " << h1->GetRMSError()
       << endl;
  cout << " SM Higgs mass peak: "       << h1->GetFunction("gaus")->GetParameter(1) 
       << " +/- " << h1->GetFunction("gaus")->GetParError(1) << endl;
  tl->DrawLatex(0.16, 0.8, Form("resolution: %4.1f GeV", h1->GetFunction("gaus")->GetParameter(2))); 


  // -------
  // -- HIPT
  // -------
  // -- Background 
  c0->cd(3);
  gPad->SetLogy(1);
  h1 = (TH1D*)fHists["pt_sherpa_hipt"]->Clone("h1");
  normHist(h1, "sherpa", LUMI); 
  fBg = h1->GetSumOfWeights(); 
  h1->SetMinimum(0.001); 
  h1->SetMaximum(10000.); 
  h1->Fit("expo", "r", "", PTLO, PTHI); 
  h1->GetFunction("expo")->SetLineColor(fDS["sherpa"]->fColor); 
  fBgTau  = h1->GetFunction("expo")->GetParameter(1); 
  fBgTauE = h1->GetFunction("expo")->GetParError(1); 
  h1->Draw("hist");
  h1->GetFunction("expo")->Draw("same");
  cout << "  SHERPA background:        " << fBg << " fitted exponential slope = " << fBgTau << "+/-" << fBgTauE << endl;
  tl->SetTextColor(fDS["sherpa"]->fColor); tl->DrawLatex(0.48, 0.80, Form("%6.1f", fBg)); 

  // -- SM Higgs (Sg0)
  h1 = (TH1D*)fHists["pt_mcatnlo_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 
  fSg0 = h1->GetSumOfWeights(); 
  h1->Fit("expo", "r", "histsame", PTLO, PTHI); 
  h1->GetFunction("expo")->SetLineColor(fDS["mcatnlo0"]->fColor); 
  fSg0Tau  = h1->GetFunction("expo")->GetParameter(1); 
  fSg0TauE = h1->GetFunction("expo")->GetParError(1); 
  h1->Draw("histsame");
  h1->GetFunction("expo")->Draw("same");
  cout << "  SM Higgs count:           " << fSg0  << " fitted exponential slope = " << fSg0Tau << "+/-" << fSg0TauE << endl;
   tl->SetTextColor(fDS["mcatnlo0"]->fColor); tl->DrawLatex(0.48, 0.75, Form("%6.1f", fSg0)); 
  // -- Contact Higgs (Sg1)
  h1 = (TH1D*)fHists["pt_mcatnlo5_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo5", LUMI); 
  fSg1 = h1->GetSumOfWeights(); 
  h1->Fit("expo", "r", "histsame", PTLO, PTHI); 
  h1->GetFunction("expo")->SetLineColor(fDS["mcatnlo5"]->fColor); 
  fSg1Tau  = h1->GetFunction("expo")->GetParameter(1); 
  fSg1TauE = h1->GetFunction("expo")->GetParError(1); 
  h1->Draw("samehist");
  h1->GetFunction("expo")->Draw("same");
  cout << "  Contact Higgs count:      " << fSg1  << " fitted exponential slope = " << fSg1Tau << "+/-" << fSg1TauE << endl;
  tl->SetTextColor(fDS["mcatnlo5"]->fColor); tl->DrawLatex(0.48, 0.70, Form("%6.1f", fSg1)); 
  cout << "  Signal/Background:      " << fSg0/fBg << " and " << fSg1/fBg << endl;
  tl->SetTextColor(fDS["mcatnlo5"]->fColor); tl->DrawLatex(0.60, 0.60, Form("S/B = %4.3f", fSg0/fBg)); 
  tl->SetTextColor(fDS["mcatnlo0"]->fColor);  tl->DrawLatex(0.60, 0.55, Form("S/B = %4.3f", fSg1/fBg)); 


  // -------
  // -- LOPT
  // -------
  h1 = (TH1D*)fHists["pt_sherpa_lopt"]->Clone("h1");
  normHist(h1, "sherpa", LUMI); 
  c0->cd(3);
  h1->Draw("histsame");
  fNormBg = h1->GetSumOfWeights(); 
  cout << "  SHERPA norm background:   " << fNormBg << endl;
  tl->SetTextColor(fDS["sherpa"]->fColor); tl->DrawLatex(0.7, 0.80, Form("%6.1f", fNormBg)); 

  h1 = (TH1D*)fHists["pt_mcatnlo_lopt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 
  h1->Draw("histsame");
  fNormSg0 = h1->GetSumOfWeights(); 
  cout << "  SM Higgs norm count:      " << fNormSg0 << endl;
  tl->SetTextColor(fDS["mcatnlo0"]->fColor); tl->DrawLatex(0.7, 0.75, Form("%6.1f", fNormSg0)); 

  h1 = (TH1D*)fHists["pt_mcatnlo5_lopt"]->Clone("h1");
  normHist(h1, "mcatnlo5", LUMI); 
  h1->Draw("histsame");
  fNormSg1 = h1->GetSumOfWeights(); 
  cout << "  Contact Higgs norm count: " << fNormSg1 << endl;
  tl->SetTextColor(fDS["mcatnlo5"]->fColor); tl->DrawLatex(0.7, 0.70, Form("%6.1f", fNormSg1)); 


  // -- overlay mass histograms: hipt
  TH1D *mhipt = (TH1D*)fHists["m_sherpa_hipt"]->Clone("mhipt");
  mhipt->SetTitle("HIPT"); 
  normHist(mhipt, "sherpa", LUMI); 

  h1 = (TH1D*)fHists["m_mcatnlo_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 

  // -- debug consistency printout
  c0->cd(3); tl->SetTextColor(fDS["sherpa"]->fColor); tl->DrawLatex(0.31, 0.35, Form("%6.1f", mhipt->Integral())); 
  c0->cd(3); tl->SetTextColor(fDS["mcatnlo0"]->fColor); tl->DrawLatex(0.31, 0.30, Form("%6.1f", h1->Integral())); 

  mhipt->Add(h1); 
  c0->cd(5);
  mhipt->SetMinimum(0.);
  mhipt->Draw();


  // -- overlay mass histograms: lopt
  TH1D *mlopt = (TH1D*)fHists["m_sherpa_lopt"]->Clone("mlopt");
  normHist(mlopt, "sherpa", LUMI); 
  mlopt->SetTitle("LOPT"); 

  h1 = (TH1D*)fHists["m_mcatnlo_lopt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 

  // -- debug consistency printout
  c0->cd(3); tl->SetTextColor(fDS["sherpa"]->fColor); tl->DrawLatex(0.31, 0.20, Form("%6.1f", mlopt->Integral())); 
  c0->cd(3); tl->SetTextColor(fDS["mcatnlo0"]->fColor); tl->DrawLatex(0.31, 0.15, Form("%6.1f", h1->Integral())); 

  mlopt->Add(h1); 
  mlopt->SetMinimum(0.);
  c0->cd(6);
  mlopt->Draw();


  // -- run toys
  if (1) {
    c0->cd(4);
    int ntoy = 200; 
    toy4(fSg0, fSg1, fBg, ntoy); 
  }
  
  c0->SaveAs(Form("%s/allNumbers-%s.pdf", fDirectory.c_str(), fSetup.c_str())); 

  fTEX.open(fTexFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << formatTex(g0pt, Form("%s:g0pt:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(g1pt, Form("%s:g1pt:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(ptlo, Form("%s:ptlo:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(pthi, Form("%s:pthi:val", fSetup.c_str()), 2) << endl;
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << formatTex(fSg0, Form("%s:sg0:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fSg1, Form("%s:sg1:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fBg, Form("%s:bg:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fSg0/fBg, Form("%s:sb0:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fSg1/fBg, Form("%s:sb1:val", fSetup.c_str()), 2) << endl;
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << formatTex(fBgTau, Form("%s:bgTau:val", fSetup.c_str()), 6) << endl;
  fTEX << formatTex(fBgTauE, Form("%s:bgTau:err", fSetup.c_str()), 6) << endl;
  fTEX << formatTex(fSg0Tau, Form("%s:sg0Tau:val", fSetup.c_str()), 6) << endl;
  fTEX << formatTex(fSg0TauE, Form("%s:sg0Tau:err", fSetup.c_str()), 6) << endl;
  fTEX << formatTex(fSg1Tau, Form("%s:sg1Tau:val", fSetup.c_str()), 6) << endl;
  fTEX << formatTex(fSg1TauE, Form("%s:sg1Tau:err", fSetup.c_str()), 6) << endl;
  fTEX.close();
  
} 


// ----------------------------------------------------------------------
void plotHpt::setCuts(string cuts) {
  cout << "==> plotHpt::setCuts: " << cuts << endl;

  istringstream ss(cuts);
  string token, name, sval;

  while (getline(ss, token, ',')) {
    
    string::size_type m1 = token.find("="); 
    name = token.substr(0, m1);
    sval = token.substr(m1+1);

    if (string::npos != name.find("PTLO")) {
      float val; 
      val = atof(sval.c_str()); 
      PTLO = val;
    }

    if (string::npos != name.find("PTHI")) {
      float val; 
      val = atof(sval.c_str()); 
      PTHI = val;
    }

    if (string::npos != name.find("PTNL")) {
      float val; 
      val = atof(sval.c_str()); 
      PTNL = val;
    }

    if (string::npos != name.find("PTNH")) {
      float val; 
      val = atof(sval.c_str()); 
      PTNH = val;
    }

    if (string::npos != name.find("G0PT")) {
      float val; 
      val = atof(sval.c_str()); 
      G0PT = val;
    }

    if (string::npos != name.find("G1PT")) {
      float val; 
      val = atof(sval.c_str()); 
      G1PT = val;
    }
  }
}
