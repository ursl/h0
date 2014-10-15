#include "plotHpt.hh"

#include <fstream>
#include <iostream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
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
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"

#include "dataset.hh"
#include "util.hh"

ClassImp(plotHpt)

using namespace std; 
using namespace RooFit; 
using namespace RooStats; 

// ----------------------------------------------------------------------
plotHpt::plotHpt(string dir,  string files, string setup): plotClass(dir, files, setup), fWorkspace("w") {
  loadFiles(files);

  fHistFile = TFile::Open(Form("%s/plotHpt.root", dir.c_str())); 

  NBINS = 25;  
  GETA = 2.5;
  G0ISO = G1ISO = 0.2; 
  G0ISOR = 0.005; 
  G1ISOR = 0.010; 

  PTNL = 150.;
  PTNH = 250.;
  PTLO = 300.;
  PTHI = 10000.;

  MGGLO = 100.;
  MGGHI = 150.;

  fLumi = 1000.;

  G0PT = 150.;
  G1PT =  70.;

}


// ----------------------------------------------------------------------
plotHpt::~plotHpt() {

}


// ----------------------------------------------------------------------
void plotHpt::readHistograms() {
  fHists.clear();

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
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 800.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g1pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1pt");
  sprintf(ahist, "%s", "p_T(#gamma_{1}) [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 300.))); 
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
  if (bitmask & 0x1) treeAnalysis();
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
void plotHpt::treeAnalysis(int nevts) {
  
  TFile *f = TFile::Open(Form("%s/plotHpt.root", fDirectory.c_str()), "RECREATE"); 
  
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

  cout << "hello world" << endl;

  // =======
  // -- mass
  // =======

  // -- Signal
  RooRealVar m("m", "m", 100., 150.); 
  RooRealVar meanM("meanM", "meanM", 125., 100, 150);  
  RooRealVar sigmaM("sigmaM", "sigmaM", 5, 0., 20.);  
  RooGaussian sgM("sgM", "signal mass", m, meanM, sigmaM); 

  // -- Background
  RooRealVar C0("C0", "coefficient #0", 0.086, -1., 1.); 
  RooRealVar C1("C1", "coefficient #1", 0.006, -1., 1.); 
  RooRealVar C2("C2", "coefficient #2", 0.008, -1., 1.); 
  RooChebychev bgM("bgM", "background mass", m, RooArgList(C0, C1, C2)); 


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
  //exponential BG:   RooAddPdf modelPt("modelPt", "model for pT", RooArgList(sgPt, bgPt), fsig) ;
  RooAddPdf modelPt("modelPt", "model for pT", RooArgList(sgHistPt, bgHistPt), fsig) ;

  RooProdPdf model("model", "complete model", RooArgSet(modelM, modelPt)); 

  

  RooDataSet *data = model.generate(RooArgSet(m, pt), 10000);  
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
void plotHpt::overlayAll() {
  // -- simple overlays
  c0->cd(1); 
  overlay("mcatnlo5", "H1pt", "sherpa", "H1pt"); 
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

  if (fb.g0pt < G0PT) fGoodCand = false; 
  if (fb.g1pt < G1PT) fGoodCand = false; 

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
    if (fb.pt > PTNL && fb.pt < PTNH) { 
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
    if (fb.pt > PTLO && fb.pt < PTHI) {
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
