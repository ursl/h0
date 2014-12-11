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
#include "RooMCStudy.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooDLLSignificanceMCSModule.h"

#include "RooStats/NumberCountingUtils.h"
#include "RooGlobalFunc.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestPlot.h"

#include "dataset.hh"
#include "util.hh"

ClassImp(plotHpt)

using namespace std; 
using namespace RooFit; 
using namespace RooStats; 

#include "plotHpt.icc"

// ----------------------------------------------------------------------
plotHpt::plotHpt(string dir,  string files, string setup): plotClass(dir, files, setup), fWorkspace("w") {
  loadFiles(files);

  fNtoy = 1000; 
  fSetup = setup;

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
  
  NPTBINS = 12; 
  GETA = 2.5;
  G0ISO = G1ISO = 0.2; 
  G0ISOR = 0.005; 
  G1ISOR = 0.010; 

  PTNL = 150.;
  PTNH = 400.;
  PTLO = 400.;
  PTHI = 10000.;
  
  NBINS = 55;  
  MGGLO = 70.;
  MGGHI = 180.;

  fLumi = 1000.;
  
  G0PT = 150.;
  G1PT =  70.;
  
  G0PTLO = 60.; 
  G1PTLO = 40.;
}


// ----------------------------------------------------------------------
plotHpt::~plotHpt() {
  if (fHistFile) fHistFile->Close();
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
  h1 = (TH1*)fHistFile->Get("cuts");
  fHists.insert(make_pair("cuts", h1));

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
  
  //  fHistFile->Close();

}

// ----------------------------------------------------------------------
void plotHpt::bookHist(string name, string cuts) {
  char hist[200], thist[200], ahist[200];

  // -- gen pT for efficiency
  sprintf(hist, "%s_%s_%s", "genpt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "genpt");
  sprintf(ahist, "%s", "p_{T}(#gamma#gamma)^{gen} [GeV]"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  // -- gen pT vs g1pt
  sprintf(hist, "genpt_g1pt_%s_%s", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "pt");
  sprintf(ahist, "%s", "p_T^{gen}(#gamma#gamma) [GeV]"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH2D(hist, thist, 100, 0, 1000., 100, 0., 300.))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "pt_{#gamma1}^{gen} [GeV]");
  }

  // -- reco overlays
  sprintf(hist, "%s_%s_%s", "mpt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "mpt");
  sprintf(ahist, "%s", "p_T(#gamma#gamma) [GeV]"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH2D(hist, thist, 100, 0., 1000., NBINS, MGGLO, MGGHI))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "m_{#gamma #gamma} [GeV]");
  }

  sprintf(hist, "%s_%s_%s", "m", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "m");
  sprintf(ahist, "%s", "m(#gamma#gamma) [GeV]"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, NBINS, MGGLO, MGGHI))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  sprintf(hist, "%s_%s_%s", "pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "pt");
  sprintf(ahist, "%s", "pT(#gamma#gamma) [GeV]"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  sprintf(hist, "%s_%s_%s", "eta", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "eta");
  sprintf(ahist, "%s", "#eta(#gamma#gamma)"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, -5., 5.))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  sprintf(hist, "%s_%s_%s", "g0pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0pt");
  sprintf(ahist, "%s", "p_{T}(#gamma_{0}) [GeV]"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 800, 0., 800.))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  sprintf(hist, "%s_%s_%s", "g1pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1pt");
  sprintf(ahist, "%s", "p_T(#gamma_{1}) [GeV]"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 300, 0., 300.))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  sprintf(hist, "%s_%s_%s", "g0iso", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0iso");
  sprintf(ahist, "%s", "I(#gamma_{0})"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 2.))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  sprintf(hist, "%s_%s_%s", "g1iso", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1iso");
  sprintf(ahist, "%s", "I(#gamma_{1})"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 2.))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  sprintf(hist, "%s_%s_%s", "g0isor", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0isor");
  sprintf(ahist, "%s", "I^{rel}(#gamma_{0})"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 0.006))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }

  sprintf(hist, "%s_%s_%s", "g1isor", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1isor");
  sprintf(ahist, "%s", "I^{rel}(#gamma_{1})"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 0.012))); 
    setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  }
}


// ----------------------------------------------------------------------
void plotHpt::makeAll(int bitmask) {
  if (bitmask & 0x1) {
    //    treeAnalysis();
    allNumbers1(); 
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

  if (bitmask & 0x20) {
    // -- different mass cuts
    // MSTW08
    bgSyst("sys132", "sys133", "sys134"); 
    // CTEQ6
    bgSyst("sys150", "sys151", "sys152"); 

    // -- different processes
    bgSyst("sys132", "sys141"); 

    // -- MSTW08 vs CTEQ6
    bgSyst("sys133", "sys152"); 
    bgSyst("sys134", "sys150"); 
    bgSyst("sys135", "sys151"); 

    
  }

}


// ----------------------------------------------------------------------
void plotHpt::bgSyst(string ds0, string ds1, string ds2) {
  
  vector<string> vds; 
  vds.push_back(ds0);
  vds.push_back(ds1);
  if (string::npos == ds2.find("nada")) vds.push_back(ds2);

  for (unsigned int i = 0; i < vds.size(); ++i) {
    fCds = vds[i]; 
    bookHist(vds[i], "nopt"); 
    bookHist(vds[i], "goodcand"); 
    bookHist(vds[i], "lopt"); 
    bookHist(vds[i], "hipt"); 
    TTree *t = getTree(vds[i]); 
    setupTree(t); 
    loopOverTree(t, 1); 
  }

  int OTYPE = LUMI;
  string what("pt");
  string sel("hipt");

  zone();
  if (3 == vds.size()) {
    overlay(fHists[Form("%s_%s_%s", what.c_str(), vds[0].c_str(), sel.c_str())], vds[0], 
	    fHists[Form("%s_%s_%s", what.c_str(), vds[1].c_str(), sel.c_str())], vds[1], 
	    fHists[Form("%s_%s_%s", what.c_str(), vds[2].c_str(), sel.c_str())], vds[2], 
	    OTYPE, false, true, 0.2); 
    c0->SaveAs(Form("%s/bgsyst-%s-%s-%s.pdf", fDirectory.c_str(), vds[0].c_str(), vds[1].c_str(), vds[2].c_str())); 
  } else {
    overlay(fHists[Form("%s_%s_%s", what.c_str(), vds[0].c_str(), sel.c_str())], vds[0], 
	    fHists[Form("%s_%s_%s", what.c_str(), vds[1].c_str(), sel.c_str())], vds[1], 
	    OTYPE, false, true, 0.2); 
    c0->SaveAs(Form("%s/bgsyst-%s-%s.pdf", fDirectory.c_str(), vds[0].c_str(), vds[1].c_str())); 
    
  }
}

// ----------------------------------------------------------------------
void plotHpt::bgShape(int nevts) {
  string ds("sherpa");
  fCds = ds; 

  char hist[200], thist[200];
  sprintf(hist, "bgshape_lopt_m");
  sprintf(thist, "bgshape lopt m");
  fHists.insert(make_pair(hist, new TH1D(hist, thist, NBINS, MGGLO, MGGHI))); 
  fHists[hist]->Sumw2();
  sprintf(hist, "bgshape_lopt_pt");
  sprintf(thist, "bgshape lopt pt");
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
  fHists[hist]->Sumw2();

  sprintf(hist, "bgshape_hipt_m");
  sprintf(thist, "bgshape hipt m");
  fHists.insert(make_pair(hist, new TH1D(hist, thist, NBINS, MGGLO, MGGHI))); 
  fHists[hist]->Sumw2();
  sprintf(hist, "bgshape_hipt_pt");
  sprintf(thist, "bgshape hipt pt");
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
  fHists[hist]->Sumw2();


  sprintf(hist, "bgshape_inpt_m");
  sprintf(thist, "bgshape inpt m");
  fHists.insert(make_pair(hist, new TH1D(hist, thist, NBINS, MGGLO, MGGHI))); 
  fHists[hist]->Sumw2();
  sprintf(hist, "bgshape_inpt_pt");
  sprintf(thist, "bgshape inpt pt");
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

  TH1D *h1 = new TH1D("h1", "", NPTBINS, 0., NPTBINS); 
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetBinLabel(11, "pT > 500 GeV"); 
  h1->GetXaxis()->SetBinLabel(10, "pT > 450 GeV"); 
  h1->GetXaxis()->SetBinLabel(9,  "pT > 400 GeV"); 
  h1->GetXaxis()->SetBinLabel(8,  "pT > 350 GeV"); 
  h1->GetXaxis()->SetBinLabel(7,  "pT > 300 GeV"); 
  h1->GetXaxis()->SetBinLabel(6,  "pT > 250 GeV"); 
  h1->GetXaxis()->SetBinLabel(5,  "pT > 200 GeV"); 
  h1->GetXaxis()->SetBinLabel(4,  "pT > 150 GeV"); 
  h1->GetXaxis()->SetBinLabel(3,  "inpt"); 
  h1->GetXaxis()->SetBinLabel(2,  "hipt"); 
  h1->GetXaxis()->SetBinLabel(1,  "lopt"); 
  
  
  tl->SetNDC(kTRUE);
  zone(1, 4); 
  fHists[Form("bgshape_lopt_m")]->Fit("pol1", "l");
  TF1 *f1 = fHists[Form("bgshape_lopt_m")]->GetFunction("pol1");
  double defx = f1->GetParameter(1)/fHists[Form("bgshape_lopt_m")]->Integral();
  double defxe = f1->GetParError(1)/fHists[Form("bgshape_lopt_m")]->Integral();
  h1->SetBinContent(1, f1->GetParameter(1)/fHists[Form("bgshape_lopt_m")]->Integral()); 
  h1->SetBinError(1, f1->GetParError(1)/fHists[Form("bgshape_lopt_m")]->Integral()); 
  tl->DrawLatex(0.25, 0.8, Form("p1 = %8.7f +/- %8.7f", 
				f1->GetParameter(1)/fHists["bgshape_lopt_m"]->Integral(),  
				f1->GetParError(1)/fHists["bgshape_lopt_m"]->Integral())
		); 
  
  c0->cd(2);
  fHists[Form("bgshape_hipt_m")]->Fit("pol1", "l");
  f1 = fHists[Form("bgshape_hipt_m")]->GetFunction("pol1");
  defx = f1->GetParameter(1)/fHists[Form("bgshape_hipt_m")]->Integral();
  defxe = f1->GetParError(1)/fHists[Form("bgshape_hipt_m")]->Integral();
  h1->SetBinContent(2, f1->GetParameter(1)/fHists[Form("bgshape_hipt_m")]->Integral()); 
  h1->SetBinError(2, f1->GetParError(1)/fHists[Form("bgshape_hipt_m")]->Integral()); 
  tl->DrawLatex(0.25, 0.8, Form("p1 = %8.7f +/- %8.7f", 
				f1->GetParameter(1)/fHists["bgshape_hipt_m"]->Integral(),  
				f1->GetParError(1)/fHists["bgshape_hipt_m"]->Integral())
		); 

  c0->cd(3); 
  fHists[Form("bgshape_inpt_m")]->Fit("pol1", "l");
  f1 = fHists[Form("bgshape_inpt_m")]->GetFunction("pol1");
  defx = f1->GetParameter(1)/fHists[Form("bgshape_inpt_m")]->Integral();
  defxe = f1->GetParError(1)/fHists[Form("bgshape_inpt_m")]->Integral();
  h1->SetBinContent(3, f1->GetParameter(1)/fHists[Form("bgshape_inpt_m")]->Integral()); 
  h1->SetBinError(3, f1->GetParError(1)/fHists[Form("bgshape_inpt_m")]->Integral()); 
  tl->DrawLatex(0.25, 0.8, Form("p1 = %8.7f +/- %8.7f", 
				f1->GetParameter(1)/fHists["bgshape_inpt_m"]->Integral(),  
				f1->GetParError(1)/fHists["bgshape_inpt_m"]->Integral())
		); 
  
  c0->SaveAs(Form("%s/bgshape-default.pdf", fDirectory.c_str())); 
  
  
  tl->SetNDC(kTRUE); 
  TH1D *h0(0); 
  for (int i = 4; i < NPTBINS; ++i) {
    zone(1,3);
    h0 = (TH1D*)fHists[Form("bgshape_ptbin%d_m", i)]->Clone("h0"); 
    h0->Scale(1./h0->Integral()); 
    cout << "XXXXXX FITTING PTBIN "
	 << " with " << fHists[Form("bgshape_ptbin%d_pt", i)]->GetBinLowEdge(fHists[Form("bgshape_ptbin%d_pt", i)]->FindFirstBinAbove(0.5))
	 << " < pT < " << fHists[Form("bgshape_ptbin%d_pt", i)]->GetBinLowEdge(fHists[Form("bgshape_ptbin%d_pt", i)]->FindLastBinAbove(0.5))
	 << endl;
    fHists[Form("bgshape_ptbin%d_m", i)]->Fit("pol1", "l"); 
    f1 = fHists[Form("bgshape_ptbin%d_m", i)]->GetFunction("pol1");
    h1->SetBinContent(i, f1->GetParameter(1)/fHists[Form("bgshape_ptbin%d_m", i)]->Integral()); 
    h1->SetBinError(i, f1->GetParError(1)/fHists[Form("bgshape_ptbin%d_m", i)]->Integral()); 
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
  pl->DrawLine(0., defx, NPTBINS, defx); 
  pl->SetLineStyle(kDashed);
  pl->DrawLine(0., defx+defxe, NPTBINS, defx+defxe); 
  pl->SetLineStyle(kDashed);
  pl->DrawLine(0., defx-defxe, NPTBINS, defx-defxe); 
  
  c0->SaveAs(Form("%s/bgshape-summary-%s.pdf", fDirectory.c_str(), fSetup.c_str())); 

}


// ----------------------------------------------------------------------
void plotHpt::treeAnalysis(int nevts) {

  cout << "treeAnalysis: open " << fHistFileName << endl;
  cout << " cuts: " << endl;
  cout << "  PT > " << PTLO << endl;
  cout << "  G0PT > " << G0PT << endl;
  cout << "  G1PT > " << G1PT << endl;
  
  TFile *f = TFile::Open(fHistFileName.c_str(), "RECREATE"); 

  TH1D *h = new TH1D("cuts", "cuts", 20, 0, 20.); 

  h->SetBinContent(1, PTNL); 
  h->SetBinContent(2, PTNH); 

  h->SetBinContent(3, PTLO); 
  h->SetBinContent(4, PTHI); 

  h->SetBinContent(10, G0PTLO); 
  h->SetBinContent(11, G1PTLO); 
  h->SetBinContent(12, G0PT); 
  h->SetBinContent(13, G1PT); 

  
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

  if (fb.pt > PTLO && fb.g0pt > G0PT && fb.g1pt > G1PT) {
    fHists["bgshape_hipt_pt"]->Fill(fb.pt); 
    fHists["bgshape_hipt_m"]->Fill(fb.m); 
  }

  if (fb.pt > PTNL && fb.pt < PTNH && fb.g0pt > G0PTLO && fb.g1pt > G1PTLO)  { 
    fHists["bgshape_lopt_pt"]->Fill(fb.pt); 
    fHists["bgshape_lopt_m"]->Fill(fb.m); 
  }

  if (fb.pt > PTNH && fb.pt < PTLO && fb.g0pt > G0PTLO && fb.g1pt > G1PTLO)  { 
    fHists["bgshape_inpt_pt"]->Fill(fb.pt); 
    fHists["bgshape_inpt_m"]->Fill(fb.m); 
  }


  int ptbin(0); 
  if (fb.pt > 500)      ptbin = 11;
  else if (fb.pt > 450) ptbin = 10;
  else if (fb.pt > 400) ptbin =  9;
  else if (fb.pt > 350) ptbin =  8;
  else if (fb.pt > 300) ptbin =  7;
  else if (fb.pt > 250) ptbin =  6;
  else if (fb.pt > 200) ptbin =  5;
  else if (fb.pt > 150) ptbin =  4;
  else if (fb.pt > 100) ptbin =  3;
  else if (fb.pt >  50) ptbin =  2;
  else if (fb.pt >   0) ptbin =  1;


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
void plotHpt::allNumbers(int ntoy) {

  tl->SetNDC(kTRUE);
  tl->SetTextColor(kBlack);
  tl->SetTextSize(0.05);

  TH1D *h1(0); 

  readHistograms();

  zone(3, 4);

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
  fNormHiggsMpeak = h1->GetFunction("gaus")->GetParameter(1);
  fNormHiggsMres = h1->GetFunction("gaus")->GetParameter(2); 

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
  TH1D *m0hipt = (TH1D*)fHists["m_sherpa_hipt"]->Clone("m0hipt");
  TH1D *m1hipt = (TH1D*)fHists["m_sherpa_hipt"]->Clone("m1hipt");
  TH1D *mshipt = (TH1D*)fHists["m_sherpa_hipt"]->Clone("mshipt");
  m0hipt->SetTitle("HIPT"); 
  normHist(m0hipt, "sherpa", LUMI); 
  normHist(m1hipt, "sherpa", LUMI); 
  normHist(mshipt, "sherpa", LUMI); 

  h1 = (TH1D*)fHists["m_mcatnlo5_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo5", LUMI); 

  // -- debug consistency printout
  c0->cd(3); tl->SetTextColor(fDS["sherpa"]->fColor); tl->DrawLatex(0.31, 0.35, Form("%6.1f", m0hipt->Integral())); 
  c0->cd(3); tl->SetTextColor(fDS["mcatnlo0"]->fColor); tl->DrawLatex(0.31, 0.30, Form("%6.1f", h1->Integral())); 

  m1hipt->Add(h1); 
  c0->cd(5);
  m1hipt->SetMinimum(0.);
  m1hipt->Draw();

  h1 = (TH1D*)fHists["m_mcatnlo_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 
  m0hipt->Add(h1); 
  c0->cd(5);
  m0hipt->SetMinimum(0.);
  m0hipt->SetMarkerColor(kBlue);
  m0hipt->SetLineColor(kBlue);
  m0hipt->Draw("same");

  mshipt->SetMarkerColor(kBlack);
  mshipt->SetLineColor(kBlack);
  mshipt->Draw("same");


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
  //  mlopt->Draw();


  // -- Fit it!
  RooRealVar m("m", "m", MGGLO, MGGHI, "GeV"); 
  RooRealVar sgP("sgP", "signal peak mass", fNormHiggsMpeak, 100, 150);  
  sgP.setConstant(kTRUE);
  RooRealVar sgS("sgS", "signal sigma mass", fNormHiggsMres, 0., 15.);  
  sgS.setConstant(kTRUE);
  RooGaussian sgM("sgM", "signal mass", m, sgP, sgS); 

  // -- Background
  double bgc0  = 0.0245329; // FIXME!!!
  RooRealVar C0("C0", "coefficient #0", bgc0, -1., 1.); 
  RooPolynomial bgM("bgM", "background gamma gamma mass", m, RooArgList(C0)); 

  RooRealVar nSg("nSg", "signal events", h1->Integral(), 0., 100*h1->Integral());
  RooRealVar nBg("nBg", "background events", mlopt->Integral(), 0., 100.*mlopt->Integral());

  RooAddPdf modelM("modelM", "model for mass", RooArgList(sgM, bgM), RooArgList(nSg, nBg));

  RooDataHist mloptRdh("mloptRdh","lopt mass distribution", RooArgList(m), mlopt);
  modelM.fitTo(mloptRdh, Range(80., 160.)); 

  RooPlot *f0M = m.frame(); 
  f0M->SetTitle("normalization region");
  mloptRdh.plotOn(f0M);  
  modelM.plotOn(f0M); 
  modelM.plotOn(f0M, Components("bgM"), LineStyle(kDashed)) ;
  //  modelM.paramOn(f0M, Layout(0.15, 0.85)) ;
  f0M->Draw();

  double f0 = fSg0/fNormSg0; 
  double f1 = fSg1/fNormSg1; 

  tl->SetTextSize(0.08);
  tl->SetTextColor(kBlack); 
  tl->DrawLatex(0.25, 0.5, Form("SG: %4.1f +/- %3.1f", nSg.getVal(), nSg.getError())); 
  tl->DrawLatex(0.25, 0.4, Form("BG: %4.1f +/- %3.1f", nBg.getVal(), nBg.getError())); 
  tl->DrawLatex(0.25, 0.3, Form("f0: %4.3f", f0));
  tl->DrawLatex(0.25, 0.2, Form("f1: %4.3f", f1)); 
  tl->SetTextSize(0.05);
  

  

  // -- run toys
  if (1) {
    c0->cd(4);
    ntoy = 1; 
    //    toy4(fSg0, fSg1, fBg, ntoy); 
    toy5(fSg0, fSg1, fBg, ntoy); 
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
void plotHpt::allNumbers1(int ntoy) {
  
  tl->SetNDC(kTRUE);
  tl->SetTextColor(kBlack);
  tl->SetTextSize(0.05);
  
  TH1D *h1(0); 

  readHistograms();

  double ptlo(0.), pthi(0.), ptnl(0.), ptnh(0.), g0pt(0.), g1pt(0.), g0ptlo(0.), g1ptlo(0.); 

  cout << "allNumbers1: " << endl;
  TH1 *cuts = (TH1D*)fHists["cuts"]; 
  cout << "cuts = " << cuts << endl;
  ptnl = cuts->GetBinContent(1); 
  ptnh = cuts->GetBinContent(2); 
  ptlo = cuts->GetBinContent(3); 
  pthi = cuts->GetBinContent(4); 

  g0ptlo = cuts->GetBinContent(10); 
  g1ptlo = cuts->GetBinContent(11); 
  g0pt   = cuts->GetBinContent(12); 
  g1pt   = cuts->GetBinContent(13); 
  
  cout << " signal region: " << ptlo << " < pT < " << pthi << endl;
  cout << " norm   region: " << ptnl << " < pT < " << ptnh << endl;
  cout << "  Gamma0 pT > " << g0pt << endl;
  cout << "  Gamma1 pT > " << g1pt << endl;

  zone(3,4);

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
  fNormHiggsMpeak = h1->GetFunction("gaus")->GetParameter(1);
  fNormHiggsMres = h1->GetFunction("gaus")->GetParameter(2); 

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
  TH1D *m0hipt = (TH1D*)fHists["m_sherpa_hipt"]->Clone("m0hipt");
  TH1D *m1hipt = (TH1D*)fHists["m_sherpa_hipt"]->Clone("m1hipt");
  TH1D *mshipt = (TH1D*)fHists["m_sherpa_hipt"]->Clone("mshipt");
  m0hipt->SetTitle("HIPT"); 
  normHist(m0hipt, "sherpa", LUMI); 
  normHist(m1hipt, "sherpa", LUMI); 
  normHist(mshipt, "sherpa", LUMI); 

  h1 = (TH1D*)fHists["m_mcatnlo5_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo5", LUMI); 

  m1hipt->Add(h1); 

  c0->cd(4);
  cuts->SetMinimum(0.);
  cuts->Draw();
  tl->SetTextSize(0.05);
  tl->SetTextColor(kBlack); 
  tl->DrawLatex(0.5, 0.80, Form("Lumi: %4.0f", fLumi)); 

  tl->DrawLatex(0.5, 0.75, Form("pthi: %4.0f", pthi)); 
  tl->DrawLatex(0.5, 0.70, Form("ptlo: %4.0f", ptlo)); 

  tl->DrawLatex(0.5, 0.65, Form("g0pt: %4.0f", g0pt)); 
  tl->DrawLatex(0.5, 0.60, Form("g1pt: %4.0f", g1pt)); 

  tl->DrawLatex(0.5, 0.55, Form("ptnh: %4.0f", ptnh)); 
  tl->DrawLatex(0.5, 0.50, Form("ptnl: %4.0f", ptnl)); 

  tl->DrawLatex(0.5, 0.45, Form("g0ptlo: %4.0f", g0ptlo)); 
  tl->DrawLatex(0.5, 0.40, Form("g1ptlo: %4.0f", g1ptlo)); 

  


  h1 = (TH1D*)fHists["m_mcatnlo_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 
  m0hipt->Add(h1); 

  c0->cd(5);
  m0hipt->SetMinimum(0.);
  m0hipt->SetMarkerColor(kBlue);
  m0hipt->SetLineColor(kBlue);
  m0hipt->Draw("same");

  mshipt->SetMarkerColor(kBlack);
  mshipt->SetLineColor(kBlack);
  mshipt->Draw("same");


  // -- overlay mass histograms: lopt
  TH1D *mlopt = (TH1D*)fHists["m_sherpa_lopt"]->Clone("mlopt");
  normHist(mlopt, "sherpa", LUMI); 
  mlopt->SetTitle("LOPT"); 

  h1 = (TH1D*)fHists["m_mcatnlo_lopt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 

  mlopt->Add(h1); 
  mlopt->SetMinimum(0.);

  c0->cd(6);

  // -- Fit normalization mass
  RooRealVar nm("nm", "m", MGGLO, MGGHI, "GeV"); 
  RooRealVar nsgP("nsgP", "signal peak mass", fNormHiggsMpeak, 100, 150);  
  nsgP.setConstant(kTRUE);
  RooRealVar nsgS("nsgS", "signal sigma mass", fNormHiggsMres, 0., 15.);  
  nsgS.setConstant(kTRUE);
  RooGaussian nsgM("nsgM", "signal mass", nm, nsgP, nsgS); 

  double bgc0  = 0.0245329; // FIXME!!!
  RooRealVar nC0("nC0", "coefficient #0", bgc0, -1., 1.); 
  RooPolynomial nbgM("nbgM", "background gamma gamma mass", nm, RooArgList(nC0)); 

  RooRealVar nnSg("nnSg", "signal events", h1->Integral(), 0., 100*h1->Integral());
  RooRealVar nnBg("nnBg", "background events", mlopt->Integral(), 0., 100.*mlopt->Integral());

  RooAddPdf nmodelM("nmodelM", "model for mass", RooArgList(nsgM, nbgM), RooArgList(nnSg, nnBg));

  RooDataHist nmloptRdh("nmloptRdh","lopt mass distribution", RooArgList(nm), mlopt);
  nmodelM.fitTo(nmloptRdh, Range(80., 160.)); 

  RooPlot *nf0M = nm.frame(); 
  nf0M->SetTitle("normalization region");
  nmloptRdh.plotOn(nf0M);  
  nmodelM.plotOn(nf0M); 
  nmodelM.plotOn(nf0M, Components("bgM"), LineStyle(kDashed)) ;
  nf0M->Draw("hist");

  double f0 = fSg0/fNormSg0; 
  double f1 = fSg1/fNormSg1; 

  tl->SetTextSize(0.08);
  tl->SetTextColor(kBlack); 
  tl->DrawLatex(0.25, 0.5, Form("SG: %4.1f +/- %3.1f", nnSg.getVal(), nnSg.getError())); 
  tl->DrawLatex(0.25, 0.4, Form("BG: %4.1f +/- %3.1f", nnBg.getVal(), nnBg.getError())); 
  tl->DrawLatex(0.25, 0.3, Form("f0: %4.3f", f0));
  tl->DrawLatex(0.25, 0.2, Form("f1: %4.3f", f1)); 
  tl->SetTextSize(0.05);
  

  // -- do the analysis
  int NBINS(11); 

  RooRealVar m("m", "m", MGGLO, MGGHI, "GeV"); 
  RooRealVar sgP("sgP", "signal peak mass", 125., MGGLO, MGGHI);  
  sgP.setConstant(kTRUE);
  RooRealVar sgS("sgS", "signal sigma mass", fHiggsMres, 0., 15.);  
  sgS.setConstant(kTRUE);
  RooGaussian sg0M("sg0M", "signal mass", m, sgP, sgS); 
  RooGaussian sg1M("sg1M", "signal mass", m, sgP, sgS); 

  // -- Background
  //FIXME  double bgc0  = 0.0245329; // FIXME!!!
  //  double bgc0e = 0.00418385;  // derived from the normalized unscaled fit! 
  RooRealVar C0("C0", "coefficient #0", bgc0, -1., 1.); 

  RooPolynomial bg0M("bg0M", "background 0 gamma gamma mass", m, RooArgList(C0)); 
  RooPolynomial bg1M("bg1M", "background 1 gamma gamma mass", m, RooArgList(C0)); 

  RooRealVar sg0N("sg0N","Number of Higgs 0 signal events", fSg0, 0., 1.e4);
  RooRealVar sg1N("sg1N","Numberof Higgs 1 signal events", fSg1, 0., 1.e4);
  RooRealVar bgN("bgN","Number of background events", fBg, 0., 1.e5);

  RooAddPdf model0("model0M", "model 0 for mass", RooArgList(sg0M, bg0M), RooArgList(sg0N, bgN));
  RooAddPdf model1("model1M", "model 1 for mass", RooArgList(sg1M, bg1M), RooArgList(sg1N, bgN));

 
  RooDataSet *data1 = model1.generate(m, fBg+fSg1); 

  RooMCStudy* mcstudy = new RooMCStudy(model1, m, Binned(), Silence(), Extended(kTRUE),
				       FitOptions(Save(kTRUE), Extended(kTRUE), PrintEvalErrors(-1))
				       );

  double nullHypo = fSg0; 
  RooDLLSignificanceMCSModule sigModule(sg1N, nullHypo);
  mcstudy->addModule(sigModule);

  // Generate and fit 1000 samples of Poisson(nExpected) events
  if (ntoy > 0) fNtoy = ntoy; 
  mcstudy->generateAndFit(fNtoy);
  
  TH1 *hdll = mcstudy->fitParDataSet().createHistogram("dll_nullhypo_sg1N",-40,-40) ;
  TH1 *hz = mcstudy->fitParDataSet().createHistogram("significance_nullhypo_sg1N",-40,-40) ;
  
  RooPlot* frame1 = mcstudy->plotParam(sg1N, Bins(40)) ;
  RooPlot* frame2 = mcstudy->plotError(sg1N, Bins(40)) ;
  RooPlot* frame3 = mcstudy->plotPull(sg1N, Bins(40), FitGauss(kTRUE)) ;

  delete mcstudy; 
  delete data1;

  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;

  c0->cd(7) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  c0->cd(8) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  c0->cd(9) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;

  c0->cd(10) ; gPad->SetLeftMargin(0.15) ; hdll->GetYaxis()->SetTitleOffset(1.6) ; hdll->Draw("") ;
  c0->cd(11) ; gPad->SetLeftMargin(0.15) ; hz->GetYaxis()->SetTitleOffset(1.6) ; hz->Draw("") ;

  double zMean = hz->GetMean();
  double zMedian = median(hz); 
  double zRMS  = hz->GetRMS(); 

  tl->SetNDC(kTRUE);
  tl->SetTextColor(kBlack);
  tl->SetTextSize(0.05);
  tl->DrawLatex(0.2, 0.8, Form("Null Hypo: %4.0f", nullHypo)); 
  tl->DrawLatex(0.2, 0.7, Form("Mean: %4.3f +/- %4.3f", zMedian, zRMS)); 
  

  
  c0->SaveAs(Form("%s/allNumbers-%s.pdf", fDirectory.c_str(), fSetup.c_str())); 
  
  fTEX.open(fTexFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << formatTex(g0pt, Form("%s:g0pt:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(g1pt, Form("%s:g1pt:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(ptlo, Form("%s:ptlo:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(pthi, Form("%s:pthi:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(ptnl, Form("%s:ptnl:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(ptnh, Form("%s:ptnh:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(g0ptlo, Form("%s:g0ptlo:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(g1ptlo, Form("%s:g1ptlo:val", fSetup.c_str()), 2) << endl;
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
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << formatTex(zMedian, Form("%s:zMedian:val", fSetup.c_str()), 6) << endl;
  fTEX << formatTex(zMean, Form("%s:zMean:val", fSetup.c_str()), 6) << endl;
  fTEX << formatTex(zRMS, Form("%s:zRMS:val", fSetup.c_str()), 6) << endl;
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


