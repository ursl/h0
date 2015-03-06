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

  PTNL = 100.;
  PTNH = 250.;
  PTLO = 300.;
  PTHI = 10000.;
  
  NBINS = 55;  
  MGGLO = 70.;
  MGGHI = 180.;

  fLumi = 1000.;
  
  G0PT =  80.;
  G1PT =  50.;
  
  G0PTLO = 60.; 
  G1PTLO = 40.;

  fMu = 1.0; 
}


// ----------------------------------------------------------------------
plotHpt::~plotHpt() {
  if (fHistFile) fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotHpt::readHistograms(vector<string> extrads) {

  TH1::SetDefaultSumw2(kTRUE);

  fHists.clear();
  cout << "readHistograms from " << fHistFileName << endl;
  fHistFile = TFile::Open(fHistFileName.c_str()); 

  vector<string> ds;
  ds.push_back("sherpa");
  ds.push_back("mcatnlo5");
  ds.push_back("mcatnlo");

  if (extrads.size() > 0) {
    copy(extrads.begin(), extrads.end(), back_inserter(ds));
  }

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
	//	if (!dataset.compare("mcatnlo")) dataset = "mcatnlo0";
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
    cout << "name ->" << name << "<-" << " hist ->" << hist << "<-" << endl;
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
    treeAnalysis();
    allNumbers5(); 
  }
  if (bitmask & 0x2) {
    allNumbers5(); 
  }
  if (bitmask & 0x4) {
    optimizeCuts(); 
  }
  if (bitmask & 0x8) {
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
void plotHpt::sgSyst(string ds0, string ds1, string ds2) {

  static int first(1); 

  vector<string> v; 
  v.push_back("mcatnlov5");
  v.push_back("mcatnlo5v5"); 
  if (first) {
    first = 0;
    readHistograms(v); 
  }

  vector<string> vds; 
  vds.push_back(ds0);
  vds.push_back(ds1);
  if (string::npos == ds2.find("nada")) vds.push_back(ds2);

//   for (unsigned int i = 0; i < vds.size(); ++i) {
//     fCds = vds[i]; 
//     bookHist(vds[i], "nopt"); 
//     bookHist(vds[i], "goodcand"); 
//     bookHist(vds[i], "lopt"); 
//     bookHist(vds[i], "hipt"); 
//     TTree *t = getTree(vds[i]); 
//     setupTree(t); 
//     loopOverTree(t, 1); 
//   }

  int OTYPE = LUMI;
  string what("pt");
  string sel("hipt");

  zone(2,2);

  TH1D *h0 = (TH1D*)fHists[Form("%s_%s_%s", what.c_str(), vds[0].c_str(), sel.c_str())];
  normHist(h0, vds[0], LUMI); 
  TH1D *h1 = (TH1D*)fHists[Form("%s_%s_%s", what.c_str(), vds[1].c_str(), sel.c_str())];
  normHist(h1, vds[1], LUMI); 
  TH1D *hr0 = (TH1D*)h1->Clone("ratio0"); hr0->Sumw2(); hr0->Reset();
  hr0->Divide(h0, h1); 
  hr0->SetMinimum(0.);
  hr0->SetMaximum(2.);
  TH1D *hr1 = (TH1D*)h1->Clone("ratio1"); hr1->Sumw2(); hr1->Reset();
  hr1->Divide(h0, h1); 
  hr1->SetMinimum(0.);
  hr1->SetMaximum(2.);

  tl->SetNDC(kTRUE); 
  c0->cd(1); 
  h0->Draw(); 

  c0->cd(2); 
  h1->Draw(); 

  c0->cd(3); 
  hr0->Fit("pol0", "R", "", 300., 700.); 
  tl->DrawLatex(0.2, 0.8, Form("offset %5.4f+/-%5.4f", hr0->GetFunction("pol0")->GetParameter(0), 
			       hr0->GetFunction("pol0")->GetParError(0))); 
  tl->DrawLatex(0.2, 0.92, Form("%s/%s", vds[0].c_str(), vds[1].c_str()));

  c0->cd(4); 
  hr1->Fit("pol1", "R", "", 300., 700.); 
  tl->DrawLatex(0.2, 0.80, Form("offset %5.4f+/-%5.4f", hr1->GetFunction("pol1")->GetParameter(0), 
			       hr1->GetFunction("pol1")->GetParError(0))); 
  tl->DrawLatex(0.2, 0.75, Form("slope  %5.4f+/-%5.4f", hr1->GetFunction("pol1")->GetParameter(1), 
			       hr1->GetFunction("pol1")->GetParError(1))); 
  tl->DrawLatex(0.2, 0.92, Form("%s/%s", vds[0].c_str(), vds[1].c_str()));
  
  c0->SaveAs(Form("%s/sgSyst-%s-%s-%s.pdf", fDirectory.c_str(), what.c_str(), vds[0].c_str(), vds[1].c_str())); 

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
void plotHpt::treeAnalysis(int mask, int nevents) {

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

  vector<string> vds; 
  if (mask & 0x1) {
    vds.push_back("sherpa");
    vds.push_back("mcatnlo5");
    vds.push_back("mcatnlo1");
    vds.push_back("mcatnlo0");
  }

  if (mask & 0x2) {
    vds.push_back("man180");
    vds.push_back("man181");
    vds.push_back("man185");

    vds.push_back("man170");
    vds.push_back("man171");
    vds.push_back("man175");

    vds.push_back("man160");
    vds.push_back("man161");
    vds.push_back("man165");

    vds.push_back("man150");
    vds.push_back("man151");
    vds.push_back("man155");
  }

  if (mask & 0x4) {
    vds.push_back("sys132");
    vds.push_back("sys133");
    vds.push_back("sys134");
    vds.push_back("sys135");
    vds.push_back("sys141");
    vds.push_back("sys150");
    vds.push_back("sys151");
    vds.push_back("sys152");
  }

  
  for (unsigned int i = 0; i < vds.size(); ++i) {
    fCds = vds[i]; 
    bookHist(vds[i], "nopt"); 
    bookHist(vds[i], "goodcand"); 
    bookHist(vds[i], "lopt"); 
    bookHist(vds[i], "hipt"); 
    TTree *t = getTree(vds[i]); 
    setupTree(t); 
    loopOverTree(t, 1, nevents); 
  }
  
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
    fHists["bgshape_hipt_pt"]->Fill(fb.pt, fb.w8); 
    fHists["bgshape_hipt_m"]->Fill(fb.m, fb.w8); 
  }

  if (fb.pt > PTNL && fb.pt < PTNH && fb.g0pt > G0PTLO && fb.g1pt > G1PTLO)  { 
    fHists["bgshape_lopt_pt"]->Fill(fb.pt, fb.w8); 
    fHists["bgshape_lopt_m"]->Fill(fb.m, fb.w8); 
  }

  if (fb.pt > PTNH && fb.pt < PTLO && fb.g0pt > G0PTLO && fb.g1pt > G1PTLO)  { 
    fHists["bgshape_inpt_pt"]->Fill(fb.pt, fb.w8); 
    fHists["bgshape_inpt_m"]->Fill(fb.m, fb.w8); 
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
fHists[Form("bgshape_ptbin%d_pt", ptbin)]->Fill(fb.pt, fb.w8); 
  fHists[Form("bgshape_ptbin%d_m", ptbin)]->Fill(fb.m, fb.w8); 
  
}

// ----------------------------------------------------------------------
void plotHpt::loopFunction1() {
  char cds[100], cut[100];

  sprintf(cds, "%s", fCds.c_str());
  
  fHists[Form("genpt_%s_goodcand", cds)]->Fill(fb.gpt, fb.w8); 
  TH2D* h2(0); 
  h2 = (TH2D*)(fHists[Form("genpt_g1pt_%s_goodcand", cds)]);
  h2->Fill(fb.gpt, fb.gg1pt, fb.w8); 

  if (fGoodCandNoPtCuts) { 
    fHists[Form("genpt_%s_nopt", cds)]->Fill(fb.gpt, fb.w8); 
    h2 = (TH2D*)(fHists[Form("genpt_g1pt_%s_nopt", cds)]);
    h2->Fill(fb.gpt, fb.gg1pt, fb.w8); 
    fHists[Form("pt_%s_nopt", cds)]->Fill(fb.pt, fb.w8); 
    fHists[Form("m_%s_nopt", cds)]->Fill(fb.m, fb.w8); 
    h2 = (TH2D*)(fHists[Form("mpt_%s_nopt", cds)]);
    h2->Fill(fb.pt, fb.m, fb.w8); 
  }

  if (fGoodCand) { 
    fHists[Form("pt_%s_goodcand", cds)]->Fill(fb.pt, fb.w8); 
    fHists[Form("m_%s_goodcand", cds)]->Fill(fb.m, fb.w8); 
    h2 = (TH2D*)(fHists[Form("mpt_%s_goodcand", cds)]);
    h2->Fill(fb.pt, fb.m, fb.w8); 
    
    if (fb.pt < PTNL) return;
    if (fb.pt > PTNL && fb.pt < PTNH && fb.g0pt > G0PTLO && fb.g1pt > G1PTLO)  { 
      sprintf(cut, "lopt"); 
      fHists[Form("m_%s_%s", cds, cut)]->Fill(fb.m, fb.w8); 
      fHists[Form("pt_%s_%s", cds, cut)]->Fill(fb.pt, fb.w8); 
      h2 = (TH2D*)(fHists[Form("mpt_%s_%s", cds, cut)]);
      h2->Fill(fb.pt, fb.m, fb.w8); 
      fHists[Form("eta_%s_%s", cds, cut)]->Fill(fb.eta, fb.w8); 
      fHists[Form("g0pt_%s_%s", cds, cut)]->Fill(fb.g0pt, fb.w8); 
      fHists[Form("g1pt_%s_%s", cds, cut)]->Fill(fb.g1pt, fb.w8); 
      fHists[Form("g0iso_%s_%s", cds, cut)]->Fill(fb.g0iso, fb.w8); 
      fHists[Form("g1iso_%s_%s", cds, cut)]->Fill(fb.g1iso, fb.w8); 
      fHists[Form("g0isor_%s_%s", cds, cut)]->Fill(fb.g0iso/fb.g0pt, fb.w8); 
      fHists[Form("g1isor_%s_%s", cds, cut)]->Fill(fb.g1iso/fb.g1pt, fb.w8); 
    } 
    if (fb.pt > PTLO && fb.pt < PTHI && fb.g0pt > G0PT && fb.g1pt > G1PT) {
      sprintf(cut, "hipt"); 
      fHists[Form("m_%s_%s", cds, cut)]->Fill(fb.m, fb.w8); 
      fHists[Form("pt_%s_%s", cds, cut)]->Fill(fb.pt, fb.w8); 
      h2 = (TH2D*)(fHists[Form("mpt_%s_%s", cds, cut)]);
      h2->Fill(fb.pt, fb.m, fb.w8); 
      fHists[Form("eta_%s_%s", cds, cut)]->Fill(fb.eta, fb.w8); 
      fHists[Form("g0pt_%s_%s", cds, cut)]->Fill(fb.g0pt, fb.w8); 
      fHists[Form("g1pt_%s_%s", cds, cut)]->Fill(fb.g1pt, fb.w8); 
      fHists[Form("g0iso_%s_%s", cds, cut)]->Fill(fb.g0iso, fb.w8); 
      fHists[Form("g1iso_%s_%s", cds, cut)]->Fill(fb.g1iso, fb.w8); 
      fHists[Form("g0isor_%s_%s", cds, cut)]->Fill(fb.g0iso/fb.g0pt, fb.w8); 
      fHists[Form("g1isor_%s_%s", cds, cut)]->Fill(fb.g1iso/fb.g1pt, fb.w8); 
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
  double xs0(36.44), xs1(2.37); 

  string h0name = h0->GetName(); 
  replaceAll(h0name, "mcatnlo0", "mcatnlo"); 
  cout << "  combMCAtNLOHist:  " << h0->GetName() << " .. " << h1->GetName() << endl;

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



// ----------------------------------------------------------------------
void plotHpt::anaOptimization() {

  


}


// ----------------------------------------------------------------------
void plotHpt::fitMass(TH1D *hs, double &peak, double &reso) {
  hs->Fit("gaus", "", "e"); 
  cout << " SM Higgs high-mass fitted resolution: " << hs->GetFunction("gaus")->GetParameter(2) 
       << " +/- " << hs->GetFunction("gaus")->GetParError(2) 
       << " RMS: " << hs->GetRMS() << " +/- " << hs->GetRMSError()
       << endl;
  cout << " SM Higgs mass peak: "       << hs->GetFunction("gaus")->GetParameter(1) 
       << " +/- " << hs->GetFunction("gaus")->GetParError(1) << endl;
  tl->DrawLatex(0.16, 0.8, Form("resolution: %4.1f GeV", hs->GetFunction("gaus")->GetParameter(2))); 
  peak  = hs->GetFunction("gaus")->GetParameter(2);
  reso = hs->GetFunction("gaus")->GetParameter(1);
}


// ----------------------------------------------------------------------
void plotHpt::ptDistributions() {
  gPad->SetLogy(1);
  TH1D *h1 = (TH1D*)fHists["pt_sherpa_hipt"]->Clone("h1");
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
}


// ----------------------------------------------------------------------
void plotHpt::displayCuts() {

  TH1 *cuts = (TH1D*)fHists["cuts"]; 
  cout << "cuts = " << cuts << endl;
  fptnl = cuts->GetBinContent(1); 
  fptnh = cuts->GetBinContent(2); 
  fptlo = cuts->GetBinContent(3); 
  fpthi = cuts->GetBinContent(4); 
  
  fg0ptlo = cuts->GetBinContent(10); 
  fg1ptlo = cuts->GetBinContent(11); 
  fg0pt   = cuts->GetBinContent(12); 
  fg1pt   = cuts->GetBinContent(13); 
  
  
  //   cuts->SetMinimum(0.);
  //   cuts->Draw();
  
  double xpos(0.2); 
  
  tl->SetTextSize(0.04);
  tl->SetTextColor(kBlack); 
  tl->DrawLatex(xpos, 0.80, Form("Lumi: %4.0f", fLumi)); 

  tl->DrawLatex(xpos, 0.75, Form("pthi: %4.0f", fpthi)); 
  tl->DrawLatex(xpos, 0.70, Form("ptlo: %4.0f", fptlo)); 

  tl->DrawLatex(xpos, 0.65, Form("g0pt: %4.0f", fg0pt)); 
  tl->DrawLatex(xpos, 0.60, Form("g1pt: %4.0f", fg1pt)); 

  tl->DrawLatex(xpos, 0.55, Form("ptnh: %4.0f", fptnh)); 
  tl->DrawLatex(xpos, 0.50, Form("ptnl: %4.0f", fptnl)); 

  tl->DrawLatex(xpos, 0.45, Form("g0ptlo: %4.0f", fg0ptlo)); 
  tl->DrawLatex(xpos, 0.40, Form("g1ptlo: %4.0f", fg1ptlo)); 

  tl->DrawLatex(xpos, 0.40, Form("g1ptlo: %4.0f", fHiggsMres)); 

  
  
  cout << " signal region: " << fptlo << " < pT < " << fpthi << endl;
  cout << " norm   region: " << fptnl << " < pT < " << fptnh << endl;
  cout << "  Gamma0 pT > " << fg0pt << endl;
  cout << "  Gamma1 pT > " << fg1pt << endl;
}


// ----------------------------------------------------------------------
void plotHpt::massHistograms(TH1D *hb, TH1D *hs0, TH1D *hs1) {
  hs0 = (TH1D*)fHists["m_sherpa_hipt"]->Clone("hs0");
  hs1 = (TH1D*)fHists["m_sherpa_hipt"]->Clone("hs1");
  hb = (TH1D*)fHists["m_sherpa_hipt"]->Clone("hb");
  hs0->SetTitle("HIPT"); 
  normHist(hs0, "sherpa", LUMI); 
  normHist(hs1, "sherpa", LUMI); 
  normHist(hb, "sherpa", LUMI); 

  TH1D *h1 = (TH1D*)fHists["m_mcatnlo5_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo5", LUMI); 

  hs1->Add(h1); 

  h1 = (TH1D*)fHists["m_mcatnlo_hipt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 
  hs0->Add(h1); 

  hs1->Draw();
  hs0->SetMinimum(0.);
  hs0->SetMarkerColor(kBlue);
  hs0->SetLineColor(kBlue);
  hs0->Draw("same");


  RooRealVar m("m", "m", 70., 180.); 
  RooRealVar a0("a0","a0", 0.02, 0., 1.) ;
  RooChebychev bgM("bgM","bgM", m, RooArgSet(a0)) ;

  RooDataHist hdata("hdata","hdata", m, hb);
  RooFitResult *r1 = bgM.fitTo(hdata, SumW2Error(kTRUE), RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));

  RooPlot *plotMh = m.frame(Title("mass"));
  hdata.plotOn(plotMh);
  bgM.plotOn(plotMh);
  bgM.paramOn(plotMh, Layout(0.5, 0.9, 0.98));
  plotMh->Draw("same");

  hb->SetMarkerColor(kBlack);
  hb->SetLineColor(kBlack);
  hb->Fit("pol1", "", "same");
  hb->GetFunction("pol1")->SetLineColor(kBlack);

  fBgMc0  = hb->GetFunction("pol1")->GetParameter(1); 
  fBgMc0E = hb->GetFunction("pol1")->GetParError(1); 
  cout << "xxx> pol1 fBgMp1 = " << fBgMc0 << " +/- " << fBgMc0E << endl;

  fBgMp1  = a0.getVal();
  fBgMp1E = a0.getError(); 
  cout << "xxx> cheb fBgMp1 = " << fBgMp1 << " +/- " << fBgMp1E << endl;
  
  // -- overlay mass histograms: lopt
  TH1D *mlopt = (TH1D*)fHists["m_sherpa_lopt"]->Clone("mlopt");
  normHist(mlopt, "sherpa", LUMI); 
  mlopt->SetTitle("LOPT"); 

  h1 = (TH1D*)fHists["m_mcatnlo_lopt"]->Clone("h1");
  normHist(h1, "mcatnlo0", LUMI); 

  mlopt->Add(h1); 
  mlopt->SetMinimum(0.);
}



// ----------------------------------------------------------------------
void  plotHpt::findMidPoint(TH1D* hq0, TH1D* hq1, double &midpoint, double &tailprob) {
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
double plotHpt::oneSidedGaussianSigma(double prob) {

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


// ----------------------------------------------------------------------
void plotHpt::iniRooVars() {

  cout << "########################################" << endl;
  fMu = 1.; 

  fRsg0N->setConstant(kFALSE);
  fRsg0N->setVal(fSg0); 

  fRsg1N->setConstant(kFALSE);
  fRsg1N->setVal(fSg1); 

  fRsg2N->setConstant(kFALSE);
  fRsg2N->setVal(fSg1); 

  fRsg3N->setConstant(kFALSE);
  fRsg3N->setVal(fSg1); 
  fRsg3Tau->setConstant(kFALSE);
  fRsg3Tau->setVal(fSg1Tau);

  fRbg3N->setConstant(kFALSE);
  fRbg3N->setVal(fBg); 
  fRbg3Slope->setConstant(kFALSE);
  fRbg3Slope->setVal(fBgMp1);
  fRbg3Tau->setConstant(kFALSE);
  fRbg3Tau->setVal(fBgTau);

  fRsg4N->setConstant(kFALSE);
  fRsg4N->setVal(fSg1); 
  fRsg4Tau->setConstant(kFALSE);
  fRsg4Tau->setVal(fSg1Tau);

  fRbg4N->setConstant(kFALSE);
  fRbg4N->setVal(fBg); 
  fRbg4Slope->setConstant(kFALSE);
  fRbg4Slope->setVal(fBgMp1);
  fRbg4Tau->setConstant(kFALSE);
  fRbg4Tau->setVal(fBgTau);

  fRbg0N->setConstant(kFALSE);
  fRbg0N->setVal(fBg); 
  cout << "# fRbg0 = " << fBg << endl;
  fRbg1N->setConstant(kFALSE);
  fRbg1N->setVal(fBg); 
  cout << "# fRbg1 = " << fBg << endl;
  fRbg2N->setConstant(kFALSE);
  fRbg2N->setVal(fBg); 
  cout << "# fRbg2 = " << fBg << endl;
  
  fRsg0Tau->setConstant(kFALSE);
  fRsg0Tau->setVal(fSg0Tau);

  fRsg1Tau->setConstant(kFALSE);
  fRsg1Tau->setVal(fSg1Tau);

  fRbg0Tau->setConstant(kFALSE);
  fRbg0Tau->setVal(fBgTau);

  fRbg1Tau->setConstant(kFALSE);
  fRbg1Tau->setVal(fBgTau);

//   fRbgSlope->setConstant(kFALSE);
//   fRbgSlope->setVal(fBgMp1);

  fRbg0Slope->setConstant(kFALSE);
  fRbg0Slope->setVal(fBgMp1);
  fRbg1Slope->setConstant(kFALSE);
  fRbg1Slope->setVal(fBgMp1);

  cout << "# fBgMp1 = " << fBgMp1 << endl;
  cout << "########################################" << endl;

  

}

// ----------------------------------------------------------------------
void plotHpt::shutup() {

  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(0).removeTopic(Tracing);
  //  RooMsgService::instance().Print() ;
}






// ----------------------------------------------------------------------
void plotHpt::allNumbers5(int ntoy) {
  
  tl->SetNDC(kTRUE);
  tl->SetTextColor(kBlack);
  tl->SetTextSize(0.05);
  
  readHistograms();

  cout << "allNumbers5: " << endl;

  fitMass((TH1D*)fHists["m_mcatnlo_hipt"], fHiggsMres, fHiggsMpeak); 
  fitMass((TH1D*)fHists["m_mcatnlo_lopt"], fNormHiggsMres, fNormHiggsMpeak); 

  zone(3, 3);
  c0->cd(1);
  ptDistributions();

  c0->cd(4);
  TH1D *hb(0), *hs0(0), *hs1(0); 
  massHistograms(hb, hs0, hs1);
  displayCuts(); 
  
  // -- SIGNAL variables
  fRm = new RooRealVar("m", "m", MGGLO, MGGHI, "GeV"); 
  fRsgP = new RooRealVar("sgP", "signal peak mass", 125., MGGLO, MGGHI);  
  fRsgP->setConstant(kTRUE);
  fRsgS = new RooRealVar("sgS", "signal sigma mass", fHiggsMres, 0., 15.);  
  fRsgS->setConstant(kTRUE);

  fRpt = new RooRealVar("pt", "pt", 200., 1000., "GeV"); 

  // -- SM Higgs 1D pdfs
  fRsg0N = new RooRealVar("sg0N","Number of Higgs 0 signal events", fSg0, 0., 1.e4);
  fRsg1N = new RooRealVar("sg1N","Number of Higgs 1 signal events", fSg1, 0., 1.e4);
  fRsg2N = new RooRealVar("sg2N","Number of Higgs 2 signal events", fSg1, 0., 1.e4);

  fRsg0Tau = new RooRealVar("sg0Tau", "signal 0 tau", fSg0Tau, -10., 10.); 
  fRsg1Tau = new RooRealVar("sg1Tau", "signal 1 tau", fSg1Tau, -10., 10.); 
  RooGaussian sg0M("sg0M", "signal 0 mass", *fRm, *fRsgP, *fRsgS); 
  RooExponential sg0Pt("sg0Pt", "signal 0 pT", *fRpt, *fRsg0Tau);   

  RooGaussian sg1M("sg1M", "signal 1 mass", *fRm, *fRsgP, *fRsgS); 
  RooExponential sg1Pt("sg1Pt", "signal 1 pT", *fRpt, *fRsg1Tau);   

  RooGaussian sg2M("sg2M", "signal 2 mass", *fRm, *fRsgP, *fRsgS); 
  fRbg2Slope = new RooRealVar("bg2Slope", "coefficient #0 for bg 2", fBgMp1, -10., 10.); 
  RooPolynomial bg2M("bg2M", "background 2 gamma gamma mass", *fRm, RooArgList(*fRbg2Slope)); 


  // -- BACKGROUND variables and 1D pdfs
  fRbg0Slope = new RooRealVar("bg0Slope", "coefficient #0 for bg 0", fBgMp1, -10., 10.); 
  fRbgSlope  = new RooRealVar("bgSlope", "coefficient #0 for bg", fBgMp1, -10., 10.); 
  fRbg0Tau   = new RooRealVar("bg0Tau", "background 0 gamma gamma tau", fBgTau, -10., 10.); 
  fRbg1Tau   = new RooRealVar("bg1Tau", "background 1 gamma gamma tau", fBgTau, -10., 10.); 
  fRbgTau    = new RooRealVar("bgTau", "background gamma gamma tau", fBgTau, -10., 10.); 

  fRbg0N   = new RooRealVar("bg0N","Number of background 0 events", fBg, 0., 1.e5);
  fRbg1N   = new RooRealVar("bg1N","Number of background 1 events", fBg, 0., 1.e5);
  fRbg2N   = new RooRealVar("bg2N","Number of background 2 events", fBg, 0., 1.e5);
  fRbgN    = new RooRealVar("bgN","Number of background  events", fBg, 0., 1.e5);
  fRbg0Tau = new RooRealVar("bg0Tau", "background 0 gamma gamma tau", fBgTau, -10., 10.); 
  fRbg1Tau = new RooRealVar("bg1Tau", "background 1 gamma gamma tau", fBgTau, -10., 10.); 
  fRbgTau  = new RooRealVar("bgTau", "background gamma gamma tau", fBgTau, -10., 10.); 

  RooChebychev bg0M("bg0M", "background 0 gamma gamma mass", *fRm, RooArgList(*fRbg0Slope)); 
  RooExponential bg0Pt("bg0Pt", "background 0 gamma gamma pT", *fRpt, *fRbg0Tau);   

  //  fRbg1Slope = new RooRealVar("bg1Slope", "coefficient #0 for bg 1", fBgMc0, -10., 10.); 
  //  RooChebychev bg1M("bg1M", "background 1 gamma gamma mass", *fRm, RooArgList(*fRbg1Slope)); 
  fRbg1Slope = new RooRealVar("bg1Slope", "coefficient #0 for bg 1", fBgMp1, -10., 10.); 
  RooPolynomial bg1M("bg1M", "background 1 gamma gamma mass", *fRm, RooArgList(*fRbg1Slope)); 
  RooExponential bg1Pt("bg1Pt", "background 1 gamma gamma pT", *fRpt, *fRbg1Tau);   


  // -- model 3 setup: manual pLLR for Sig1
  fRsg3N       = new RooRealVar("sg3N","Number of Higgs 3 signal events", fSg1, 0., 1.e4);
  fRsg3Tau     = new RooRealVar("sg3Tau", "signal 3 tau", fSg1Tau, -10., 10.); 
  fRbg3N       = new RooRealVar("bg3N","Number of background 3  events", fBg, 0., 1.e4);
  fRbg3Slope   = new RooRealVar("bg3Slope", "coefficient #0 for bg 3", fBgMp1, -10., 10.); 
  fRbg3Tau     = new RooRealVar("bg3Tau", "background 3 gamma gamma tau", fBgTau, -10., 10.); 
  RooPolynomial  bg3M("bg3M", "background 3 gamma gamma mass", *fRm, RooArgList(*fRbg3Slope)); 
  RooExponential bg3Pt("bg3Pt", "background 3 gamma gamma pT", *fRpt, *fRbg3Tau);   
  RooGaussian    sg3M("sg3M", "signal 3 mass", *fRm, *fRsgP, *fRsgS); 
  RooExponential sg3Pt("sg3Pt", "signal 3 pT", *fRpt, *fRsg3Tau);   
  RooProdPdf     sg3Pdf("sg3Pdf", "sg3Pdf", RooArgSet(sg3M, sg3Pt));
  RooProdPdf     bg3Pdf("bg3Pdf", "bg3Pdf", RooArgSet(bg3M, bg3Pt));
  RooAddPdf      model3("model3", "model3", RooArgList(sg3Pdf, bg3Pdf), RooArgList(*fRsg3N, *fRbg3N));

  // -- model 4 setup: manual pLLR for Sig0
  fRsg4N       = new RooRealVar("sg4N","Number of Higgs 4 signal events", fSg1, 0., 1.e4);
  fRsg4Tau     = new RooRealVar("sg4Tau", "signal 4 tau", fSg1Tau, -10., 10.); 
  fRbg4N       = new RooRealVar("bg4N","Number of background 4  events", fBg, 0., 1.e4);
  fRbg4Slope   = new RooRealVar("bg4Slope", "coefficient #0 for bg 4", fBgMp1, -10., 10.); 
  fRbg4Tau     = new RooRealVar("bg4Tau", "background 4 gamma gamma tau", fBgTau, -10., 10.); 
  RooPolynomial  bg4M("bg4M", "background 4 gamma gamma mass", *fRm, RooArgList(*fRbg4Slope)); 
  RooExponential bg4Pt("bg4Pt", "background 4 gamma gamma pT", *fRpt, *fRbg4Tau);   
  RooGaussian    sg4M("sg4M", "signal 4 mass", *fRm, *fRsgP, *fRsgS); 
  RooExponential sg4Pt("sg4Pt", "signal 4 pT", *fRpt, *fRsg4Tau);   
  RooProdPdf     sg4Pdf("sg4Pdf", "sg4Pdf", RooArgSet(sg4M, sg4Pt));
  RooProdPdf     bg4Pdf("bg4Pdf", "bg4Pdf", RooArgSet(bg4M, bg4Pt));
  RooAddPdf      model4("model4", "model4", RooArgList(sg4Pdf, bg4Pdf), RooArgList(*fRsg4N, *fRbg4N));
  

//   RooChebychev bgM("bgM", "background gamma gamma mass", *fRm, RooArgList(*fRbgSlope)); 
//   RooExponential bgPt("bgPt", "background gamma gamma pT", *fRpt, *fRbgTau);   

  // -- parametrization in terms of mu
//   fRmu   = new RooRealVar("mu", "mu", fMu, -1., 10.);
//   fRsgN   = new RooFormulaVar("sgN", "sgN", "@0 + @1*(@2-@0)", RooArgList(*fRsg0N, *fRmu, *fRsg1N)); 
//   fRsgTau = new RooFormulaVar("sgTau", "sgTau", "@0 + @1*(@2-@0)", RooArgList(*fRsg0Tau, *fRmu, *fRsg1Tau)); 
//   RooGaussian sgM("sgM", "signal mass", *fRm, *fRsgP, *fRsgS); 
//   RooExponential sgPt("sgPt", "signal pT", *fRpt, *fRsgTau);   

  // -- signal and background 2D pdfs
  RooProdPdf sg0Pdf = RooProdPdf("sg0Pdf", "sg0Pdf", RooArgSet(sg0M, sg0Pt));
  RooProdPdf sg1Pdf = RooProdPdf("sg1Pdf", "sg1Pdf", RooArgSet(sg1M, sg1Pt));
//   RooProdPdf sgPdf  = RooProdPdf("sgPdf",  "sgPdf",  RooArgSet(sgM,  sgPt));

  RooProdPdf bg0Pdf = RooProdPdf("bg0Pdf", "bg0Pdf", RooArgSet(bg0M, bg0Pt));
  RooProdPdf bg1Pdf = RooProdPdf("bg1Pdf", "bg1Pdf", RooArgSet(bg1M, bg1Pt));
//   RooProdPdf bgPdf  = RooProdPdf("bgPdf",  "bgPdf",  RooArgSet(bgM,  bgPt));

//   RooAddPdf model("model", "model", RooArgList(sgPdf, bgPdf), RooArgList(*fRsgN, *fRbgN));
  RooAddPdf model0("model0", "model0", RooArgList(sg0Pdf, bg0Pdf), RooArgList(*fRsg0N, *fRbg0N));
  RooAddPdf model1("model1", "model1", RooArgList(sg1Pdf, bg1Pdf), RooArgList(*fRsg1N, *fRbg1N));
  
  // -- 1D mass fit only
  RooAddPdf model2("model2", "model2", RooArgList(sg2M, bg2M), RooArgList(*fRsg2N, *fRbg2N));



  // -- run toys 
  if (ntoy > 0) fNtoy = ntoy; 
  RooDataSet *bgData(0), *sg0Data(0), *sg1Data(0), *data(0), *data1(0), *data0(0);
  RooFitResult *fr1(0), *fr0(0);
  
  RooCmdArg fitargs; 
  if (1) {
    fitargs.addArg(Minos(kTRUE)); 
    fitargs.addArg(PrintLevel(-1));
    fitargs.addArg(PrintEvalErrors(-1));
    fitargs.addArg(Verbose(kFALSE));
    fitargs.addArg(Save()); 
    fitargs.addArg(Minimizer("Minuit2")); 
  }

  RooWorkspace *w(0); 

  RooRandom::randomGenerator()->SetSeed(111);

  TH1D *hsig1 = new TH1D("hsig1", "", 50, 0., 5.); hsig1->SetLineColor(kBlack); 
  TH1D *hsig2 = new TH1D("hsig2", "", 50, 0., 5.); hsig2->SetLineColor(kRed); 
  TH1D *hsig3 = new TH1D("hsig3", "", 50, 0., 5.); hsig3->SetLineColor(kMagenta); 

  TH1D *ht0 = new TH1D("ht0", "prof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 400, -40., 40.); ht0->SetLineColor(kBlue); 
  TH1D *ht1 = new TH1D("ht1", "prof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 400, -40., 40.); ht1->SetLineColor(kRed); 

  TH1D *hT0 = new TH1D("hT0", "unprof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 400, -40., 40.); hT0->SetLineColor(kBlue); 
  TH1D *hT1 = new TH1D("hT1", "unprof t = -2ln(lambda), lambda = L(t, v^^)/L(t^, v^)", 400, -40., 40.); hT1->SetLineColor(kRed); 


  for (int i = 0; i < fNtoy; ++i) {
    cout << "==================> TOY " << i << endl;
    w = new RooWorkspace("w"); 

    iniRooVars(); 
    bgData = bg0Pdf.generate(RooArgSet(*fRm, *fRpt), fBg); 
    sg0Data = sg0Pdf.generate(RooArgSet(*fRm, *fRpt), fSg0);
    sg1Data = sg0Pdf.generate(RooArgSet(*fRm, *fRpt), fSg1);

    data = new RooDataSet(*bgData);
    data->append(*sg1Data);

    data1 = new RooDataSet(*bgData);
    data1->append(*sg1Data);

    data0 = new RooDataSet(*bgData);
    data0->append(*sg0Data);

    // -----------------------
    // -- mu parametrization?!
    // -----------------------
    //     fRsg0N->setConstant(kTRUE);
    //     fRsg1N->setConstant(kTRUE);
    //     fRsg0Tau->setConstant(kTRUE);
    //     fRsg1Tau->setConstant(kTRUE);
    
    //     RooAbsPdf * pdf = &model;
    //     RooFitResult *r0 = pdf->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
    //     cout << "### RooFitResult r0 ###################################################################" << endl;
    //     r0->Print("v");
    //     cout << "#######################################################################################" << endl;
    
    //     // -- plot 0
    //     c0->cd(6);
    //     RooPlot *plotM = fRm->frame(Title("mass"), Bins(NBINS));
    //     data->plotOn(plotM);
    //     pdf->plotOn(plotM);
    //     pdf->paramOn(plotM, Layout(0.5,0.9,0.85));
    //     plotM->Draw();
    
    
    //     c0->cd(7);
    //     RooPlot *plotPt = fRpt->frame(Title("pt"));
    //     data->plotOn(plotPt);
    //     pdf->plotOn(plotPt);
    //     pdf->paramOn(plotPt, Layout(0.5,0.9,0.85));
    //     plotPt->Draw();
    

    // ----------------------------
    // -- "2D" mass-lifetime
    // ----------------------------
    iniRooVars(); 
    RooFitResult *r1 = model1.fitTo(*data1, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
    cout << "### RooFitResult r1 ###################################################################" << endl;
    r1->Print("v");
    cout << "#######################################################################################" << endl;
    c0->cd(2);
    if (0 == i) {
      RooPlot *plotM1 = fRm->frame(Title("mass"), Bins(NBINS));
      plotM1->remove("", kTRUE);
      data1->plotOn(plotM1);
      model1.plotOn(plotM1);
      model1.paramOn(plotM1, Layout(0.5, 0.8, 0.45));
      plotM1->Draw();
    }
    
    c0->cd(3);
    if (0 == i) {
      RooPlot *plotPt1 = fRpt->frame(Title("pt"));
      data1->plotOn(plotPt1);
      model1.plotOn(plotPt1);
      model1.paramOn(plotPt1, Layout(0.4,0.8,0.85));
      plotPt1->Draw();
    }

    RooStats::ModelConfig mc1("mc1", w);
    mc1.SetPdf(model1);
    mc1.SetObservables(RooArgSet(*fRm, *fRpt));
    mc1.SetParametersOfInterest(*fRsg1N);
    RooArgList nuisParams(*fRbg1Tau, *fRbg1N, *fRbg1Slope); 
    mc1.SetNuisanceParameters(nuisParams);


    //     // -- plot 1: PLL
    //     if (0) {
    //       ProfileLikelihoodCalculator pl(*data1, mc1);
    //       pl.SetConfidenceLevel(0.683); // 68% interval
    //       LikelihoodInterval* interval0 = pl.GetInterval();
    
    //       // find the iterval on the first Parameter of Interest
    //       RooRealVar* firstPOI = (RooRealVar*) mc1.GetParametersOfInterest()->first();
    
    //       double lowerLimit = interval0->LowerLimit(*firstPOI);
    //       double upperLimit = interval0->UpperLimit(*firstPOI);
    
    
    //       cout << "### 68% interval on " << firstPOI->GetName()<< " ###################################" << endl
    // 	   << "### [" << lowerLimit << ", " << upperLimit << "] "
    // 	   << endl;
    
    //       c0->cd(5);
    //       LikelihoodIntervalPlot *liplot = new LikelihoodIntervalPlot(interval0);
    //       liplot->SetRange(50., 200.); 
    //       liplot->Draw("");  
    //     }

    // ----------------------------------------------------------------------
    // -- Asymptotic significance: works for one POI only!
    // ----------------------------------------------------------------------
    ModelConfig *sbModel = (ModelConfig*)mc1.Clone();
    sbModel->SetName("S+B Model");      
    RooRealVar *poi = (RooRealVar*)sbModel->GetParametersOfInterest()->first();
    poi->setVal(fSg1);  // set POI snapshot in S+B model for expected significance
    sbModel->SetSnapshot(*poi);
    
    ModelConfig *bModel = (ModelConfig*)mc1.Clone();
    bModel->SetName("B Model");      
    poi->setVal(fSg0);
    bModel->SetSnapshot(*poi);

    AsymptoticCalculator  ac(*data1, *sbModel, *bModel);
    //T3    ac.SetOneSidedDiscovery(true);  // for one-side discovery test
    
    HypoTestResult *asResult = ac.GetHypoTest();
    cout << "-------------------------------------------------" << endl;
    asResult->Print();
    cout << "-------------------------------------------------" << endl;
    hsig1->Fill(asResult->Significance()); 
    

    // ----------------------------------------------------------------------
    // -- profile likelihood, from tutorials/roostats/rs102_hypotestwithshapes.C
    //    https://root.cern.ch/root/html/RooStats__ProfileLikelihoodCalculator.html
    // ----------------------------------------------------------------------
    ModelConfig mcpll; 
    mcpll.SetWorkspace(*w);
    mcpll.SetPdf(model1);
    ProfileLikelihoodCalculator plc; 
    plc.SetData(*data1); 
    
    RooArgSet ras_poi(*fRsg1N);
    //    RooArgSet ras_poi(*fRsg1N, *fRsg1Tau);
    RooArgSet *nullParams = (RooArgSet*) ras_poi.snapshot(); 
    nullParams->setRealValue("sg1N", fSg0); 
    //    nullParams->setRealValue("sg1Tau", fSg0Tau); 
    
    plc.SetModel(mcpll);
    plc.SetNullParameters(*nullParams);

    HypoTestResult* htr = plc.GetHypoTest();
    cout << "-------------------------------------------------" << endl;
    cout << "The p-value for the null is " << htr->NullPValue() << endl;
    cout << "Corresponding to a significance of " << htr->Significance() << endl;
    cout << "-------------------------------------------------\n\n" << endl;
    hsig2->Fill(htr->Significance()); 


    // ----------------------------------------------------------------------
    // -- profile likelihood, attempt with 2POI, but results in low significance?!
    //    https://root.cern.ch/root/html/RooStats__ProfileLikelihoodCalculator.html
    // ----------------------------------------------------------------------
    ModelConfig mcpll2d; 
    mcpll2d.SetWorkspace(*w);
    mcpll2d.SetPdf(model1);
    ProfileLikelihoodCalculator plc2d; 
    plc2d.SetData(*data1); 
    
    //    RooArgSet ras_poi(*fRsg1N);
    RooArgSet poi2d(*fRsg1N, *fRsg1Tau);
    RooArgSet *nullParams2d = (RooArgSet*) poi2d.snapshot(); 
    nullParams2d->setRealValue("sg1N", fSg0); 
    nullParams2d->setRealValue("sg1Tau", fSg0Tau); 
    
    plc2d.SetModel(mcpll2d);
    plc2d.SetNullParameters(*nullParams2d);

    HypoTestResult* htr2d = plc2d.GetHypoTest();
    cout << "-------------------------------------------------" << endl;
    cout << "The p-value for the null is " << htr2d->NullPValue() << endl;
    cout << "Corresponding to a significance of " << htr2d->Significance() << endl;
    cout << "-------------------------------------------------\n\n" << endl;
    hsig3->Fill(htr2d->Significance()); 


    // ----------------------------------------------------------------------
    // -- manual implementation of pLLR
    // ----------------------------------------------------------------------
    // http://www.hep.manchester.ac.uk/u/ukyang/chiggs/codes/Limit.cc ??


    // ----------------------------------------------------------------------
    // -- and now for SIG0
    // ----------------------------------------------------------------------
    RooArgSet m4poi(*fRsg4N, *fRsg4Tau); 
    
    iniRooVars(); 
    fRsg4N->setConstant(kFALSE);
    fRsg4Tau->setConstant(kFALSE);
    cout << "xxxxx > create nll" << endl;
    RooAbsReal *nll0 = model4.createNLL(*data0, Extended(kTRUE), Verbose(kFALSE));
    RooMinuit tminuit0(*nll0);
    tminuit0.setPrintLevel(-1);
    tminuit0.setNoWarn();
    tminuit0.migrad();
    RooFitResult* tfitres0 = tminuit0.save();  

    double fitRsg4N   = fRsg4N->getVal();
    double fitRsg4Tau = fRsg4Tau->getVal();

    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "XXXXXXXXXXXXXX tfitres0" << endl;
    tfitres0->Print(); 
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

    double valNll0 = nll0->getVal(); 

    cout << endl;
    cout << "xxxxx > valNll0 = " << valNll0 << endl;

    RooProfileLL* pll0 = dynamic_cast<RooProfileLL*>(nll0->createProfile(m4poi));
    
    // -- const-set SG0-hypothesis
    fRsg4N->setVal(fSg0); 
    fRsg4N->setConstant(kTRUE);
    fRsg4Tau->setVal(fSg0Tau); 
    fRsg4Tau->setConstant(kTRUE);
    double val4Nll0   = nll0->getVal();
    double val4Pll0   = pll0->getVal();

    double fit2Rsg4N = fRsg4N->getVal();

    cout << endl;
    cout << "xxxxx > val4Nll0 = " << val4Nll0 << endl;
    cout << "xxxxx > val4Pll0 = " << val4Pll0 << endl;

    // -- const-set SG1-hypothesis
    fRsg4N->setVal(fSg1); 
    fRsg4N->setConstant(kTRUE);
    fRsg4Tau->setVal(fSg1Tau); 
    fRsg4Tau->setConstant(kTRUE);
    double val4Nll1   = nll0->getVal();
    double val4Pll1   = pll0->getVal();

    cout << endl;
    cout << "xxxxx > val4Nll1 = " << val4Nll1 << endl;
    cout << "xxxxx > val4Pll1 = " << val4Pll1 << endl;

    ht0->Fill(2.*(val4Pll1 - val4Pll0));    
    hT0->Fill(2.*(val4Nll1 - val4Nll0));

    cout << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << " val4Nll0 = " << val4Nll0 << " val4Pll0 = " << val4Pll0 << " before const-setting null-hypo parameters Nsig = " << fitRsg4N << endl
	 << " val4Nll1 = " << val4Nll1 << " val4Pll1 = " << val4Pll1 << " after const-setting null-hypo parameters Nsig = " << fit2Rsg4N << endl
	 << "note: fSg0 = " << fSg0 
	 << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << endl;
    
    // -- reset
    fRsg4N->setConstant(kFALSE);
    fRsg4Tau->setConstant(kFALSE);


    // ----------------------------------------------------------------------
    // -- and now for SIG1
    // ----------------------------------------------------------------------
    RooArgSet m3poi(*fRsg3N, *fRsg3Tau); 
    
    cout << "xxxxx > create nll" << endl;
    iniRooVars(); 
    fRsg3N->setConstant(kFALSE);
    fRsg3Tau->setConstant(kFALSE);
    RooAbsReal *nll1 = model3.createNLL(*data1, Extended(kTRUE), Verbose(kFALSE));
    RooMinuit tminuit(*nll1);
    tminuit.setPrintLevel(-1);
    tminuit.setNoWarn();
    tminuit.migrad();
    RooFitResult* tfitres = tminuit.save();  

    double fitRsg3N   = fRsg3N->getVal();
    double fitRsg3Tau = fRsg3Tau->getVal();

    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "XXXXXXXXXXXXXX tfitres" << endl;
    tfitres->Print(); 
    cout << "fRsg3N->getVal() = " << fRsg3N->getVal()  << endl;
    cout << "fRsg3Tau->getVal() = " << fRsg3Tau->getVal()  << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

    double valNll1 = nll1->getVal(); 

    cout << endl;
    cout << "xxxxx > valNll1 = " << valNll1 << endl;
    // https://root.cern.ch/root/html/RooProfileLL.html:
    // The value return by RooProfileLL is the input likelihood nll minimized w.r.t all nuisance parameters 
    // (which are all parameters except for those listed in the constructor) minus the -log(L) of the best fit.
    // -> i.e. the return value is 0.5 * test statistic t
    RooProfileLL* pll1 = dynamic_cast<RooProfileLL*>(nll1->createProfile(m3poi));
    
    double val3Nll = nll1->getVal(); 
    double val3Pll = pll1->getVal(); 
    cout << endl;
    cout << "xxxxx > val3Nll = " << val3Nll << endl;
    cout << "xxxxx > val3Pll = " << val3Pll << endl;

    // -- const-set SG0-hypothesis
    fRsg3N->setVal(fSg0); 
    fRsg3N->setConstant(kTRUE);
    fRsg3Tau->setVal(fSg0Tau); 
    fRsg3Tau->setConstant(kTRUE);
    double val3Nll0   = nll1->getVal(); 
    double val3Pll0   = pll1->getVal();

    cout << endl;
    cout << "xxxxx > val3Nll0 = " << val3Nll0 << endl;
    cout << "xxxxx > val3Pll0 = " << val3Pll0 << endl;


    // -- const-set SG1-hypothesis
    fRsg3N->setVal(fSg1); 
    fRsg3N->setConstant(kTRUE);
    fRsg3Tau->setVal(fSg1Tau); 
    fRsg3Tau->setConstant(kTRUE);
    double val3Nll1   = nll1->getVal(); 
    double val3Pll1   = pll1->getVal();

    cout << endl;
    cout << "xxxxx > val3Nll1 = " << val3Nll1 << endl;
    cout << "xxxxx > val3Pll1 = " << val3Pll1 << endl;

    double fit2Rsg3N = fRsg3N->getVal();
    ht1->Fill(2.*(val3Pll1 - val3Pll0));
    hT1->Fill(2.*(val3Nll1 - val3Nll0));
    c0->cd(7);
    RooPlot *frame2 = fRsg3N->frame(Name("plot2"), Range(0., 200.));
    nll1->plotOn(frame2, LineColor(kBlack), ShiftToZero());
    pll1->plotOn(frame2, LineColor(kBlue), ShiftToZero());
    frame2->Draw();
    
    // -- reset
    fRsg3N->setConstant(kFALSE);
    fRsg3Tau->setConstant(kFALSE);

    cout << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << " val3Nll0 = " << val3Nll0 << " val3Pll0 = " << val3Pll0 << " before const-setting null-hypo parameters Nsig = " << fitRsg3N << endl
	 << " val3Nll1 = " << val3Nll1 << " val3Pll1 = " << val3Pll1 << " after const-setting null-hypo parameters Nsig = " << fit2Rsg3N  << endl
	 << "note: fSg0 = " << fSg0 
	 << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << endl;



    // ----------------------------------------------------------------------
    // -- this is consistent with the other tests 
    //    (AsymptoticCalculator, ProfileLikelihoodCalculator),
    //    but takes much more time
    // ----------------------------------------------------------------------
    if (0) {
      FrequentistCalculator   fc(*data1, *sbModel, *bModel);
      fc.SetToys(2000,500);    // 2000 for null (B) and 500 for alt (S+B) 
      
      // create the test statistics
      ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
      // use one-sided profile likelihood
      profll.SetOneSidedDiscovery(true);
      
      // configure  ToyMCSampler and set the test statistics
      ToyMCSampler *toymcs = (ToyMCSampler*)fc.GetTestStatSampler();
      toymcs->SetTestStatistic(&profll);
      
      if (!sbModel->GetPdf()->canBeExtended())
	toymcs->SetNEventsPerToy(1);
      
      // run the test
      HypoTestResult * fqResult = fc.GetHypoTest();
      fqResult->Print();
      
      // plot test statistic distributions
      c0->cd(3);
      HypoTestPlot *hplot = new HypoTestPlot(*fqResult);
      hplot->SetLogYaxis(true);
      hplot->Draw();
    }



    delete w;
    delete r1; 

    delete pll1; 
    delete pll0; 

    delete data; 
    delete data1; 
    delete data0; 

    delete bgData;
    delete sg0Data;
    delete sg1Data; 

    delete nullParams; 
    delete bModel; 
    delete sbModel; 

    delete asResult; 
    delete htr; 

    c0->cd(5);
    gStyle->SetOptStat(0); 
    hsig2->Draw();
    hsig1->Draw("same");
    
    tl->SetTextSize(0.05);
    tl->SetTextColor(kBlack); 
    tl->DrawLatex(0.10, 0.92, Form("Asymp: %3.2f", median(hsig1))); 
    tl->SetTextColor(kRed); 
    tl->DrawLatex(0.4, 0.92, Form("ProfL: %3.2f", median(hsig2))); 
    tl->SetTextColor(kBlack); 
    tl->DrawLatex(0.67, 0.92, Form("ntoys: %4d", i+1)); 

    hsig3->Draw("same");
    tl->SetTextColor(kMagenta); 
    tl->DrawLatex(0.55, 0.82, Form("ProfL2d: %3.2f", median(hsig3))); 

    c0->cd(6); 
    gStyle->SetOptStat(0); 
    if (ht0->GetMaximum() > ht1->GetMaximum()) {
      ht0->Draw();
      ht1->Draw("same");
    } else {
      ht1->Draw();
      ht0->Draw("same");
    }
    double midpoint, tailprob; 
    findMidPoint(ht0, ht1, midpoint, tailprob); 
    double sigma = 2.*oneSidedGaussianSigma(tailprob);
    tl->SetTextSize(0.04);
    tl->DrawLatex(0.6, 0.80, Form("midpoint: %4.3f", midpoint)); 
    tl->DrawLatex(0.6, 0.70, Form("tailprob: %4.3f", tailprob)); 
    tl->DrawLatex(0.6, 0.60, Form("Sep: %4.3f", sigma)); 


    tl->SetTextSize(0.05);
    tl->SetTextColor(kRed); 
    tl->DrawLatex(0.10, 0.92, Form("pLLR1: %3.2f", median(ht1))); 

    tl->SetTextSize(0.05);
    tl->SetTextColor(kBlue); 
    tl->DrawLatex(0.50, 0.92, Form("pLLR0: %3.2f", median(ht0))); 


    c0->cd(9); 
    gStyle->SetOptStat(0); 
    if (hT0->GetMaximum() > hT1->GetMaximum()) {
      hT0->Draw();
      hT1->Draw("same");
    } else {
      hT1->Draw();
      hT0->Draw("same");
    }
    findMidPoint(hT0, hT1, midpoint, tailprob); 
    sigma = 2.*oneSidedGaussianSigma(tailprob);
    tl->SetTextSize(0.04);
    tl->DrawLatex(0.6, 0.80, Form("midpoint: %4.3f", midpoint)); 
    tl->DrawLatex(0.6, 0.70, Form("tailprob: %4.3f", tailprob)); 
    tl->DrawLatex(0.6, 0.60, Form("Sep: %4.3f", sigma)); 


    tl->SetTextSize(0.05);
    tl->SetTextColor(kRed); 
    tl->DrawLatex(0.10, 0.92, Form("LLR1: %3.2f", median(hT1))); 

    tl->SetTextSize(0.05);
    tl->SetTextColor(kBlue); 
    tl->DrawLatex(0.50, 0.92, Form("LLR0: %3.2f", median(hT0))); 


    c0->SaveAs(Form("%s/allNumbers5-toy-%d.pdf", fDirectory.c_str(), i)); 

  }


  // -- 1d hypothesis test with RooDLLSignificanceMCSModule
  if (0) {
    c0->cd(7);
    iniRooVars();
    RooDataSet *data2 = model2.generate(*fRm, fBg + fSg1); 
    RooMCStudy* mcstudy = new RooMCStudy(model2, *fRm, Binned(), Silence(), Extended(kTRUE),
					 FitOptions(Save(kTRUE), Extended(kTRUE), PrintEvalErrors(-1))
					 );
    double nullHypo = fSg0; 
    RooDLLSignificanceMCSModule sigModule(*fRsg2N, nullHypo);
    mcstudy->addModule(sigModule);
    mcstudy->generateAndFit(1000);
    
    TH1 *hz = mcstudy->fitParDataSet().createHistogram("significance_nullhypo_sg2N") ;
    hz->Draw("hist");
    tl->SetTextColor(kBlack); 
    tl->DrawLatex(0.2, 0.92, Form("RooDLLSig: %3.2f", median(hz))); 
  }




  c0->SaveAs(Form("%s/allNumbers5-%s.pdf", fDirectory.c_str(), fSetup.c_str())); 
  
  fTEX.open(fTexFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << formatTex(fg0pt, Form("%s:g0pt:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fg1pt, Form("%s:g1pt:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fptlo, Form("%s:ptlo:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fpthi, Form("%s:pthi:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fptnl, Form("%s:ptnl:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fptnh, Form("%s:ptnh:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fg0ptlo, Form("%s:g0ptlo:val", fSetup.c_str()), 2) << endl;
  fTEX << formatTex(fg1ptlo, Form("%s:g1ptlo:val", fSetup.c_str()), 2) << endl;
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
  fTEX.close();

} 



