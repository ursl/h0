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

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooDataSet.h"

#include "RooAbsPdf.h"
#include "RooRandom.h"
#include "RooFitResult.h"

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

  NBINS = 25;  
  GETA = 2.5;
  G0ISO = G1ISO = 0.2; 
  G0ISOR = 0.005; 
  G1ISOR = 0.010; 
  G0PT = 100.;
  G1PT = 40.;
  PTLO = 250.;
  PTHI = 999.;
  MGGLO = 100.;
  MGGHI = 150.;

  fLumi = 1000.;

  // test
  G0PT = 140.;
  G1PT = 100.;

  G0PT = 0.;
  G1PT = 40.;

}


// ----------------------------------------------------------------------
plotHpt::~plotHpt() {

}


// ----------------------------------------------------------------------
void plotHpt::bookHist(string name, string cuts) {
  char hist[200], thist[200], ahist[200];

  // -- gen pT for efficiency
  sprintf(hist, "gen%s_%s_%s", "pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "pt");
  sprintf(ahist, "%s", "pt_{#gamma#gamma} [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  // -- reco overlays
  sprintf(hist, "%s_%s_%s", "m", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "m");
  sprintf(ahist, "%s", "m_{#gamma#gamma} [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, NBINS, MGGLO, MGGHI))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "pt");
  sprintf(ahist, "%s", "pt_{#gamma#gamma} [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0, 1000.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  
  sprintf(hist, "%s_%s_%s", "eta", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "eta");
  sprintf(ahist, "%s", "eta_{#gamma#gamma}"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, -5., 5.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");
  
  sprintf(hist, "%s_%s_%s", "g0pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0pt");
  sprintf(ahist, "%s", "pt_{#gamma0} [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 800.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g1pt", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1pt");
  sprintf(ahist, "%s", "pt_{#gamma1} [GeV]"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 300.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g0iso", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0iso");
  sprintf(ahist, "%s", "I_{#gamma0}"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 2.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g1iso", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1iso");
  sprintf(ahist, "%s", "I_{#gamma1}"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 2.))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g0isor", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g0isor");
  sprintf(ahist, "%s", "I_{#gamma0}^{rel}"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 0.006))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

  sprintf(hist, "%s_%s_%s", "g1isor", name.c_str(), cuts.c_str());
  sprintf(thist, "%s", "g1isor");
  sprintf(ahist, "%s", "I_{#gamma1}^{rel}"); 
  fHists.insert(make_pair(hist, new TH1D(hist, thist, 100, 0., 0.012))); 
  setHistTitles(fHists[hist], fDS[name], ahist, "Entries/bin");

}


// ----------------------------------------------------------------------
void plotHpt::makeAll(int bitmask) {
  if (bitmask & 0x1) treeAnalysis();
  if (bitmask & 0x2) toy1(); 
}


// ----------------------------------------------------------------------
void plotHpt::treeAnalysis(int nevts) {
  string ds("sherpa");
  fCds = ds; 
  bookHist(ds, "goodcand"); 
  bookHist(ds, "lopt"); 
  bookHist(ds, "hipt"); 
  TTree *t = getTree(ds); 
  setupTree(t); 
  loopOverTree(t, 1, nevts); 
  

  ds = "mcatnlo5"; 
  fCds = ds; 
  bookHist(ds, "goodcand"); 
  bookHist(ds, "lopt"); 
  bookHist(ds, "hipt"); 
  t = getTree(ds); 
  setupTree(t);
  loopOverTree(t, 1, nevts); 


  ds = "mcatnlo0"; 
  fCds = ds; 
  bookHist(ds, "goodcand"); 
  bookHist(ds, "lopt"); 
  bookHist(ds, "hipt"); 
  t = getTree(ds); 
  setupTree(t);
  loopOverTree(t, 1, nevts); 

  ds = "mcatnlo1"; 
  fCds = ds; 
  bookHist(ds, "goodcand"); 
  bookHist(ds, "lopt"); 
  bookHist(ds, "hipt"); 
  t = getTree(ds); 
  setupTree(t);
  loopOverTree(t, 1, nevts); 

  map<string, TH1*>::iterator hit = fHists.begin();
  map<string, TH1*>::iterator hite = fHists.end();
  string h0name(""), h1name("");
  TH1D *h1(0); 
  for (; hit != hite; ++hit) {
    h0name = hit->second->GetName(); 
    if (string::npos == h0name.compare("mcatnlo0")) continue;
    h1name = h0name;
    replaceAll(h1name, "mcatnlo0", "mcatnlo1"); 
    h1 = (TH1D*)fHists[h1name];
    combMCAtNLOHist((TH1D*)hit->second, h1); 
  }

  makeCanvas(8);
  zone(2,2, c4);
  TH1D *h0 = (TH1D*)fHists["genpt_mcatnlo5_goodcand"]->DrawCopy();
  c4->cd(2);
  h1 = (TH1D*)fHists["pt_mcatnlo5_goodcand"]->DrawCopy();
  c4->cd(3);
  TH1D *heff = (TH1D*)h0->Clone("heff");
  heff->Reset();
  heff->Divide(h1, h0); 
  heff->Draw();


  gStyle->SetOptStat(0); 
  string what("pt"); 
  string sel("lopt");
  int ipad(0); 
  int OTYPE = LUMI;


  zone(4,2);

  ++ipad;
  c0->cd(ipad);
  what = "m"; 
  sel = "hipt";
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
	  OTYPE, false, true); 

  ++ipad;
  c0->cd(ipad);
  what = "pt"; 
  sel = "goodcand";
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
	  OTYPE, true, true); 


  ++ipad; 
  c0->cd(ipad);
  what = "pt"; 
  sel = "lopt";
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
	  OTYPE, false, true); 

  ++ipad;
  c0->cd(ipad);
  what = "pt"; 
  sel = "hipt";
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
	  OTYPE, true, true); 

  ++ipad;
  c0->cd(ipad);
  what = "g0pt"; 
  sel = "lopt";
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
	  OTYPE, true, true); 


  ++ipad;
  c0->cd(ipad);
  what = "g0pt"; 
  sel = "hipt";
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
	  OTYPE, true, true); 

  ++ipad;
  c0->cd(ipad);
  what = "g1pt"; 
  sel = "lopt";
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
	  OTYPE, true, true); 

  ++ipad;
  c0->cd(ipad);
  what = "g1pt"; 
  sel = "hipt";
  overlay(fHists[Form("%s_sherpa_%s", what.c_str(), sel.c_str())], "sherpa", 
	  fHists[Form("%s_mcatnlo5_%s", what.c_str(), sel.c_str())], "mcatnlo5", 
	  fHists[Form("%s_mcatnlo_%s", what.c_str(), sel.c_str())], "mcatnlo0", 
	  OTYPE, true, true); 



}

// ----------------------------------------------------------------------
void plotHpt::toy1(int nsg, int nbg) {

  cout << "hello world" << endl;

  RooRandom::randomGenerator()->SetSeed(1);

  RooWorkspace *wspace = new RooWorkspace(); 
  wspace->factory("Gaussian::normal(x[-10., 10.], mu[-1., 1.], sigma[1])"); 
  wspace->defineSet("poi", "mu"); 
  wspace->defineSet("obs", "x"); 
  
  ModelConfig *modelConfig = new ModelConfig("G(x|mu,1)"); 
  modelConfig->SetWorkspace(*wspace); 
  modelConfig->SetPdf(*wspace->pdf("normal"));
  modelConfig->SetParametersOfInterest(*wspace->set("poi"));
  modelConfig->SetObservables(*wspace->set("obs")); 

  RooDataSet *data = wspace->pdf("normal")->generate(*wspace->set("obs"), 100); 

  double confidenceLevel(0.95); 

  // for convenience later on 
  RooRealVar* x = wspace->var("x");
  RooRealVar* mu = wspace->var("mu");

  ProfileLikelihoodCalculator plc(*data, *modelConfig); 
  plc.SetConfidenceLevel(confidenceLevel); 
  LikelihoodInterval *plInt = plc.GetInterval(); 

  cout << "plInt: " << plInt->LowerLimit(*mu) << " .. " << plInt->UpperLimit(*mu) << endl;
  mu->setVal(0.); 
  cout << "mu = 0? " << plInt->IsInInterval(*mu) << endl;

  RooPlot *plot = x->frame(Title("blah"));
  data->plotOn(plot);
  plot->Draw();
 

  
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
	  sdecay = "#gamma #gamma";
	  ds->fColor = kBlue; 
	  ds->fLcolor = kBlue; 
	  ds->fFcolor = kBlue; 
	  ds->fSymbol = 24; 
	} else if (string::npos != stype.find("1")) {
	  sname = "mcatnlo1";
	  sdecay = "#gamma #gamma";
	  ds->fColor = kBlue+2; 
	  ds->fLcolor = kBlue+2; 
	  ds->fFcolor = kBlue+2; 
	  ds->fSymbol = 25; 
	} else if (string::npos != stype.find("5")) {
	  sname = "mcatnlo5";
	  sdecay = "#gamma #gamma";
	  ds->fColor = kBlack; 
	  ds->fLcolor = kBlack; 
	  ds->fFcolor = kBlack; 
	  ds->fSymbol = 26; 
	} 
	ds->fF = pF; 
	ds->fXsec = atof(sxsec.c_str());
	ds->fBf   = 2.28E-03;
	ds->fLumi = nevt/ds->fXsec/ds->fBf/1000.;
	ds->fName = "MC@NLO " + sdecay; 
	ds->fFillStyle = 3356; 
	ds->fSize = 1; 
	ds->fWidth = 2; 
	fDS.insert(make_pair(sname, ds)); 
      }


      if (string::npos != stype.find("sherpa")) {
	dataset *ds = new dataset(); 
	sname = "sherpa";
	sdecay = "#gamma #gamma";
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
      cout << "opened MC file "  << sfile  << " as " << sname << " (" << stype << ") with xsec = " << sxsec
	   << " equivalent lumi = " << fDS[sname]->fLumi << "/fb"
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
  if (fb.g0pt < G0PT) fGoodCand = false; 
  if (fb.g1pt < G1PT) fGoodCand = false; 

  //   if (fb.g0iso/fb.g0pt > G0ISOR) fGoodCand = false; 
  //   if (fb.g1iso/fb.g1pt > G1ISOR) fGoodCand = false; 

  if (fb.g0iso > G0ISO) fGoodCand = false; 
  if (fb.g1iso > G1ISO) fGoodCand = false; 

}

// ----------------------------------------------------------------------
void plotHpt::loopFunction2() {
}

// ----------------------------------------------------------------------
void plotHpt::loopFunction1() {
  char cds[100], cut[100];

  sprintf(cds, "%s", fCds.c_str());
  
  fHists[Form("genpt_%s_goodcand", cds)]->Fill(fb.gpt); 

  if (fGoodCand) { 
    fHists[Form("pt_%s_goodcand", cds)]->Fill(fb.pt); 
    
    if (fb.pt < 250) return;
    if (fb.pt > 250 && fb.pt < 400) { 
      sprintf(cut, "lopt"); 
      fHists[Form("m_%s_%s", cds, cut)]->Fill(fb.m); 
      fHists[Form("pt_%s_%s", cds, cut)]->Fill(fb.pt); 
      fHists[Form("eta_%s_%s", cds, cut)]->Fill(fb.eta); 
      fHists[Form("g0pt_%s_%s", cds, cut)]->Fill(fb.g0pt); 
      fHists[Form("g1pt_%s_%s", cds, cut)]->Fill(fb.g1pt); 
      fHists[Form("g0iso_%s_%s", cds, cut)]->Fill(fb.g0iso); 
      fHists[Form("g1iso_%s_%s", cds, cut)]->Fill(fb.g1iso); 
      fHists[Form("g0isor_%s_%s", cds, cut)]->Fill(fb.g0iso/fb.g0pt); 
      fHists[Form("g1isor_%s_%s", cds, cut)]->Fill(fb.g1iso/fb.g1pt); 
    } else {
      sprintf(cut, "hipt"); 
      fHists[Form("m_%s_%s", cds, cut)]->Fill(fb.m); 
      fHists[Form("pt_%s_%s", cds, cut)]->Fill(fb.pt); 
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
void plotHpt::combMCAtNLOHist(TH1D *h0, TH1D *h1) {
  double xs0(36.44), xs1(-2.37); 

  string h1name = h0->GetName(); 
  replaceAll(h1name, "mcatnlo0", "mcatnlo"); 

  TH1D *h = (TH1D*)h0->Clone(h1name.c_str()); 
  h->Reset();
  h->Add(h0, h1, 1., xs1/xs0); 

  fHists.insert(make_pair(h->GetName(), h)); 

}
