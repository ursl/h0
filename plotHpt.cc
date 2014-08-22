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

  NBINS = 50;  
  GETA = 2.5;
  G0ISO = G1ISO = 0.2; 
  G0PT = 100.;
  G1PT = 40.;
  PTLO = 300.;
  PTHI = 999.;
  MGGLO = 100.;
  MGGHI = 150.;

}


// ----------------------------------------------------------------------
plotHpt::~plotHpt() {

}


// ----------------------------------------------------------------------
void plotHpt::bookHist(string name) {
  fHists.insert(make_pair(Form("m_%s", name.c_str()), 
			  new TH1D(Form("m_%s", name.c_str()), Form("m_%s", name.c_str()), NBINS, MGGLO, MGGHI))); 

  fHists.insert(make_pair(Form("pt_%s", name.c_str()), 
			  new TH1D(Form("pt_%s", name.c_str()), Form("pt_%s", name.c_str()), 100, 0, 1000.))); 

  fHists.insert(make_pair(Form("g0pt_%s", name.c_str()), 
			  new TH1D(Form("g0pt_%s", name.c_str()), Form("g0pt_%s", name.c_str()), 100, 0., 300.))); 

  fHists.insert(make_pair(Form("g1pt_%s", name.c_str()), 
			  new TH1D(Form("g1pt_%s", name.c_str()), Form("g1pt_%s", name.c_str()), 100, 0., 300.))); 

  fHists.insert(make_pair(Form("g0iso_%s", name.c_str()), 
			  new TH1D(Form("g0iso_%s", name.c_str()), Form("g0iso_%s", name.c_str()), 100, 0., 1.))); 

  fHists.insert(make_pair(Form("g1iso_%s", name.c_str()), 
			  new TH1D(Form("g1iso_%s", name.c_str()), Form("g1iso_%s", name.c_str()), 100, 0., 1.))); 


}


// ----------------------------------------------------------------------
void plotHpt::makeAll(int bitmask) {
  if (bitmask & 0x1) treeAnalysis();
  if (bitmask & 0x2) toy1(); 
}


// ----------------------------------------------------------------------
void plotHpt::treeAnalysis() {

  string ds("sherpa");
  fCds = ds; 
  bookHist(ds); 
  TTree *t = getTree(ds); 
  setupTree(t); 
  loopOverTree(t); 

  ds = "mcatnlo5"; 
  fCds = ds; 
  bookHist(ds); 
  t = getTree(ds); 
  setupTree(t);
  loopOverTree(t); 


  fHists["m_sherpa"]->Draw(); 
  fHists["m_mcatnlo5"]->Draw("same"); 

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
	ds->fLumi = nevt/ds->fXsec/ds->fBf;
	ds->fName = "MC@NLO " + sdecay; 
	ds->fFillStyle = 3356; 
	ds->fSize = 1; 
	ds->fWidth = 2; 
	fDS.insert(make_pair(sname, ds)); 
      }


      if (string::npos != stype.find("sherpa")) {
	dataset *ds = new dataset(); 
	if (string::npos != stype.find("1")) {
	  sname = "sherpa";
	  sdecay = "#gamma #gamma";
	  ds->fColor = kRed; 
	  ds->fLcolor = kRed; 
	  ds->fFcolor = kRed; 
	  ds->fSymbol = 27; 
	} 
	ds->fF = pF; 
	ds->fXsec = atof(sxsec.c_str());
	ds->fBf   = 1.;
	ds->fLumi = nevt/ds->fXsec/ds->fBf;
	ds->fName = "SHERPA " + sdecay; 
	ds->fFillStyle = 3365; 
	ds->fSize = 1; 
	ds->fWidth = 2; 
	fDS.insert(make_pair(sname, ds)); 
      } 

      // mb ub nb pb fb 
      cout << "opened MC file "  << sfile  << " as " << sname << " (" << stype << ") with xsec = " << sxsec
	   << " equivalent lumi = " << fDS[sname]->fLumi/1000. << "/fb"
	   << endl;

    }
  }
  

}



// // ----------------------------------------------------------------------
// void plotHpt::normHist(TH1D *h, double integral, string type) {
//   double scale(1.); 
//   if (TMath::Abs(integral - 1.) < fEpsilon) {
//     // -- normalize to 1
//     scale = (h->GetSumOfWeights() > 0 ? integral/h->GetSumOfWeights() : 1.); 
//   } else if (TMath::Abs(integral + 1.) < fEpsilon) {
//     // -- normalize to xsec*bf
//     //    n = xsec * L
//     //    "integral" over histogram should be xsec

//     scale = (h->GetSumOfWeights() > 0 ? fDS[type]->fXsec*fDS[type]->fBf/h->Integral() : 1.); 
//     setTitles(h, h->GetXaxis()->GetTitle(), "pb");
//   } else {
//     scale = 1.;
//   }
//   double c(0.), e(0.); 
//   for (int i = 0; i <= h->GetNbinsX(); ++i) {
//     c = h->GetBinContent(i); 
//     e = h->GetBinError(i); 
//     h->SetBinContent(i, c*scale);
//     h->SetBinError(i, e*scale);
//   }
  
// }


// ----------------------------------------------------------------------
void plotHpt::overlayAll() {

  // -- simple overlays
  c0->cd(1); 
  overlay("mcatnlo5", "H1pt", "sherpa", "H1pt"); 

}

// ----------------------------------------------------------------------
void plotHpt::candAnalysis() {
  fGoodCand = true; 
  if (fb.m < 100 && fb.m > 150) fGoodCand = false; 
  if (fb.pt < 250) fGoodCand = false; 
  if (TMath::Abs(fb.eta) > GETA) fGoodCand = false; 
  if (fb.g0pt < G0PT) fGoodCand = false; 
  if (fb.g1pt < G1PT) fGoodCand = false; 
  if (fb.g0iso > G0ISO) fGoodCand = false; 
  if (fb.g1iso > G1ISO) fGoodCand = false; 

}

// ----------------------------------------------------------------------
void plotHpt::loopFunction() {
  char cds[100];
  sprintf(cds, "%s", fCds.c_str());
  if (fGoodCand) { 
    fHists[Form("m_%s", cds)]->Fill(fb.m); 
    fHists[Form("pt_%s", cds)]->Fill(fb.pt); 
    fHists[Form("g0pt_%s", cds)]->Fill(fb.g0pt); 
    fHists[Form("g1pt_%s", cds)]->Fill(fb.g1pt); 
    fHists[Form("g0iso_%s", cds)]->Fill(fb.g0iso); 
    fHists[Form("g1iso_%s", cds)]->Fill(fb.g1iso); 
  }
}


// ----------------------------------------------------------------------
void plotHpt::loopOverTree(TTree *t, int nevts, int nstart) {
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

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
    
    candAnalysis();
    loopFunction();
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

