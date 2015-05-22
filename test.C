#ifndef __CINT__
#include "RooCFunction1Binding.h" 
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h" 

using namespace RooFit;

// ----------------------------------------------------------------------
double fExp(double *x, double *par) {
  return TMath::Exp(par[0]+par[1]*x[0]);
}

// ----------------------------------------------------------------------
double fPow(double *x, double *par) {
  return par[2]*TMath::Power(par[0]*x[0], par[1]);
}

// ----------------------------------------------------------------------
double fPareto(double *x, double *par) {
  double xm(300.); 
  if (x[0] > xm)  return par[0]*(par[1]*TMath::Power(xm, par[1]))/TMath::Power(x[0], par[1]+1); 
}

// ----------------------------------------------------------------------
double rPareto(double *x, double *par) {
  double xm(300.); 
  if (x[0] > xm)  return (par[1]*TMath::Power(xm, par[1]))/TMath::Power(x[0], par[1]+1); 
}


// ----------------------------------------------------------------------
double fLogNormal(double *x, double *par) {
  return par[0]*TMath::LogNormal(x[0], par[1], par[2], par[3]); 
}

// ----------------------------------------------------------------------
double fWeibull(double *x, double *par) {
  double k = par[1]; 
  double l = par[2]; 
  return par[0]*(k/l)*TMath::Power(x[0]/l, k-1)*TMath::Exp(TMath::Power(-(x[0]/l),k)); 
}



// ----------------------------------------------------------------------
void noExp(string name = "sherpa", string fname = "power") {

  if (name == "all") {
    noExp("sherpa", fname);
    noExp("manx", fname);
    noExp("man5", fname);
    return;
  }

  TFile *f = TFile::Open("hpt0/plotHpt-m.root");

  TH1D *hb = (TH1D*)(f->Get(Form("pt_%s_goodcand", name.c_str())));
  hb->Scale(1.e6);
  hb->SetAxisRange(300., 1000., "X");

  gPad->SetLogy(1);

  if (fname == "expo") {
    TF1 *f0 = new TF1("fExp", fExp, 300., 1000., 2); 
    f0->SetParameters(1., -0.01); 
  }

  if (fname == "power") {
    TF1 *f0 = new TF1("fPow", fPow, 300., 1000., 3); 
    f0->SetParameters(1., -3., 3.); 
  }

  if (fname == "pareto") {
    TF1 *f0 = new TF1("fPareto", fPareto, 300., 1000., 2); 
    f0->SetParameters(1., 3.); 
  }

  if (fname == "lognormal") {
    TF1 *f0 = new TF1("fLogNormal", fLogNormal, 300., 1000., 4); 
    f0->SetParameters(1., 200., 10., 100.); 
  }

  if (fname == "weibull") {
    TF1 *f0 = new TF1("fWeibull", fWeibull, 300., 1000., 3); 
    f0->SetParameters(1., 1., 0.1); 
  }
  
  hb->Fit(f0); 

  c0.SaveAs(Form("noExp-%s-%s.pdf", name.c_str(), f0->GetName())); 
}


// ----------------------------------------------------------------------
vector<double> vk, vm; 
vector<string> samples;
void RooLognormal(string name = "sherpa", double m0Val = 200.) {
  if (name == "all") {
    samples.clear(); 
    samples.push_back("manx");
    samples.push_back("man5");
    samples.push_back("sherpa");

    vk.clear();
    vm.clear();

    zone(2,2); 
    RooLognormal("manx"); 
    c0.cd(2); 
    RooLognormal("man5"); 
    c0.cd(3); 
    RooLognormal("sherpa"); 
    c0.cd(4); 
    gPad->SetLogy(1);

    RooRealVar pt("pt", "pt", 300., 1000.); 
    RooRealVar m0("m0","m0", 0., 0., 1000.);
    RooRealVar k("k","k", 0., 0., 1000.) ;
    RooLognormal *p(0); 
    RooPlot *plotPt = pt.frame(Title(Form("%s", name.c_str())), Bins(70));

    for (unsigned int i = 0; i < vk.size(); ++i) {
      m0.setVal(vm[i]); 
      k.setVal(vk[i]); 
      p = new RooLognormal(Form("bgPt%d", i), "bgPt", pt, m0, k); 

      int col = kBlue;
      if (samples[i] == "sherpa") col = kRed; 
      if (samples[i] == "man5") col = kGreen+2; 

      p->plotOn(plotPt, LineColor(col)); 
    }
    
    plotPt->Draw();

    
    c0.SaveAs(Form("hpt0/RooLognormal-%.0f.pdf", m0Val)); 
    return;
  }

  
  TFile *f = TFile::Open("hpt0/plotHpt-m.root");
  
  TH1D *hb = (TH1D*)(f->Get(Form("pt_%s_goodcand", name.c_str())));
  hb->SetAxisRange(300., 1000., "X");
  hb->SetLineColor(kRed); 
  double scale = 1.e6; // 1000/fb*1000pb/fb
  hb->Scale(scale); 

  TH1D *hbl = new TH1D("hbl", "", 70, 300., 1000.); 
  int hbbin; 
  for (int i = 1; i <= 71; ++i) {
    hbbin = hb->FindBin(299.) + i;
    hbl->SetBinContent(i, hb->GetBinContent(hbbin)); 
    hbl->SetBinError(i, hb->GetBinError(hbbin)); 
  }

  // -- RooFit part
  RooRealVar m("m", "m", 70, 110, "GeV"); 
  RooRealVar mslope("mslope", "mslope", 0.071276, 0., 1.); 
  
  RooRealVar pt("pt", "pt", 300., 1000.); 

  // -- #1: Lognormal
  RooRealVar m0("m0","m0", m0Val, 0., 1000.);
  //  m0.setConstant(kTRUE); 
  RooRealVar k("k","k", 3.0, 0., 1000.) ;
  RooLognormal bgPt("bgPt", "bgPt", pt, m0, k); 

  // -- #2: Landau
  //   RooRealVar mean("mean", "mean", m0Val, 0., 1000.);
  //   RooRealVar sigma("sigma", "sigma", 50.0, 0., 1000.) ;
  //   RooLandau bgPt("bgPt", "bgPt", pt, mean, sigma); 

  RooDataHist hdata("hdata","hdata", pt, Import(*hbl));
  
 
  RooFitResult *r1 = bgPt.fitTo(hdata);
  vk.push_back(k->getVal()); 
  vm.push_back(m0->getVal()); 

  gPad->SetLogy(1);
  RooPlot *plotPt = pt.frame(Title(Form("%s", name.c_str())), Bins(70));
  
  if (0) {
    hdata.plotOn(plotPt, DataError(RooAbsData::SumW2));
    bgPt.paramOn(plotPt, Layout(0.4, 0.8, 0.85));

    //    RooDataSet *d0 = bgPt.generate(RooArgSet(pt), 1363); 
    //    d0->plotOn(plotPt, MarkerColor(kRed)); 
    //     RooHist *h = plotPt->getHist("h_bgPtData");
    //     cout << "h = " << h << endl;
    //     if (h) {
    //       for (int i = 0; i < h->GetN(); ++i) {
    // 	cout << h->GetX()[i] << " " << h->GetY()[i] 
    // 	     << " " << h->GetErrorYlow(i)	     << " " << h->GetErrorYhigh(i)
    // 	     << endl;
    // 	if (h->GetY()[i] < 1.e-6) {
    // 	  h->SetPointError(i, 0., 0., 0., 0.);
    // 	  h->GetY()[i] = 1.e-4;
    // 	}
    //       }
    //     }
    
    int col = kBlue;
    if (name == "sherpa") col = kRed; 
    if (name == "man5") col = kGreen+2; 

    bgPt.plotOn(plotPt, LineColor(col)); 
    plotPt->SetMinimum(1.e-2);
    plotPt->Draw("");
  
  } else {


    RooDataSet *data  = bgPt->generate(pt, 1360); 

    RooFitResult *r2 = bgPt.fitTo(*data);

    data->plotOn(plotPt); 


    plotPt->Print();
    RooHist *h = plotPt->getHist("h_bgPtData");
    if (h) {
      for (int i = 0; i < h->GetN(); ++i) {
 	if (h->GetY()[i] < 1.e-6) {
 	  h->SetPointError(i, 0., 0., 0., 0.);
 	  h->GetY()[i] = 1.e-4;
 	}
      }
    }

    bgPt.plotOn(plotPt, LineColor(kRed)); 

    plotPt->SetMinimum(0.2);
    plotPt->Draw("");
    

  }
  //  hbl->Draw("hist");
}


// ----------------------------------------------------------------------
void expTest() {

  TF1 *fExp = new TF1("fExp",  "expo", 0., 1000.);
  fExp->SetLineColor(kBlack);
  // -- this is harder
  fExp->SetParameters(1., -0.012);

  TF1 *fAxp = new TF1("fAxp",  "expo", 0., 1000.);
  fAxp->SetLineColor(kRed);
  // -- this is softer
  fAxp->SetParameters(1., -0.013);

  fExp->Draw();
  fAxp->Draw("same");
}

// ----------------------------------------------------------------------
void expModification(int evt = 10000) {

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

  zone(2,2);
  c0.cd(1);
  gPad->SetLogy(1);
  h0->Fit("expo");

  c0.cd(2);
  gPad->SetLogy(1);
  h1->Fit("expo");

  c0.cd(3);
  gPad->SetLogy(1);
  h2->Fit("expo");

} 



// ----------------------------------------------------------------------
void modLognormal(int evt = 258300) {

  TF1 *f1 = new TF1("f1",  "1.0 + 0.01*(x-40)*0.01*(x-40)*0.015", 0., 1000.);
  TF1 *f2 = new TF1("f2",  "1.0 - 0.01*(x-40)*0.01*(x-40)*0.015", 0., 1000.); // w/o modifications the pdf is no longer > 0 everywhere.

  RooRealVar pt("pt", "pt", 0., 1000.); 
  RooRealVar m0("m0","m0", 40., 0., 1000.);
  RooRealVar k("k","k", 2.2, 0., 1000.);

  RooGenericPdf f1Pdf("f1Pdf", "f1Pdf","1.0 + 0.01*(pt - 40)*0.01*(pt - 40)*0.015", pt);
  RooGenericPdf f2Pdf("f2Pdf", "f2Pdf","1.0 - 0.008*(pt - 40)*0.01*(pt - 40)*0.012", pt);

  RooLognormal sg0Pdf("sg0Pt", "sg0Pt", pt, m0, k); 
  RooLognormal sg1Pdf("sg1Pt", "sg1Pt", pt, m0, k); 
  RooLognormal sg2Pdf("sg2Pt", "sg2Pt", pt, m0, k); 

  RooProdPdf sg1f1Pdf("sg1f1Pt", "sg1Pt", RooArgSet(sg1Pdf, f1Pdf)); 
  RooProdPdf sg2f2Pdf("sg2f2Pt", "sg2Pt", RooArgSet(sg2Pdf, f2Pdf)); 

  k.setVal(2.2);
  m0.setVal(40.0);

  RooDataSet *data  = sg0Pdf.generate(pt, evt); 
  RooDataSet *data1 = sg1f1Pdf.generate(pt, evt); 
  RooDataSet *data2 = sg2f2Pdf.generate(pt, evt); 
  RooDataSet *data0 = sg0Pdf.generate(pt, evt); 


  cout << "f1(300)  = " << f1->Eval(300.) << endl;
  cout << "f1(600)  = " << f1->Eval(600.) << endl;
  cout << "f1(1000) = " << f1->Eval(1000.) << endl;

  pt.setVal(300.);
  cout << "300: " << sg1f1Pdf.getVal() << " / " << sg0Pdf.getVal() << " = " << sg1f1Pdf.getVal()/sg0Pdf.getVal() << endl;
  pt.setVal(22.);
  cout << "22: " << sg1f1Pdf.getVal() << " / " << sg0Pdf.getVal() << " = " << sg1f1Pdf.getVal()/sg0Pdf.getVal() << endl;
  pt.setVal(400.);
  cout << "400: " << sg1f1Pdf.getVal() << " / " << sg0Pdf.getVal() << " = " << sg1f1Pdf.getVal()/sg0Pdf.getVal() << endl;
  pt.setVal(600.);
  cout << "600: " << sg1f1Pdf.getVal() << " / " << sg0Pdf.getVal() << " = " << sg1f1Pdf.getVal()/sg0Pdf.getVal() << endl;
  pt.setVal(800.);
  cout << "800: " << sg1f1Pdf.getVal() << " / " << sg0Pdf.getVal() << " = " << sg1f1Pdf.getVal()/sg0Pdf.getVal() << endl;
  pt.setVal(1000.);
  cout << "1000 " << sg1f1Pdf.getVal() << " / " << sg0Pdf.getVal() << " = " << sg1f1Pdf.getVal()/sg0Pdf.getVal() << endl;

  pt.setRange("signal", 300., 1000.); 
  RooAbsReal *integral = sg0Pdf.createIntegral(pt, Range("signal")); 
  cout << "sg0Pdf integral = " << integral->getVal()*evt << endl;


  zone(2,2);
  RooPlot *plotPt = pt.frame(Title("  "));
  data->plotOn(plotPt); 
  sg0Pdf->plotOn(plotPt); 
  sg0Pdf->paramOn(plotPt, Layout(0.3, 0.8, 0.85));
  sg1f1Pdf->plotOn(plotPt, LineColor(kRed)); 
  sg2f2Pdf->plotOn(plotPt, LineColor(kCyan)); 

  plotPt->SetMinimum(0.05);
  plotPt->Draw();
  gPad->SetLogy(1); 

  c0->cd(2); 
  RooFitResult *r1  = sg1f1Pdf.fitTo(*data1);
  
  RooPlot *plot1 = pt.frame(Title("  "));
  data1->plotOn(plot1); 
  sg1f1Pdf->plotOn(plot1, LineColor(kRed)); 
  sg1f1Pdf->paramOn(plot1, Layout(0.3, 0.8, 0.85));
  plot1->SetMinimum(0.05);
  plot1->Draw();
  gPad->SetLogy(1); 


  c0->cd(3); 
  RooFitResult *r2  = sg2f2Pdf.fitTo(*data1);
  
  RooPlot *plot2 = pt.frame(Title("  "));
  data2->plotOn(plot2); 
  sg2f2Pdf->plotOn(plot2, LineColor(kCyan)); 
  sg2f2Pdf->paramOn(plot2, Layout(0.3, 0.8, 0.85));
  plot2->SetMinimum(0.05);
  plot2->Draw();
  gPad->SetLogy(1); 

  c0->cd(4); 
  RooFitResult *r2  = sg0Pdf.fitTo(*data0);
  
  RooPlot *plot0 = pt.frame(Title("  "));
  data0->plotOn(plot0); 
  sg0Pdf->plotOn(plot0, LineColor(kBlue)); 
  sg0Pdf->paramOn(plot0, Layout(0.3, 0.8, 0.85));
  plot0->SetMinimum(0.05);
  plot0->Draw();
  gPad->SetLogy(1); 
  return;
  
  TH1D *h0 = new TH1D("h0", "", 80, 200., 1000.); 
  h0->FillRandom("fExp", evt); 

  TH1D *h1 = (TH1D*)h0->Clone("h1");
  h1->Multiply(f1); 
  h1->Scale(evt/h1->Integral());

  TH1D *h2 = (TH1D*)h0->Clone("h2");
  h2->Multiply(f2); 
  h2->Scale(evt/h2->Integral());

  zone(2,2);
  c0.cd(1);
  gPad->SetLogy(1);
  h0->Fit("expo");

  c0.cd(2);
  gPad->SetLogy(1);
  h1->Fit("expo");

  c0.cd(3);
  gPad->SetLogy(1);
  h2->Fit("expo");

} 


// ----------------------------------------------------------------------
void test(int ntoys) {
  double fSeparationE(0.); 
  if (ntoys < 2200) fSeparationE = 0.02; 
  if (ntoys < 1200) fSeparationE = 0.02; 
  if (ntoys < 700) fSeparationE = 0.03; 
  if (ntoys < 200) fSeparationE = 0.08; 
  cout << "ntoys = " << ntoys << " -> fSeparationE = " << fSeparationE << endl;
}


// -- hstat-s2m12n1000.log
void test1() {
  TH1D *h1 = new TH1D("h1", "h1", 400, 1., 5.); 
  h1->Fill(3.93635); 
  h1->Fill(3.97173);
  h1->Fill(3.90296);
  h1->Fill(4.00931);
  h1->Fill(4.00888);
  h1->Fill(4.00845);
  h1->Fill(4.0089);
  h1->Fill(4.02818);
  h1->Fill(4.00931);
  h1->Fill(3.99079);
  h1->Fill(3.88584);
  h1->Fill(4.0089);
  h1->Fill(3.9195);
  h1->Fill(4.24014);
  h1->Draw();
  c0->SaveAs("hstat-sigStudies-12-1000.pdf");
}

// -- hstat-s2m12n2000.log
void test2() {
  TH1D *h1 = new TH1D("h1", "h1", 400, 1., 5.); 
  h1->Fill(3.92765); 
  h1->Fill(3.93718);
  h1->Fill(4.01806);
  h1->Fill(4.02797);
  h1->Fill(3.9726);
  h1->Fill(3.9455);
  h1->Fill(4.11768);
  h1->Draw();
  c0->SaveAs("hstat-sigStudies-12-2000.pdf");
}
