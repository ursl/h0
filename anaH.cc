#include "anaH.hh"

#include <TProfile.h>

#include "util.hh"

using namespace std;

// ----------------------------------------------------------------------
// Run with: ./runH -c chains/bg-test -D root 
//           ./runH -f test.root 
// ----------------------------------------------------------------------

#include "anaH.icc"


bool sortPtP(GenParticle *a, GenParticle *b) {
  return ( a->PT > b->PT );
}

bool sortPtE(Electron *a, Electron *b) {
    return ( a->PT > b->PT );
}

bool sortPtM(Muon *a, Muon *b) {
    return ( a->PT > b->PT );
}

// ----------------------------------------------------------------------
void anaH::startAnalysis() {
  cout << "anaH: startAnalysis: ..." << endl;
}

// ----------------------------------------------------------------------
void anaH::endAnalysis() {
  cout << "anaH: endAnalysis: ..." << endl;
}


// ----------------------------------------------------------------------
void anaH::eventProcessing() {

  initVariables();

//   cout << "w8: " <<  getEvent(0)->Weight 
//        << " processid: " << getEvent(0)->ProcessID
//        << " MPI: " << getEvent(0)->MPI
//        << " X1: " << getEvent(0)->X1 
//        << " X2: " << getEvent(0)->X2
//        << " qcd: " << getEvent(0)->AlphaQCD
//        << " nevents = " << fbEvent->GetEntries() 
//        << endl;

  genLevelAnalysis(); 
  ggAnalysis();
  fillHist(); 
  

  if (0) {
    dumpDaughters(fH0); 
    
    cout << "reco electrons" << endl;
    GenParticle *pGen(0); 
    Electron *pE(0); 
    for (int i = 0; i < fbElectrons->GetEntries(); ++i) {
      pE   = getElectron(i); 
      pGen = 0; 
      pGen = (GenParticle*)pE->Particle.GetObject();
      cout << Form("%2d reco PT: %9.3f %9.3f %9.3f ", i, pE->PT, pE->Eta, pE->Phi);
      if (pGen > 0) cout << Form(" gen PT: %9.3f %9.3f %9.3f %x", pGen->PT, pGen->Eta, pGen->Phi, pGen); 
      cout << endl;
    }

    cout << "reco muons" << endl;
    Muon *pM(0); 
    for (int i = 0; i < fbMuons->GetEntries(); ++i) {
      pM   = getMuon(i); 
      pGen = 0; 
      pGen = (GenParticle*)pM->Particle.GetObject();
      cout << Form("%2d reco PT: %9.3f %9.3f %9.3f ", i, pM->PT, pM->Eta, pM->Phi);
      if (pGen > 0) cout << Form(" gen PT: %9.3f %9.3f %9.3f %x", pGen->PT, pGen->Eta, pGen->Phi, pGen); 
      cout << endl;
    }
  }


}


// ----------------------------------------------------------------------
void anaH::genLevelAnalysis() {

  fH0 = fH1 = 0; 
  fG0 = fG1 = 0; 
  // -- find NLO and SMC Higgs bosons
  GenParticle *pGen(0);
  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    pGen = getParticle(i); 
    if (25 == pGen->PID) {
      if (0 == fH0) {
	fH0 = pGen;
	fp4genH0 = fH0->P4();
      }
      fH1 = pGen;
      fp4genH1 = fH1->P4();
    }
  }
  
  if (0 != fH1) {
    // -- fish out two photons from Higgs
    int dauID = getParticle(fH1->D1)->PID; 
    if (22 == dauID) { 
      fClass   = 0; 
      fG0      =  getParticle(fH1->D1);
      fG1      =  getParticle(fH1->D2);
      fp4genG0 = fG0->P4();
      fp4genG1 = fG1->P4();
    } else {    
      cout << "Higgs decay into non-photons?!" << endl;
      return;
    }
  } else {
    // -- fish out two highest-pT photons
    fClass = 1; 
    for (int i = 0; i < fbParticles->GetEntries(); ++i) {
      pGen = getParticle(i); 
      if (22 == pGen->PID) {
	if (fG0 && pGen->PT > fG0->PT) {
	  //	  cout << "replacing 0: " << fG0->PT << " with " << pGen->PT << endl;
	  fG1 = fG0; 
 	  fG0 = pGen; 
	} else if (fG1 && pGen->PT > fG1->PT) {
	  //	  cout << "replacing 1: " << fG1->PT << " with " << pGen->PT << endl;
	  fG1 = pGen;
	} else if (fG0 == 0) {
	  //	  cout << "setting 0:  with " << pGen->PT << endl;
	  fG0 = pGen;
	} else if (fG1 == 0) {
	  //	  cout << "setting 1:  with " << pGen->PT << endl;
	  fG1 = pGen;
	}
      }
    }

    if (fG0 != 0 && fG1 != 0) {
      //      cout << "0: " << fG0->PT << " 1: " << fG1->PT << endl;
      fp4genG0 = fG0->P4();
      fp4genG1 = fG1->P4();
      fp4genH0 = fG0->P4() + fG1->P4();
      fp4genH1 = fG0->P4() + fG1->P4();
    }
  }
}


// ----------------------------------------------------------------------
void anaH::ggAnalysis() {

  Photon *pG(0), *pG0(0), *pG1(0); 
  GenParticle *pGen(0); 
  int nphotons(0); 
  for (int i = 0; i < fbPhotons->GetEntries(); ++i) {
    pG = getPhoton(i); 
    if (pG->Particles.GetEntriesFast() != 1) {
      //      cout << "photon with references to multiple particles" << endl;
      continue;
    }
    pGen = (GenParticle*)pG->Particles.At(0);
    //    cout << "pGen: " << pGen << " G0: " << fG0 << " G1: " << fG1 << " "; 
    //    cout << pG->PT << " GenParticle pT = " << pGen->PT << endl;
    if (fG0 == pGen) { 
      pG0 = pG;
      ++nphotons;
    } 

    if (fG1 == pGen) { 
      pG1 = pG;
      ++nphotons;
    } 

  }

  fNRecoPhotons = nphotons;

  if (2 == nphotons) {
    fp4G0.SetPtEtaPhiM(pG0->PT, pG0->Eta, pG0->Phi, 0);
    fp4G1.SetPtEtaPhiM(pG1->PT, pG1->Eta, pG1->Phi, 0);
    fp4H = fp4G0 + fp4G1; 

    fG0Iso = iso(pG0, 0.5, 0.5);
    fG1Iso = iso(pG1, 0.5, 0.5);
    //     cout << "==============> built H candidate with m = " << fp4H.M() 
    // 	 << " photons: " << fG0->PT << " " << fp4G0.Pt() << " .. " << fG1->PT << " " << fp4G1.Pt() 
    // 	 << endl;
  }
}


// ----------------------------------------------------------------------
void anaH::llAnalysis() {

}




// ----------------------------------------------------------------------
void anaH::llReco() {

}


// ----------------------------------------------------------------------
void anaH::ggReco() {

}


// ----------------------------------------------------------------------
void anaH::fillHist(int cat) {

  fillRedTreeData(); 
  fTree->Fill();

  string histdir("");
  //  histdir = Form("class%d", cat); 
  //  cout << histdir << endl;
  
  if (fp4genH0.Pt() > 0.) {
    ((TH1D*)fpHistFile->Get("H0pt"))->Fill(fp4genH0.Pt());
    //    printParticle(fH0); 
  }
  if (fp4genH1.Pt() > 0.) {
    ((TH1D*)fpHistFile->Get("H1pt"))->Fill(fp4genH1.Pt());
    ((TH1D*)fpHistFile->Get("Hm"))->Fill(fp4genH1.M());
    //    printParticle(fH0); 
  }

  if (fp4G0.Pt() > 0.) {
    ((TH1D*)fpHistFile->Get("gpt"))->Fill(fp4G0.Pt());
    ((TH1D*)fpHistFile->Get("geta"))->Fill(fp4G0.Eta());
  }    

  if (fp4G1.Pt() > 0.) {
    ((TH1D*)fpHistFile->Get("gpt"))->Fill(fp4G1.Pt());
    ((TH1D*)fpHistFile->Get("geta"))->Fill(fp4G1.Eta());
  }    

  if (fp4H.Pt() > 10.) {
    ((TH1D*)fpHistFile->Get("Hrm"))->Fill(fp4H.M());
    ((TH1D*)fpHistFile->Get("Hrpt"))->Fill(fp4H.Pt());
    if (fp4H.Pt() > 0.) ((TH1D*)fpHistFile->Get("Hdpt"))->Fill(fp4H.Pt() - fp4genH1.Pt());
    if (fp4H.Pt() > 0.) ((TH1D*)fpHistFile->Get("Hd2pt"))->Fill((fp4H.Pt() - fp4genH1.Pt())/fp4genH1.Pt());

    if (fp4H.Pt() > 0.) ((TH1D*)fpHistFile->Get("Hdpt"))->Fill(fp4H.Pt() - fp4genH1.Pt());
    if (fp4H.Pt() > 0.) ((TH1D*)fpHistFile->Get("Hd2pt"))->Fill((fp4H.Pt() - fp4genH1.Pt())/fp4genH1.Pt());
  }

}



// ----------------------------------------------------------------------
void anaH::bookHist() {
  fpHistFile->cd();			   
  cout << "==> anaH: bookHist"  << endl;

  TH1D *h1(0); 
  TH2D *h2(0); 
  TProfile *hp(0);

  //  fpHistFile->mkdir(Form("class%d", i)); 
  //  fpHistFile->cd(Form("class%d", i));
  
  h1 = new TH1D("H0pt", "", PTN, 0.0, PTMAX); 
  setTitles(h1, "p_{T}^{H} [GeV]", Form("Candidates / %3.0f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  h1 = new TH1D("H1pt", "", PTN, 0.0, PTMAX);
  setTitles(h1, "p_{T}^{H} [GeV]", Form("Candidates / %3.0f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  h1 = new TH1D("gpt",  "", 100, 0.0, 200.0);
  setTitles(h1, "p_{T}^{#gamma} [GeV]", Form("Candidates / %3.0f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  h1 = new TH1D("geta", "", 100, -5.0, 5.0);
  setTitles(h1, "#eta^{#gamma} ", Form("Candidates / %3.2f", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  
  h1 = new TH1D("Hpt", "", PTN, 0.0, PTMAX);
  setTitles(h1, "p_{T}^{H} [GeV]", Form("Candidates / %3.0f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  h1 = new TH1D("Hdau", "", 100, 0, 100.0);
  setTitles(h1, "H daughters", Form("Candidates / Bin"), 0.05, 1.1, 2.0);
  h1 = new TH1D("Hm", "", 200, 124.5, 125.5);
  setTitles(h1, "m^{H} [GeV]", Form("Candidates / %3.0f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  
  h1 = new TH1D("Hrm", "", 200, 100, 140);
  setTitles(h1, "m^{H} [GeV]", Form("Candidates / %3.2f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  h1 = new TH1D("Hrpt", "", PTN, 0, PTMAX);
  setTitles(h1, "p_{T}^{H} [GeV]", Form("Candidates / %3.0f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  h1 = new TH1D("Hdpt", "", 200, -100, 100);
  setTitles(h1, "#Delta(p_{T}^{H}) [GeV]", Form("Candidates / %3.0f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  h1 = new TH1D("Hd2pt", "", 100, -0.2, 0.2);
  setTitles(h1, "#Delta(p_{T}^{H})/p_{T}^{H gen}", Form("Candidates / %3.3f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  
  hp = new TProfile("Pdpt", "", PTN, 0., PTMAX, -100, 100);
  setTitles(h1, "#Delta(p_{T}^{H}) [GeV]", Form("Candidates / %3.0f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  hp = new TProfile("Pd2pt", "", PTN, 0., PTMAX, -0.2, 0.2);
  setTitles(h1, "#Delta(p_{T}^{H})/p_{T}^{H gen}", Form("Candidates / %3.3f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 2.0);
  
  h1 = new TH1D("recophotons", "", 10, 0, 10);
  setTitles(h1, "N^{#gamma} ", Form("Events"), 0.05, 1.2, 2.0);
  
  
  h2 = new TH2D("mall", "muon (all) pt eta", 50, -2.5, 2.5, 50, 0., 100.);
  h2 = new TH2D("mpass", "muon (pass) pt eta", 50, -2.5, 2.5, 50, 0., 100.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupReducedTree();
  
  
}


// ----------------------------------------------------------------------
void anaH::initVariables() {
  fW8 = 0.;
  fH0 = fH1 = 0; 
  fG0 = fG1 = 0; 
  
  fp4genH0 = fp4genH1 = fp4genG0 = fp4genG1
    = fp4H = fp4G0 = fp4G1
    = TLorentzVector(0,0,0,0);

  fClass = -1; 
  fNRecoPhotons = 0; 

  fRtd.type = -1; 
  
  fRtd.m   = -99.; 
  fRtd.pt  = -99.; 
  fRtd.eta = -99.; 
  fRtd.phi = -99.; 

  fRtd.g0pt  = -99.;
  fRtd.g0eta = -99.; 
  fRtd.g0phi = -99.; 
  fRtd.g0iso = -99.;

  fRtd.g1pt  = -99.;
  fRtd.g1eta = -99.;
  fRtd.g1phi = -99.;
  fRtd.g1iso = -99.;

  fRtd.gm   = -99.;
  fRtd.gpt  = -99.;
  fRtd.geta = -99.;
  fRtd.gphi = -99.;

  fRtd.gg0pt  = -99.;
  fRtd.gg0eta = -99.;
  fRtd.gg0phi = -99.;
  fRtd.gg0iso = -99.;

  fRtd.gg1pt  = -99.;
  fRtd.gg1eta = -99.;
  fRtd.gg1phi = -99.;
  fRtd.gg1iso = -99.;

}

// ----------------------------------------------------------------------
void anaH::printParticle(GenParticle *p) {
  if (p) {
    for (int i = 0; i < fbParticles->GetEntries(); ++i) {
      if (p == getParticle(i)) {
	cout << Form("%4d %+6d  M: %4d %4d D: %4d %4d P: %+9.3f %+9.3f %+9.3f PT: %+9.3f %+9.3f %+9.3f %x" , 
		     i, p->PID, p->M1, p->M2, p->D1, p->D2, 
		     p->Px, p->Py, p->Pz,
		     p->PT, p->Eta, p->Phi,
		     p
		     ) 
	     << endl;
	break;
      }  
    }
  }
}


// ----------------------------------------------------------------------
bool anaH::isAncestor(GenParticle *mo, GenParticle *dau) {

  if (mo == dau)  return true;

  GenParticle *pGen(dau); 
  while ((pGen = getParticle(pGen->M1))) {
    if (mo == pGen) {
      return true;
    }
    if (pGen->M1 < 0) break;
    if (pGen->M1 > fbParticles->GetEntries()) break;
  }
  return false;
}

// ----------------------------------------------------------------------
void anaH::dumpGenBlock(bool withGluons) {

  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    if (21 == getParticle(i)->PID) {
      if (withGluons) printParticle(getParticle(i)); 
    } else {
      printParticle(getParticle(i)); 
    }
  }

}


// ----------------------------------------------------------------------
void anaH::muonEfficiency() {

  string histdir = Form("class%d", fClass); 
  
  GenParticle *gm(0); 
  vector <GenParticle *> vM;
  // create unique list of muons throwing away radiative descendants
  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    gm = getParticle(i);
    if (13 != TMath::Abs(gm->PID)) continue;
    bool fresh(true);
    for (unsigned int j = 0; j < vM.size(); ++j) {
      if (isAncestor(vM[j], gm)) {
	fresh = false; 
	break;
      }
    }
    if (fresh) {
      vM.push_back(gm); 
    }
  }

  Muon *pM;
  for (unsigned int i = 0; i < vM.size(); ++i) {
    ((TH1D*)fpHistFile->Get(Form("%s/mall", histdir.c_str())))->Fill(vM[i]->Eta, vM[i]->PT);
    for (int j = 0; j < fbMuons->GetEntries(); ++j) {
      pM = getMuon(j); 
      // cone matching ...
    }
  }

    

}


// ----------------------------------------------------------------------
void anaH::dumpDaughters(GenParticle *pMom) {
  GenParticle *pGen(0); 
  for (int i = pMom->D1; i <= pMom->D2; ++i) {
    pGen = getParticle(i); 
    printParticle(pGen); 
    if (pGen->D1 > -1) dumpDaughters(pGen); 
  }
      
  
}

// ----------------------------------------------------------------------
void anaH::setupReducedTree() {
  fTree->Branch("type",    &fRtd.type,            "type/I");
  fTree->Branch("m",       &fRtd.m,               "m/D");
  fTree->Branch("w8",      &fRtd.w8,              "w8/D");
  fTree->Branch("pt",      &fRtd.pt,              "pt/D");
  fTree->Branch("eta",     &fRtd.eta,             "eta/D");
  fTree->Branch("phi",     &fRtd.phi,             "phi/D");

  fTree->Branch("gm",      &fRtd.gm,              "gm/D");
  fTree->Branch("gpt",     &fRtd.gpt,             "gpt/D");
  fTree->Branch("geta",    &fRtd.geta,            "geta/D");
  fTree->Branch("gphi",    &fRtd.gphi,            "gphi/D");

  fTree->Branch("g0pt",    &fRtd.g0pt,            "g0pt/D");
  fTree->Branch("g0eta",   &fRtd.g0eta,           "g0eta/D");
  fTree->Branch("g0phi",   &fRtd.g0phi,           "g0phi/D");
  fTree->Branch("g0iso",   &fRtd.g0iso,           "g0iso/D");

  fTree->Branch("g1pt",    &fRtd.g1pt,            "g1pt/D");
  fTree->Branch("g1eta",   &fRtd.g1eta,           "g1eta/D");
  fTree->Branch("g1phi",   &fRtd.g1phi,           "g1phi/D");
  fTree->Branch("g1iso",   &fRtd.g1iso,           "g1iso/D");

  fTree->Branch("gg0pt",   &fRtd.gg0pt,           "gg0pt/D");
  fTree->Branch("gg0eta",  &fRtd.gg0eta,          "gg0eta/D");
  fTree->Branch("gg0phi",  &fRtd.gg0phi,          "gg0phi/D");
  fTree->Branch("gg0iso",  &fRtd.gg0iso,          "gg0iso/D");

  fTree->Branch("gg1pt",   &fRtd.gg1pt,           "gg1pt/D");
  fTree->Branch("gg1eta",  &fRtd.gg1eta,          "gg1eta/D");
  fTree->Branch("gg1phi",  &fRtd.gg1phi,          "gg1phi/D");
  fTree->Branch("gg1iso",  &fRtd.gg1iso,          "gg1iso/D");


}


// ----------------------------------------------------------------------
void anaH::fillRedTreeData(int type) {

  fRtd.type = fClass; 
  fRtd.w8   = fW8;

  //  cout << "  fillRedTreeData: " << fW8 << endl;

  if (fp4H.Pt() > 0.) {
    fRtd.m   = fp4H.M(); 
    fRtd.pt  = fp4H.Pt(); 
    fRtd.eta = fp4H.Eta(); 
    fRtd.phi = fp4H.Phi(); 
  }

  if (fp4G0.Pt() > 0.) {
    fRtd.g0pt  = fp4G0.Pt(); 
    fRtd.g0eta = fp4G0.Eta(); 
    fRtd.g0phi = fp4G0.Phi(); 
    fRtd.g0iso = fG0Iso; 
  }

  if (fp4G1.Pt() > 0.) {
    fRtd.g1pt  = fp4G1.Pt(); 
    fRtd.g1eta = fp4G1.Eta(); 
    fRtd.g1phi = fp4G1.Phi(); 
    fRtd.g1iso = fG1Iso; 
  }

  if (fp4genH1.Pt() > 0.) {
    fRtd.gm   = fp4genH1.M(); 
    fRtd.gpt  = fp4genH1.Pt(); 
    fRtd.geta = fp4genH1.Eta(); 
    fRtd.gphi = fp4genH1.Phi(); 
  }
  
  if (fp4genG0.Pt() > 0.) {
    fRtd.gg0pt  = fp4genG0.Pt(); 
    fRtd.gg0eta = fp4genG0.Eta(); 
    fRtd.gg0phi = fp4genG0.Phi(); 
  }

  if (fp4genG1.Pt() > 0.) {
    fRtd.gg1pt  = fp4genG1.Pt(); 
    fRtd.gg1eta = fp4genG1.Eta(); 
    fRtd.gg1phi = fp4genG1.Phi(); 
  }

}

// ----------------------------------------------------------------------
double anaH::iso(Photon *gamma, double radius, double ptmin) {

  double et(0.); 
  //  cout << "gamma: " << gamma->PT << " " << gamma->Eta << endl;
  Tower *t; 
  TLorentzVector p4; 
  TLorentzVector g4 = gamma->P4(); 
  for (int i = 0; i < fbPFphotons->GetEntries(); ++i) {
    t = (Tower*)fbPFphotons->At(i);
    p4 = t->P4(); 
    if (p4.Pt() < ptmin) continue;
    if (t->Particles.At(0) == gamma->Particles.At(0)) {
      //      cout << " tower: " << t->ET << " " << t->Eta << endl;
    } else {
      if (g4.DeltaR(p4) < radius) {
	et += t->ET; 
      }
    }
  }

  Track *tr;
  for (int i = 0; i < fbPFtracks->GetEntries(); ++i) {
    tr = (Track*)fbPFtracks->At(i);
    p4 = tr->P4(); 
    if (p4.Pt() < ptmin) continue;
    if (tr->Particle == gamma->Particles.At(0)) {
      //      cout << " track: " << tr->PT << " " << tr->Eta << endl;
    } else {
      if (g4.DeltaR(p4) < radius) {
	et += tr->PT; 
      }
    }
  }

  return et; 

}
