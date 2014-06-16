// C++
#include <iostream>
#include <stdio.h>  
#include <vector>
#include <math.h>
#include <iomanip>

// ROOT
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLegend.h"
#include "TRandom1.h"
#include "TPad.h"
#include "TH2F.h"
#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"

// CMS2
#include "CMS2.h"
#include "muonSelections.h"
#include "susySelections.h"
#include "ssSelections.h"
#include "eventSelections.h"
#include "MT2/MT2.h"
#include "MT2/MT2Utility.h"


// Good run list

#include "/home/users/jgran/CMSSW_5_3_2_patch4_V05-03-23/src/CMS2/NtupleMacros/Tools/goodrun.cc"

// My includes
#include "myheader.h"
#include "NSBABY.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVec;

using namespace std;
//using namespace tas;
using namespace ROOT::Math;
using namespace nsbaby;


int ScanChain( TChain* chain, char* suffix = "", bool ismc = true, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  const int FIRTBIN = 0;
  const int LASTBIN = 180;
  const int BINNUM = 60;
  const int INVMLASTBIN = 300;
  const int METLASTBIN = 300;
  const int MTLASTBIN = 200;

  TFile* f1 = new TFile("./hists/hists_odata_50met0btehlt_7.root");
  TH1F* RWP = (TH1F*)f1->Get("wp"); 
  // f1->Close();

////////////////////////////////////////////////////////////////////////////////////////////
  TH1F* MT2_hist  = new TH1F("MT2","MT2 distribution", BINNUM, FIRTBIN, LASTBIN+1);
  TH1F* MT2_el    = new TH1F("MT2e","MT2 distribution for e  events", BINNUM, FIRTBIN, LASTBIN+1);
  TH1F* MT2_mu    = new TH1F("MT2m","MT2 distribution for mu events", BINNUM, FIRTBIN, LASTBIN+1);
  TH1F* h2_nozero = new TH1F("hnz","MT2 without 0 bin", BINNUM, 1, LASTBIN+1);
  TH1F* MT2_Cut50 = new TH1F("h50","MT2 without 0 bin with met cut 50", BINNUM, 1, LASTBIN+1);
        
  TH1F* MT_hist  = new TH1F("MT",  "MT distribution", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MTc_hist = new TH1F("MTc", "MT distribution", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MT_el    = new TH1F("MTe", "MT distribution for e  events", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MT_mu    = new TH1F("MTm", "MT distribution for mu events", BINNUM, FIRTBIN, MTLASTBIN+1);
        
  TH1F* JetMult_a = new TH1F("jeta", "Jet mutiplicity all", 13, 0, 13);
  TH1F* Jeta_el   = new TH1F("jae",  "Jet mutiplicity for e  events", 13, 0, 13);
  TH1F* Jeta_mu   = new TH1F("jam",  "Jet mutiplicity for mu events", 13, 0, 13);
  TH1F* JetMult_b = new TH1F("jetb", "b-Jet mutiplicity", 13, 0, 13);
  TH1F* Jetb_el   = new TH1F("jbe",  "b-Jet mutiplicity for e  events", 13, 0, 13);
  TH1F* Jetb_mu   = new TH1F("jbm",  "b-Jet mutiplicity for mu events", 13, 0, 13);
        
  TH1F* Met_all = new TH1F("cmet",  "MET for all events", 70, 0, METLASTBIN+1);
  TH1F* Met_el  = new TH1F("mete",  "MET for e  events",  70, 0, METLASTBIN+1);
  TH1F* Met_mu  = new TH1F("metm",  "MET for mu events",  70, 0, METLASTBIN+1);
  TH1F* FMet    = new TH1F("fmet",  "Real Met + Fake \nu",70, 0, METLASTBIN+1);
        
  TH1F* MetPhi    = new TH1F("metphi",  "MetPhi ",		   90, -3.3, 3.3);
  TH1F* MetPhiCor = new TH1F("cmetphi", "Corrected MetPhi ",       90, -3.3, 3.3);
  TH1F* FMetPhi  = new TH1F("fmetphi", "Metphi + fake nu",         90, -3.3, 3.3);
  TH1F* Metx    = new TH1F("metx", "Met x component",           120, -240, 240);
  TH1F* Mety    = new TH1F("mety", "Met y component",           120, -240, 240);
  TH1F* MetCorx = new TH1F("cmetx","Corrected Met x component", 120, -240, 240);
  TH1F* MetCory = new TH1F("cmety","Corrected Met y component", 120, -240, 240);
  TH1F* DiffPhi  = new TH1F("dphi",  "Difference in phi for lep and met", 70, 0, 3.3);
  TH1F* DiffPhic = new TH1F("dphic", "Difference in cphi for lep and met", 70, 0, 3.3);
  TH1F* DiffPhi_e  = new TH1F("dpm",  "Difference in phi for lep and met", 70, 0, 3.3);
  TH1F* DiffPhi_m  = new TH1F("dpe",  "Difference in phi for lep and met", 70, 0, 3.3);
        
  TH1F* InvM_lep = new TH1F("invm","Dileptons InvMass all",   70, 0, INVMLASTBIN+1);
  TH1F* W_P      = new TH1F("rwp",  "W_P value by (lpt + met)/cos\theta ",  80, 0, 1601);
  TH1F* W_P_el   = new TH1F("rwpe",  "W_P value by (lpt + met)/cos\theta mu",   80, 0, 1601);
  TH1F* W_P_mu   = new TH1F("rwpm",  "W_P value by (lpt + met)/cos\theta ele",  80, 0, 1601);
  TH1F* W_PT     = new TH1F("rwpt","Pt value of lpt + met ",  70, 0, 351);
  TH1F* W_PT_el  = new TH1F("wpte","Pt value of lpt + met for ele",   70, 0, 351);
  TH1F* W_PT_mu  = new TH1F("wptm","Pt value of lpt + met for muon",  70, 0, 351);
  TH1F* Lep_PT = new TH1F("lpt","Pt value of leptons ",    70, 0, 351);
  TH1F* Mu_PT  = new TH1F("mpt","Pt value of muons ",      70, 0, 351);
  TH1F* TEST  = new TH1F("test","Testing histo for an interesting variable",    70, 0, 181);
  TH1F* TEST3 = new TH1F("test3","Testing histo 3 for an interesting variable", 70, 0, 351);
  TH1F* TEST2 = new TH1F("test2","Testing histo 2 for an interesting variable", 70, -3, 3);
  TH1F* LepEta = new TH1F("lepEta","Eta Distribution of real leps",  90, -2.7, 2.7);
  TH1F* ElEta  = new TH1F("elEta" ,"Eta Distribution of electrons",  90, -2.7, 2.7);
  TH1F* MuEta  = new TH1F("muEta" ,"Eta Distribution of muons",      90, -2.7, 2.7);

  MT_hist   ->Sumw2(); 
  MTc_hist  ->Sumw2(); 
  MT2_hist  ->Sumw2(); 
  h2_nozero ->Sumw2();
  MT2_Cut50 ->Sumw2();

  InvM_lep  ->Sumw2();
  JetMult_a ->Sumw2();
  JetMult_b ->Sumw2();
  Met_all   ->Sumw2();
  FMet      ->Sumw2();


///////////////////////////////////////////////////////////////////////////////////////////

  int file_count = 0;

  // int less_jets = 0;
  // int nGoodEvents = 0;
  // float bTagDiscriminant = 0.244; 

  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // Set Good Run List
  // set_goodrun_file("/home/users/jgran/analysis/sswh/fakes/json/final_19p49fb_cms2.txt");

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    
    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("tree");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    // cms2.Init(tree);
    baby.Init(tree);  
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      // cms2.GetEntry(event);
      baby.GetEntry(event);  
      ++nEventsTotal;
    
      // Progress	
      CMS2::progress( nEventsTotal, nEventsChain );

      // Select Good Runs

      Int_t n_jets = 0;
      Int_t n_bTag = 0;
	
      // Applying cuts
      if( lr_p4().pt() < 30 )			 continue;
      if( fabs(lr_p4().eta()) > 2.1 )		 continue;
      if( isRealData() && !isHLT1Lep() )         continue;
      if( isRealData() && !isHLT2Lep() )         continue;
      //if( met() < 20 )                           continue;
      // if( njets() < 5)				 continue;

    //// To be done: ////
      //// Later ///    GenRandomW:   //////			
      TRandom1 rand; 

      float M_W = 80.2;

      ///////////////// fake w initialization  ///////////////////
      float fakeW_eta = -4;	
      float fakeW_theta = -4;
      float fakeW_phi = -4;
      float fakeW_p  = -1;
      Float_t fakeW_px = 0, fakeW_py = 0 , fakeW_pz = 0;

      ///////////////// fake lepton initialization  ///////////////////
      Float_t fakeLep_px = 0, fakeLep_py = 0 , fakeLep_pz = 0;
      float fakeLep_p  = M_W/2;		// Later can apply a width here
      float fakeLep_eta = -4;
      float fakeLep_theta = -4;
      float fakeLep_phi = -4;

      LorentzVec lf_p4;   // lf_p4 == fakeLep_p4;
      
      float ww = 0; 

      // do{
      // 	ww = rand.Gaus(60,20);
      // }while(ww < 10);

      //fakeW_p = 30; // WP->GetRandom();
      fakeW_p = RWP->GetRandom();

      // float rdmWpt = RWPT->GetRandom();
      // Fill1F(TEST, rdmWpt, scale_1fb());
      // rdmWpt = RWPT->GetRandom();
      
      // do{
      	float w1 = rand.Rndm();
      	float w2 = rand.Rndm();
      	fakeW_theta = acos(2*w2-1);
      	fakeW_eta = - log(0.5*fakeW_theta);
      	fakeW_phi = 2*3.14156265359*(w1-0.5);
      // } while(fabs(fakeW_eta) > 2.4); 	

      fakeW_px = fakeW_p* sin(fakeW_theta)* cos(fakeW_phi);
      fakeW_py = fakeW_p* sin(fakeW_theta)* sin(fakeW_phi);
      fakeW_pz = fakeW_p* cos(fakeW_theta);

      fakeW_px = 10;
      fakeW_py = 30;
      fakeW_pz = 10;

      TLorentzVector fakeLepP4, fakeNeuP4;

     // Generate rest decay 
      do{
	float r1 = rand.Rndm();
	float r2 = rand.Rndm();
	fakeLep_theta = acos(2*r2-1);
	fakeLep_eta = - log(0.5*fakeLep_theta);
	fakeLep_phi = 2*3.14156265359*(r1-0.5);
	//} while(fabs(fakeLep_eta) > 2.4); 	// Want those would pass selection

	fakeLep_px = fakeLep_p* sin(fakeLep_theta)* cos(fakeLep_phi);
	fakeLep_py = fakeLep_p* sin(fakeLep_theta)* sin(fakeLep_phi);
	fakeLep_pz = fakeLep_p* cos(fakeLep_theta);

	// float v2 = 1/((M_W/(fakeW_p))*(M_W/(fakeW_p))+1);
	// float gamma = 1/(sqrt(1-v2));
	float gammaM = sqrt(fakeW_p*fakeW_p + M_W*M_W); // gammaM == gamma times Mass of W

	// Initialize lep and nu at W rest frame
	fakeLepP4.SetPxPyPzE( fakeLep_px,  fakeLep_py,  fakeLep_pz, fakeLep_p);
	fakeNeuP4.SetPxPyPzE(-fakeLep_px, -fakeLep_py, -fakeLep_pz, fakeLep_p);

	// Boost from W rest frame to Lab frame as W_p / gamma*Mass == v 
	fakeLepP4.Boost(fakeW_px / gammaM, fakeW_py / gammaM, fakeW_pz / gammaM );
	fakeNeuP4.Boost(fakeW_px / gammaM, fakeW_py / gammaM, fakeW_pz / gammaM );

      } while(fakeLepP4.Eta() > 2.1);
      
      // cout << "Sanity check: M_W = (fakeLepP4+fakeNeuP4).M() = " << (fakeLepP4+fakeNeuP4).M() << "  p_W = (fakeLepP4 + fakeNeuP4).P() = " << (fakeLepP4+fakeNeuP4).P() << endl;
      // cout << "original vs boosted: " 
      // 	   << "vector: " << fakeW_px / M_W << "  " << fakeW_py / M_W << "  " <<  fakeW_pz / M_W << endl
      //      << "org fakeLep: " << fakeLep_px << "  " <<  fakeLep_py << "  " << fakeLep_pz << endl
      // 	   << "boosted fakeLep: " << fakeLepP4.Px() << "  " << fakeLepP4.Py() << "  " << fakeLepP4.Pz() << endl << endl;

      lf_p4.SetPxPyPzE(fakeLepP4.Px(), fakeLepP4.Py(), fakeLepP4.Pz(), fakeLepP4.P());
      // lf_p4.SetPxPyPzE(fakeLep_px, fakeLep_py, fakeLep_pz, fakeLep_p);

      float nux = - fakeLep_px;   // fakeNeuP4.Px();
      float nuy = - fakeLep_py;   // fakeNeuP4.Py();

      // float M_ll = (lr_p4() + lf_p4).M();
     
      float scale1fb = ( isRealData() )? 1 
	: (abs(lr_id()) == 13) ? getTriggerEfficiency_HLT_IsoMu24(lr_p4().pt(), lr_p4().eta()) * scale_1fb()
	: getTriggerEfficiency_HLT_Ele27_WP80(lr_p4().pt(), lr_p4().eta()) * scale_1fb();


      Fill1F(Met_all, met(), scale1fb);
      Fill1F(MetPhi, metPhi(), scale1fb);
      
      float diffPhi = fabs(lr_p4().Phi() - metPhi());
      if(diffPhi >  3.1415926) diffPhi = 6.2831853 - diffPhi;
      
      Fill1F(DiffPhi, diffPhi , scale1fb);

      // for metphi correction
      float metx = met() * cos( metPhi() );
      float mety = met() * sin( metPhi() );

      float wr_px = metx + lr_p4().Px();
      float wr_py = mety + lr_p4().Py();
      float wr_pt = sqrt(wr_px* wr_px  + wr_py*wr_py);

      // Plot PT of real W
      wr_px = metx + lr_p4().Px();
      wr_py = mety + lr_p4().Py();
      wr_pt = sqrt(wr_px* wr_px  + wr_py*wr_py);
      
      Fill1F(W_PT, wr_pt, scale1fb);

      float tmetx = metx + nux;	                 // fmet == met after adding fake neutrinos
      float tmety = mety + nuy;

      float tmet = sqrt( tmetx*tmetx + tmety*tmety );
      float tmetphi = atan2( tmety , tmetx );

      float shiftx = 0.;
      float shifty = 0.;

      shiftx = (! isRealData()) ? (0.1166 + 0.0200*nvtxs()) : (0.2661 + 0.3217*nvtxs());
      shifty = (! isRealData()) ? (0.2764 - 0.1280*nvtxs()) : (-0.2251 - 0.1747*nvtxs());
      
      //Fill1F(TEST, wr_p, scale_1fb());
    
      metx -= shiftx;
      mety -= shifty;

      float cmet = sqrt( metx*metx + mety*mety ); // cmet == corrected met
      float cmetphi = atan2( mety , metx );

      // if( cmet < 50 ) continue;
      float diffPhic = fabs(lr_p4().phi() - cmetphi);      
      if(diffPhic >  3.1415926) diffPhic = 6.2831853 - diffPhic;
      Fill1F(DiffPhic, diffPhic , scale1fb);

      float fmetx = metx + nux;	                 // fmet == met after adding fake neutrinos
      float fmety = mety + nuy;

      float fmet = sqrt( fmetx*fmetx + fmety*fmety );
      float fmetphi = atan2( fmety , fmetx );

      Fill1F(FMet, fmet, scale1fb);
      Fill1F(FMetPhi, fmetphi, scale1fb);


      Fill1F(Lep_PT, lr_p4().pt(), scale1fb);


      if(njets() == 0) cerr << "nJets Branch 0 !! \n";
      // if(njets() > 9)
      // 	cout << "nJets = " << njets() << " at eventNumber: " << eventNumber() << ", runNumber: "
      // 	     << runNumber()  << ", lumiBlock: " << lumiBlock() << endl;
 
      // Fill the jet mutiplicity to the histogram
      Fill1F(JetMult_a, njets(), scale1fb);
      Fill1F(JetMult_b, nbTag(), scale1fb);

      float mt2 = MT2(fmet, fmetphi, lr_p4() , lf_p4);

      Fill1F(MT2_hist, mt2, scale1fb); 
      h2_nozero->Fill(mt2, scale1fb); 

      float mt2t = MT2(tmet, tmetphi, lr_p4() , lf_p4);
      Fill1F(TEST, mt2t, scale1fb);

      float mt  = sqrt(2*lr_p4().Pt()* met()*(1 - cos(lr_p4().Phi() - metPhi())));
      float mtc = sqrt(2*lr_p4().Pt()* cmet *(1 - cos(lr_p4().Phi() - cmetphi )));

      Fill1F(MT_hist,  mt,  scale1fb); 
      Fill1F(MTc_hist, mtc, scale1fb); 

      float wr_cosTheta = mtc/80.2;			// assume W mass at 80, calculate cosTheta from M_T/80.2
      float wr_p = wr_pt / wr_cosTheta;			// calculate W momentum from W_pt by cosTheta

      Fill1F(W_P, wr_p, scale1fb);

      Fill1F(LepEta, lr_p4().eta(), scale1fb);

      if(abs(lr_id()) == 11 ){
	Fill1F(Jeta_el, njets(), scale1fb); 
	Fill1F(Jetb_el, nbTag(), scale1fb); 
	Fill1F(MT2_el, mt2, scale1fb); 
	Fill1F(MT_el, mt, scale1fb); 
	Fill1F(Met_el, met(), scale1fb); 
	Fill1F(ElEta, lr_p4().eta(), scale1fb);
	Fill1F(W_PT_el, wr_pt, scale1fb);
	Fill1F(W_P_el, wr_p, scale1fb);
	Fill1F(DiffPhi_e, diffPhi , scale1fb);

      }
      if(abs(lr_id()) == 13 ){
	Fill1F(Jeta_mu, njets(), scale1fb);
	Fill1F(Jetb_mu, nbTag(), scale1fb);
	Fill1F(MT2_mu, mt2, scale1fb); 
	Fill1F(MT_mu, mt, scale1fb); 
	Fill1F(Met_mu, met(), scale1fb); 
	Fill1F(MuEta, lr_p4().eta(), scale1fb);
	Fill1F(Mu_PT, lr_p4().pt(), scale1fb);
	Fill1F(W_PT_mu, wr_pt, scale1fb);
	Fill1F(W_P_mu, wr_p, scale1fb);
	Fill1F(DiffPhi_m, diffPhi , scale1fb);
	
      }

    } //loop over events in the current file

    file_count++;

    // Clean Up
    delete tree;
    file->Close();
    delete file;
  } //file loop
  
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  cout << "\nNumber of total Events: " << nEventsTotal 
       //<< " at 20GeV cut, " << nttbarEvents-metCut30 
       <<endl <<endl;
  cout << "For the data samples: "
       //<< "\n# notGoodRun: " << notGoodRun 
       //<< "\n# goodRun: " << goodRun
       //<< "\n# noGoodVtx: " << noGoodVtx
       << endl;


  TFile* fout = new TFile(Form("./hists/hists%s.root",suffix),"RECREATE");

  MT2_hist ->Write();
  MT2_el   ->Write();
  MT2_mu   ->Write();
  MT_hist  ->Write();
  MTc_hist ->Write();
  MT_el    ->Write();
  MT_mu    ->Write();
  MT2_Cut50->Write();
  h2_nozero->Write();
  Met_all  ->Write();
  Met_el   ->Write();
  Met_mu   ->Write();
  FMet     ->Write();
  InvM_lep ->Write();
  JetMult_a->Write();
  JetMult_b->Write();
  Jeta_el  ->Write();
  Jeta_mu  ->Write();
  Jetb_el  ->Write();
  Jetb_mu  ->Write();
  W_P      ->Write();
  W_PT     ->Write();
  Lep_PT   ->Write();
  Mu_PT    ->Write();
  W_PT_el  ->Write();
  W_PT_mu  ->Write();
  W_P_el   ->Write();
  W_P_mu   ->Write();
  LepEta   ->Write();
  MuEta    ->Write();
  ElEta    ->Write();
  DiffPhi  ->Write();
  DiffPhi_e->Write();
  DiffPhi_m->Write();
  DiffPhic ->Write();
  TEST	   ->Write();
  TEST2    ->Write();
  TEST3    ->Write();

  fout->Close();

  //TCanvas* c1 = new TCanvas;
  // TLegend* l1 = new TLegend(0.4,0.1,1,0.4,"Legend");
  // l1->AddEntry(MT2_hist,"MT2","f");

  // MT2_ttbar->Draw();
  // // c1->BuildLegend(0.2, 0.1, 0.3, 0.2);
  // c1->SaveAs(Form("./hists/mt2%s.png", suffix));
  // h2_nozero->Draw();
  // c1->SaveAs(Form("./hists/mt2_nonzero%s.png", suffix));

  // cout << "Entries FMet vs Met: " << FMet->GetEntries() << "\t" << Met_all->GetEntries() << endl;


  // TCanvas* c3 = new TCanvas;		   
  // W_P->Draw();
  //MT2_hist->SetLineStyle(1);			   
  //MT2_hist->Draw("hist");			   

  // TCanvas* c4 = new TCanvas;		   
  //RWPT->Draw();

  //MT2_Cut50->Draw("same");
  //MT2_Cut50->SetLineColor(kRed+1);
  // c3->SaveAs(Form("./hists/mt2_emuonly%s.png", suffix));
					   
  //  c1->SetLogy();			   
  // c1->Close();			   
  					   

  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << "Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:   " << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:   " << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;

  return 0;
 }

