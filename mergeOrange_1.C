#include "TStyle.h"
#include "TFile.h"
#include "dataMCplotMaker.h"

//int mergeHist(char* datafile, char* mcfile)
int mergeOrange_1()
{
  // Get File 
  TFile* data = new TFile("/home/users/sicheng/build/hists/hists_odata_40met0b30lpthlt_8.root");
  TFile* mcwj = new TFile("/home/users/sicheng/build/hists/hists_omcwj_40met0b30lpthlt_8.root");
  TFile* mctt = new TFile("/home/users/sicheng/build/hists/hists_omctt_40met0b30lpthlt_8.root");

  // Param Defining
  float lumi = 5.2;		// for AB files

  // Hists Getting
  
  TH1F* data_MT2   = (TH1F*) data->Get("MT2"); 
  TH1F* mcwj_MT2   = (TH1F*) mcwj->Get("MT2"); 
  TH1F* mctt_MT2   = (TH1F*) mctt->Get("MT2"); 

  TH1F* data_MT2e  = (TH1F*) data->Get("MT2e"); 
  TH1F* mcwj_MT2e  = (TH1F*) mcwj->Get("MT2e"); 
  TH1F* mctt_MT2e  = (TH1F*) mctt->Get("MT2e"); 

  TH1F* data_MT2m  = (TH1F*) data->Get("MT2m"); 
  TH1F* mcwj_MT2m  = (TH1F*) mcwj->Get("MT2m"); 
  TH1F* mctt_MT2m  = (TH1F*) mctt->Get("MT2m"); 

  
  TH1F* data_hnz   = (TH1F*) data->Get("hnz"); 
  TH1F* mcwj_hnz   = (TH1F*) mcwj->Get("hnz"); 
  TH1F* mctt_hnz   = (TH1F*) mctt->Get("hnz"); 

  TH1F* data_h50   = (TH1F*) data->Get("h50"); 
  TH1F* mcwj_h50   = (TH1F*) mcwj->Get("h50"); 
  TH1F* mctt_h50   = (TH1F*) mctt->Get("h50"); 

  TH1F* data_jeta  = (TH1F*) data->Get("jeta"); 
  TH1F* mcwj_jeta  = (TH1F*) mcwj->Get("jeta"); 
  TH1F* mctt_jeta  = (TH1F*) mctt->Get("jeta"); 
		   
  TH1F* data_jetb  = (TH1F*) data->Get("jetb"); 
  TH1F* mcwj_jetb  = (TH1F*) mcwj->Get("jetb"); 
  TH1F* mctt_jetb  = (TH1F*) mctt->Get("jetb"); 

  TH1F* data_jae   = (TH1F*) data->Get("jae"); 
  TH1F* mcwj_jae   = (TH1F*) mcwj->Get("jae"); 
  TH1F* mctt_jae   = (TH1F*) mctt->Get("jae"); 
  TH1F* data_jbe   = (TH1F*) data->Get("jbe"); 
  TH1F* mcwj_jbe   = (TH1F*) mcwj->Get("jbe"); 
  TH1F* mctt_jbe   = (TH1F*) mctt->Get("jbe"); 
  TH1F* data_jam   = (TH1F*) data->Get("jam"); 
  TH1F* mcwj_jam   = (TH1F*) mcwj->Get("jam"); 
  TH1F* mctt_jam   = (TH1F*) mctt->Get("jam"); 
  TH1F* data_jbm   = (TH1F*) data->Get("jbm"); 
  TH1F* mcwj_jbm   = (TH1F*) mcwj->Get("jbm"); 
  TH1F* mctt_jbm   = (TH1F*) mctt->Get("jbm"); 


  TH1F* data_Mll   = (TH1F*) data->Get("invm"); 
  TH1F* mcwj_Mll   = (TH1F*) mcwj->Get("invm"); 
  TH1F* mctt_Mll   = (TH1F*) mctt->Get("invm"); 

  TH1F* data_cmet  = (TH1F*) data->Get("cmet"); 
  TH1F* mcwj_cmet  = (TH1F*) mcwj->Get("cmet"); 
  TH1F* mctt_cmet  = (TH1F*) mctt->Get("cmet"); 

  TH1F* data_mete  = (TH1F*) data->Get("mete"); 
  TH1F* mcwj_mete  = (TH1F*) mcwj->Get("mete"); 
  TH1F* mctt_mete  = (TH1F*) mctt->Get("mete"); 
  TH1F* data_metm  = (TH1F*) data->Get("metm"); 
  TH1F* mcwj_metm  = (TH1F*) mcwj->Get("metm"); 
  TH1F* mctt_metm  = (TH1F*) mctt->Get("metm"); 


  TH1F* data_fmet  = (TH1F*) data->Get("fmet"); 
  TH1F* mcwj_fmet  = (TH1F*) mcwj->Get("fmet"); 
  TH1F* mctt_fmet  = (TH1F*) mctt->Get("fmet"); 

  TH1F* data_mt    = (TH1F*) data->Get("MT"); 
  TH1F* mcwj_mt    = (TH1F*) mcwj->Get("MT"); 
  TH1F* mctt_mt    = (TH1F*) mctt->Get("MT"); 
  TH1F* data_mtc   = (TH1F*) data->Get("MTc"); 
  TH1F* mcwj_mtc   = (TH1F*) mcwj->Get("MTc"); 
  TH1F* mctt_mtc   = (TH1F*) mctt->Get("MTc"); 

  TH1F* data_mte   = (TH1F*) data->Get("MTe"); 
  TH1F* mcwj_mte   = (TH1F*) mcwj->Get("MTe"); 
  TH1F* mctt_mte   = (TH1F*) mctt->Get("MTe"); 
  TH1F* data_mtm   = (TH1F*) data->Get("MTm"); 
  TH1F* mcwj_mtm   = (TH1F*) mcwj->Get("MTm"); 
  TH1F* mctt_mtm   = (TH1F*) mctt->Get("MTm"); 

  TH1F* data_wpt   = (TH1F*) data->Get("rwpt"); 
  TH1F* mcwj_wpt   = (TH1F*) mcwj->Get("rwpt"); 
  TH1F* mctt_wpt   = (TH1F*) mctt->Get("rwpt"); 
  TH1F* data_wpte  = (TH1F*) data->Get("wpte"); 
  TH1F* mcwj_wpte  = (TH1F*) mcwj->Get("wpte"); 
  TH1F* mctt_wpte  = (TH1F*) mctt->Get("wpte"); 
  TH1F* data_wptm  = (TH1F*) data->Get("wptm"); 
  TH1F* mcwj_wptm  = (TH1F*) mcwj->Get("wptm"); 
  TH1F* mctt_wptm  = (TH1F*) mctt->Get("wptm"); 

  TH1F* data_rwp   = (TH1F*) data->Get("rwp"); 
  TH1F* mcwj_rwp   = (TH1F*) mcwj->Get("rwp"); 
  TH1F* mctt_rwp   = (TH1F*) mctt->Get("rwp"); 
  TH1F* data_rwpe  = (TH1F*) data->Get("rwpe"); 
  TH1F* mcwj_rwpe  = (TH1F*) mcwj->Get("rwpe"); 
  TH1F* mctt_rwpe  = (TH1F*) mctt->Get("rwpe"); 
  TH1F* data_rwpm  = (TH1F*) data->Get("rwpm"); 
  TH1F* mcwj_rwpm  = (TH1F*) mcwj->Get("rwpm"); 
  TH1F* mctt_rwpm  = (TH1F*) mctt->Get("rwpm"); 

  TH1F* data_lpt   = (TH1F*) data->Get("lpt"); 
  TH1F* mcwj_lpt   = (TH1F*) mcwj->Get("lpt"); 
  TH1F* mctt_lpt   = (TH1F*) mctt->Get("lpt"); 
  TH1F* data_mpt   = (TH1F*) data->Get("mpt"); 
  TH1F* mcwj_mpt   = (TH1F*) mcwj->Get("mpt"); 
  TH1F* mctt_mpt   = (TH1F*) mctt->Get("mpt"); 

  TH1F* data_dphi  = (TH1F*) data->Get("dphi"); 
  TH1F* mcwj_dphi  = (TH1F*) mcwj->Get("dphi"); 
  TH1F* mctt_dphi  = (TH1F*) mctt->Get("dphi"); 
  TH1F* data_dpc   = (TH1F*) data->Get("dphic"); 
  TH1F* mcwj_dpc   = (TH1F*) mcwj->Get("dphic"); 
  TH1F* mctt_dpc   = (TH1F*) mctt->Get("dphic"); 
  TH1F* data_dpe   = (TH1F*) data->Get("dpe"); 
  TH1F* mcwj_dpe   = (TH1F*) mcwj->Get("dpe"); 
  TH1F* mctt_dpe   = (TH1F*) mctt->Get("dpe"); 
  TH1F* data_dpm   = (TH1F*) data->Get("dpm"); 
  TH1F* mcwj_dpm   = (TH1F*) mcwj->Get("dpm"); 
  TH1F* mctt_dpm   = (TH1F*) mctt->Get("dpm"); 


  TH1F* data_test  = (TH1F*) data->Get("test"); 
  TH1F* mcwj_test  = (TH1F*) mcwj->Get("test"); 
  TH1F* mctt_test  = (TH1F*) mctt->Get("test"); 

  TH1F* data_etam  = (TH1F*) data->Get("muEta"); 
  TH1F* data_etae  = (TH1F*) data->Get("elEta"); 
  TH1F* mcwj_etam  = (TH1F*) mcwj->Get("muEta"); 
  TH1F* mcwj_etae  = (TH1F*) mcwj->Get("elEta"); 
  TH1F* mctt_etam  = (TH1F*) mctt->Get("muEta"); 
  TH1F* mctt_etae  = (TH1F*) mctt->Get("elEta"); 

  
  // MC hists Scaling to lumi
  mctt_MT2 ->Scale(lumi);        mcwj_MT2 ->Scale(lumi);
  mctt_MT2e->Scale(lumi);        mcwj_MT2e->Scale(lumi);
  mctt_MT2m->Scale(lumi);        mcwj_MT2m->Scale(lumi);
  mctt_hnz ->Scale(lumi);        mcwj_hnz ->Scale(lumi);
  mctt_h50 ->Scale(lumi);        mcwj_h50 ->Scale(lumi);
  mctt_jeta->Scale(lumi);        mcwj_jeta->Scale(lumi);
  mctt_jetb->Scale(lumi);        mcwj_jetb->Scale(lumi);
  mctt_jae ->Scale(lumi);        mcwj_jae ->Scale(lumi);
  mctt_jam ->Scale(lumi);        mcwj_jam ->Scale(lumi);
  mctt_jbe ->Scale(lumi);        mcwj_jbe ->Scale(lumi);
  mctt_jbm ->Scale(lumi);        mcwj_jbm ->Scale(lumi);
  mctt_Mll ->Scale(lumi);        mcwj_Mll ->Scale(lumi);
  mctt_cmet->Scale(lumi);        mcwj_cmet->Scale(lumi);
  mctt_mete->Scale(lumi);        mcwj_mete->Scale(lumi);
  mctt_metm->Scale(lumi);        mcwj_metm->Scale(lumi);
  mctt_fmet->Scale(lumi);        mcwj_fmet->Scale(lumi);
  mctt_mt  ->Scale(lumi);        mcwj_mt  ->Scale(lumi);
  mctt_mtc ->Scale(lumi);        mcwj_mtc ->Scale(lumi);
  mctt_mte ->Scale(lumi);        mcwj_mte ->Scale(lumi);
  mctt_mtm ->Scale(lumi);        mcwj_mtm ->Scale(lumi);
  mctt_rwp ->Scale(lumi);        mcwj_rwp ->Scale(lumi);
  mctt_rwpe->Scale(lumi);        mcwj_rwpe->Scale(lumi);
  mctt_rwpm->Scale(lumi);        mcwj_rwpm->Scale(lumi);
  mctt_dphi->Scale(lumi);        mcwj_dphi->Scale(lumi);
  mctt_dpc ->Scale(lumi);        mcwj_dpc ->Scale(lumi);
  mctt_dpe ->Scale(lumi);        mcwj_dpe ->Scale(lumi);
  mctt_dpm ->Scale(lumi);        mcwj_dpm ->Scale(lumi);
  mctt_wpt ->Scale(lumi);        mcwj_wpt ->Scale(lumi);
  mctt_wpte->Scale(lumi);        mcwj_wpte->Scale(lumi);
  mctt_wptm->Scale(lumi);        mcwj_wptm->Scale(lumi);
  mctt_lpt ->Scale(lumi);        mcwj_lpt ->Scale(lumi);
  mctt_mpt ->Scale(lumi);        mcwj_mpt ->Scale(lumi);
  mctt_etam->Scale(lumi);        mcwj_etam->Scale(lumi);
  mctt_etae->Scale(lumi);        mcwj_etae->Scale(lumi);
  mctt_test->Scale(lumi);        mcwj_test->Scale(lumi);
    
  // MC Color Setting

  // mctt_MT2 ->SetFillColor(kRed+1);         mcwj_MT2 ->SetFillColor(kGreen+1);
  // mctt_hnz ->SetFillColor(kRed+1);         mcwj_hnz ->SetFillColor(kGreen+1);
  // mctt_h50 ->SetFillColor(kRed+1);         mcwj_h50 ->SetFillColor(kGreen+1);
  // mctt_jeta->SetFillColor(kRed+1);         mcwj_jeta->SetFillColor(kGreen+1);
  // mctt_jetb->SetFillColor(kRed+1);         mcwj_jetb->SetFillColor(kGreen+1);
  // mctt_Mll ->SetFillColor(kRed+1);         mcwj_Mll ->SetFillColor(kGreen+1);
  // mctt_cmet->SetFillColor(kRed+1);         mcwj_cmet->SetFillColor(kGreen+1);
  // mctt_fmet->SetFillColor(kRed+1);         mcwj_fmet->SetFillColor(kGreen+1);


  // MC Stacking
  
  std::vector<TH1F*> hmcs_MT2 ;
  std::vector<TH1F*> hmcs_MT2e;
  std::vector<TH1F*> hmcs_MT2m;
  std::vector<TH1F*> hmcs_hnz ;
  std::vector<TH1F*> hmcs_h50 ;
  std::vector<TH1F*> hmcs_jeta;
  std::vector<TH1F*> hmcs_jae ;
  std::vector<TH1F*> hmcs_jam ;
  std::vector<TH1F*> hmcs_jetb;
  std::vector<TH1F*> hmcs_jbe ;
  std::vector<TH1F*> hmcs_jbm ;
  std::vector<TH1F*> hmcs_Mll ;
  std::vector<TH1F*> hmcs_cmet;
  std::vector<TH1F*> hmcs_mete;
  std::vector<TH1F*> hmcs_metm;
  std::vector<TH1F*> hmcs_fmet;
  std::vector<TH1F*> hmcs_mt  ;
  std::vector<TH1F*> hmcs_mtc ;
  std::vector<TH1F*> hmcs_mte ;
  std::vector<TH1F*> hmcs_mtm ;
  std::vector<TH1F*> hmcs_rwp ;
  std::vector<TH1F*> hmcs_rwpe;
  std::vector<TH1F*> hmcs_rwpm;
  std::vector<TH1F*> hmcs_dphi;
  std::vector<TH1F*> hmcs_dpc ;
  std::vector<TH1F*> hmcs_dpe ;
  std::vector<TH1F*> hmcs_dpm ;
  std::vector<TH1F*> hmcs_wpt ;
  std::vector<TH1F*> hmcs_wpte;
  std::vector<TH1F*> hmcs_wptm;
  std::vector<TH1F*> hmcs_lpt ;
  std::vector<TH1F*> hmcs_mpt ;
  std::vector<TH1F*> hmcs_etam;
  std::vector<TH1F*> hmcs_etae;
  std::vector<TH1F*> hmcs_test;
  
  hmcs_MT2 .push_back(mctt_MT2 );	  hmcs_MT2 .push_back(mcwj_MT2 );
  hmcs_MT2e.push_back(mctt_MT2e);	  hmcs_MT2e.push_back(mcwj_MT2e);
  hmcs_MT2m.push_back(mctt_MT2m);	  hmcs_MT2m.push_back(mcwj_MT2m);
  hmcs_hnz .push_back(mctt_hnz );	  hmcs_hnz .push_back(mcwj_hnz );
  hmcs_h50 .push_back(mctt_h50 );	  hmcs_h50 .push_back(mcwj_h50 );
  hmcs_jeta.push_back(mctt_jeta);	  hmcs_jeta.push_back(mcwj_jeta);
  hmcs_jae .push_back(mctt_jae );	  hmcs_jae .push_back(mcwj_jae );
  hmcs_jam .push_back(mctt_jam );	  hmcs_jam .push_back(mcwj_jam );
  hmcs_jetb.push_back(mctt_jetb);	  hmcs_jetb.push_back(mcwj_jetb);
  hmcs_jbe .push_back(mctt_jbe );	  hmcs_jbe .push_back(mcwj_jbe );
  hmcs_jbm .push_back(mctt_jbm );	  hmcs_jbm .push_back(mcwj_jbm );
  hmcs_Mll .push_back(mctt_Mll );	  hmcs_Mll .push_back(mcwj_Mll );
  hmcs_cmet.push_back(mctt_cmet);	  hmcs_cmet.push_back(mcwj_cmet);
  hmcs_mete.push_back(mctt_mete);	  hmcs_mete.push_back(mcwj_mete);
  hmcs_metm.push_back(mctt_metm);	  hmcs_metm.push_back(mcwj_metm);
  hmcs_fmet.push_back(mctt_fmet);	  hmcs_fmet.push_back(mcwj_fmet);
  hmcs_mt  .push_back(mctt_mt  );	  hmcs_mt  .push_back(mcwj_mt  );
  hmcs_mtc .push_back(mctt_mtc );	  hmcs_mtc .push_back(mcwj_mtc );
  hmcs_mte .push_back(mctt_mte );	  hmcs_mte .push_back(mcwj_mte );
  hmcs_mtm .push_back(mctt_mtm );	  hmcs_mtm .push_back(mcwj_mtm );
  hmcs_rwp .push_back(mctt_rwp );	  hmcs_rwp .push_back(mcwj_rwp );
  hmcs_rwpe.push_back(mctt_rwpe);	  hmcs_rwpe.push_back(mcwj_rwpe);
  hmcs_rwpm.push_back(mctt_rwpm);	  hmcs_rwpm.push_back(mcwj_rwpm);
  hmcs_wpt .push_back(mctt_wpt );	  hmcs_wpt .push_back(mcwj_wpt );
  hmcs_wpte.push_back(mctt_wpte);	  hmcs_wpte.push_back(mcwj_wpte);
  hmcs_wptm.push_back(mctt_wptm);	  hmcs_wptm.push_back(mcwj_wptm);
  hmcs_dphi.push_back(mctt_dphi);	  hmcs_dphi.push_back(mcwj_dphi);
  hmcs_dpc .push_back(mctt_dpc );	  hmcs_dpc .push_back(mcwj_dpc );
  hmcs_dpe .push_back(mctt_dpe );	  hmcs_dpe .push_back(mcwj_dpe );
  hmcs_dpm .push_back(mctt_dpm );	  hmcs_dpm .push_back(mcwj_dpm );
  hmcs_lpt .push_back(mctt_lpt );	  hmcs_lpt .push_back(mcwj_lpt );
  hmcs_mpt .push_back(mctt_mpt );	  hmcs_mpt .push_back(mcwj_mpt );
  hmcs_etam.push_back(mctt_etam);	  hmcs_etam.push_back(mcwj_etam);
  hmcs_etae.push_back(mctt_etae);	  hmcs_etae.push_back(mcwj_etae);
  hmcs_test.push_back(mctt_test);	  hmcs_test.push_back(mcwj_test);

  
  // Titles setting
  std::vector<char*> titmcs;
  titmcs.push_back("ttbar");	  
  titmcs.push_back("Wjets");


  // Hists Drawing
  // TCanvas* c1 = new TCanvas;

  char* suffix = "40met0b30lpt_8";
  char* opts = "--noOverflow --noDivisionLabel --outputName ";
  
  // dataMCplotMaker(data_MT2 , hmcs_MT2,  titmcs, (char*)data_MT2 ->GetTitle(), suffix, Form("--xAxisLabel MT2 %s MT2_%s",    opts, suffix) );
  // dataMCplotMaker(data_MT2e, hmcs_MT2e, titmcs, (char*)data_MT2e->GetTitle(), suffix, Form("--xAxisLabel MT2 %s MT2_el_%s", opts, suffix) );
  // dataMCplotMaker(data_MT2m, hmcs_MT2m, titmcs, (char*)data_MT2m->GetTitle(), suffix, Form("--xAxisLabel MT2 %s MT2_mu_%s", opts, suffix) );
  // dataMCplotMaker(data_hnz , hmcs_hnz,  titmcs, (char*)data_hnz ->GetTitle(), suffix, Form("--xAxisLabel MT2 %s hnz_%s",    opts, suffix) );
  // dataMCplotMaker(data_mt  , hmcs_mt ,  titmcs, (char*)data_mt  ->GetTitle(), suffix, Form("--xAxisLabel M_T %s MT_%s", opts, suffix) );
  // dataMCplotMaker(data_mtc , hmcs_mtc,  titmcs, (char*)data_mtc ->GetTitle(), suffix, Form("--xAxisLabel M_T %s MTc_%s", opts, suffix) );
  // dataMCplotMaker(data_mte , hmcs_mte,  titmcs, (char*)data_mte ->GetTitle(), suffix, Form("--xAxisLabel M_T %s MTe_%s", opts, suffix) );
  dataMCplotMaker(data_mtm , hmcs_mtm,  titmcs, (char*)data_mtm ->GetTitle(), suffix, Form("--xAxisLabel M_T %s MTm_%s", opts, suffix) );
  // dataMCplotMaker(data_jeta, hmcs_jeta, titmcs, (char*)data_jeta->GetTitle(), suffix, Form("--xAxisLabel nJets --noXaxisUnit %s jeta_%s", opts, suffix) );
  // dataMCplotMaker(data_jae , hmcs_jae,  titmcs, (char*)data_jae ->GetTitle(), suffix, Form("--xAxisLabel nJets --noXaxisUnit %s jae_%s", opts, suffix) );
  dataMCplotMaker(data_jam , hmcs_jam,  titmcs, (char*)data_jam ->GetTitle(), suffix, Form("--xAxisLabel nJets --noXaxisUnit %s jam_%s", opts, suffix) );
  // dataMCplotMaker(data_cmet, hmcs_cmet, titmcs, (char*)data_cmet->GetTitle(), suffix, Form("-xAxisLabel MET %s cmet_%s", opts, suffix) );
  // dataMCplotMaker(data_mete, hmcs_mete, titmcs, (char*)data_mete->GetTitle(), suffix, Form("-xAxisLabel MET %s mete_%s", opts, suffix) );
  dataMCplotMaker(data_metm, hmcs_metm, titmcs, (char*)data_metm->GetTitle(), suffix, Form("-xAxisLabel MET %s metm_%s", opts, suffix) );
  // dataMCplotMaker(data_rwp,  hmcs_rwp,  titmcs, (char*)data_rwp ->GetTitle(), suffix, Form("-xAxisLabel W_{p} %s wp_%s", opts, suffix) );
  // dataMCplotMaker(data_rwpe, hmcs_rwpe, titmcs, (char*)data_rwpe->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{p} %s rwpe_%s", opts, suffix) );
  // dataMCplotMaker(data_rwpm, hmcs_rwpm, titmcs, (char*)data_rwpm->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{p} %s rwpm_%s", opts, suffix) );
  //dataMCplotMaker(data_wpt , hmcs_wpt,  titmcs, (char*)data_wpt ->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{pt} %s wpt_l_%s", opts, suffix) );
  //dataMCplotMaker(data_wpte, hmcs_wpte, titmcs, (char*)data_wpte->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{pt} %s wpte_%s", opts, suffix) );
  // dataMCplotMaker(data_wptm, hmcs_wptm, titmcs, (char*)data_wptm->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{pt} %s wptm_%s", opts, suffix) );
  //dataMCplotMaker(data_dphi, hmcs_dphi, titmcs, (char*)data_dphi->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{pt} %s dphi_%s", opts, suffix) );
  // dataMCplotMaker(data_dpc , hmcs_dpc , titmcs, (char*)data_dpc ->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{pt} %s dphic_%s", opts, suffix) );
  //dataMCplotMaker(data_dpe , hmcs_dpe , titmcs, (char*)data_dpe ->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{pt} %s dpe_%s", opts, suffix) );
  // dataMCplotMaker(data_dpm , hmcs_dpm , titmcs, (char*)data_dpm ->GetTitle(), suffix, Form("--isLinear -xAxisLabel W_{pt} %s dpm_%s", opts, suffix) );
  //dataMCplotMaker(data_lpt , hmcs_lpt , titmcs, (char*)data_lpt ->GetTitle(), suffix, Form("-xAxisLabel l_{pt} %s lpt_%s", opts, suffix) );
  // dataMCplotMaker(data_mpt , hmcs_mpt , titmcs, (char*)data_mpt ->GetTitle(), suffix, Form("-xAxisLabel l_{pt} %s mpt_%s", opts, suffix) );
  // dataMCplotMaker(data_etam, hmcs_etam, titmcs, (char*)data_etam->GetTitle(), suffix, Form("-xAxisLabel Eta --isLinear --noXaxisUnit %s etam_l_%s", opts, suffix) );
  // dataMCplotMaker(data_etae, hmcs_etae, titmcs, (char*)data_etae->GetTitle(), suffix, Form("-xAxisLabel Eta --isLinear --noXaxisUnit %s etae_l_%s", opts, suffix) );
  //dataMCplotMaker(data_test, hmcs_test, titmcs, "Electron P_T", suffix, Form("-xAxisLabel l_{pt} %s ept_%s", opts, suffix) );

  cout << "MT2-:"  << data_MT2->Integral() << "  " << (hmcs_MT2[0]->Integral() + hmcs_MT2[1]->Integral()) << "  ("
       << data_MT2->Integral()/(hmcs_MT2[0]->Integral() + hmcs_MT2[1]->Integral()) << ") "<< endl;
  cout << "MT2m:"  << data_MT2m->Integral() << "  " << (hmcs_MT2m[0]->Integral() + hmcs_MT2m[1]->Integral()) << "  ("
       << data_MT2m->Integral()/(hmcs_MT2m[0]->Integral() + hmcs_MT2m[1]->Integral()) << ") "<< endl;
  cout << "mtc-:"  << data_mtc->Integral() << "  " << (hmcs_mtc[0]->Integral() + hmcs_mtc[1]->Integral()) << endl;
  cout << "mtm-:"  << data_mtm->Integral() << "  " << (hmcs_mtm[0]->Integral() + hmcs_mtm[1]->Integral()) << endl;
  cout << "metm:"  << data_metm->Integral() << "  " << (hmcs_metm[0]->Integral() + hmcs_metm[1]->Integral()) << endl;
  cout << "jeta:"  << data_jeta->Integral() << "  " << (hmcs_jeta[0]->Integral() + hmcs_jeta[1]->Integral()) << endl;
  cout << "wpt-:"  << data_wpt->Integral() << "  " << (hmcs_wpt[0]->Integral() + hmcs_wpt[1]->Integral()) << endl;
  cout << "lpt-:"  << data_lpt->Integral() << "  " << (hmcs_lpt[0]->Integral() + hmcs_lpt[1]->Integral()) << endl;
  cout << "mpt-:"  << data_mpt->Integral() << "  " << (hmcs_mpt[0]->Integral() + hmcs_mpt[1]->Integral()) << endl;
  cout << "etam:"  << data_etam->Integral() << "  " << (hmcs_etam[0]->Integral() + hmcs_etam[1]->Integral()) << endl;
  cout << "etae:"  << data_etae->Integral() << "  " << (hmcs_etae[0]->Integral() + hmcs_etae[1]->Integral()) << endl;
  //cout << data_test->Integral() << "  " << (hmcs_test[0]->Integral() + hmcs_test[1]->Integral()) << endl;


  // dataMCplotMaker(data_jetb, hmcs_jetb, titmcs, "# bTag distribution", suffix,		Form("%s  jetb_40met0b30w_4" );
  // dataMCplotMaker(data_jbe, hmcs_jbe, titmcs, "# bTag distribution for e events ", suffix,   Form("%s  jbe_40met0b30w_4" );
  // dataMCplotMaker(data_jbm, hmcs_jbm, titmcs, "# bTag distribution for mu events", suffix,   Form("%s  jbm_40met0b30w_4" );
  
  return 0;
}


  // h7_da->Rebin(2);
  // h7_mc->Rebin(2);
  // gStyle->SetErrorX(0);
  // h7_da->SetMarkerStyle(20);

