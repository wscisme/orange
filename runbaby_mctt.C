{
  gSystem->AddIncludePath(Form("-I%s/CORE", gSystem->Getenv("HOME")));
  gSystem->AddIncludePath(Form("-I%s/CORE/MT2", gSystem->Getenv("HOME")));
  gSystem->Load(Form("%s/CORE/libCMS2NtupleMacrosCORE.so", gSystem->Getenv("HOME")));  
  gROOT->ProcessLine(".L NSBABY.cc+");
  gROOT->ProcessLine(".L orange_2.C+");
    // gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  //gROOT->ProcessLine(".L /home/users/jgran/CMSSW_5_3_2_patch4_V05-03-23/src/CMS2/NtupleMacros/Tools/goodrun.cc");

  TChain *ch = new TChain("tree"); 
  // ch->Add("/hadoop/cms/store/user/sicheng/babies/ttbar_sl_4.root");
  // ch->Add("/home/users/sicheng/makebaby/bob/babies/ttbar_g2jg1b_4.root");
  ch->Add("/home/users/sicheng/makebaby/bob/babies/ttbar_g1j0b_8.root");

  ScanChain(ch, "_omctt_");

}
