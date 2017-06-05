#ifndef LostMuon_H
#define LostMuon_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH3.h"
#include <TH2Poly.h>
#include <TProfile.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

class LostMuon : public NtupleVariables{

 public:
  LostMuon(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~LostMuon();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  TLorentzVector getBestPhoton();
  bool check_eMatchedtoGamma();
  double getGendRLepPho();
  void print(Long64_t);
  void findObjMatchedtoG(TLorentzVector);
  //Variables defined
  int bestPhotonIndxAmongPhotons=-100;
  double gendRLepPho=1000.;
  TLorentzVector bestPhoton;//(0.,0.,0.,0.);
    
  float HT_PtCut=30;
  float MHT_PtCut=30;//keep MHT_PtCut <= HT_PtCut and <= Njets_PtCut
  float Njets_PtCut=30;

  float HT_EtaCut=2.4;
  float MHT_EtaCut=5;
  float Njets_EtaCut=2.4;
  double wt=0,lumiInfb=35.86;//36.814;
  bool isSignal=false;
  vector<double> METBinLowEdge={0,20,40,60,80,100,125,160,200,270,350,500,600};
  vector<double> ptBin={51,73,129,163,230,299,365,435,566,1000,10000};
  vector<double> etaBin={0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
  vector<string> ptR={"51","73","129","163","230","299","365","435","566","1000","Inf"};
  vector<string> etaR={"0p000", "0p261", "0p522", "0p783", "1p044", "1p305", "1p479", "1p653", "1p930", "2p172", "2p322", "2p500", "2p650", "2p853", "2p964", "3p139", "3p489", "3p839", "5p191"};
  static const int pts=10,etas=18;
  //histograms
  TH1I *h_RunNum;
  TH1D *h_intLumi;
  TH1D *h_TrigPS;
  TH1D *h_TrigDec;

  TH1D *h_JetPt_JetID1;
  TH1D *h_JetEta_JetID1;
  TH1D *h_NEMF_JetID1;
  TH1D *h_CEMF_JetID1;
  TH1D *h_NHF_JetID1;
  TH1D *h_CHF_JetID1;

  TH1D *h_CM_JetID1;
  TH1D *h_CEleM_JetID1;
  TH1D *h_CHadM_JetID1;

  TH1D *h_NHM_JetID1;
  TH1D *h_NM_JetID1;

  TH1D *h_PhoMult_JetID1;
  TH1D *h_PhoFrac_JetID1;
  //--------
  TH1D *h_JetPt_JetID0;
  TH1D *h_JetEta_JetID0;
  TH1D *h_NEMF_JetID0;
  TH1D *h_CEMF_JetID0;
  TH1D *h_NHF_JetID0;
  TH1D *h_CHF_JetID0;

  TH1D *h_CM_JetID0;
  TH1D *h_CEleM_JetID0;
  TH1D *h_CHadM_JetID0;

  TH1D *h_NHM_JetID0;
  TH1D *h_NM_JetID0;

  TH1D *h_PhoMult_JetID0;
  TH1D *h_PhoFrac_JetID0;
  //-------------
  TH1D *h_JetPt;
  TH1D *h_JetEta;
  TH1D *h_NEMF;
  TH1D *h_CEMF;
  TH1D *h_NHF;
  TH1D *h_CHF;

  TH1D *h_CM;
  TH1D *h_CEleM;
  TH1D *h_CHadM;

  TH1D *h_NHM;
  TH1D *h_NM;

  TH1D *h_PhoMult;
  TH1D *h_PhoFrac;
  //-----------
  TH1D *h_JetPt_R[pts][etas];
  TH1D *h_JetEta_R[pts][etas];
  TH1D *h_NEMF_R[pts][etas];
  TH1D *h_CEMF_R[pts][etas];
  TH1D *h_NHF_R[pts][etas];
  TH1D *h_CHF_R[pts][etas];

  TH1D *h_CM_R[pts][etas];
  TH1D *h_CEleM_R[pts][etas];
  TH1D *h_CHadM_R[pts][etas];

  TH1D *h_NHM_R[pts][etas];
  TH1D *h_NM_R[pts][etas];

  TH1D *h_PhoMult_R[pts][etas];
  TH1D *h_PhoFrac_R[pts][etas];
  //-----------------
  TFile *oFile;
 
};
#endif

#ifdef LostMuon_cxx

void LostMuon::BookHistogram(const char *outFileName) {

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;
  char name[100],title[100];
 
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  TH2::SetDefaultSumw2(1);
  h_RunNum=new TH1I("runs","Run nos.",300000,0,300000);
  h_intLumi=new TH1D("intLumi","integrated luminosity in /fb",2500,25,50);
  h_TrigPS=new TH1D("TrigPS","trigger PS",1000,0,1000);
  h_TrigDec=new TH1D("TrigDec","trig dec 0:fail,1:pass,-1:not run",10,-5,5);

  h_JetPt_JetID1=new TH1D("JetPt_JetID1","JetPt_JetID1",200,0,2000);
  h_JetEta_JetID1=new TH1D("JetEta_JetID1","JetEta_JetID1",120,-6,6);
  h_NEMF_JetID1=new TH1D("NEMF_JetID1","NEMF_JetID1",220,0,1.1);
  h_CEMF_JetID1=new TH1D("CEMF_JetID1","CEMF_JetID1",220,0,1.1);
  h_NHF_JetID1=new TH1D("NHF_JetID1","NHF_JetID1",220,0,1.1);
  h_CHF_JetID1=new TH1D("CHF_JetID1","CHF_JetID1",220,0,1.1);

  h_CM_JetID1=new TH1D("CM_JetID1","CM_JetID1",50,0,50);
  h_CEleM_JetID1=new TH1D("CEleM_JetID1","CEleM_JetID1",50,0,50);
  h_CHadM_JetID1=new TH1D("CHadM_JetID1","CHadM_JetID1",50,0,50);

  h_NHM_JetID1=new TH1D("NHM_JetID1","NHM_JetID1",50,0,50);
  h_NM_JetID1=new TH1D("NM_JetID1","NM_JetID1",50,0,50);

  h_PhoMult_JetID1=new TH1D("PhoMult_JetID1","PhoMult_JetID1",50,0,50);
  h_PhoFrac_JetID1=new TH1D("PhoFrac_JetID1","PhoFrac_JetID1",220,0,1.1);
  //--------------
  h_JetPt_JetID0=new TH1D("JetPt_JetID0","JetPt_JetID0",200,0,2000);
  h_JetEta_JetID0=new TH1D("JetEta_JetID0","JetEta_JetID0",120,-6,6);
  h_NEMF_JetID0=new TH1D("NEMF_JetID0","NEMF_JetID0",220,0,1.1);
  h_CEMF_JetID0=new TH1D("CEMF_JetID0","CEMF_JetID0",220,0,1.1);
  h_NHF_JetID0=new TH1D("NHF_JetID0","NHF_JetID0",220,0,1.1);
  h_CHF_JetID0=new TH1D("CHF_JetID0","CHF_JetID0",220,0,1.1);

  h_CM_JetID0=new TH1D("CM_JetID0","CM_JetID0",50,0,50);
  h_CEleM_JetID0=new TH1D("CEleM_JetID0","CEleM_JetID0",50,0,50);
  h_CHadM_JetID0=new TH1D("CHadM_JetID0","CHadM_JetID0",50,0,50);

  h_NHM_JetID0=new TH1D("NHM_JetID0","NHM_JetID0",50,0,50);
  h_NM_JetID0=new TH1D("NM_JetID0","NM_JetID0",50,0,50);

  h_PhoMult_JetID0=new TH1D("PhoMult_JetID0","PhoMult_JetID0",50,0,50);
  h_PhoFrac_JetID0=new TH1D("PhoFrac_JetID0","PhoFrac_JetID0",220,0,1.1);
  //-------------
  h_JetPt=new TH1D("JetPt","JetPt",200,0,2000);
  h_JetEta=new TH1D("JetEta","JetEta",120,-6,6);
  h_NEMF=new TH1D("NEMF","NEMF",220,0,1.1);
  h_CEMF=new TH1D("CEMF","CEMF",220,0,1.1);
  h_NHF=new TH1D("NHF","NHF",220,0,1.1);
  h_CHF=new TH1D("CHF","CHF",220,0,1.1);

  h_CM=new TH1D("CM","CM",50,0,50);
  h_CEleM=new TH1D("CEleM","CEleM",50,0,50);
  h_CHadM=new TH1D("CHadM","CHadM",50,0,50);

  h_NHM=new TH1D("NHM","NHM",50,0,50);
  h_NM=new TH1D("NM","NM",50,0,50);

  h_PhoMult=new TH1D("PhoMult","PhoMult",50,0,50);
  h_PhoFrac=new TH1D("PhoFrac","PhoFrac",220,0,1.1);
  //----------------
  for(int i=0;i<pts;i++){
    for(int j=0;j<etas;j++){
      TString st1;
      st1="JetPt_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_JetPt_R[i][j]=new TH1D(st1,st1,200,0,2000);
      
      st1="JetEta_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_JetEta_R[i][j]=new TH1D(st1,st1,120,-6,6);

      st1="NEMF_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_NEMF_R[i][j]=new TH1D(st1,st1,220,0,1.1);
      //      cout<<st1<<" ";

      st1="CEMF_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_CEMF_R[i][j]=new TH1D(st1,st1,220,0,1.1);
      //      cout<<st1<<" ";

      st1="NHF_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_NHF_R[i][j]=new TH1D(st1,st1,220,0,1.1);
      //      cout<<st1<<" ";

      st1="CHF_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_CHF_R[i][j]=new TH1D(st1,st1,220,0,1.1);
      //      cout<<st1<<" ";
      
      st1="CM_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_CM_R[i][j]=new TH1D(st1,st1,50,0,50);

      st1="CEleM_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_CEleM_R[i][j]=new TH1D(st1,st1,50,0,50);

      st1="CHadM_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_CHadM_R[i][j]=new TH1D(st1,st1,50,0,50);
      
      st1="NHM_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_NHM_R[i][j]=new TH1D(st1,st1,50,0,50);

      st1="NM_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_NM_R[i][j]=new TH1D(st1,st1,50,0,50);

      st1="PhoMult_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];      
      h_PhoMult_R[i][j]=new TH1D(st1,st1,50,0,50);

      st1="PhoFrac_Pt"+ptR[i]+"To"+ptR[i+1]+"_Eta"+etaR[j]+"To"+etaR[j+1];
      h_PhoFrac_R[i][j]=new TH1D(st1,st1,220,0,1.1);
      
      //      sprintf(name,"JetPt_Pt%sTo%s_Eta%sTo%s",ptR[i],ptR[i+1],etaR[j],etaR[j+1]);
      //      cout<<st1<<endl;
    }
  }

}


LostMuon::LostMuon(const TString &inputFileList, const char *outFileName, const char* dataset) {
  string nameData=dataset;//vvv
  TChain *tree = new TChain("PreSelection");
  tree = new TChain("TreeMaker2/PreSelection");//vvv
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  if(nameData!="signalH") nameData="BG";
  if(nameData=="signalH") nameData="signal";
  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree,nameData);

  BookHistogram(outFileName);
  
}

Bool_t LostMuon::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t LostMuon::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

LostMuon::~LostMuon() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif

