#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TTree.h"
#include"TBrowser.h"
#include"TF1.h"
#include<iomanip>
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"
#include"TGraph.h"
#include"TCanvas.h"
#include"TPaveStats.h"
#include"TStyle.h"
#include"TLegend.h"

//#include"Variables.h"

using namespace std;

void nameLegend2(const char*);
void decorate(TH1D*,int);
void compHists(TString);
void compHists();
// double getYmax(TH1D*);
// double getYmin(TH1D*);
TList *FileList;
vector<TString> legName;
vector<TString> histName1,histName2;
TLatex textOnTop,intLumiE;
bool drawNorm=1;
TLegend *legend1=new TLegend(0.85,0.55,0.95,0.9);
TString xAxisLabel;
char name[100];
int col[7]={kBlue,kRed,kOrange-3,kTeal+9,kMagenta,kOrange,kCyan};

int rebin=1;
void compHists(){
  TString histName="NHF";
  compHists(histName);
}
void compHists(TString histName){
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  FileList = new TList();
  int nfiles=7;
  TFile *f[nfiles];

  xAxisLabel="";

  f[0] = new TFile("Run2016B_03Feb2017_HTMHT.root");
  f[1] = new TFile("Run2016C_03Feb2017_HTMHT.root");
  f[2] = new TFile("Run2016D_03Feb2017_HTMHT.root");
  f[3] = new TFile("Run2016E_03Feb2017_HTMHT.root");
  f[4] = new TFile("Run2016F_03Feb2017_HTMHT.root");
  f[5] = new TFile("Run2016G_03Feb2017_HTMHT.root");
  f[6] = new TFile("Run2016H_03Feb2017_HTMHT.root");

  //  f[0] = new TFile("Run2016C_03Feb2017_HTMHT.root");
  
  // histName1.push_back("PhotonPt_PhoPt");
  // histName2.push_back("PhotonPt_HT");
  // histName1.push_back("HT_HT");
  // histName2.push_back("HT_PhoPt");
  for(int i=0;i<nfiles;i++){
    //    histName1.push_back("JetPt_JetID1");
    //    histName1.push_back("NHF");
    //histName1.push_back("NEMF");
    //    histName1.push_back("NM_Pt435To566_Eta1p044To1p305");
    histName1.push_back(histName);
  }

  // histName1.push_back("MET_PhoPt"); 
  // histName2.push_back("MET_HT"); 
  rebin=10;
  //  legName.push_back("2016C to 2016H (35.1 fb^{-1})");
  //  legName.push_back("2016 (35.9 fb^{-1})");
  legName.push_back("2016B");
  legName.push_back("2016C");
  legName.push_back("2016D");
  legName.push_back("2016E");
  legName.push_back("2016F");
  legName.push_back("2016G");
  legName.push_back("2016H");

  TCanvas *cA = new TCanvas("his","hist",1500,850);
  cA->SetTickx();cA->SetTicky();
  TH1D *h1[nfiles],*h2[nfiles];

  double ymax=0.0,ymin=0.0;
  
  for(int i=0;i<nfiles;i++){
    h1[i] = (TH1D*)f[i]->FindObjectAny(histName1[i]);
    h1[i]->Rebin(rebin);
    
    cA->cd();
    cA->SetGridx();
    cA->SetGridy();
    decorate(h1[i],i);
    if(i==0) {
      h1[i]->SetTitle(0);
      h1[i]->GetXaxis()->SetTitle(histName1[i]);
      //      h1[i]->GetYaxis()->SetTitle("Events");
      h1[i]->GetYaxis()->SetTitleOffset(0.5);
      h1[i]->GetYaxis()->SetTitleSize(0.06);
      if(drawNorm)  h1[i]->Scale(1.0/h1[i]->Integral());
      h1[i]->Draw();
    }
    else{
      if(drawNorm)  h1[i]->Scale(1.0/h1[i]->Integral());
      h1[i]->Draw("sames");
    }
    //++++++++++++++++++++++++++ get minimum and maximum of histograms +++++++++++++++++++++++++++++++
    if(ymax < h1[i]->GetBinContent(h1[i]->GetMaximumBin())+1.2*h1[i]->GetBinError(h1[i]->GetMaximumBin()) ) 
      ymax = h1[i]->GetBinContent(h1[i]->GetMaximumBin())+1.2*h1[i]->GetBinError(h1[i]->GetMaximumBin());
   
    if(ymin > h1[i]->GetBinContent(h1[i]->GetMinimumBin())-1.2*h1[i]->GetBinError(h1[i]->GetMinimumBin()) ) 
      ymin = h1[i]->GetBinContent(h1[i]->GetMinimumBin())-1.2*h1[i]->GetBinError(h1[i]->GetMinimumBin());
    if(i==nfiles-1) h1[0]->SetMaximum(ymax);
    if(i==nfiles-1) h1[0]->SetMinimum(ymin);
    //-------------------------- get minimum and maximum of histograms -------------------------------
    
    legend1->AddEntry(h1[i],legName[i],"lpe");
    gPad->Update();
  }
  legend1->Draw();
  cA->cd();

  char name2[100];  
  textOnTop.SetTextSize(0.035);
  intLumiE.SetTextSize(0.035);
  textOnTop.DrawLatexNDC(0.1,0.91,"CMS #it{#bf{Preliminary}}");
  sprintf(name2,"#bf{35.9 fb^{-1}(13TeV)}");
  intLumiE.DrawLatexNDC(0.73,0.91,name2);
  cA->SaveAs("plots_rebin10/"+histName1[0]+".pdf");
}

  //  FileList->Add(TFile::Open("b.root") );  histName.push_back(name); legName.push_back("HE_RBX_plus0p10");
   
void decorate(TH1D* h,int i){
  h->SetLineColor(col[i]);
  h->SetLineWidth(2);
  //  h->SetFillColor(h->GetLineColor());
  h->SetMarkerStyle(20+i);
  h->SetMarkerColor(col[i]);
}

// double getYmax(TH1D* h1){
//   double ymax=0.0;
//   for(int i=0;i<h1->GetNbinsX();i++){
//     if(ymax<h1->GetBinContent
//   }
// }
