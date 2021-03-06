#define LostMuon_cxx
#include "LostMuon.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  LostMuon ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void LostMuon::EventLoop(const char *data,const char *inputFileList) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  string s_data=data;
  
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  
  int evtSurvived=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" <<endl;
    decade = k;
    // cout<<"j:"<<jentry<<" fcurrent:"<<fCurrent<<endl;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    bool process=true;
    if(!(CSCTightHaloFilter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0)) continue;
    bool passHTMHTtrig = false;
    int trigPS = -1;
    for(int i=0;i<TriggerNames->size();i++){
      string trgName=(*TriggerNames)[i];
      trgName.pop_back();
      if( trgName=="HLT_PFHT300_PFMET110_v"){
	trigPS=(*TriggerPrescales)[i];
	h_TrigPS->Fill(trigPS);
	h_TrigDec->Fill((*TriggerPass)[i]);
	if( (*TriggerPass)[i]==1 ) passHTMHTtrig = true;
      }
    }
    if(passHTMHTtrig){
      evtSurvived++;
      double wt=1.0;
      h_RunNum->Fill(RunNum);
      h_intLumi->Fill(lumiInfb);
      //      print(jentry);
      // h_MET->Fill(MET,wt);
      // h_nHadJets->Fill(nHadJets,wt);
      // h_BTags->Fill(BTags,wt);
      // h_HT->Fill(HT,wt);
      // h_MHT->Fill(MHT,wt);
      // h_nJets->Fill(NJets,wt);
      // h_METvBin->Fill(MET,wt);
      // h_nVtx->Fill(NVtx,wt);

      for(int i=0;i<Jets->size();i++){
	//	if( ((*Jets)[i].Pt() > 200) && (abs((*Jets)[i].Eta()) >= 2.7) && (abs((*Jets)[i].Eta()) <= 3.0) ){
	if( ((*Jets)[i].Pt() > 30) ){
	  if((*Jets_ID)[i]){
	    h_JetPt_JetID1->Fill((*Jets)[i].Pt(),wt);
	    h_JetEta_JetID1->Fill((*Jets)[i].Eta(),wt);
	    h_NEMF_JetID1->Fill((*Jets_neutralEmEnergyFraction)[i],wt);
	    h_CEMF_JetID1->Fill((*Jets_chargedEmEnergyFraction)[i],wt);
	    h_NHF_JetID1->Fill((*Jets_neutralHadronEnergyFraction)[i],wt);
	    h_CHF_JetID1->Fill((*Jets_chargedHadronEnergyFraction)[i],wt);

	    h_CM_JetID1->Fill((*Jets_chargedMultiplicity)[i],wt);
	    h_CEleM_JetID1->Fill((*Jets_electronMultiplicity)[i],wt);
	    h_CHadM_JetID1->Fill((*Jets_chargedHadronMultiplicity)[i],wt);
	  
	    h_NHM_JetID1->Fill((*Jets_neutralHadronMultiplicity)[i],wt);
	    h_NM_JetID1->Fill((*Jets_neutralMultiplicity)[i],wt);
	  
	    h_PhoMult_JetID1->Fill((*Jets_photonMultiplicity)[i],wt);
	    h_PhoFrac_JetID1->Fill((*Jets_photonEnergyFraction)[i],wt);
	  }
	  else{
	    h_JetPt_JetID0->Fill((*Jets)[i].Pt(),wt);
	    h_JetEta_JetID0->Fill((*Jets)[i].Eta(),wt);
	    h_NEMF_JetID0->Fill((*Jets_neutralEmEnergyFraction)[i],wt);
	    h_CEMF_JetID0->Fill((*Jets_chargedEmEnergyFraction)[i],wt);
	    h_NHF_JetID0->Fill((*Jets_neutralHadronEnergyFraction)[i],wt);
	    h_CHF_JetID0->Fill((*Jets_chargedHadronEnergyFraction)[i],wt);

	    h_CM_JetID0->Fill((*Jets_chargedMultiplicity)[i],wt);
	    h_CEleM_JetID0->Fill((*Jets_electronMultiplicity)[i],wt);
	    h_CHadM_JetID0->Fill((*Jets_chargedHadronMultiplicity)[i],wt);
	  
	    h_NHM_JetID0->Fill((*Jets_neutralHadronMultiplicity)[i],wt);
	    h_NM_JetID0->Fill((*Jets_neutralMultiplicity)[i],wt);
	  
	    h_PhoMult_JetID0->Fill((*Jets_photonMultiplicity)[i],wt);
	    h_PhoFrac_JetID0->Fill((*Jets_photonEnergyFraction)[i],wt);
	  }
	  h_JetPt->Fill((*Jets)[i].Pt(),wt);
	  h_JetEta->Fill((*Jets)[i].Eta(),wt);
	  h_NEMF->Fill((*Jets_neutralEmEnergyFraction)[i],wt);
	  h_CEMF->Fill((*Jets_chargedEmEnergyFraction)[i],wt);
	  h_NHF->Fill((*Jets_neutralHadronEnergyFraction)[i],wt);
	  h_CHF->Fill((*Jets_chargedHadronEnergyFraction)[i],wt);

	  h_CM->Fill((*Jets_chargedMultiplicity)[i],wt);
	  h_CEleM->Fill((*Jets_electronMultiplicity)[i],wt);
	  h_CHadM->Fill((*Jets_chargedHadronMultiplicity)[i],wt);
	  
	  h_NHM->Fill((*Jets_neutralHadronMultiplicity)[i],wt);
	  h_NM->Fill((*Jets_neutralMultiplicity)[i],wt);
	  
	  h_PhoMult->Fill((*Jets_photonMultiplicity)[i],wt);
	  h_PhoFrac->Fill((*Jets_photonEnergyFraction)[i],wt);

	  for(int j=0;j<pts;j++){
	    for(int k=0;k<etas;k++){
	      if( ((*Jets)[i].Pt() >= ptBin[j]) && ((*Jets)[i].Pt() < ptBin[j+1]) &&
		  (abs((*Jets)[i].Eta()) >= etaBin[k]) && (abs((*Jets)[i].Eta()) < etaBin[k+1]) ){
		h_JetPt_R[j][k]->Fill((*Jets)[i].Pt(),wt);
		h_JetEta_R[j][k]->Fill((*Jets)[i].Eta(),wt);
		h_NEMF_R[j][k]->Fill((*Jets_neutralEmEnergyFraction)[i],wt);
		h_CEMF_R[j][k]->Fill((*Jets_chargedEmEnergyFraction)[i],wt);
		h_NHF_R[j][k]->Fill((*Jets_neutralHadronEnergyFraction)[i],wt);
		h_CHF_R[j][k]->Fill((*Jets_chargedHadronEnergyFraction)[i],wt);

		h_CM_R[j][k]->Fill((*Jets_chargedMultiplicity)[i],wt);
		h_CEleM_R[j][k]->Fill((*Jets_electronMultiplicity)[i],wt);
		h_CHadM_R[j][k]->Fill((*Jets_chargedHadronMultiplicity)[i],wt);

		h_NHM_R[j][k]->Fill((*Jets_neutralHadronMultiplicity)[i],wt);
		h_NM_R[j][k]->Fill((*Jets_neutralMultiplicity)[i],wt);

		h_PhoMult_R[j][k]->Fill((*Jets_photonMultiplicity)[i],wt);
		h_PhoFrac_R[j][k]->Fill((*Jets_photonEnergyFraction)[i],wt);
	      }
	    }//for pt
	  }//for eta
	}//if jet.pt>30
      }//jets
    }//trigger
  }// loop over entries
  cout<<"Events Survied:"<<evtSurvived<<endl;
}


TLorentzVector LostMuon::getBestPhoton(){
  vector<TLorentzVector> goodPho;
  vector<int> goodPhoIndx;

  for(int iPho=0;iPho<Photons->size();iPho++){
    if( ((*Photons_fullID)[iPho]) && ((*Photons_hasPixelSeed)[iPho]<0.001) ) {
      goodPho.push_back( (*Photons)[iPho] );
      goodPhoIndx.push_back(iPho);
    }
  }
  int highPtIndx=-100;
  for(int i=0;i<goodPho.size();i++){
    if(i==0) highPtIndx=0;
    else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
  }

  if(highPtIndx>=0){
    bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
    return goodPho[highPtIndx];
  }
  else{
    bestPhotonIndxAmongPhotons = -100;
    TLorentzVector v0;return v0;
  }
}


bool LostMuon::check_eMatchedtoGamma(){
  if(bestPhotonIndxAmongPhotons>=0){
    for(int i=0;i<Electrons->size();i++){
      if( (*Photons)[bestPhotonIndxAmongPhotons].DeltaR( (*Electrons)[i] ) < 0.1){
	//	cout<<(*Electrons)[i].Pt()<<" "<<(*Electrons)[i].Eta()<<" "<<(*Electrons)[i].Phi()<<" "<<(*Photons)[bestPhotonIndxAmongPhotons].Pt()<<" "<<(*Photons)[bestPhotonIndxAmongPhotons].Eta()<<" "<<(*Photons)[bestPhotonIndxAmongPhotons].Phi()<<" dR:"<<(*Photons)[bestPhotonIndxAmongPhotons].DeltaR( (*Electrons)[i])<<endl;
	return true;
      }
    }
  }
  else
    return false;
}

double LostMuon::getGendRLepPho(){//MC only
  TLorentzVector genPho1,genLep1;
  /*  int leadGenPhoIdx=-100;
      for(int i=0;i<GenParticles->size();i++){
      if((*GenParticles)[i].Pt()!=0){
      if((abs((*GenParticles_PdgId)[i])==22) && ((abs((*GenParticles_ParentId)[i])<=25) || ((*GenParticles_ParentId)[i]==2212) ) && (*GenParticles_Status)[i]==1 ){
      if(genPho1.Pt() < (*GenParticles)[i].Pt()){
      leadGenPhoIdx = i;
      genPho1 = ((*GenParticles)[i]);
      }
      }
      if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
      if(genLep1.Pt() < ((*GenParticles)[i]).Pt()) genLep1 = ((*GenParticles)[i]);
      }
      }
      }*/ //for
  if(genPho1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(genPho1);
  else return 1000.0;
}

void  LostMuon::findObjMatchedtoG(TLorentzVector bestPhoton){//MC only
  /*  double dR=100;
      int match=-100;
      for(int i=0;i<GenParticles->size();i++){
      if((*GenParticles)[i].Pt()!=0){
      if(i==0){dR=DeltaR( bestPhoton.Eta(),bestPhoton.Phi(),(*GenParticles)[i].Eta(),(*GenParticles)[i].Phi() );}
      else if(dR > (DeltaR( bestPhoton.Eta(),bestPhoton.Phi(),(*GenParticles)[i].Eta(),(*GenParticles)[i].Phi())) ){
      dR=(DeltaR( bestPhoton.Eta(),bestPhoton.Phi(),(*GenParticles)[i].Eta(),(*GenParticles)[i].Phi()));
      match=i;
      }
      }
      }
      if(dR<0.1){
      h_GmatchedObj->Fill(abs((*GenParticles_PdgId)[match]),wt);
      if(abs((*GenParticles_PdgId)[match])==22)  h_PdgIdPhoParent->Fill(abs((*GenParticles_ParentId)[match]),wt);
      }
      else{
      h_GmatchedObj->Fill(0.0,wt);
      h_PdgIdPhoParent->Fill(0.0,wt);
      }
      //find obj matched to muon
      dR=100;match=-100;
      if(Muons->size()==1){
      for(int i=0;i<GenParticles->size();i++){
      if((*GenParticles)[i].Pt()!=0){
      if(i==0){dR=(*GenParticles)[i].DeltaR((*Muons)[0]);}
      else if(dR > ((*GenParticles)[i].DeltaR( (*Muons)[0]) ) ){
      dR= (*GenParticles)[i].DeltaR((*Muons)[0]);
      match=i;
      }
      }
      }
      if(dR<0.1){
      h_MuMatchedObj->Fill(abs((*GenParticles_PdgId)[match]),wt);
      if(abs((*GenParticles_PdgId)[match])==13){
      h_PdgIdMuParent->Fill(abs((*GenParticles_ParentId)[match]),wt);
      }
      }
      else{
      h_MuMatchedObj->Fill(0.0,wt);
      h_PdgIdMuParent->Fill(0.0,wt);
      }
      }
  */  
}



void LostMuon::print(Long64_t jentry){
  cout<<"*********************************************************************************"<<endl;
  cout<<"Photons:"<<endl;
  for(int i=0;i<Photons->size();i++){
    double dR=0;//DeltaR( bestPhoton.Eta(),bestPhoton.Phi(),(*Photons)[i].Eta(),(*Photons)[i].Phi() );
    //    cout<<jentry<<" i:"<<i<<" phoSize:"<<Photons->size()<<" Pt:"<<bestPhoton.Pt()<<" eta:"<<bestPhoton.Eta()<<" phi:"<<bestPhoton.Phi()<<" otherP:"<<(*Photons)[i].Pt()<<" eta:"<<(*Photons)[i].Eta()<<" phi:"<<(*Photons)[i].Phi()<<" dR:"<<dR<<endl;
    cout<<"PhotonsPt:"<<(*Photons)[i].Pt()<<" eta:"<<(*Photons)[i].Eta()<<" phi:"<<(*Photons)[i].Phi()<<endl;
  }
  cout<<"bestPhoton Pt: "<<bestPhoton.Pt()<<" eta:"<<bestPhoton.Eta()<<" phi:"<<bestPhoton.Phi()<<" E: "<<bestPhoton.Energy()<<endl;
  cout<<"Muons:"<<endl;
  for(int i=0;i<Muons->size();i++){
    cout<<"MuonPt: "<<(*Muons)[i].Pt()<<" Eta: "<<(*Muons)[i].Eta()<<" Phi: "<<(*Muons)[i].Phi()<<" M: "<<(*Muons)[i].M()<<endl;
  }
  cout<<"Electrons:"<<endl;
  for(int i=0;i<Electrons->size();i++){
    cout<<"ElePt: "<<(*Electrons)[i].Pt()<<" Eta: "<<(*Electrons)[i].Eta()<<" Phi: "<<(*Electrons)[i].Phi()<<" M: "<<(*Electrons)[i].M()<<endl;
  }
  cout<<"Jets:"<<endl; 
  for(int i=0;i<Jets->size();i++){
    cout<<"JetPt:"<<(*Jets)[i].Pt()<<" JetEta:"<<(*Jets)[i].Eta()<<" JetPhi:"<<(*Jets)[i].Phi()<<endl;
  }
  //------------------------- MC only -------------------------------------------------
  /*  for(int i=0;i<GenJets->size();i++){
      cout<<"GenJetPt:"<<(*GenJets)[i].Pt()<<" JetEta:"<<(*GenJets)[i].Eta()<<" JetPhi:"<<(*GenJets)[i].Phi()<<endl;
      }
  
      for(int i=0;i<GenParticles->size();i++){
      // cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" Status:"<<(*GenParticles_Status)[i]<<"\tPx:"<<(*GenParticles)[i].Px()<<" Py:"<<(*GenParticles)[i].Py()<<" Pz:"<<(*GenParticles)[i].Pz()<<" E:"<<(*GenParticles)[i].Energy()<<endl;
      cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<<"\tPx:"<<(*GenParticles)[i].Px()<<" Py:"<<(*GenParticles)[i].Py()<<" Pz:"<<(*GenParticles)[i].Pz()<<"\tPt:"<<(*GenParticles)[i].Pt()<<" Eta:"<<(*GenParticles)[i].Eta()<<" Phi:"<<(*GenParticles)[i].Phi()<<" E:"<<(*GenParticles)[i].Energy()<<endl;
      }*/
  //-------------------------------------------------------------------------
  cout<<"^^^^^^^^^^^^^^^^^^ Event ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl<<endl;
}
