#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvWeights.h"
#include "UserCode/ICHiggsTauTau/interface/TriggerPath.hh"
#include "UserCode/ICHiggsTauTau/interface/TriggerObject.hh"
#include "UserCode/ICHiggsTauTau/interface/EventInfo.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/SimpleParamParser.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "iostream"
#include "fstream"


namespace ic {//namespace

  HinvWeights::HinvWeights(std::string const& name) : ModuleBase(name),
    mc_(mc::summer12_53X),
    era_(era::data_2012_moriond) {
    save_weights_ = true;
    do_trg_weights_     = false;
    do_3dtrg_weights_   = false;
    trg_applied_in_mc_  = false;
    do_idiso_tight_weights_   = false;
    do_idiso_veto_weights_   = false;
    do_w_soup_          = false;
    do_w_reweighting_   = false;
    do_dy_soup_         = false;
    do_dy_reweighting_   = false;
    do_idiso_err_       = false;
    do_idiso_errupordown_ = true;
    do_idiso_errmuore_ = true;
    do_lumixs_weights_ = false;
    save_lumixs_weights_ = true;
    input_params_ = "filelists/Apr04/ParamsApr04.dat";
    sample_name_= "test";
    input_met_ = "metNoMuons";
    hist_trigSF_METL1 = 0;
    hist_trigSF_METHLT = 0;
    hist_trigSF_MjjHLT = 0;
    hist_trigSF_JetHLT = 0;
    hist_trigSF_h3DHLT_Dijet40_A = 0;
    hist_trigSF_h3DHLT_Dijet35_BC = 0;
    hist_trigSF_h3DHLT_Dijet30_D = 0;
    trg_weight_file_="data/scale_factors/MergedHistos.root";
  }

  HinvWeights::~HinvWeights() {
    ;
  }

  int HinvWeights::PreAnalysis() {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "PreAnalysis Info for HinvWeights" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Era: " << Era2String(era_) << std::endl;
    std::cout << "MC: " << MC2String(mc_) << std::endl;
    std::cout << "Save weights: " << save_weights_ << std::endl;
    std::cout << "Do Trg Weights?: \t\t" << do_trg_weights_ << std::endl;
    if(trg_applied_in_mc_) std::cout << "Trg Sel Is Applied In MC make sure you use DataMC scale factors for trigger weightst" << std::endl;
    else std::cout << "Trg Sel Not Applied In MC make sure you use raw Data efficiencies for trigger weightst" << std::endl;
    std::cout << "Trg Sel Applied?: \t\t" << trg_applied_in_mc_ << std::endl;
    std::cout << "Do ID & iso weights for Tight leptons ?: \t\t" << do_idiso_tight_weights_ << std::endl;
    std::cout << "Do ID & iso weights for veto leptons ?: \t\t" << do_idiso_veto_weights_ << std::endl;
    std::cout << "Do ID & iso weight errors ?: \t\t" << do_idiso_err_ <<std::endl;
    std::cout << "Do ID & iso weight error for muore ?: \t\t" << do_idiso_errmuore_ <<std::endl;
 std::cout << "Do ID & iso weight errors up or down?: \t\t" << do_idiso_errupordown_ <<std::endl;
    
    std::cout << "Input MET for MET HLT:  \t\t" << input_met_ << std::endl;
    std::cout << "Note: Input MET for MET L1 is always metNoMuons." << std::endl;

    //Making output histograms
    TFileDirectory const& dir = fs_->mkdir("leptoneffweights");
    tighteleweight = dir.make<TH1F>("tighteleweight","tighteleweight",2000,0.,2.);
    tightmuweight  = dir.make<TH1F>("tightmuweight","tightmuweight",2000,0.,2.);
    vetoeleweight = dir.make<TH1F>("vetoeleweight","vetoeleweight",2000,0.,2.);
    vetomuweight  = dir.make<TH1F>("vetomuweight","vetomuweight",2000,0.,2.);

    if(do_lumixs_weights_){
      //Get Sample name file
      std::string suffix=".root";
      std::size_t found = sample_name_.find(suffix);
      if(found==std::string::npos){
	lumixsweight=1;
	std::cout<<"Non-standard sample name format not doing lumixs weight"<<std::endl;
      }
      else{
	sample_name_.erase(found,5);
	std::cout << "Sample Name: "<<sample_name_<<std::endl;

	//Get lumi xs and events from params file
	SimpleParamParser parser;
	std::cout << "** Parsing parameter file... **" << input_params_ << std::endl;
	parser.ParseFile(input_params_);
	std::cout<<"parsed"<<std::endl;
	double xs=parser.GetParam<double>("XS_"+sample_name_);
	std::cout<<"got xs"<<std::endl;
	double events=parser.GetParam<double>("EVT_"+sample_name_);
	std::cout<<"got events"<<std::endl;
	double lumi=parser.GetParam<double>("LUMI_DATA");
	std::cout<<"got lumi"<<std::endl;

	std::cout<<"XS is: "<<xs<<"pb"<<std::endl;
	std::cout<<"EVT is: "<<events<<std::endl;
	std::cout<<"LUMI is: "<<lumi<<"pb^-1"<<std::endl;

	lumixsweight=xs*lumi/events;
	std::cout<<"LUMIXSWEIGHT is: "<<lumixsweight<<std::endl;
      }
    }

    if (do_w_soup_) {
      std::cout << "Making W Soup:" << std::endl;
      std::cout << "nInc = " << n_inc_ << std::endl;
      w1_ = (n_inc_*f1_) / ( (n_inc_*f1_) + n1_ );
      w2_ = (n_inc_*f2_) / ( (n_inc_*f2_) + n2_ );
      w3_ = (n_inc_*f3_) / ( (n_inc_*f3_) + n3_ );
      w4_ = (n_inc_*f4_) / ( (n_inc_*f4_) + n4_ );
      std::cout << "f1 = " << f1_ << "\t" << "n1 = " << n1_ << "\t" << "w1 = " << w1_ << std::endl;
      std::cout << "f2 = " << f2_ << "\t" << "n2 = " << n2_ << "\t" << "w2 = " << w2_ << std::endl;
      std::cout << "f3 = " << f3_ << "\t" << "n3 = " << n3_ << "\t" << "w3 = " << w3_ << std::endl;
      std::cout << "f4 = " << f4_ << "\t" << "n4 = " << n4_ << "\t" << "w4 = " << w4_ << std::endl;
    }
    if (do_dy_soup_) {
      std::cout << "Making DY Soup:" << std::endl;
      std::cout << "nInc = " << zn_inc_ << std::endl;
      zw1_ = (zn_inc_*zf1_) / ( (zn_inc_*zf1_) + zn1_ );
      zw2_ = (zn_inc_*zf2_) / ( (zn_inc_*zf2_) + zn2_ );
      zw3_ = (zn_inc_*zf3_) / ( (zn_inc_*zf3_) + zn3_ );
      zw4_ = (zn_inc_*zf4_) / ( (zn_inc_*zf4_) + zn4_ );
      std::cout << "f1 = " << zf1_ << "\t" << "n1 = " << zn1_ << "\t" << "w1 = " << zw1_ << std::endl;
      std::cout << "f2 = " << zf2_ << "\t" << "n2 = " << zn2_ << "\t" << "w2 = " << zw2_ << std::endl;
      std::cout << "f3 = " << zf3_ << "\t" << "n3 = " << zn3_ << "\t" << "w3 = " << zw3_ << std::endl;
      std::cout << "f4 = " << zf4_ << "\t" << "n4 = " << zn4_ << "\t" << "w4 = " << zw4_ << std::endl;
    }

    if (do_w_reweighting_) {
      std::cout << "Applying reweighting of W events to NLO MCFM." << std::endl;
    }
    if (do_dy_reweighting_) {
      std::cout << "Applying reweighting of DY events to NLO MCFM." << std::endl;
    }



    if (save_weights_){
    // do weights even if not applied, to fill histo with weight for comparison !
    

    //get trigger scale factor histograms from file
      triggerSF_ = new TFile(trg_weight_file_.c_str());
    hist_trigSF_METL1 = (TH1F*)gDirectory->Get("METL1");
    hist_trigSF_METHLT = (TH1F*)gDirectory->Get("METHLT");
    hist_trigSF_MjjHLT = (TH1F*)gDirectory->Get("MjjHLT");
    hist_trigSF_JetHLT = (TH1F*)gDirectory->Get("JetHLT");
    hist_trigSF_h3DHLT_Dijet40_A = (TH3F*)gDirectory->Get("h3DHLT_Dijet40_A");
    hist_trigSF_h3DHLT_Dijet35_BC = (TH3F*)gDirectory->Get("h3DHLT_Dijet35_BC");
    hist_trigSF_h3DHLT_Dijet30_D = (TH3F*)gDirectory->Get("h3DHLT_Dijet30_D");    
  

//    std::cout << " -- Content of histogram DijetA : " << std::endl;
//     for (int i(0); i< hist_trigSF_h3DHLT_Dijet40_A->GetNbinsX()+2; ++i){
//       std::cout << " -- bin " << i << " ["
// 		<< hist_trigSF_h3DHLT_Dijet40_A->GetXaxis()->GetBinLowEdge(i) << "-"
//                 << hist_trigSF_h3DHLT_Dijet40_A->GetXaxis()->GetBinUpEdge(i) << "] : "
//                 <<  hist_trigSF_h3DHLT_Dijet40_A->GetBinContent(i)
//                 << std::endl;
//     }

//     std::cout << " -- Content of histogram DijetBC : " << std::endl;
//     for (int i(0); i< hist_trigSF_h3DHLT_Dijet35_BC->GetNbinsX()+2; ++i){
//       std::cout << " -- bin " << i << " ["
//                 << hist_trigSF_h3DHLT_Dijet35_BC->GetXaxis()->GetBinLowEdge(i) << "-"
//                 << hist_trigSF_h3DHLT_Dijet35_BC->GetXaxis()->GetBinUpEdge(i) << "] : "
//                 <<  hist_trigSF_h3DHLT_Dijet35_BC->GetBinContent(i)
//                 << std::endl;
//     }

//     std::cout << " -- Content of histogram DijetD : " << std::endl;
//     for (int i(0); i< hist_trigSF_h3DHLT_Dijet30_D ->GetNbinsX()+2; ++i){
//       std::cout << " -- bin " << i << " ["
//                 << hist_trigSF_h3DHLT_Dijet30_D->GetXaxis()->GetBinLowEdge(i) << "-"
//                 << hist_trigSF_h3DHLT_Dijet30_D->GetXaxis()->GetBinUpEdge(i) << "] : "
//                 <<  hist_trigSF_h3DHLT_Dijet30_D->GetBinContent(i)
//                 << std::endl;
//     }

  std::cout << " -- Content of histogram METL1 : " << std::endl;
    for (int i(0); i< hist_trigSF_METL1->GetNbinsX()+2; ++i){
      std::cout << " -- bin " << i << " [" 
		<< hist_trigSF_METL1->GetXaxis()->GetBinLowEdge(i) << "-" 
		<< hist_trigSF_METL1->GetXaxis()->GetBinUpEdge(i) << "] : "
		<<  hist_trigSF_METL1->GetBinContent(i)
		<< std::endl;
    }
    
    std::cout << " -- Content of histogram METHLT : " << std::endl;
    for (int i(0); i< hist_trigSF_METHLT->GetNbinsX()+2; ++i){
      std::cout << " -- bin " << i << " [" 
		<< hist_trigSF_METHLT->GetXaxis()->GetBinLowEdge(i) << "-" 
		<< hist_trigSF_METHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
		<<  hist_trigSF_METHLT->GetBinContent(i);
      
      //change first bins to "1": no data/MC scale factor applied...
      if (i>0 && hist_trigSF_METHLT->GetBinContent(i) == 0) {
	hist_trigSF_METHLT->SetBinContent(i,1);
	std::cout << " -> changed to " <<  hist_trigSF_METHLT->GetBinContent(i);
      }
      std::cout	<< std::endl;
      
    }
    
    std::cout << " -- Content of histogram MjjHLT : " << std::endl;
    for (int i(0); i< hist_trigSF_MjjHLT->GetNbinsX()+2; ++i){
      std::cout << " -- bin " << i << " [" 
		<< hist_trigSF_MjjHLT->GetXaxis()->GetBinLowEdge(i) << "-" 
		<< hist_trigSF_MjjHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
		<<  hist_trigSF_MjjHLT->GetBinContent(i);
      std::cout	<< std::endl;
    }
    
    std::cout << " -- Content of histogram JetHLT : " << std::endl;
    for (int i(0); i< hist_trigSF_JetHLT->GetNbinsX()+2; ++i){
      std::cout << " -- bin " << i << " [" 
		<< hist_trigSF_JetHLT->GetXaxis()->GetBinLowEdge(i) << "-" 
		<< hist_trigSF_JetHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
		<<  hist_trigSF_JetHLT->GetBinContent(i);
      
      //change first bins to "1": no data/MC scale factor applied...
      if (i>0 && hist_trigSF_JetHLT->GetBinContent(i) == 0) {
	hist_trigSF_JetHLT->SetBinContent(i,1);
	std::cout << " -> changed to " <<  hist_trigSF_JetHLT->GetBinContent(i);
      }
      std::cout	<< std::endl;
    }

    if(!do_idiso_err_){
      fillVector("data/scale_factors/ele_tight_id_with_syst.txt",eTight_idisoSF_);
      fillVector("data/scale_factors/ele_veto_id_data_eff_with_syst.txt",eVeto_idisoDataEff_);
      fillVector("data/scale_factors/ele_veto_id_mc_eff.txt",eVeto_idisoMCEff_);
      
      fillVector("data/scale_factors/mu_tight_id_SF.txt",muTight_idSF_);
      fillVector("data/scale_factors/mu_tight_iso_SF.txt",muTight_isoSF_);
      fillVector("data/scale_factors/mu_loose_id_data_eff.txt",muVeto_idDataEff_);
      fillVector("data/scale_factors/mu_loose_iso_data_eff.txt",muVeto_isoDataEff_);
      fillVector("data/scale_factors/mu_loose_id_mc_eff.txt",muVeto_idMCEff_);
      fillVector("data/scale_factors/mu_loose_iso_mc_eff.txt",muVeto_isoMCEff_);
    }
    else if(do_idiso_errmuore_){
      fillVector("data/scale_factors/ele_tight_id_with_syst.txt",eTight_idisoSF_);
      fillVector("data/scale_factors/ele_veto_id_data_eff_with_syst.txt",eVeto_idisoDataEff_);
      fillVector("data/scale_factors/ele_veto_id_mc_eff.txt",eVeto_idisoMCEff_);
      
      fillVectorError("data/scale_factors/mu_tight_id_SF.txt",muTight_idSF_,do_idiso_errupordown_);
      fillVectorError("data/scale_factors/mu_tight_iso_SF.txt",muTight_isoSF_,do_idiso_errupordown_);
      fillVectorError("data/scale_factors/mu_loose_id_data_eff.txt",muVeto_idDataEff_,do_idiso_errupordown_);
      fillVectorError("data/scale_factors/mu_loose_iso_data_eff.txt",muVeto_isoDataEff_,do_idiso_errupordown_);
      fillVectorError("data/scale_factors/mu_loose_id_mc_eff.txt",muVeto_idMCEff_,do_idiso_errupordown_);
      fillVectorError("data/scale_factors/mu_loose_iso_mc_eff.txt",muVeto_isoMCEff_,do_idiso_errupordown_);    
    }
    else{
      fillVectorError("data/scale_factors/ele_tight_id_with_syst.txt",eTight_idisoSF_,do_idiso_errupordown_);
      fillVectorError("data/scale_factors/ele_veto_id_data_eff_with_syst.txt",eVeto_idisoDataEff_,do_idiso_errupordown_);
      fillVectorError("data/scale_factors/ele_veto_id_mc_eff.txt",eVeto_idisoMCEff_,do_idiso_errupordown_);
      
      fillVector("data/scale_factors/mu_tight_id_SF.txt",muTight_idSF_);
      fillVector("data/scale_factors/mu_tight_iso_SF.txt",muTight_isoSF_);
      fillVector("data/scale_factors/mu_loose_id_data_eff.txt",muVeto_idDataEff_);
      fillVector("data/scale_factors/mu_loose_iso_data_eff.txt",muVeto_isoDataEff_);
      fillVector("data/scale_factors/mu_loose_id_mc_eff.txt",muVeto_idMCEff_);
      fillVector("data/scale_factors/mu_loose_iso_mc_eff.txt",muVeto_isoMCEff_);    
    }

    
    for (unsigned iBin(0); iBin<muTight_idSF_.size();++iBin){
      muTight_idisoSF_.push_back(muTight_idSF_[iBin]*muTight_isoSF_[iBin]);
      if(muVeto_idDataEff_[iBin]>1)muVeto_idDataEff_[iBin]=1;
      if(muVeto_isoDataEff_[iBin]>1)muVeto_isoDataEff_[iBin]=1;
      if(muVeto_idMCEff_[iBin]>1)muVeto_idMCEff_[iBin]=1;
      if(muVeto_isoMCEff_[iBin]>1)muVeto_isoMCEff_[iBin]=1;
      if(muVeto_idDataEff_[iBin]<0)muVeto_idDataEff_[iBin]=0;
      if(muVeto_isoDataEff_[iBin]<0)muVeto_isoDataEff_[iBin]=0;
      if(muVeto_idMCEff_[iBin]<0)muVeto_idMCEff_[iBin]=0;
      if(muVeto_isoMCEff_[iBin]<0)muVeto_isoMCEff_[iBin]=0;
      muVeto_idisoDataEff_.push_back(muVeto_idDataEff_[iBin]*muVeto_isoDataEff_[iBin]);
      muVeto_idisoMCEff_.push_back(muVeto_idMCEff_[iBin]*muVeto_isoMCEff_[iBin]);
    }

    for (unsigned iBin(0); iBin<eVeto_idisoDataEff_.size();++iBin){
      if(eVeto_idisoDataEff_[iBin]>1)eVeto_idisoDataEff_[iBin]=1;
      if(eVeto_idisoMCEff_[iBin]>1)eVeto_idisoMCEff_[iBin]=1;
      if(eVeto_idisoDataEff_[iBin]<0)eVeto_idisoDataEff_[iBin]=0;
      if(eVeto_idisoMCEff_[iBin]<0)eVeto_idisoMCEff_[iBin]=0;
    }
    

    eventsWithGenElectron_ = 0;
    eventsWithGenElectronFromTau_ = 0;
    eventsWithGenMuon_ = 0;
    eventsWithGenMuonFromTau_ = 0;
    eventsWithGenElectronInAcc_ = 0;
    eventsWithGenElectronFromTauInAcc_ = 0;
    eventsWithGenMuonInAcc_ = 0;
    eventsWithGenMuonFromTauInAcc_ = 0;
    eventsWithGenTau_ = 0;
    }

    return 0;
  }
  //    std::cout<<"!!!!!!!!!!!!!!!!! IMPORTANT BIT !!!!!!!!!!!!!!!!!!" << std::endl;    
  int HinvWeights::Execute(TreeEvent *event) {

    EventInfo * eventInfo = event->GetPtr<EventInfo>("eventInfo");

    if(do_lumixs_weights_){
      //std::cout<<"weight before lumixs: "<<eventInfo->total_weight()<<std::endl;
      //std::cout<<"intended lumixsweight: "<<lumixsweight<<std::endl;
      eventInfo->set_weight("lumixs",lumixsweight);
      //std::cout<<"weight after lumixs: "<<eventInfo->total_weight()<<std::endl;
      //std::cout<<"lumixsweight info "<<eventInfo->weight_is_enabled("lumixs") << " " << eventInfo->weight_defined("lumixs") << " " << eventInfo->weight("lumixs") << std::endl;
    }
    //else eventInfo->set_weight("!lumixs",lumixsweight);
    
    if (save_weights_){

          //double weight = 1.0;

 //    if (do_btag_weight_) {
//       std::vector<PFJet*> jets = event->GetPtrVec<PFJet>("pfJetsPFlow"); // Make a copy of the jet collection
//       ic::erase_if(jets,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));
//       //double no_btag_weight = btag_weight.GetWeight(jets, "CSVM", 0, 0, is_2012_);
//       //double inclusive_btag_weight = btag_weight.GetWeight(jets, "CSVM", 1, 99, is_2012_);
//       double no_btag_weight = 1.0;
//       double inclusive_btag_weight = 1.0;
//       btag_weight.ReTag(jets, mc_ == mc::summer12_53X);
//       event->Add("no_btag_weight", no_btag_weight);
//       event->Add("inclusive_btag_weight", inclusive_btag_weight);
//     }

    //get METnoMuons:
    Met const* metHLT = event->GetPtr<Met>(input_met_);
    Met const* metL1 = event->GetPtr<Met>("metNoMuons");

    int lBinmetHLTA=0;
    int lBinjet2ptA=0;
    int lBinMjjA=0;
    int lBinmetHLTBC=0;
    int lBinjet2ptBC=0;
    int lBinMjjBC=0;
    int lBinmetHLTD=0;
    int lBinjet2ptD=0;
    int lBinMjjD=0;

    double l1met = metL1->pt();
    double lMax = hist_trigSF_METL1->GetXaxis()->GetBinCenter(hist_trigSF_METL1->GetNbinsX());
    double lMin = hist_trigSF_METL1->GetXaxis()->GetBinCenter(1);
    if (l1met > lMax)  l1met = lMax;
    if (l1met < lMin)  l1met = lMin;
    int lBin = hist_trigSF_METL1->GetXaxis()->FindFixBin(l1met);
    double metl1weight = hist_trigSF_METL1->GetBinContent(lBin);
    if (do_trg_weights_ && !do_3dtrg_weights_) eventInfo->set_weight("trig_metL1",metl1weight);
    else eventInfo->set_weight("!trig_metL1",metl1weight);
    //std::cout << " -- MET L1 " << l1met << " " << metl1weight << std::endl;

    double hltmet = metHLT->pt();
    if (do_trg_weights_ && !do_3dtrg_weights_){
    lMax = hist_trigSF_METHLT->GetXaxis()->GetBinCenter(hist_trigSF_METHLT->GetNbinsX());
    lMin = hist_trigSF_METHLT->GetXaxis()->GetBinCenter(1);
    if (hltmet > lMax)  hltmet = lMax;
    if (hltmet < lMin)  hltmet = lMin;
    }
    lBin = hist_trigSF_METHLT->GetXaxis()->FindFixBin(hltmet);
    double methltweight = hist_trigSF_METHLT->GetBinContent(lBin);
    if (do_trg_weights_ && !do_3dtrg_weights_) eventInfo->set_weight("trig_metHLT",methltweight);
    else eventInfo->set_weight("!trig_metHLT",methltweight);
    //std::cout << " -- MET HLT " << hltmet << " " << methltweight << std::endl;

    double mjjhltweight = 1.0;
    double jet1hltweight = 1.0;
    double jet2hltweight = 1.0;

    double mjj=0.;
    double jet1pt=0.;
    double jet2pt=0.;

    //get 2 leading jets
    std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>("jjLeadingCandidates");
    if (dijet_vec.size() > 0) {
      
      CompositeCandidate const* dijet = dijet_vec.at(0);
      
      Candidate const* jet1 = dijet->GetCandidate("jet1");
      Candidate const* jet2 = dijet->GetCandidate("jet2");
      
      mjj = dijet->M();
      if (do_trg_weights_ && !do_3dtrg_weights_){
      lMax = hist_trigSF_MjjHLT->GetXaxis()->GetBinCenter(hist_trigSF_MjjHLT->GetNbinsX());
      lMin = hist_trigSF_MjjHLT->GetXaxis()->GetBinCenter(1);
      if (mjj > lMax)  mjj = lMax;
      if (mjj < lMin)  mjj = lMin;
      }
      lBin = hist_trigSF_MjjHLT->GetXaxis()->FindFixBin(mjj);
      mjjhltweight = hist_trigSF_MjjHLT->GetBinContent(lBin);
      if (do_trg_weights_ && !do_3dtrg_weights_) eventInfo->set_weight("trig_mjjHLT",mjjhltweight);
      else eventInfo->set_weight("!trig_mjjHLT",mjjhltweight);
      //      std::cout << " -- Mjj HLT " << mjj << " " << mjjhltweight << std::endl;
      
      lMax = hist_trigSF_JetHLT->GetXaxis()->GetBinCenter(hist_trigSF_JetHLT->GetNbinsX());
      lMin = hist_trigSF_JetHLT->GetXaxis()->GetBinCenter(1);
      jet1pt = jet1->pt();
      if (jet1pt < lMin)  jet1pt = lMin;
      if (jet1pt > lMax)  jet1pt = lMax;
      lBin = hist_trigSF_JetHLT->GetXaxis()->FindFixBin(jet1pt);
      jet1hltweight = hist_trigSF_JetHLT->GetBinContent(lBin);
      if (do_trg_weights_ && !do_3dtrg_weights_) eventInfo->set_weight("trig_jet1HLT",jet1hltweight);
      else eventInfo->set_weight("!trig_jet1HLT",jet1hltweight);
      //std::cout << " -- Jet1 HLT " << jet1pt << " " << jet1hltweight << std::endl;
      
      jet2pt = jet2->pt();
      if (do_trg_weights_ && !do_3dtrg_weights_){
      if (jet2pt > lMax)  jet2pt = lMax;
      if (jet2pt < lMin)  jet2pt = lMin;
      }
      lBin = hist_trigSF_JetHLT->GetXaxis()->FindFixBin(jet2pt);
      jet2hltweight = hist_trigSF_JetHLT->GetBinContent(lBin);
      if (do_trg_weights_ && !do_3dtrg_weights_) eventInfo->set_weight("trig_jet2HLT",jet2hltweight);
      else eventInfo->set_weight("!trig_jet2HLT",jet2hltweight);
      //std::cout << " -- Jet2pt " << jet2pt << " " << jet2hltweight << std::endl;
      lBinjet2ptA = hist_trigSF_h3DHLT_Dijet40_A->GetXaxis()->FindFixBin(jet2pt);
      //      std::cout << " Weight Ax " << lBinjet2ptA << std::endl;
      lBinmetHLTA = hist_trigSF_h3DHLT_Dijet40_A->GetYaxis()->FindFixBin(hltmet);
      //      std::cout << " Weight Ay " << lBinmetHLTA << std::endl;
      lBinMjjA = hist_trigSF_h3DHLT_Dijet40_A->GetZaxis()->FindFixBin(mjj);
      //      std::cout << " Weight Az " << lBinMjjA << std::endl;
      lBinjet2ptBC = hist_trigSF_h3DHLT_Dijet35_BC->GetXaxis()->FindFixBin(jet2pt);
      lBinmetHLTBC = hist_trigSF_h3DHLT_Dijet35_BC->GetYaxis()->FindFixBin(hltmet);     
      lBinMjjBC = hist_trigSF_h3DHLT_Dijet35_BC->GetZaxis()->FindFixBin(mjj);
      lBinjet2ptD = hist_trigSF_h3DHLT_Dijet30_D->GetXaxis()->FindFixBin(jet2pt);
      lBinmetHLTD = hist_trigSF_h3DHLT_Dijet30_D->GetYaxis()->FindFixBin(hltmet);
      lBinMjjD = hist_trigSF_h3DHLT_Dijet30_D->GetZaxis()->FindFixBin(mjj);	  
      double binweightA = hist_trigSF_h3DHLT_Dijet40_A->GetBinContent(lBinjet2ptA,lBinmetHLTA,lBinMjjA);
      //      std::cout << " Weight A " << binweightA << std::endl;
      double binweightBC = hist_trigSF_h3DHLT_Dijet35_BC->GetBinContent(lBinjet2ptBC,lBinmetHLTBC,lBinMjjBC);
      double binweightD = hist_trigSF_h3DHLT_Dijet30_D->GetBinContent(lBinjet2ptD,lBinmetHLTD,lBinMjjD);
      double binweight =((0.0449*binweightA)+(0.5853*binweightBC)+(0.3698*binweightD));
      //      std::cout << " Weight Total " << binweight << std::endl;
      
//       std::cout << " Event # " << eventInfo->event() << std::endl;
//       std::cout << " jet2pt " << jet2pt << " hltmet " << hltmet << " mjj " << mjj << std::endl;
//       std::cout << " Bin Jet2pt0 " << lBinjet2ptA << " Bin metHLT0 " << lBinmetHLTA << " BinMjj0 " << lBinMjjA << std::endl;
//       std::cout << " Bin Jet2pt1 " << lBinjet2ptBC << "Bin metHLT1 " << lBinmetHLTBC << " BinMjj1 " << lBinMjjBC << std::endl;
//       std::cout << " Bin Jet2pt2 " << lBinjet2ptD << "Bin metHLT2 " << lBinmetHLTD << " BinMjj2 " << lBinMjjD << std::endl;
//       std::cout << " Weight 0 " << binweightA << " Weight 1 " << binweightBC << " Weight 2 " << binweightD << std::endl;
//       std::cout << " Total Weight " << binweight << std::endl;

      //Calculate trg weights for 3D case
	if (do_trg_weights_ && do_3dtrg_weights_) eventInfo->set_weight("trig_binweight",binweight);	  
	else eventInfo->set_weight("!trig_binweight", binweight);
	  
      //std::cout << " -- Jet2 HLT " << jet2pt << " " << jet2hlt << std::endl;

      //weight *= (ele_trg * tau_trg);
      //event->Add("trigweight_1", ele_trg);
      //event->Add("trigweight_2", tau_trg);
    }

    //eventInfo->set_weight("!trigger",metl1weight*methltweight*mjjhltweight*jet1hltweight*jet2hltweight);

    //eventInfo->set_weight("lepton", weight);

    //ID+iso tight leptons
    std::vector<Electron*> const& elecs = event->GetPtrVec<Electron>("selElectrons");
    double ele_weight = 1.0;
    for (unsigned iEle(0); iEle<elecs.size();++iEle){
      ele_weight *= eTight_idisoSF_[findElectronPtEtaBin(elecs[iEle]->pt(),elecs[iEle]->eta())];
    }
    eventInfo->set_weight("!eleTight_idisoSF",ele_weight);
    tighteleweight->Fill(ele_weight);

   std::vector<Muon*> const& mus = event->GetPtrVec<Muon>("selMuons");
    double mu_weight = 1.0;
    for (unsigned iEle(0); iEle<mus.size();++iEle){
      unsigned lBin = findMuonPtEtaBin(mus[iEle]->pt(),mus[iEle]->eta());
      mu_weight *= muTight_idisoSF_[lBin];
    }
    eventInfo->set_weight("!muTight_idisoSF",mu_weight);
    tightmuweight->Fill(mu_weight);

    if(do_idiso_tight_weights_){
      eventInfo->set_weight("idisoTight",ele_weight*mu_weight);
    }
    else{
      eventInfo->set_weight("!idisoTight",ele_weight*mu_weight);
    }

    //TO DO: id+iso veto leptons
    //first try: take leptons from W in pT,eta acceptance
    std::vector<GenParticle*> const& genParts = event->GetPtrVec<GenParticle>("genParticles");
    double ele_veto_weight = 1.0;
    double mu_veto_weight = 1.0;

    for (unsigned iEle(0); iEle<genParts.size(); ++iEle){//loop on genparticles

      if (genParts[iEle]->status()!=3) continue;
 
      //do Electrons
      if (abs(genParts[iEle]->pdgid())==11){
	eventsWithGenElectron_++;
	if (genParts[iEle]->pt() > 10 && fabs(genParts[iEle]->eta()) < 2.4) {
	  unsigned lBin = findElectronPtEtaBin(genParts[iEle]->pt(),genParts[iEle]->eta());
	  ele_veto_weight *= (1-eVeto_idisoDataEff_[lBin])/(1-eVeto_idisoMCEff_[lBin]);
	  eventsWithGenElectronInAcc_++;
	}
      }

      //doMuons
      if (abs(genParts[iEle]->pdgid())==13){
	eventsWithGenMuon_++;
	if (genParts[iEle]->pt() > 10 && fabs(genParts[iEle]->eta()) < 2.1) {
	  unsigned lBin = findMuonPtEtaBin(genParts[iEle]->pt(),genParts[iEle]->eta());
	  mu_veto_weight *= (1-muVeto_idisoDataEff_[lBin])/(1-muVeto_idisoMCEff_[lBin]);
	  eventsWithGenMuonInAcc_++;
	}
      }

      //doTaus
      if (abs(genParts[iEle]->pdgid())==15){
	eventsWithGenTau_++;
	//get the specific taus collection with daughters filled
	std::vector<GenParticle*> const& taus = event->GetPtrVec<GenParticle>("genParticlesTaus");
	for (unsigned j = 0; j < taus.size(); ++j) {
	  unsigned idDau = abs(taus[j]->pdgid());
	  if (idDau==11) {
	    eventsWithGenElectronFromTau_++;
	    if (taus[j]->pt() > 10 && fabs(taus[j]->eta()) < 2.4){
	      unsigned lBin = findElectronPtEtaBin(taus[j]->pt(),taus[j]->eta());
	      ele_veto_weight *= (1-eVeto_idisoDataEff_[lBin])/(1-eVeto_idisoMCEff_[lBin]);
	      eventsWithGenElectronFromTauInAcc_++;
	    }
	  }
	  if (idDau==13) {
	    eventsWithGenMuonFromTau_++;
	    if (taus[j]->pt() > 10 && fabs(taus[j]->eta()) < 2.1){
	      unsigned lBin = findMuonPtEtaBin(taus[j]->pt(),taus[j]->eta());
	      mu_veto_weight *= (1-muVeto_idisoDataEff_[lBin])/(1-muVeto_idisoMCEff_[lBin]);
	      eventsWithGenMuonFromTauInAcc_++;
	    }
	  }
	}
      }

    }//loop on genparticles
    vetoeleweight->Fill(ele_veto_weight);
    vetomuweight->Fill(mu_veto_weight);
    if (do_idiso_veto_weights_) eventInfo->set_weight("idisoVeto",ele_veto_weight*mu_veto_weight);
    else eventInfo->set_weight("!idisoVeto",ele_veto_weight*mu_veto_weight);

    


    }


    bool zeroParton = false;

    if (do_w_soup_) {
      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      bool count_jets = false;
      unsigned partons = 0;
      for (unsigned i = 0; i < parts.size(); ++i) {
        if (parts[i]->status() != 3) continue;
        unsigned id = abs(parts[i]->pdgid());
        if (count_jets) { 
          if (id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6 || id == 21) partons++;
        }
        if (id == 24) count_jets = true; 
      }
      if (partons > 4) {
        std::cerr << "Error making soup, event has " << partons << " partons!" << std::endl;
        throw;
      }
      if (partons == 1) eventInfo->set_weight("wsoup", w1_);
      if (partons == 2) eventInfo->set_weight("wsoup", w2_);
      if (partons == 3) eventInfo->set_weight("wsoup", w3_);
      if (partons == 4) eventInfo->set_weight("wsoup", w4_);

      if (partons == 0) zeroParton = true;

    }

    if (do_w_reweighting_ || do_dy_reweighting_) {
      double vReweight = 1.0;

      std::vector<GenJet *> genjets = event->GetPtrVec<GenJet>("genJets");
      if (genjets.size() > 1) {
	double lMjj = (genjets[0]->vector()+genjets[1]->vector()).M();

	GenParticle* lBoson = 0;
	std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
	
	for (unsigned i = 0; i < parts.size(); ++i) {
	  if (parts[i]->status() != 3) continue;
	  unsigned id = abs(parts[i]->pdgid());
	  if (id==24 || id==23) lBoson = parts[i];
	}
	if (lBoson){
	  double y_star = fabs(lBoson->vector().Rapidity() - (genjets[0]->vector().Rapidity()+genjets[1]->vector().Rapidity())/2.);
	  vReweight = nloReweighting(lMjj,y_star);
	}
      }
      
      eventInfo->set_weight("vReweighting", vReweight);

    }

    if (do_dy_soup_) {
      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      bool count_jets = false;
      unsigned partons = 0;
      for (unsigned i = 0; i < parts.size(); ++i) {
        // std::cout << i << "\t" << parts[i]->status() << "\t" << parts[i]->pdgid() << "\t" << parts[i]->vector() << std::endl;
        if (parts[i]->status() != 3) continue;
        unsigned id = abs(parts[i]->pdgid());
        if (count_jets) { 
          if (id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6 || id == 21) partons++;
        }
        if (id == 23) count_jets = true; 
      }
      if (partons > 4) {
        std::cerr << "Error making soup, event has " << partons << " partons!" << std::endl;
        throw;
      }
      if (partons == 1) eventInfo->set_weight("dysoup", zw1_);
      if (partons == 2) eventInfo->set_weight("dysoup", zw2_);
      if (partons == 3) eventInfo->set_weight("dysoup", zw3_);
      if (partons == 4) eventInfo->set_weight("dysoup", zw4_);
      if (partons == 0) zeroParton = true;
    }

    if (!save_weights_) event->Add("NoParton",zeroParton);
    //std::cout<<"Final weight: "<<eventInfo->total_weight()<<std::endl;
    return 0;
  }

  int HinvWeights::PostAnalysis() {
    if (save_weights_) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "PostAnalysis Info for HinvWeights" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << " -- eventsWithGenElectron_ = " << eventsWithGenElectron_ << std::endl;
    std::cout << " -- eventsWithGenElectronInAcc_ = " << eventsWithGenElectronInAcc_ << std::endl;
    std::cout << " -- eventsWithGenElectronFromTau_ = " << eventsWithGenElectronFromTau_ << std::endl;
    std::cout << " -- eventsWithGenElectronFromTauInAcc_ = " << eventsWithGenElectronFromTauInAcc_ << std::endl;
    std::cout << " -- eventsWithGenMuon_ = " << eventsWithGenMuon_ << std::endl;
    std::cout << " -- eventsWithGenMuonInAcc_ = " << eventsWithGenMuonInAcc_ << std::endl;
    std::cout << " -- eventsWithGenMuonFromTau_ = " << eventsWithGenMuonFromTau_ << std::endl;
    std::cout << " -- eventsWithGenMuonFromTauInAcc_ = " << eventsWithGenMuonFromTauInAcc_ << std::endl;
    std::cout << " -- eventsWithGenTau_ = " << eventsWithGenTau_ << std::endl;
    }

    return 0;
  }

  void HinvWeights::PrintInfo() {
    ;
  }

  void HinvWeights::SetWTargetFractions(double f0, double f1, double f2, double f3, double f4) {
    f0_ = f0;
    f1_ = f1;
    f2_ = f2;
    f3_ = f3;
    f4_ = f4;

  }
  void HinvWeights::SetWInputYields(double n_inc, double n1, double n2, double n3, double n4) {
    n_inc_ = n_inc;
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
  }

  void HinvWeights::SetDYTargetFractions(double zf0, double zf1, double zf2, double zf3, double zf4) {
    zf0_ = zf0;
    zf1_ = zf1;
    zf2_ = zf2;
    zf3_ = zf3;
    zf4_ = zf4;

  }
  void HinvWeights::SetDYInputYields(double zn_inc, double zn1, double zn2, double zn3, double zn4) {
    zn_inc_ = zn_inc;
    zn1_ = zn1;
    zn2_ = zn2;
    zn3_ = zn3;
    zn4_ = zn4;
  }

  double HinvWeights::Efficiency(double m, double m0, double sigma, double alpha,
    double n, double norm)
  {
    const double sqrtPiOver2 = 1.2533141373;
    const double sqrt2 = 1.4142135624;
    double sig = fabs((double) sigma);
    double t = (m - m0)/sig;
    if(alpha < 0)
      t = -t;
    double absAlpha = fabs(alpha/sig);
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = absAlpha - n/absAlpha;
    double ApproxErf;
    double arg = absAlpha / sqrt2;
    if (arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    double leftArea = (1 + ApproxErf) * sqrtPiOver2;
    double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
    double area = leftArea + rightArea;
    if( t <= absAlpha ){
      arg = t / sqrt2;
      if(arg > 5.) ApproxErf = 1;
      else if (arg < -5.) ApproxErf = -1;
      else ApproxErf = TMath::Erf(arg);
      return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
    }
    else{
      return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) -
        1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
    }
  }

  unsigned HinvWeights::findElectronPtEtaBin(double pt, double eta){

    unsigned etaBin = 0;
    if (fabs(eta) < 0.8) etaBin=0;
    else if (fabs(eta) < 1.442) etaBin=1;
    else if (fabs(eta) < 1.556) etaBin=2;
    else if (fabs(eta) < 2.0) etaBin=3;
    else etaBin=4;
    unsigned ptBin = 0;
    if (pt<15) ptBin=0;
    else if (pt<20) ptBin=1;
    else if (pt<30) ptBin=2;
    else if (pt<40) ptBin=3;
    else if (pt<50) ptBin=4;
    else ptBin=5;
    
    return 6*etaBin+ptBin;

  }

  unsigned HinvWeights::findMuonPtEtaBin(double pt, double eta){
    unsigned etaBin = 0;
    if (fabs(eta) < 0.9) etaBin=0;
    else if (fabs(eta) < 1.2) etaBin=1;
    else etaBin=2;

    unsigned ptBin = 0;
    if (pt<20) ptBin=0;
    else if (pt<25) ptBin=1;
    else if (pt<30) ptBin=2;
    else if (pt<35) ptBin=3;
    else if (pt<40) ptBin=4;
    else if (pt<50) ptBin=5;
    else if (pt<60) ptBin=6;
    else if (pt<90) ptBin=7;
    else if (pt<140) ptBin=8;
    else if (pt<300) ptBin=9;
    else ptBin=10;

    return 11*etaBin+ptBin;
  }

  void HinvWeights::fillVector(const std::string & aFileName, std::vector<double> & aVector){
    std::ifstream lInput;
    lInput.open(aFileName);
    if(!lInput.is_open()) {
      std::cerr << "Unable to open file: " << aFileName << ". Setting vector content to 1." << std::endl;
      //max expected size for e and mu is 33...
      aVector.resize(33,1);
      return;
    }
    while(1){
      double pTmin = 0;
      double pTmax = 0;
      double etaMin = 0;
      double etaMax = 0;
      double SF = 0;
      double SFerrPlus = 0;
      double SFerrMinus = 0;
      lInput>>pTmin>>pTmax>>etaMin>>etaMax>>SF>>SFerrMinus>>SFerrPlus;
      //protect against blank line at the end of the file
      if (pTmin > 1) aVector.push_back(SF);
      if(lInput.eof()){
	break; 
      }
    }

    std::cout << " ---- Size of vector for file " << aFileName << " = " << aVector.size() << std::endl;
    lInput.close();

  }


  void HinvWeights::fillVectorError(const std::string & aFileName, std::vector<double> & aVector, bool upordown){

    std::ifstream lInput;
    lInput.open(aFileName);
    if(!lInput.is_open()){
      std::cerr << "Unable to open file: " << aFileName << ". Setting vector content to 1." << std::endl;
      //max expected size for e and mu is 33...
      aVector.resize(33,1);
      return;
    }
    while(1){
      double pTmin = 0;
      double pTmax = 0;
      double etaMin = 0;
      double etaMax = 0;
      double SF = 0;
      double SFerrPlus = 0;
      double SFerrMinus = 0;
      lInput>>pTmin>>pTmax>>etaMin>>etaMax>>SF>>SFerrMinus>>SFerrPlus;
      //protect against blank line at the end of the file
      if(upordown){
	if (pTmin > 1) aVector.push_back(SF+SFerrPlus);
      }
      else{
	if (pTmin > 1) aVector.push_back(SF-SFerrMinus);
      }

      if(lInput.eof()){
	break; 
      }
    }

    std::cout << " ---- Size of vector for file " << aFileName << " = " << aVector.size() << std::endl;
    lInput.close();

  }

  double HinvWeights::nloReweighting(const double & aMjj, const double & aYstar){

    double weight = 1.0;

    //double y_par0 = 8.49667e-01;
    //double y_par1 = 1.49687e-01;
    TF1 *f0c = new TF1("f0c","[0]+[1]*x "); // |Ystar|
    f0c->SetParameters(8.49667e-01, 1.49687e-01 );  //NLOmcfm/MadgraphGenJ 80M
    
    TF1 *f1l = new TF1("f1l","[0]+[1]*log(x)+[2]*x "); // mjj (in GeV)
    f1l->SetParameters( 3.92568e-01, 1.20734e-01, -2.55622e-04  ); //NLOmcfm/MadgraphGenJ 80M 
    
    weight = f0c->Eval(aYstar)*f1l->Eval(aMjj);
    
    return weight;
  };


}//namespace
