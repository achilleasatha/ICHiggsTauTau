#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataZNormShape.h"
#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"
#include <map>

namespace ic{

  DataZNormShape::DataZNormShape(std::string name) : LTModule(name){
    sigmcweight_="total_weight_lepveto";
    contmcweight_="total_weight_leptight";
    contdataweight_="weight_nolep";
    contdataextrasel_="";
    sigcontextrafactor_=1;
    std::vector<std::string> shapes;
    shapes.push_back("jet2_pt(200,0.,1000.)");
    shape_=shapes;
    dirname_="";
  };

  DataZNormShape::~DataZNormShape(){ ;};

  int DataZNormShape::Init(TFile* fs){
    fs_=fs;
    std::cout<<"Initialisation info for "<<module_name_<<":"<<std::endl;
    std::cout<<"Signal MC ewk set is: "<<sigmcewkset_<<std::endl;
    std::cout<<"Signal MC qcd set is: "<<sigmcqcdset_<<std::endl;
    std::cout<<"Control MC ewk set  is: "<<contmcewkset_<<std::endl;
    std::cout<<"Control MC qcd set is: "<<contmcqcdset_<<std::endl;
    std::cout<<"Control data set is: "<<contdataset_<<std::endl;
    std::cout<<"Base selection is: "<<basesel_<<std::endl;
    std::cout<<"Signal extra selection is: "<<sigcat_<<std::endl;
    std::cout<<"Control extra selection is: "<<contcat_<<std::endl;
    
    std::cout<<"Inclusive control ewk cross-section: "<<sigmainccontewk_<<std::endl;
    std::cout<<"Inclusive control qcd cross-section: "<<sigmainccontqcd_<<std::endl;
    std::cout<<"Inclusive signal ewk cross-section: "<<sigmaincsigewk_<<std::endl;
    std::cout<<"Inclusive signal qcd cross-section: "<<sigmaincsigqcd_<<std::endl;
    std::cout<<"Number of generated ewk Z events: "<<ngenincewk_<<std::endl;
    std::cout<<"Number of generated qcd Z events: "<<ngenincqcd_<<std::endl;
    std::cout<<"Number of generated ewk Z events after mass filtering: "<<ngenmassfilteredewk_<<std::endl;
    std::cout<<"Number of generated qcd Z events after mass filtering: "<<ngenmassfilteredqcd_<<std::endl;
    return 0;
  };

  int DataZNormShape::Run(LTFiles* filemanager){
    std::cout<<module_name_<<":"<<std::endl;
    TFile *file=fs_;
    TDirectory* dir;
    if(dirname_==""){
      dir=file->mkdir("zmumu");
    }
    else if(!fs_->GetDirectory(dirname_.c_str())){
      dir=file->mkdir(dirname_.c_str());
    }
    else{
      dir=fs_->GetDirectory(dirname_.c_str());
    }
    dir->cd();
    //Get Shapes for NSMC, NCMC, NCData and NCBkg
    std::cout<<"  Getting control MC shape"<<std::endl;
    TH1F  contmcewkshape = filemanager->GetSetShape(contmcewkset_,"jet2_pt(200,0.,1000.)",basesel_,contcat_,"total_weight_leptight",false);
    TH1F  contmcqcdshape = filemanager->GetSetShape(contmcqcdset_,"jet2_pt(200,0.,1000.)",basesel_,contcat_,"total_weight_leptight",false);
    std::cout<<"  Getting control MC Backgrounds shape"<<std::endl;
    TH1F  contbkgshape = filemanager->GetSetsShape(contbkgset_,"jet2_pt(200,0.,1000.)",basesel_,contcat_,"total_weight_leptight",false);
    std::cout<<"  Getting control Data shape"<<std::endl;
    TH1F  contdatashape = filemanager->GetSetShape(contdataset_,"jet2_pt(200,0.,1000.)",basesel_,contcat_+contdataextrasel_,"weight_nolep",false);
    
    //Integrate over shape to get number in each region
    double ncmcewk = Integral(&contmcewkshape);
    double ncmcqcd = Integral(&contmcqcdshape);
    double ncdata = Integral(&contdatashape);
    double ncbkg = Integral(&contbkgshape);

    std::cout<<"  ncmcewk: "<<ncmcewk<<", ncmcqcd: "<<ncmcqcd<<", ncdata: "<<ncdata<<", ncbkg: "<<ncbkg<<std::endl;

    //Do Weighting
    double effcvbfewk=ncmcewk/ngenincewk_;
    double effcvbfqcd=ncmcqcd/ngenincqcd_;
    double baseweight=(ncdata-ncbkg)/(sigmainccontewk_*effcvbfewk+sigmainccontqcd_*effcvbfqcd);
    for(unsigned iShape=0;iShape<shape_.size();iShape++){
      std::string histname;
      if(shapename_.size()==0){
	std::vector<std::string> strs;
	boost::split(strs, shape_[iShape], boost::is_any_of("("));
        histname=strs[0];
      }
      else{
        histname=shapename_[iShape];
      }     
      TH1F  sigmcewkshape = filemanager->GetSetShape(sigmcewkset_,shape_[iShape],basesel_,sigcat_,"weight_nolep",false);
      TH1F  sigmcqcdshape = filemanager->GetSetShape(sigmcqcdset_,shape_[iShape],basesel_,sigcat_,"weight_nolep",false);
      dir->cd();
      sigmcewkshape.Scale(sigcontextrafactor_*baseweight*sigmaincsigewk_/ngenmassfilteredewk_);
      sigmcqcdshape.Scale(sigcontextrafactor_*baseweight*sigmaincsigqcd_/ngenmassfilteredqcd_);
     
      //!!MAKE OVERALL SIGMCSHAPE HISTO
      //Get binning info from shape
      std::vector<std::string> strs;
      boost::split(strs, shape_[iShape], boost::is_any_of("(,)"));
      TH1F* sigmcshape=new TH1F(histname.c_str(),histname.c_str(),boost::lexical_cast<int>(strs[1]),boost::lexical_cast<double>(strs[2]),boost::lexical_cast<double>(strs[3]));
      sigmcshape->Add(&sigmcewkshape,&sigmcqcdshape);
      sigmcshape->SetName(histname.c_str());
      sigmcshape->Write();
    }


    return 0;
  };

}
