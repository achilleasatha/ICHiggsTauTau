#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HiggsNuNuAnalysisTools.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeAnalyser.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeModule.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeFiles.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataWEst.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataQCDEst.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/NormPlots.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/MVATrain.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/MVAApply.h"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"

namespace po=boost::program_options;
using namespace ic;

int main(int argc, char* argv[]){
  /*##########################################
  #                                          #
  #            SET UP OPTIONS                #
  #                                          #
  ##########################################*/

  //Input output and config options
  std::string cfg;
  std::string outputname;
  std::string inputfolder;
  std::string inputparams;
  std::string filelist;
  
  //Base selection cuts
  double jet1ptcut;
  double jet2ptcut;
  double metcut;
  double mjjcut;
  double cjvcut;
  double dphicut;
  double detacut;
  double etaprodcut;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options
    ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))
    ("input_folder,i",           po::value<std::string>(&inputfolder)->default_value("../output_lighttree_withrle2"))//"root://xrootd.grid.hep.ph.ic.ac.uk//store/user/pdunne/lighttree/"))
    ("input_params,p",           po::value<std::string>(&inputparams)->default_value("../filelists/Dec18/ParamsDec18test.dat"))
    ("filelist,f",               po::value<std::string>(&filelist)->default_value("filelists/filelist.dat"))
    //Base selection cuts
    ("jet1ptcut",                po::value<double>(&jet1ptcut)->default_value(50))
    ("jet2ptcut",                po::value<double>(&jet2ptcut)->default_value(50))
    ("metcut",                   po::value<double>(&metcut)->default_value(130))
    ("mjjcut",                   po::value<double>(&mjjcut)->default_value(1100))
    ("cjvcut",                   po::value<double>(&cjvcut)->default_value(1))
    ("dphicut",                  po::value<double>(&dphicut)->default_value(1))
    ("detacut",                  po::value<double>(&detacut)->default_value(4.2))
    ("etaprodcut",               po::value<double>(&etaprodcut)->default_value(0));
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  /*##########################################
  #                                          #
  #          INSTANTIATE ANALYSER            #
  #                                          #
  ##########################################*/

  LTAnalyser* analysis = new LTAnalyser(outputname);

  analysis->AddFiles(filelist);

  analysis->SetInFolder(inputfolder);
  analysis->SetInputParams(inputparams);

  //Set selection step common to all categories
  analysis->set_baseselection("");//passtrigger==1&& jet1_pt>"+boost::lexical_cast<std::string>(jet1ptcut)+"&& jet2_pt>"+boost::lexical_cast<std::string>(jet2ptcut)+" && dijet_M >"+boost::lexical_cast<std::string>(mjjcut)+"&& jet1_eta*jet2_eta<="+boost::lexical_cast<std::string>(etaprodcut)+"&& dijet_dphi<"+boost::lexical_cast<std::string>(dphicut)+"&& dijet_deta >"+boost::lexical_cast<std::string>(detacut));

  /*##########################################
  #                                          #
  #            DEFINE MODULES                #
  #                                          #
  ##########################################*/

  //MVA TRAIN
  std::vector<std::string> sigsets;
  sigsets.push_back("sig125");
  std::vector<std::string> bkgsets;
  bkgsets.push_back("VBF-QCD");
  bkgsets.push_back("VV");
  bkgsets.push_back("Top");
  bkgsets.push_back("ZJets_ll");
  bkgsets.push_back("ZJets_ll_vbf");
  bkgsets.push_back("ZJets_nunu");
  bkgsets.push_back("WJets_enu");
  bkgsets.push_back("WJets_munu");
  bkgsets.push_back("WJets_taunu");
  std::vector<std::string> variables;
  variables.push_back("jetmet_mindphi");
  variables.push_back("dijetmet_ptfraction");
  variables.push_back("dijetmet_vectorialSum_pt");
  variables.push_back("jet1met_scalarprod_frac := jet1met_scalarprod/met");
  variables.push_back("jet2met_scalarprod_frac := jet2met_scalarprod/met");
  variables.push_back("n_jets_cjv_30");
  variables.push_back("jet2_pt");
  variables.push_back("dijet_M");
  variables.push_back("dijet_deta");
  variables.push_back("dijet_dphi");
  variables.push_back("met");
  variables.push_back("met_significance");
  variables.push_back("sqrt(ht)");


//   variables.push_back("jet2_eta");
//   variables.push_back("met");
//   variables.push_back("met_significance");
//   variables.push_back("mht");
//   variables.push_back("jet2met_dphi");
//   variables.push_back("jetmet_mindphi");
//   variables.push_back("jetunclet_mindphi");
//   variables.push_back("metunclet_dphi");
//   variables.push_back("dijetmet_scalarSum_pt");
//   variables.push_back("dijetmet_vectorialSum_pt");
//   variables.push_back("dijetmet_ptfraction");
//   variables.push_back("jet1met_scalarprod");
//   variables.push_back("jet2met_scalarprod");
//   variables.push_back("jet1met_scalarprod_frac := jet1met_scalarprod/met");
//   variables.push_back("jet2met_scalarprod_frac := jet2met_scalarprod/met");
  std::vector<std::string> specvariables;

  MVATrain mvatrainer("mvatrainer");
  mvatrainer.set_sigsets(sigsets)
    .set_bkgsets(bkgsets)
    .set_variables(variables)
    .set_specvariables(specvariables)
    .set_basesel("nvetomuons==0&&nvetoelectrons==0&&jet1_pt>50&&jet2_pt>50&&jet1_eta<4.7&&jet2_eta<4.7&&metnomu_significance>3&&dijet_deta>3.6&&jetmetnomu_mindphi>1.5")
    .set_sigcat("")
    .set_bkgcat("");

  DataWEst premvaest("premvaest");
    premvaest.set_sigmcset("sig125")
    .set_contmcset("sig125")
    .set_contdataset("sig125")
    .set_contbkgset(bkgsets)
    .set_basesel("nvetomuons==0&&nvetoelectrons==0&&jet1_pt>50&&jet2_pt>50&&jet1_eta<4.7&&jet2_eta<4.7&&metnomu_significance>3&&dijet_deta>3.6&&jetmetnomu_mindphi>1.5")
      .set_sigcat("")
      .set_contcat("");
  

  /*##########################################
  #                                          #
  #   SET UP ANALYSIS SEQUENCE AND RUN       #
  #                                          #
  ##########################################*/

  analysis->AddModule(&premvaest);
  analysis->AddModule(&mvatrainer);

  analysis->RunAnalysis();

  return 0;

}
