#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HiggsNuNuAnalysisTools.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeAnalyser.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeModule.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeFiles.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataWEst.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/NormPlots.h"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"

namespace po=boost::program_options;
using namespace ic;

int main(int argc, char* argv[]){
  std::string cfg;
  std::string outputname;
  std::string inputfolder;
  std::string inputparams;
  std::string filelist;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))
    ("input_folder,i",           po::value<std::string>(&inputfolder)->default_value("../output_lighttree_newcpptrial/"))
    ("input_params,p",           po::value<std::string>(&inputparams)->default_value("../filelists/Dec18/ParamsDec18test.dat"))
    ("filelist,f",               po::value<std::string>(&filelist)->default_value("filelist.dat"));
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  LTAnalyser* analysis = new LTAnalyser(outputname);

  analysis->AddFiles(filelist);

  analysis->SetInFolder(inputfolder);
  analysis->SetInputParams(inputparams);

  //Set selection step common to all categories
  double jet1ptcut=50;
  double jet2ptcut=50;
  double metcut=130;
  double mjjcut=1100;
  double cjvcut=1;
  double dphicut=1;
  double detacut=4.2;
  double etaprodcut=0;
  analysis->set_baseselection("passtrigger==1&& jet1_pt>"+boost::lexical_cast<std::string>(jet1ptcut)+"&& jet2_pt>"+boost::lexical_cast<std::string>(jet2ptcut)+" && dijet_M >"+boost::lexical_cast<std::string>(mjjcut)+"&& jet1_eta*jet2_eta<="+boost::lexical_cast<std::string>(etaprodcut)+" && n_jets_cjv_30<"+boost::lexical_cast<std::string>(cjvcut)+"&& dijet_dphi<"+boost::lexical_cast<std::string>(dphicut)+"&& dijet_deta >"+boost::lexical_cast<std::string>(detacut));

  //Define Modules

  //WBKG
  std::vector<std::string> contbkgsets; //List of sets for ncbkg
  contbkgsets.push_back("VV");
  contbkgsets.push_back("Top");
  contbkgsets.push_back("ZJets_ll");
  contbkgsets.push_back("ZJets_ll_vbf");
  contbkgsets.push_back("ZJets_nunu");

  DataWEst wmunu("wmunu");
  wmunu.set_sigmcset("WJets_munu")
    .set_contmcset("WJets_munu")
    .set_contdataset("MET")
    .set_contbkgset(contbkgsets)
    .set_basesel(analysis->baseselection())
    .set_sigcat("nvetoelectrons==0 && nvetomuons==0&& met>"+boost::lexical_cast<std::string>(metcut))
    .set_contcat("nvetoelectrons==0 && nvetomuons==1 && nselmuons==1&& metnomuons>"+boost::lexical_cast<std::string>(metcut));

  DataWEst wenu("wenu");
  wenu.set_sigmcset("WJets_enu")
    .set_contmcset("WJets_enu")
    .set_contdataset("MET")
    .set_contbkgset(contbkgsets)
    .set_basesel(analysis->baseselection())
    .set_sigcat("nvetoelectrons==0 && nvetomuons==0&& met>"+boost::lexical_cast<std::string>(metcut))
    .set_contcat("nselelectrons==1 && nvetoelectrons ==1 && nvetomuons==0&& met>"+boost::lexical_cast<std::string>(metcut));

  //NORMALISED PLOTS FOR REFEREE
  std::vector<std::string> ewksets; //List of sets for ewk
  ewksets.push_back("VV");
  ewksets.push_back("Top");
  ewksets.push_back("ZJets_ll");
  ewksets.push_back("ZJets_ll_vbf");
  ewksets.push_back("ZJets_nunu");
  ewksets.push_back("WJets_enu");
  ewksets.push_back("WJets_munu");
  ewksets.push_back("WJets_taunu");

  std::vector<std::string> shapes; //List of shapes to draw
  shapes.push_back("dijet_M(370,150.,2000.)");
  shapes.push_back("dijet_deta(160,0.,8.)");
  shapes.push_back("dijet_dphi(310,0.,3.1)");
  // shapes.push_back("jet1_eta(200,-5.,5.)");
  //   shapes.push_back("jet1_pt(200,0.,1000.)");
  //   shapes.push_back("jet2_eta(200,-5.,5.)");
  //   shapes.push_back("jet2_pt(200,0.,1000.)");
  shapes.push_back("met(80,0.,400.)");
  shapes.push_back("cjv_30_jet3pt(100,0.,100.)");
  // shapes.push_back("n_jets_cjv_30(15,0.,15.)");
  //   shapes.push_back("n_jets_cjv_20EB_30EE(15,0.,15.)");
  //   shapes.push_back("met_significance(400,0.,20.)");
  //   shapes.push_back("dijet_sumeta(200,-5.,5.)");


  NormPlots normplots("normplots");
  normplots.set_qcdset("QCD")
    .set_sigset("sig125")
    .set_ewkset(ewksets)
    .set_cat("")
    .set_basesel("jet1_eta<4.7&&jet2_eta<4.7&&jet1_pt>30&&jet2_pt>30&&nvetoelectrons==0 && nvetomuons==0&&dijet_M>150&&met>40")
    .set_shapes(shapes);

  //analysis->AddModule(&normplots);
  analysis->AddModule(&wmunu);
  analysis->AddModule(&wenu);

  analysis->RunAnalysis();

  return 0;

}