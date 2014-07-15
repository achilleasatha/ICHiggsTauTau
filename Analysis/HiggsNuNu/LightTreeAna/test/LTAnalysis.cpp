#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HiggsNuNuAnalysisTools.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeAnalyser.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeModule.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeFiles.h"
<<<<<<< HEAD
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataWEst.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataQCDEst.h"
=======
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataNormShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataZNormShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataZEst.h"
>>>>>>> 80cf1d29fd51d2d7b2beec8256ad805433385814
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/NormPlots.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/APlot.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/SimplePlots.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/MVATrain.h"
<<<<<<< HEAD
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/ShapeSignal.h"
#include <boost/algorithm/string.hpp>
=======
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/AddFriends.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/Plotter.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/HistPlotter.h"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "TColor.h"
>>>>>>> 80cf1d29fd51d2d7b2beec8256ad805433385814

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
  std::string shape_variable;
  
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
    ("input_folder,i",           po::value<std::string>(&inputfolder)->default_value("../output_lighttree_withrle2/"))
    ("input_params,p",           po::value<std::string>(&inputparams)->default_value("../filelists/Dec18/ParamsDec18test.dat"))
    ("filelist,f",               po::value<std::string>(&filelist)->default_value("filelists/filelist.dat"))
    ("shape_variable",           po::value<std::string>(&shape_variable)->default_value("jet2_pt(200,0.,1000.)"))
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
  std::cout<<"Taking input from: "<<inputfolder<<std::endl;
  analysis->SetInFolder(inputfolder);
  analysis->SetInputParams(inputparams);

  //Set selection step common to all categories
<<<<<<< HEAD
  analysis->set_baseselection("passtrigger==1&& jet1_pt>"+boost::lexical_cast<std::string>(jet1ptcut)+"&& jet2_pt>"+boost::lexical_cast<std::string>(jet2ptcut)+" && dijet_M >"+boost::lexical_cast<std::string>(mjjcut)+"&& jet1_eta*jet2_eta<="+boost::lexical_cast<std::string>(etaprodcut)+"&& dijet_dphi<"+boost::lexical_cast<std::string>(dphicut)+"&& dijet_deta >"+boost::lexical_cast<std::string>(detacut));

=======
  //analysis->set_baseselection("passtrigger==1&&jet1_eta<4.7&&jet2_eta<4.7&& jet1_pt>"+boost::lexical_cast<std::string>(jet1ptcut)+"&& jet2_pt>"+boost::lexical_cast<std::string>(jet2ptcut)+" && dijet_M >"+boost::lexical_cast<std::string>(mjjcut)+"&& jet1_eta*jet2_eta<"+boost::lexical_cast<std::string>(etaprodcut)+"&& dijet_dphi<"+boost::lexical_cast<std::string>(dphicut)+"&& dijet_deta >"+boost::lexical_cast<std::string>(detacut));
  //analysis->set_baseselection("jet1_pt>50&&jet2_pt>50&&dijet_deta>3.6&&metnomu_significance>3&&jetmetnomu_mindphi>1.5");
  analysis->set_baseselection("jet1_pt>50&&jet2_pt>50&&dijet_deta>3.6&&metnomu_significance>3&&jetmetnomu_mindphi>1.5");
>>>>>>> 80cf1d29fd51d2d7b2beec8256ad805433385814
  /*##########################################
  #                                          #
  #            DEFINE MODULES                #
  #                                          #
  ##########################################*/

  //ADD FRIENDS TO GET EXTRA VARIABLES
  std::vector<std::string> setswithfriends;
  setswithfriends.push_back("WJets_taunu");
  setswithfriends.push_back("WJets_munu");
  setswithfriends.push_back("WJets_enu");
  setswithfriends.push_back("VBF-QCD");
  setswithfriends.push_back("MET");
  setswithfriends.push_back("PARKED");
  setswithfriends.push_back("PARKEDPLUSA");
  setswithfriends.push_back("VV");
  setswithfriends.push_back("Top");
  setswithfriends.push_back("ZJets_ll");
  setswithfriends.push_back("ZJets_ll_vbf");
  setswithfriends.push_back("ZJets_nunu");
  setswithfriends.push_back("sig125");
  
  AddFriends addfriends("addfriends");
  addfriends.set_frienddir("output/")
    .set_friendtreename("mvafriend")
    .set_sets(setswithfriends);

  std::vector<std::string> shape;
  //shape.push_back("BDT(12,-1.,0.2)");
  shape.push_back("jet2_pt(27,30.,300.)");
  shape.push_back("jet1_pt(27,30.,300.)");
  shape.push_back("metnomuons(30,50.,200.)");
  shape.push_back("dijet_M(35,600.,2000.)");
  shape.push_back("jetmetnomu_mindphi(18,1.4,3.2)");
  shape.push_back("metnomu_significance(25,3.,8.)");
  shape.push_back("dijet_dphi(30,0.,3.)");
  shape.push_back("dijet_deta(17,3.6,7.)");
  shape.push_back("dijetmetnomu_ptfraction(20,0.,1.)");
  std::string dataset="PARKEDPLUSA";
  std::string sigcat="nvetomuons==0&&nvetoelectrons==0";
  std::string zextrasigcat="";//"&&nvetomuons==3";

  //DATA SHAPE GENERATION
  DataShape data("data");
  data.set_dataset(dataset)
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+"&&(passtrigger==1||passparkedtrigger1==1||passparkedtrigger2==1)");

  DataShape signal("signal");
  signal.set_dataset("sig125")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);
  
  DataShape vv("vv");
  vv.set_dataset("VV")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape top("top");
  top.set_dataset("Top")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);


  //WBKG SHAPE GENERATION
  std::vector<std::string> Wcontbkgsets; //List of sets for ncbkg
  Wcontbkgsets.push_back("VV");
  Wcontbkgsets.push_back("Top");
  Wcontbkgsets.push_back("ZJets_ll");
  Wcontbkgsets.push_back("ZJets_ll_vbf");
  Wcontbkgsets.push_back("ZJets_nunu");

  DataNormShape wmunu("wmunu");
  wmunu.set_sigmcset("WJets_munu")
    .set_shape(shape)
    .set_contmcset("WJets_munu")
    .set_contdataset(dataset)
    .set_contbkgset(Wcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel("&&(passtrigger==1||passparkedtrigger1==1||passparkedtrigger2==1)")
    .set_sigcat(sigcat)
    .set_contcat("nvetoelectrons==0 && nvetomuons==1 && nselmuons==1");

  DataNormShape wenu("wenu");
  wenu.set_sigmcset("WJets_enu")
    .set_shape(shape)
    .set_contmcset("WJets_enu")
    .set_contdataset(dataset)
    .set_contbkgset(Wcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel("&&(passtrigger==1||passparkedtrigger1==1||passparkedtrigger2==1)")
    .set_sigcat(sigcat)
    .set_contcat("nselelectrons==1 && nvetoelectrons ==1 && nvetomuons==0");

<<<<<<< HEAD
=======
  DataNormShape wtaunu("wtaunu");
  wtaunu.set_sigmcset("WJets_taunu")
    .set_shape(shape)
    .set_contmcset("WJets_taunu")
    .set_contdataset(dataset)
    .set_contbkgset(Wcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel("&&(passtrigger==1||passparkedtrigger1==1||passparkedtrigger2==1)")
    .set_sigcat(sigcat)
    .set_contcat("ntaus>=1&&nvetoelectrons ==0 && nvetomuons==0")
     .set_sigmcweight("total_weight_lepveto")
     .set_contmcweight("total_weight_lepveto")
     .set_contdataweight("weight_nolep");

  //QCDBKG
  std::vector<std::string> QCDcontbkgsets; //list of sets for ncbkg
  QCDcontbkgsets.push_back("VV");
  QCDcontbkgsets.push_back("Top");
  QCDcontbkgsets.push_back("ZJets_ll");
  QCDcontbkgsets.push_back("ZJets_ll_vbf");
  QCDcontbkgsets.push_back("ZJets_nunu");
  QCDcontbkgsets.push_back("WJets_enu");
  QCDcontbkgsets.push_back("WJets_munu");
  QCDcontbkgsets.push_back("WJets_taunu");
  
  DataNormShape QCD("QCD");
  QCD.set_sigmcset("VBF-QCD")
    .set_shape(shape)
    .set_contmcset("VBF-QCD")
    .set_contdataset(dataset)
    .set_contbkgset(QCDcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel("&&(passtrigger==1||passparkedtrigger1==1||passparkedtrigger2==1)")
    .set_sigcat(sigcat)
    .set_contcat("nvetoelectrons==0&&nvetomuons==0&&dijetmet_ptfraction>0.6");

  //ZBKG
  std::vector<std::string> Zcontbkgsets;//List of sets for ncbkg !!FILL THIS
  Zcontbkgsets.push_back("VV");
  Zcontbkgsets.push_back("Top");
  Zcontbkgsets.push_back("WJets_enu");
  Zcontbkgsets.push_back("WJets_munu");
  Zcontbkgsets.push_back("WJets_taunu");
  
  DataZNormShape zmumu("zmumu");
  zmumu.set_sigmcewkset("ZJets_ll_vbf")
    .set_shape(shape)
    .set_sigmcqcdset("ZJets_ll")
    .set_contmcewkset("ZJets_ll_vbf")
    .set_contmcqcdset("ZJets_ll")
   .set_contbkgset(Zcontbkgsets)
    .set_contdataset(dataset)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel("&&(passtrigger==1||passparkedtrigger1==1||passparkedtrigger2==1)")
    .set_sigcat("m_mumu_gen>80&&m_mumu_gen<100"+zextrasigcat)//!!MAKE GEN WINDOW CONFIGURABLE WITH WARNING ABOUT GEN MASS FILTERED NUMBER BELOW NOT CHANGING
    .set_contcat("nvetoelectrons==0 && nvetomuons==2 && nselmuons==2&&m_mumu>60&&m_mumu<120")//!!MAKE GEN WINDOW CONFIGURABLE
    .set_sigmainccontewk(303)
    .set_sigmainccontqcd(3503700./3)
    .set_sigmaincsigewk(460*3)
    .set_sigmaincsigqcd(3503700./3*5.651)
    .set_ngenincewk(5781.91)
    .set_ngenincqcd(22789300)
    .set_ngenmassfilteredewk(4226.53)
    .set_ngenmassfilteredqcd(20334900);
    
    
>>>>>>> 80cf1d29fd51d2d7b2beec8256ad805433385814


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
  shapes.push_back("met(80,0.,400.)");
  shapes.push_back("cjvjetpt(100,0.,100.)");

  NormPlots normplots("normplots");
  normplots.set_qcdset("QCD")
    .set_sigset("sig125")
    .set_ewkset(ewksets)
    .set_cat("")
    .set_basesel("jet1_eta<4.7&&jet2_eta<4.7&&jet1_pt>50&&jet2_pt>50&&nvetoelectrons==0 && nvetomuons==0&&dijet_M>150&&met>130")
    .set_shapes(shapes);

  std::vector<std::vector<std::string> > sets;
  std::vector<std::string> otherbckgs;
  std::vector<std::string> zjets;
  std::vector<std::string> wjets;
  otherbckgs.push_back("WJets_munu");
  otherbckgs.push_back("WJets_enu");
  // zjets.push_back("ZJets_ll");
  //zjets.push_back("ZJets_ll_vbf");
  // zjets.push_back("ZJets_nunu");
//   wjets.push_back("WJets_enu");
//   wjets.push_back("WJets_munu");
//   wjets.push_back("WJets_taunu");
  sets.push_back(otherbckgs); 
//   sets.push_back(zjets);
//   sets.push_back(wjets);

  std::vector<std::string> setsnames;
  std::string otherbckg;
  std::string zjet;
  std::string wjet;
  setsnames.push_back(otherbckg);
  setsnames.push_back(zjet);
  setsnames.push_back(wjet);

  std::vector<std::string> shape;
//   shape.push_back("run(100,0.,3.)");
//   shape.push_back("lumi,(100,0.,350.)");
//   shape.push_back("event(100,0.,0.1)");
//   shape.push_back("weight_nolep(100,0.,0.65)");
//   shape.push_back("total_weight_lepveto(100,0.,0.65)");
//   shape.push_back("total_weight_leptight(100,0.,0.65)");
     shape.push_back("jet1_pt(200,0.,2000.)");
//      shape.push_back("jet2_pt(200,0.,2000.)");
//      shape.push_back("jet1_eta(200,0.,5.)");
//      shape.push_back("jet2_eta(200,0.,5.)");
 
     std::vector<std::string> drawopt;
     drawopt.push_back("");

 APlot aplot("aplot");
  aplot.set_sets(sets)
    .set_basesel("jet1_eta<4.7&&jet2_eta<4.7&&jet1_pt>50&&jet2_pt>50&&nvetoelectrons==0 && nve\
tomuons==0&&dijet_M>150&&met>130")
    .set_cat("")
    .set_shapes(shape)
    .set_drawopts(drawopt);

  //MVA TRAIN
  std::vector<std::string> sigsets;
  sigsets.push_back("sig125");
  std::vector<std::string> bkgsets;
  bkgsets.push_back("QCD");
  std::vector<std::string> variables;
  variables.push_back("jet2_eta");
  variables.push_back("met");
  variables.push_back("met_significance");
  variables.push_back("mht");
  variables.push_back("jet2met_dphi");
  variables.push_back("jetmet_mindphi");
  variables.push_back("jetunclet_mindphi");
  variables.push_back("metunclet_dphi");
  variables.push_back("dijetmet_scalarSum_pt");
  variables.push_back("dijetmet_vectorialSum_pt");
  variables.push_back("dijetmet_ptfraction");
  variables.push_back("jet1met_scalarprod");
  variables.push_back("jet2met_scalarprod");
  variables.push_back("jet1met_scalarprod_frac := jet1met_scalarprod/met");
  variables.push_back("jet2met_scalarprod_frac := jet2met_scalarprod/met");
  std::vector<std::string> specvariables;

  MVATrain mvatrainer("mvatrainer");
  mvatrainer.set_sigsets(sigsets)
    .set_bkgsets(bkgsets)
    .set_variables(variables)
    .set_specvariables(specvariables)
    .set_basesel("passtrigger==1&&nvetomuons==0&&nvetoelectrons==0&&jet1_pt>50&&jet2_pt>50&&jet1_eta<4.7&&jet2_eta<4.7")
    .set_sigcat("")
    .set_bkgcat("");

<<<<<<< HEAD
  //ShapeSignal_munu
  ShapeSignal shapesignal_munu("shapesignal_munu");
  shapesignal_munu.set_sigmcset("WJets_munu")
    .set_contmcset("WJets_munu")
    .set_contdataset("MET")
    .set_contbkgset(Wcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_sigcat("nvetoelectrons==0 && nvetomuons==0&& met>"+boost::lexical_cast<std::string>(metcut)+" && n_jets_cjv_30<"+boost::lexical_cast<std::string>(cjvcut))
    .set_contcat("nvetoelectrons==0 && nvetomuons==1 && nselmuons==1&& metnomuons>"+boost::lexical_cast<std::string>(metcut)+" && n_jets_cjv_30<"+boost::lexical_cast<std::string>(cjvcut));  
  //ShapeSignal_enu
  ShapeSignal shapesignal_enu("shapesignal_enu");
  shapesignal_enu.set_sigmcset("WJets_enu")
    .set_contmcset("WJets_enu")
    .set_contdataset("MET")
    .set_contbkgset(Wcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_sigcat("nvetoelectrons==0 && nvetomuons==0&& met>"+boost::lexical_cast<std::string>(metcut)+" && n_jets_cjv_30<"+boost::lexical_cast<std::string>(cjvcut))
    .set_contcat("nselelectrons==1 && nvetoelectrons ==1 && nvetomuons==0&& met>"+boost::lexical_cast<std::string>(metcut)+" && n_jets_cjv_30<"+boost::lexical_cast<std::string>(cjvcut));
=======
  //PLOTTER
//   std::vector<std::vector<std::string> > plottersets;
//   std::vector<std::string> sets1;
//   sets1.push_back("sig125");
//   std::vector<std::string> sets2;
//   sets2.push_back("WJets_enu");
//   plottersets.push_back(sets1);
//   plottersets.push_back(sets2);
 
//   std::vector<std::string> plottervariables;
//   plottervariables.push_back("BDT(40,-1.,1.)");

//   Plotter plotter("plotter");
//   plotter.set_sets(plottersets)
//     .set_shapes(plottervariables)
//     .set_basesel("passtrigger==1&&nvetomuons==0&&nvetoelectrons==0&&jet1_pt>50&&jet2_pt>50&&jet1_eta<4.7&&jet2_eta<4.7&&met_significance>3&&dijet_deta>3.6&&jetmet_mindphi>1.5")
//     .set_cat("");
  
  //HISTPLOTTER
  std::vector<std::string> shapevec;
  for(unsigned ishape=0;ishape<shape.size();ishape++){
    std::vector<std::string> strs;
    boost::split(strs, shape[ishape], boost::is_any_of("("));
    shapevec.push_back(strs[0]);
  }

  std::vector<LTPlotElement> elementvec;
  LTPlotElement dataele;
  dataele.set_is_data(true)
    .set_scale(1)
    .set_legname("Data")
    .set_sample("data");

  LTPlotElement wmunuele;
  wmunuele.set_is_data(false)
    .set_scale(1)
    .set_color(kOrange-4)
    .set_in_stack(true)
    .set_legname("W#rightarrow#mu#nu")
    .set_sample("WJets_munu");

  LTPlotElement wenuele;
  wenuele.set_is_data(false)
    .set_scale(1)
    .set_color(kOrange  + 2)
    .set_in_stack(true)
    .set_legname("W#rightarrow e#nu")
    .set_sample("WJets_enu");

  LTPlotElement wtaunuele;
  wtaunuele.set_is_data(false)
    .set_scale(1)
    .set_color(kOrange + 4)
    .set_in_stack(true)
    .set_legname("W#rightarrow#tau#nu")
    .set_sample("WJets_taunu");

  LTPlotElement zmumuele;
  zmumuele.set_is_data(false)
    .set_scale(1)
    .set_color(kAzure  + 2)
    .set_in_stack(true)
    .set_legname("Z#rightarrow#nu#nu")
    .set_sample("zmumu");

  LTPlotElement qcdele;
  qcdele.set_is_data(false)
    .set_scale(1)
    .set_color(kMagenta-10)
    .set_in_stack(true)
    .set_legname("QCD")
    .set_sample("VBF-QCD");

  LTPlotElement vvele;
  vvele.set_is_data(false)
    .set_scale(1)
    .set_color(kGreen-5)
    .set_in_stack(true)
    .set_legname("VV")
    .set_sample("vv");

  LTPlotElement topele;
  topele.set_is_data(false)
    .set_scale(1)
    .set_color(kBlue-8)
    .set_in_stack(true)
    .set_legname("Top")
    .set_sample("top");

  LTPlotElement sigele;
  sigele.set_is_data(false)
    .set_scale(20)
    .set_color(kRed)
    .set_in_stack(false)
    .set_legname("Signalx20")
    .set_sample("signal");

  elementvec.push_back(dataele);
  elementvec.push_back(wmunuele);
  elementvec.push_back(wenuele);
  elementvec.push_back(wtaunuele);
  elementvec.push_back(zmumuele);
  elementvec.push_back(qcdele);
  elementvec.push_back(vvele);
  elementvec.push_back(topele);
  elementvec.push_back(sigele);

  HistPlotter plotter("plotter");
  plotter.set_dirname("ControlPlots")
    .set_elements(elementvec)
    .set_shapes(shapevec);
  
  
>>>>>>> 80cf1d29fd51d2d7b2beec8256ad805433385814

  /*##########################################
  #                                          #
  #   SET UP ANALYSIS SEQUENCE AND RUN       #
  #                                          #
  ##########################################*/
<<<<<<< HEAD
=======
  
  //analysis->AddModule(&addfriends);
  //analysis->AddModule(&mvatrainer);
  //analysis->AddModule(&normplots);
  analysis->AddModule(&wmunu);
  analysis->AddModule(&wenu);
  analysis->AddModule(&wtaunu);  
  analysis->AddModule(&QCD);
  analysis->AddModule(&zmumu);
  analysis->AddModule(&vv);
  analysis->AddModule(&top);
  analysis->AddModule(&data);
  analysis->AddModule(&signal);
  analysis->AddModule(&plotter);
>>>>>>> 80cf1d29fd51d2d7b2beec8256ad805433385814

  //   analysis->AddModule(&mvatrainer);
  //   analysis->AddModule(&normplots);
  //   analysis->AddModule(&wmunu);
  //   analysis->AddModule(&aplot);
  //   analysis->AddModule(&wenu);
  //   analysis->AddModule(&QCD);
  analysis->AddModule(&shapesignal_munu);
  analysis->AddModule(&shapesignal_enu);
  
  analysis->RunAnalysis();

  return 0;

}
