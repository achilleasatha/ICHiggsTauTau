#ifndef ICHiggsTauTau_Module_PUJets_h
#define ICHiggsTauTau_Module_PUJets_h

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HistoSet.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"


#include <string>

namespace ic {

class PUJets : public ModuleBase {
 private:
  CLASS_MEMBER(PUJets, fwlite::TFileService*, fs)
 private:
  TH1F * real_jet_beta_;
  TH1F * fake_jet_beta_;
  TH1F * frac_pu_before_num_;
  TH1F * frac_pu_after_num_;
  TH1F * frac_pu_before_den_;
  TH1F * frac_pu_after_den_;
  TH1F * mva1_re_num_;
  TH1F * mva2_re_num_;
  TH1F * mva3_re_num_;
  TH1F * mva4_re_num_;
  TH1F * mva1_re_den_;
  TH1F * mva2_re_den_;
  TH1F * mva3_re_den_;
  TH1F * mva4_re_den_;
  TH1F * mva1_pu_num_;
  TH1F * mva2_pu_num_;
  TH1F * mva3_pu_num_;
  TH1F * mva4_pu_num_;
  TH1F * mva1_pu_den_;
  TH1F * mva2_pu_den_;
  TH1F * mva3_pu_den_;
  TH1F * mva4_pu_den_;

 public:
  PUJets(std::string const& name);
  virtual ~PUJets();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

};

}

#endif
