#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/PUJets.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include <boost/functional/hash.hpp>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/format.hpp"
#include "TH1.h"

namespace ic {

  PUJets::PUJets(std::string const& name) : ModuleBase(name) {
    fs_ = nullptr;
  }

  PUJets::~PUJets() {
    ;
  }

  int PUJets::PreAnalysis() {
    if (fs_) {
      double bins[5] = { 0., 2.5, 2.75, 3.0, 4.7};
      real_jet_beta_ = fs_->make<TH1F>("real_jet_beta", "real_jet_beta",100,0,1);
      fake_jet_beta_ = fs_->make<TH1F>("fake_jet_beta", "fake_jet_beta",100,0,1);
      frac_pu_before_num_ = fs_->make<TH1F>("frac_pu_before_num", "frac_pu_before_num", 4, bins);
      frac_pu_after_num_  = fs_->make<TH1F>("frac_pu_after_num", "frac_pu_after_num", 4, bins);
      frac_pu_before_den_ = fs_->make<TH1F>("frac_pu_before_den", "frac_pu_before_den", 4, bins);
      frac_pu_after_den_  = fs_->make<TH1F>("frac_pu_after_den", "frac_pu_after_den", 4, bins);
      mva1_re_num_ = fs_->make<TH1F>("mva1_re_num", "mva1_re_num", 200, -1, 1);
      mva2_re_num_ = fs_->make<TH1F>("mva2_re_num", "mva2_re_num", 200, -1, 1);
      mva3_re_num_ = fs_->make<TH1F>("mva3_re_num", "mva3_re_num", 200, -1, 1);
      mva4_re_num_ = fs_->make<TH1F>("mva4_re_num", "mva4_re_num", 200, -1, 1);
      mva1_re_den_ = fs_->make<TH1F>("mva1_re_den", "mva1_re_den", 200, -1, 1);
      mva2_re_den_ = fs_->make<TH1F>("mva2_re_den", "mva2_re_den", 200, -1, 1);
      mva3_re_den_ = fs_->make<TH1F>("mva3_re_den", "mva3_re_den", 200, -1, 1);
      mva4_re_den_ = fs_->make<TH1F>("mva4_re_den", "mva4_re_den", 200, -1, 1);
      mva1_pu_num_ = fs_->make<TH1F>("mva1_pu_num", "mva1_pu_num", 200, -1, 1);
      mva2_pu_num_ = fs_->make<TH1F>("mva2_pu_num", "mva2_pu_num", 200, -1, 1);
      mva3_pu_num_ = fs_->make<TH1F>("mva3_pu_num", "mva3_pu_num", 200, -1, 1);
      mva4_pu_num_ = fs_->make<TH1F>("mva4_pu_num", "mva4_pu_num", 200, -1, 1);
      mva1_pu_den_ = fs_->make<TH1F>("mva1_pu_den", "mva1_pu_den", 200, -1, 1);
      mva2_pu_den_ = fs_->make<TH1F>("mva2_pu_den", "mva2_pu_den", 200, -1, 1);
      mva3_pu_den_ = fs_->make<TH1F>("mva3_pu_den", "mva3_pu_den", 200, -1, 1);
      mva4_pu_den_ = fs_->make<TH1F>("mva4_pu_den", "mva4_pu_den", 200, -1, 1);
    }
    return 0;
  }

  int PUJets::Execute(TreeEvent *event) {
    EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
    // std::vector<CompositeCandidate *> const& dilepton = event->GetPtrVec<CompositeCandidate>("emtauCandidates");
    double wt = eventInfo->total_weight();
    std::vector<PFJet *> jets = event->GetPtrVec<PFJet>("pfJetsPFlow");
    std::vector<GenJet *> genjets = event->GetPtrVec<GenJet>("genJets");
    ic::erase_if(jets, !boost::bind(MinPtMaxEta, _1, 30., 4.7));
    // ic::erase_if(jets, !boost::bind(MinDRToCollection<Candidate *>, _1, dilepton[0]->AsVector(), 0.5));
    std::vector<bool> is_pu(jets.size(), false);
    for (unsigned i = 0; i < jets.size(); ++i) {
      double eta = fabs(jets[i]->eta());
      double mva = jets[i]->pu_id_mva_value();
      if (MinDRToCollection(jets[i], genjets, 0.5)) is_pu[i] = true;
      frac_pu_before_den_->Fill(eta, wt);
      if (PileupJetID(jets[i], 2)) frac_pu_after_den_->Fill(eta , wt);
      if (is_pu[i]) {
        frac_pu_before_num_->Fill(eta, wt);
        if (PileupJetID(jets[i], 2)) frac_pu_after_num_->Fill(eta, wt);
      }
      if (is_pu[i] && eta < 2.5) {
        fake_jet_beta_->Fill(jets[i]->beta(), wt);
      }
      if (!is_pu[i] && eta < 2.5) {
        real_jet_beta_->Fill(jets[i]->beta(), wt);
      }
      if (eta < 2.5 && is_pu[i]) {
        for (int i = 1; i <= mva1_pu_num_->GetNbinsX(); ++i) {
          mva1_pu_den_->Fill(mva1_pu_den_->GetXaxis()->GetBinCenter(i), wt);
          if (mva > mva1_pu_num_->GetXaxis()->GetBinLowEdge(i)) {
            mva1_pu_num_->Fill(mva1_pu_num_->GetXaxis()->GetBinCenter(i), wt);
          }
        }
      }
      if (eta < 2.5 && !is_pu[i]) {
        for (int i = 1; i <= mva1_re_num_->GetNbinsX(); ++i) {
          mva1_re_den_->Fill(mva1_re_den_->GetXaxis()->GetBinCenter(i), wt);
          if (mva > mva1_re_num_->GetXaxis()->GetBinLowEdge(i)) {
            mva1_re_num_->Fill(mva1_re_num_->GetXaxis()->GetBinCenter(i), wt);
          }
        }
      }

      if (eta >= 2.5 && eta < 2.75 && is_pu[i]) {
        for (int i = 1; i <= mva2_pu_num_->GetNbinsX(); ++i) {
          mva2_pu_den_->Fill(mva2_pu_den_->GetXaxis()->GetBinCenter(i), wt);
          if (mva > mva2_pu_num_->GetXaxis()->GetBinLowEdge(i)) {
            mva2_pu_num_->Fill(mva2_pu_num_->GetXaxis()->GetBinCenter(i), wt);
          }
        }
      }
      if (eta >= 2.5 && eta < 2.75 && !is_pu[i]) {
        for (int i = 1; i <= mva2_re_num_->GetNbinsX(); ++i) {
          mva2_re_den_->Fill(mva2_re_den_->GetXaxis()->GetBinCenter(i), wt);
          if (mva > mva2_re_num_->GetXaxis()->GetBinLowEdge(i)) {
            mva2_re_num_->Fill(mva2_re_num_->GetXaxis()->GetBinCenter(i), wt);
          }
        }
      }

      if (eta >= 2.75 && eta < 3.00 && is_pu[i]) {
        for (int i = 1; i <= mva3_pu_num_->GetNbinsX(); ++i) {
          mva3_pu_den_->Fill(mva3_pu_den_->GetXaxis()->GetBinCenter(i), wt);
          if (mva > mva3_pu_num_->GetXaxis()->GetBinLowEdge(i)) {
            mva3_pu_num_->Fill(mva3_pu_num_->GetXaxis()->GetBinCenter(i), wt);
          }
        }
      }
      if (eta >= 2.75 && eta < 3.00 && !is_pu[i]) {
        for (int i = 1; i <= mva3_re_num_->GetNbinsX(); ++i) {
          mva3_re_den_->Fill(mva3_re_den_->GetXaxis()->GetBinCenter(i), wt);
          if (mva > mva3_re_num_->GetXaxis()->GetBinLowEdge(i)) {
            mva3_re_num_->Fill(mva3_re_num_->GetXaxis()->GetBinCenter(i), wt);
          }
        }
      }

      if (eta >= 3.00 && eta < 4.79 && is_pu[i]) {
        for (int i = 1; i <= mva4_pu_num_->GetNbinsX(); ++i) {
          mva4_pu_den_->Fill(mva4_pu_den_->GetXaxis()->GetBinCenter(i), wt);
          if (mva > mva4_pu_num_->GetXaxis()->GetBinLowEdge(i)) {
            mva4_pu_num_->Fill(mva4_pu_num_->GetXaxis()->GetBinCenter(i), wt);
          }
        }
      }
      if (eta >= 3.00 && eta < 4.79 && !is_pu[i]) {
        for (int i = 1; i <= mva4_re_num_->GetNbinsX(); ++i) {
          mva4_re_den_->Fill(mva4_re_den_->GetXaxis()->GetBinCenter(i), wt);
          if (mva > mva4_re_num_->GetXaxis()->GetBinLowEdge(i)) {
            mva4_re_num_->Fill(mva4_re_num_->GetXaxis()->GetBinCenter(i), wt);
          }
        }
      }
    }
    return 0;
  }
  int PUJets::PostAnalysis() {
    return 0;
  }

  void PUJets::PrintInfo() {
    ;
  }
}
