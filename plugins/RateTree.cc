/*
 * =====================================================================================
 *
 *       Filename:  RateTree.cc
 *
 *    Description:  Produce a tree for computing rates.
 *
 *         Author:  Evan Friis, evan.friis@cern.ch
 *        Company:  UW Madison
 *
 * =====================================================================================
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "TTree.h"

typedef std::vector<edm::InputTag> VInputTag;

class RateTree : public edm::EDAnalyzer {
  public:
    RateTree(const edm::ParameterSet& pset);
    virtual ~RateTree();
    void analyze(const edm::Event& evt, const edm::EventSetup& es);
  private:
    VInputTag src_;
    TTree* tree;
    std::vector<Float_t>* pts_;
    std::vector<Float_t>* etas_;
    std::vector<Float_t>* phis_;
    UInt_t run_;
    UInt_t lumi_;
    ULong64_t event_;
    edm::InputTag scalerSrc_;
    Float_t instLumi_;
    // Add UCT-only variables
    bool isUCT_;
    std::vector<Float_t>* jetPt_;
    std::vector<Float_t>* regionPt_;

    // EM versions
    std::vector<Float_t>* jetPtEM_;
    std::vector<Float_t>* regionPtEM_;

    std::vector<Float_t>* emClusterEt_;
    std::vector<Float_t>* emClusterStripEt_;
    std::vector<Float_t>* emClusterCenterEt_;
    std::vector<Float_t>* emCluster2x1Et_;
    std::vector<Int_t>* emClusterCenterFG_;
    std::vector<Int_t>* emCluster2x1FG_;

    std::vector<Float_t>* highestCenter2x1Et_;
    std::vector<Float_t>* highestNeighbor2x1Et_;

    std::vector<Int_t>* type_;
    std::vector<Int_t>* ellIso_;
    std::vector<Float_t>* pu_;
    std::vector<Float_t>* puUIC_;
    std::vector<Float_t>* puEM_;
    std::vector<Float_t>* puUICEM_;
    std::vector<Float_t>* effArea_;
    std::vector<bool>* taus_;
    std::vector<bool>* mips_;

};

RateTree::RateTree(const edm::ParameterSet& pset) {
  // Initialize the ntuple builder
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("Ntuple", "Ntuple");
  pts_ = new std::vector<Float_t>();
  etas_ = new std::vector<Float_t>();
  phis_ = new std::vector<Float_t>();
  taus_ = new std::vector<bool>();
  mips_ = new std::vector<bool>();


  tree->Branch("pt", "std::vector<float>", &pts_);
  tree->Branch("eta", "std::vector<float>", &etas_);
  tree->Branch("phi", "std::vector<float>", &phis_);
  tree->Branch("run", &run_, "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt", &event_, "evt/l");
  tree->Branch("instlumi", &instLumi_, "instlumi/F");

  type_ = new std::vector<Int_t>();
  tree->Branch("type", "std::vector<int>", &type_);

    tree->SetAlias("l1Pt", "pt");
    tree->SetAlias("l1Eta", "eta");
    tree->SetAlias("l1Phi", "phi");
    tree->SetAlias("l1Type", "type");

  src_ = pset.getParameter<VInputTag>("src");
  scalerSrc_ = pset.exists("scalerSrc") ?
    pset.getParameter<edm::InputTag>("scalerSrc") : edm::InputTag("scalersRawToDigi");
}

RateTree::~RateTree() {
  delete pts_;
  delete etas_;
  delete phis_;

  delete jetPt_;
  delete regionPt_;

  delete jetPtEM_;
  delete regionPtEM_;

  delete emClusterEt_;
  delete emClusterCenterEt_;
  delete emCluster2x1Et_;
  delete emClusterCenterFG_;
  delete emCluster2x1FG_;
  delete emClusterStripEt_;

  delete type_;
  delete ellIso_;
  delete pu_;
  delete puUIC_;
  delete puEM_;
  delete puUICEM_;
  delete effArea_;
  delete mips_;
  delete taus_;
}


namespace {

  // Predicate to sort candidates by descending pt
  class CandPtSorter {
    public:
      bool operator()(const reco::Candidate* candA, const reco::Candidate* candB)
        const {
          return candA->pt() > candB->pt();
        }
  };

  // Turn a set of InputTags into a colleciton of candidate pointers.
  std::vector<const reco::Candidate*> getCollections(
      const edm::Event& evt, const VInputTag& collections) {
    std::vector<const reco::Candidate*> output;
    // Loop over collections
    for (size_t i = 0; i < collections.size(); ++i) {
      edm::Handle<edm::View<reco::Candidate> > handle;
      evt.getByLabel(collections[i], handle);
      // Loop over objects in current collection
      for (size_t j = 0; j < handle->size(); ++j) {
        const reco::Candidate& object = handle->at(j);
        output.push_back(&object);
      }
    }
    return output;
  }

}

void RateTree::analyze(const edm::Event& evt, const edm::EventSetup& es) {

  // Get the objects.
  std::vector<const reco::Candidate*> objects = getCollections(
      evt, src_);

  std::sort(objects.begin(), objects.end(), CandPtSorter());


  // Clear previous event's objects
  pts_->clear();
  etas_->clear();
  phis_->clear();
  type_->clear();
  // Setup meta info
  run_ = evt.id().run();
  lumi_ = evt.id().luminosityBlock();
  event_ = evt.id().event();

  // Get instantaneous lumi from the scalers
  // thx to Carlo Battilana
  edm::Handle<LumiScalersCollection> lumiScalers;
  evt.getByLabel(scalerSrc_, lumiScalers);
  instLumi_ = -1;
  if (lumiScalers->size())
    instLumi_ = lumiScalers->begin()->instantLumi();

  for (size_t i = 0; i < objects.size(); ++i) {
    pts_->push_back(objects[i]->pt());
    etas_->push_back(objects[i]->eta());
    phis_->push_back(objects[i]->phi());
      // For L1 we need to get the type
      const l1extra::L1JetParticle* jetParticle =
        dynamic_cast<const l1extra::L1JetParticle*>(objects[i]);
      const l1extra::L1EmParticle* emParticle =
        dynamic_cast<const l1extra::L1EmParticle*>(objects[i]);
      if (jetParticle) {
        type_->push_back(jetParticle->type());
      } else if (emParticle) {
        type_->push_back(emParticle->type());
      } else {
        throw cms::Exception("bad input") << "Can't case L1 candidate to "
          << "either Jet or EmParticle" << std::endl;
      }

  }

  // pad everything, to work around the MaxIf bug.
  pts_->push_back(-5);
  etas_->push_back(-5);
  phis_->push_back(-5);

  tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RateTree);
