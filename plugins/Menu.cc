/*
 * =====================================================================================
 *
 *       Filename:  Menu.cc
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
#include "L1Trigger/UCT2015/interface/UCTCandidate.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "TTree.h"

typedef std::vector<edm::InputTag> VInputTag;

class Menu : public edm::EDAnalyzer {
        public:
                Menu(const edm::ParameterSet& pset);
                virtual ~Menu();
                void analyze(const edm::Event& evt, const edm::EventSetup& es);
        private:
                VInputTag srcEG_;
                VInputTag srcIsoEG_;
                VInputTag srcJET_;
                VInputTag srcTAU_;
                edm::InputTag srcMET_;
                edm::InputTag srcMHT_;
                edm::InputTag srcSHT_;
                edm::InputTag srcSET_;

                TTree* tree;
                std::vector<Float_t>* ptsEG_;
                std::vector<Float_t>* etasEG_;
                std::vector<Float_t>* phisEG_;
                std::vector<Float_t>* ptsJetEG_;
                std::vector<Float_t>* ptsRegionEG_;
                std::vector<Float_t>* MIPEG_;
                std::vector<Float_t>* tauVetoEG_;
                std::vector<Float_t>* LIsoEG_;

                std::vector<Float_t>* ptsIsoEG_;
                std::vector<Float_t>* etasIsoEG_;
                std::vector<Float_t>* phisIsoEG_;
                std::vector<Float_t>* ptsJetIsoEG_;
                std::vector<Float_t>* ptsRegionIsoEG_;
                std::vector<Float_t>* MIPIsoEG_;
                std::vector<Float_t>* tauVetoIsoEG_;
                std::vector<Float_t>* LIsoIsoEG_;

                std::vector<Float_t>* ptsJET_;
                std::vector<Float_t>* etasJET_;
                std::vector<Float_t>* phisJET_;

                std::vector<Float_t>* ptsTAU_;
                std::vector<Float_t>* etasTAU_;
                std::vector<Float_t>* phisTAU_;
                std::vector<Float_t>* ptsJetTAU_;
                std::vector<Float_t>* ptsRegionTAU_;

                Float_t MET_;
                Float_t METPHI_;
                Float_t MHT_;
                Float_t MHTPHI_;
                Float_t SET_;
                Float_t SHT_;     

                UInt_t run_;
                UInt_t lumi_;
                ULong64_t event_;
                edm::InputTag scalerSrc_;
                Float_t instLumi_;
                Float_t PULevel_;
                // Add UCT-only variables
                bool isUCT_;
};

Menu::Menu(const edm::ParameterSet& pset) {
        // Initialize the ntuple builder
        edm::Service<TFileService> fs;
        tree = fs->make<TTree>("Ntuple", "Ntuple");

        ptsEG_ = new std::vector<Float_t>();
        etasEG_ = new std::vector<Float_t>();
        phisEG_ = new std::vector<Float_t>();
        ptsJetEG_ = new std::vector<Float_t>();
        ptsRegionEG_ = new std::vector<Float_t>();
        MIPEG_ = new std::vector<Float_t>();
        tauVetoEG_ = new std::vector<Float_t>();
        LIsoEG_ = new std::vector<Float_t>();

        ptsIsoEG_ = new std::vector<Float_t>();
        etasIsoEG_ = new std::vector<Float_t>();
        phisIsoEG_ = new std::vector<Float_t>();
        ptsJetIsoEG_ = new std::vector<Float_t>();
        ptsRegionIsoEG_ = new std::vector<Float_t>();
        MIPIsoEG_ = new std::vector<Float_t>();
        tauVetoIsoEG_ = new std::vector<Float_t>();
        LIsoIsoEG_ = new std::vector<Float_t>();

        ptsJET_ = new std::vector<Float_t>();
        etasJET_ = new std::vector<Float_t>();
        phisJET_ = new std::vector<Float_t>();

        ptsTAU_ = new std::vector<Float_t>();
        etasTAU_ = new std::vector<Float_t>();
        phisTAU_ = new std::vector<Float_t>();
        ptsJetTAU_ = new std::vector<Float_t>();
        ptsRegionTAU_ = new std::vector<Float_t>();

        tree->Branch("ptEG", "std::vector<float>", &ptsEG_);
        tree->Branch("etaEG", "std::vector<float>", &etasEG_);
        tree->Branch("phiEG", "std::vector<float>", &phisEG_);
        tree->Branch("ptJetEG", "std::vector<float>", &ptsJetEG_);
        tree->Branch("ptRegionEG", "std::vector<float>", &ptsRegionEG_);
        tree->Branch("MIPEG", "std::vector<float>", &MIPEG_);
        tree->Branch("tauVetoEG", "std::vector<float>", &tauVetoEG_);
        tree->Branch("LIsoEG", "std::vector<float>", &LIsoEG_);

        tree->Branch("ptIsoEG", "std::vector<float>", &ptsIsoEG_);
        tree->Branch("etaIsoEG", "std::vector<float>", &etasIsoEG_);
        tree->Branch("phiIsoEG", "std::vector<float>", &phisIsoEG_);
        tree->Branch("ptJetIsoEG", "std::vector<float>", &ptsJetIsoEG_);
        tree->Branch("ptRegionIsoEG", "std::vector<float>", &ptsRegionIsoEG_);
        tree->Branch("MIPIsoEG", "std::vector<float>", &MIPIsoEG_);
        tree->Branch("tauVetoIsoEG", "std::vector<float>", &tauVetoIsoEG_);
        tree->Branch("LIsoIsoEG", "std::vector<float>", &LIsoIsoEG_);

        tree->Branch("ptJET", "std::vector<float>", &ptsJET_);
        tree->Branch("etaJET", "std::vector<float>", &etasJET_);
        tree->Branch("phiJET", "std::vector<float>", &phisJET_);

        tree->Branch("ptTAU", "std::vector<float>", &ptsTAU_);
        tree->Branch("etaTAU", "std::vector<float>", &etasTAU_);
        tree->Branch("phiTAU", "std::vector<float>", &phisTAU_);
        tree->Branch("ptJetTAU", "std::vector<float>", &ptsJetTAU_);
        tree->Branch("ptRegionTAU", "std::vector<float>", &ptsRegionTAU_);

        tree->Branch("MET", &MET_, "MET/F");
        tree->Branch("METPHI", &METPHI_, "METPHI/F");
        tree->Branch("MHT", &MHT_, "MHT/F");
        tree->Branch("MHTPHI", &MHTPHI_, "MHTPHI/F");
        tree->Branch("SET", &SET_, "SET/F");
        tree->Branch("SHT", &SHT_, "SHT/F");

        tree->Branch("run", &run_, "run/i");
        tree->Branch("lumi", &lumi_, "lumi/i");
        tree->Branch("evt", &event_, "evt/l");
        tree->Branch("instlumi", &instLumi_, "instlumi/F");
        tree->Branch("PULevel", &PULevel_, "PULevel/F");

        // UCT variables
        isUCT_ = pset.getParameter<bool>("isUCT");

        srcEG_ = pset.getParameter<VInputTag>("srcEG");
        srcIsoEG_ = pset.getParameter<VInputTag>("srcIsoEG");
        srcJET_ = pset.getParameter<VInputTag>("srcJET");
        srcTAU_ = pset.getParameter<VInputTag>("srcTAU");
        srcMET_ = pset.getParameter<edm::InputTag>("srcMET");
        srcMHT_ = pset.getParameter<edm::InputTag>("srcMHT");
        srcSET_ = pset.getParameter<edm::InputTag>("srcSET");
        srcSHT_ = pset.getParameter<edm::InputTag>("srcSHT");

        scalerSrc_ = pset.exists("scalerSrc") ?
                pset.getParameter<edm::InputTag>("scalerSrc") : edm::InputTag("scalersRawToDigi");
}

Menu::~Menu() {
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

        void getValue(const edm::Event& evt, const edm::InputTag& tag, Float_t& et, Float_t& phi) {
                edm::Handle<edm::View<reco::Candidate> > handle;
                evt.getByLabel(tag, handle);
                et = handle->at(0).pt();
                phi = handle->at(0).phi();
        }

        void getSumEtL1(const edm::Event& evt, const edm::InputTag& tag, Float_t& sumet,bool upgrade) {
                if(!upgrade) {
                        edm::Handle<l1extra::L1EtMissParticleCollection> handle;
                        evt.getByLabel(tag, handle);
                        sumet = handle->at(0).etTotal();
                } else{
                        edm::Handle<edm::View<reco::Candidate> > handle;
                        evt.getByLabel(tag, handle);
                        sumet = handle->at(0).pt();
                }
        }
}

void Menu::analyze(const edm::Event& evt, const edm::EventSetup& es) {

        // Get the objects.
        std::vector<const reco::Candidate*> EGobjects = getCollections(
                        evt, srcEG_);
        std::vector<const reco::Candidate*> IsoEGobjects = getCollections(
                        evt, srcIsoEG_);
        std::vector<const reco::Candidate*> JETobjects = getCollections(
                        evt, srcJET_);
        std::vector<const reco::Candidate*> TAUobjects = getCollections(
                        evt, srcTAU_);

        //std::sort(EGobjects.begin(), EGobjects.end(), CandPtSorter());


        // Clear previous event's objects
        ptsEG_->clear();
        etasEG_->clear();
        phisEG_->clear();
        ptsJetEG_->clear();
        ptsRegionEG_->clear();
        MIPEG_->clear();
        tauVetoEG_->clear();
        LIsoEG_->clear();

        ptsIsoEG_->clear();
        etasIsoEG_->clear();
        phisIsoEG_->clear();
        ptsJetIsoEG_->clear();
        ptsRegionIsoEG_->clear();
        MIPIsoEG_->clear();
        tauVetoIsoEG_->clear();
        LIsoIsoEG_->clear();

        ptsJET_->clear();
        etasJET_->clear();
        phisJET_->clear();

        ptsTAU_->clear();
        etasTAU_->clear();
        phisTAU_->clear();
        ptsJetTAU_->clear();
        ptsRegionTAU_->clear();

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

        for (size_t i = 0; i < EGobjects.size(); ++i) {
                ptsEG_->push_back(EGobjects[i]->pt());
                etasEG_->push_back(EGobjects[i]->eta());
                phisEG_->push_back(EGobjects[i]->phi());
                if (isUCT_) {
                        const UCTCandidate* uct = dynamic_cast<const UCTCandidate*>(EGobjects[i]);
                        if (!uct) {
                                throw cms::Exception("bad input")
                                        << "Can't convert input into UCT format!" << std::endl;
                        }
                        ptsJetEG_->push_back(uct->getFloat("associatedJetPt", -4));
                        ptsRegionEG_->push_back(uct->getFloat("associatedRegionEt", -4));
                        MIPEG_->push_back(uct->getInt("mipBit", -4));
                        tauVetoEG_->push_back(uct->getInt("tauVeto", -4));
                        LIsoEG_->push_back(uct->getInt("ellIsolation", -4));
                } 
        }

        for (size_t i = 0; i < IsoEGobjects.size(); ++i) {
                ptsIsoEG_->push_back(IsoEGobjects[i]->pt());
                etasIsoEG_->push_back(IsoEGobjects[i]->eta());
                phisIsoEG_->push_back(IsoEGobjects[i]->phi());
                if (isUCT_) {
                        const UCTCandidate* uct = dynamic_cast<const UCTCandidate*>(IsoEGobjects[i]);
                        if (!uct) {
                                throw cms::Exception("bad input")
                                        << "Can't convert input into UCT format!" << std::endl;
                        }
                        ptsJetIsoEG_->push_back(uct->getFloat("associatedJetPt", -4));
                        ptsRegionIsoEG_->push_back(uct->getFloat("associatedRegionEt", -4));
                        MIPIsoEG_->push_back(uct->getInt("mipBit", -4));
                        tauVetoIsoEG_->push_back(uct->getInt("tauVeto", -4));
                        LIsoIsoEG_->push_back(uct->getInt("ellIsolation", -4));
                }
        }

        for (size_t i = 0; i < TAUobjects.size(); ++i) {
                ptsTAU_->push_back(TAUobjects[i]->pt());
                etasTAU_->push_back(TAUobjects[i]->eta());
                phisTAU_->push_back(TAUobjects[i]->phi());
                if (isUCT_) {
                        const UCTCandidate* uct = dynamic_cast<const UCTCandidate*>(TAUobjects[i]);
                        if (!uct) {
                                throw cms::Exception("bad input")
                                        << "Can't convert input into UCT format!" << std::endl;
                        }
                        ptsJetTAU_->push_back(uct->getFloat("associatedJetPt", -4));
                        ptsRegionTAU_->push_back(uct->getFloat("associatedRegionEt", -4));
                }
        }

        for (size_t i = 0; i < JETobjects.size(); ++i) {
                ptsJET_->push_back(JETobjects[i]->pt());
                etasJET_->push_back(JETobjects[i]->eta());
                phisJET_->push_back(JETobjects[i]->phi());
        }

        MET_=-1, METPHI_=100, SET_=-1;
        getValue(evt,srcMET_, MET_, METPHI_);
        getSumEtL1(evt,srcSET_,SET_,isUCT_);
        MHT_=-1, MHTPHI_=100, SHT_=-1;
        getValue(evt,srcMHT_, MHT_, MHTPHI_);
        getSumEtL1(evt,srcSHT_,SHT_,isUCT_);







        tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Menu);
