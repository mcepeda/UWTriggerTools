#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

class ECALMasker : public edm::EDProducer {
    public:

        ECALMasker(const edm::ParameterSet& pset);
        virtual ~ECALMasker(){}
        void produce(edm::Event& evt, const edm::EventSetup& es);
    private:
        edm::InputTag src_;
        std::vector<int> iEtas_;
        std::vector<int> iPhis_;
        int minRank_;
};

using namespace edm;
using namespace std;

ECALMasker::ECALMasker(
        const edm::ParameterSet& pset):
        src_(pset.getParameter<InputTag>("src")),
        iEtas_(pset.getParameter<vector<int>>("iEtasMasked")),
        iPhis_(pset.getParameter<vector<int>>("iPhisMasked")),
        minRank_(pset.getUntrackedParameter<int>("minRankForMasking",255))
            {
                if(iEtas_.size()!=iPhis_.size()) {
                  throw cms::Exception("ECALMaskError")
                     << "iEta and iPhi vectors should have the same length\n"<<endl;
                }

                produces<EcalTrigPrimDigiCollection>();
            }

void ECALMasker::produce(edm::Event& evt, const edm::EventSetup& es) {

   edm::Handle<EcalTrigPrimDigiCollection> ecalDigis;
   evt.getByLabel(src_, ecalDigis);
   const EcalTrigPrimDigiCollection & pGenPart = *(ecalDigis.product());

   auto_ptr<EcalTrigPrimDigiCollection> out(new EcalTrigPrimDigiCollection);


   for (size_t i = 0; i < ecalDigis->size(); ++i) {
     int ieta = (*ecalDigis)[i].id().ieta();
     int iphi = (*ecalDigis)[i].id().iphi();

     bool skip=false;
   
     if((*ecalDigis)[i].compressedEt()>=minRank_){
      for (size_t j=0; j<iEtas_.size(); j++){
          if(ieta==iEtas_.at(j) && iphi==iPhis_.at(j) ) skip=true;
      }
     } 
     if (skip) continue;
     
     out->push_back((*ecalDigis)[i]);
    
   }

   evt.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ECALMasker);
