#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"

class EWKLeptonProducer : public edm::EDProducer {
    public:

        EWKLeptonProducer(const edm::ParameterSet& pset);
        virtual ~EWKLeptonProducer(){}
        void produce(edm::Event& evt, const edm::EventSetup& es);
    private:
        edm::InputTag src_;
        int genPDGID_;
        int genSTATUS_;	
        double minPT_;
        double maxETA_;
};

using namespace edm;
using namespace std;
using namespace reco;

EWKLeptonProducer::EWKLeptonProducer(
        const edm::ParameterSet& pset):
        src_(pset.getParameter<edm::InputTag>("src")),
        genPDGID_(pset.getUntrackedParameter<int>("genPDGID",0)),
        genSTATUS_(pset.getUntrackedParameter<int>("genSTATUS",1)),
        minPT_(pset.getUntrackedParameter<double>("minPT",0)),
        maxETA_(pset.getUntrackedParameter<double>("maxETA",100))

            {
                produces<GenParticleCollection>();
            }

void EWKLeptonProducer::produce(edm::Event& evt, const edm::EventSetup& es) {

   Handle<GenParticleCollection> GPCollection;
   evt.getByLabel("genParticles", GPCollection);
   const GenParticleCollection & pGenPart = *(GPCollection.product());

   auto_ptr<GenParticleCollection> out(new GenParticleCollection);


   //cout<<"=======================  EVENT"<<endl;

          for( size_t i = 0; i < pGenPart.size(); ++ i ) {

                GenParticle genpart(pGenPart[i]);

                if(genPDGID_!=0 && fabs(genpart.pdgId())!=genPDGID_) continue;
                else if (fabs(genpart.pdgId()) != 11 && fabs(genpart.pdgId())!=13) continue;
                if(genpart.status()!=genSTATUS_)continue;

                int particleID=genpart.pdgId();
                if(genpart.pt()<minPT_ || fabs(genpart.eta())>maxETA_) continue;

//                std::cout<<"Status 1 Electron!  "<<genpart.pt()<<"     "<<genpart.eta()<<"     "<<genpart.phi()<<std::endl;


                bool tauDecay=false, bDecay=false, cDecay=false, mesonDecay=false;
                int type=0;
                int momPDGID=0;
                if(genpart.numberOfMothers()>0){
                reco::Candidate *mother=(reco::Candidate*)genpart.mother();

                momPDGID=mother->pdgId();


                if(momPDGID==particleID && mother->numberOfMothers()!=0) 
                do{
                      reco::Candidate* motherB=(reco::Candidate*)mother->mother();  
                      mother=motherB;  
                      momPDGID=mother->pdgId();
               
                }
                while(momPDGID==particleID &&   mother->numberOfMothers()!=0);
               
                if(abs(momPDGID)==15) { tauDecay=true; type=6;}
                else if(int(abs(momPDGID)/1000)==5 || int(abs(momPDGID)/100)==5) { bDecay=true; type=5;}
                else if(int(abs(momPDGID)/1000)==4 || int(abs(momPDGID)/100)==4) {cDecay=true; type=4;}
                else if(int(abs(momPDGID)/1000)==3 || int(abs(momPDGID)/1000)==2 || int(abs(momPDGID)/1000)==1  ) { type=3;}
                else if(int(abs(momPDGID)/100)==1 || int(abs(momPDGID)/100)==2 || int(abs(momPDGID)/100)==3) {mesonDecay=true; type=2;}
                else { type=1; 
                       if(abs(momPDGID)!=24 && abs(momPDGID)!=23 && abs(momPDGID)!=22){
                                        cout<<"WARNING"<<endl;
                                        cout<<"This should have been prompt but it isnt!!! We should not be here"<<endl; 
                                        cout<<"\t D: "<<genpart.pdgId()<<"   "<<genpart.status()<<endl;
                                        cout<<"\t M: "<<mother->pdgId()<<"  "<<mother->numberOfMothers()<<"   "<<mother->status()<<endl;
                        }

                }
/*
                cout<<"\t D: "<<genpart.pdgId()<<"   "<<genpart.status()<<endl;
                cout<<"\t M: "<<mother->pdgId()<<"  "<<mother->numberOfMothers()<<"   "<<mother->status()<<endl;
                cout<<"      -->"<<type<<endl;
*/


                }

                if(type!=1) continue; 
                std::cout<<"Status 1 Electron!  "<<genpart.pt()<<"     "<<genpart.eta()<<"     "<<genpart.phi()<<std::endl;
                out->push_back(genpart);
        }

   evt.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EWKLeptonProducer);
