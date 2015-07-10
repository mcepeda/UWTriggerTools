#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

class FilterWeNu : public edm::EDProducer {

public:
  FilterWeNu (const edm::ParameterSet &);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
private:
  edm::InputTag GenParticleTag_;
  double minPtThreshold_;
  double maxEtaThreshold_;
  int genLevelSelect_;
  int genLevelStatus_;
  double maxIsolation_;
  double isolationCone_;

  double nall;
  double nsel;
  double ntype;
  double nselfin;

};
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <map>
#include <memory>

using namespace edm;
using namespace std;
using namespace reco;


FilterWeNu::FilterWeNu( const ParameterSet & cfg ) :
      GenParticleTag_(cfg.getUntrackedParameter<edm::InputTag> ("GenTag", edm::InputTag("genParticles"))),
      minPtThreshold_(cfg.getUntrackedParameter<double> ("MinPtThreshold",5)),
      maxEtaThreshold_(cfg.getUntrackedParameter<double> ("MaxEtaThreshold",5)),
      genLevelSelect_(cfg.getUntrackedParameter<int> ("GenLevelSelect",11)),
      genLevelStatus_(cfg.getUntrackedParameter<int> ("GenLevelStatus",1)),
      maxIsolation_(cfg.getUntrackedParameter<double> ("MaxIsolation",5)),
      isolationCone_(cfg.getUntrackedParameter<double> ("IsolationCone",0.4))
{
          produces<GenParticleCollection>();
}

void FilterWeNu::beginJob() {
      nall=0;
      nsel=0;	
      ntype=0;
      nselfin=0;
}

void FilterWeNu::endJob() {
     cout<<"********************************************************************"<<endl;
     cout<<"GEN LEVEL FILTERING"<<endl<<endl;
     cout<<"Total Analyzed =   "<<nall<<endl;
     cout<<"GEN Selection  =   "<<nsel<<endl;
     cout<<"Type Selection  =   "<<ntype<<endl;
     cout<<"Type Selection By Event =   "<<nselfin<<endl;

     cout<<"********************************************************************"<<endl;
}

void FilterWeNu::produce (Event & ev, const EventSetup &) {
  nall++;

  bool found=false;

  std::auto_ptr<GenParticleCollection> cleanPart(new GenParticleCollection);


  edm::Handle< vector<reco::GenParticle> >pGenPart;
  if(ev.getByLabel(GenParticleTag_, pGenPart)){

  for( size_t i = 0; i < pGenPart->size(); ++ i ) {
        const reco::GenParticle& genpart = (*pGenPart)[i];
        if ( fabs(genpart.pdgId())!=genLevelSelect_) continue;
        //        cout<<" \t"<<genpart.status()<<" "<<genpart.pt()<<"   "<<fabs(genpart.eta())<<"   "<<genpart.pdgId()<<endl;
        if ( genpart.status()!=genLevelStatus_) continue;
        if ( genpart.pt()<minPtThreshold_) continue;
        if ( fabs(genpart.eta())>maxEtaThreshold_) continue;

                bool tauDecay=false, bDecay=false, cDecay=false, mesonDecay=false;
                int type=0;
                int momPDGID=0;
                int momStatus=0;
                if(genpart.numberOfMothers()>0){
                reco::Candidate *mother=(reco::Candidate*)genpart.mother();

                momPDGID=mother->pdgId();
                momStatus=mother->status();
                if(abs(momPDGID)==genLevelSelect_ && mother->status()!=3)
                do{
                      reco::Candidate* motherB=(reco::Candidate*)mother->mother();
                      mother=motherB;
                      momPDGID=mother->pdgId();
                      momStatus=mother->status();  
                }
                while(abs(momPDGID)==genLevelSelect_ &&   mother->status()!=3);

                if(abs(momPDGID)==15) { tauDecay=true; type=3;}
                else if(int(abs(momPDGID)/1000)==5 || int(abs(momPDGID)/100)==5) { bDecay=true; type=5;}
                else if(int(abs(momPDGID)/1000)==4 || int(abs(momPDGID)/100)==4) {cDecay=true; type=4;}
                else if(int(abs(momPDGID)/100)==1 || int(abs(momPDGID)/100)==2 || int(abs(momPDGID)/100)==3) {mesonDecay=true; type=2;}
                else if(abs(momPDGID)==24 || abs(momPDGID)==23 || (abs(momPDGID)==genLevelSelect_&&momStatus==3)) { type=1;}
                else {
                                        cout<<"WARNING"<<endl;
                                        cout<<"This should have been prompt but it isnt!!! We should not be here"<<endl;
                                        cout<<"\t D: "<<genpart.pdgId()<<"   "<<genpart.status()<<endl;
                                        cout<<"\t M: "<<mother->pdgId()<<"  "<<mother->numberOfMothers()<<"   "<<mother->status()<<endl;
                        }
                }
        
        nsel++;

        if(type!=1) continue;

        ntype++;
        cleanPart->push_back(genpart);
        found=true;
  }
  }

  ev.put(cleanPart); 
  if (found) nselfin++;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FilterWeNu);
