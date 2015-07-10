#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

class FilterWeNuFromW : public edm::EDProducer {

 public:
  FilterWeNuFromW (const edm::ParameterSet &);
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


FilterWeNuFromW::FilterWeNuFromW( const ParameterSet & cfg ) :
 GenParticleTag_(cfg.getUntrackedParameter<edm::InputTag> ("GenTag", edm::InputTag("genParticles"))),
 minPtThreshold_(cfg.getUntrackedParameter<double> ("MinPtThreshold",5)),
 maxEtaThreshold_(cfg.getUntrackedParameter<double> ("MaxEtaThreshold",5)),
 genLevelSelect_(cfg.getUntrackedParameter<int> ("GenLevelSelect",11)),
 genLevelStatus_(cfg.getUntrackedParameter<int> ("GenLevelStatus",1)),
 maxIsolation_(cfg.getUntrackedParameter<double> ("MaxIsolation",5)),
 isolationCone_(cfg.getUntrackedParameter<double> ("IsolationCone",0.4))
{
 produces<CandidateCollection>();
}

void FilterWeNuFromW::beginJob() {
 nall=0;
 nsel=0;	
 ntype=0;
 nselfin=0;
}

void FilterWeNuFromW::endJob() {
 cout<<"********************************************************************"<<endl;
 cout<<"GEN LEVEL FILTERING"<<endl<<endl;
 cout<<"Total Analyzed =   "<<nall<<endl;
 cout<<"GEN Selection  =   "<<nsel<<endl;
 cout<<"Type Selection  =   "<<ntype<<endl;
 cout<<"Type Selection By Event =   "<<nselfin<<endl;

 cout<<"********************************************************************"<<endl;
}

void FilterWeNuFromW::produce (Event & ev, const EventSetup &) {
 nall++;

 std::auto_ptr<CandidateCollection> cleanPart(new CandidateCollection);

  std::cout<<nall<<std::endl;


 edm::Handle< vector<reco::GenParticle> >pGenPart;
 if(ev.getByLabel(GenParticleTag_, pGenPart)){


  std::cout<<"Read!!!"<<std::endl;
  for( size_t i = 0; i < pGenPart->size(); ++ i ) {
   const reco::GenParticle& genpart = (*pGenPart)[i];
   if ( fabs(genpart.pdgId())!=24) continue;
   nsel++;
   std::cout<<nsel<<"  W!  "<<genpart.numberOfDaughters()<<std::endl;
   for (unsigned int j=0; j<genpart.numberOfDaughters(); j++){
    reco::Candidate *dau=(reco::Candidate*)genpart.daughter(j);
    std::cout<<"\t  j"<<"   "<<dau->pt()<<"  "<<dau->pdgId()<<"   "<<fabs(dau->eta())<<std::endl;
    if(fabs(dau->pdgId())!=genLevelSelect_) continue;
    if(dau->pt()<minPtThreshold_) continue;
    if(fabs(dau->eta())>maxEtaThreshold_) continue;
    std::cout<<"Passed"<<std::endl;
    ntype++;        
    cleanPart->push_back(dau);
    std::cout<<"Saved!"<<std::endl;  
   }
  }

 }
 ev.put(cleanPart); 
 std::cout<<"Next Event"<<std::endl;  
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FilterWeNuFromW);
