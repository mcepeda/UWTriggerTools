#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

class FilterGoodElectrons : public edm::EDProducer {

public:
  FilterGoodElectrons (const edm::ParameterSet &);
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


FilterGoodElectrons::FilterGoodElectrons( const ParameterSet & cfg ) :
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

void FilterGoodElectrons::beginJob() {
      nall=0;
      nsel=0;	

}

void FilterGoodElectrons::endJob() {
     cout<<"********************************************************************"<<endl;
     cout<<"GEN LEVEL FILTERING"<<endl<<endl;
     cout<<"Total Analyzed =   "<<nall<<endl;
     cout<<"GEN Selection  =   "<<nsel<<endl;
     cout<<"********************************************************************"<<endl;
}

void FilterGoodElectrons::produce (Event & ev, const EventSetup &) {
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









        cleanPart->push_back(genpart);
        found=true;
        // cout<<" \t"<<genpart.status()<<" "<<genpart.pt()<<"   "<<fabs(genpart.eta())<<"   "<<genpart.pdgId()<<endl;
  }
  }

  ev.put(cleanPart); 
  if (found) nsel++;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FilterGoodElectrons);
