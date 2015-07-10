#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

class FilterTT : public edm::EDFilter {

public:
  FilterTT (const edm::ParameterSet &);
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
private:
  int nMuons_;
  int nElectrons_;

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


FilterTT::FilterTT( const ParameterSet & cfg ) :
      nMuons_(cfg.getUntrackedParameter<int> ("nMuons",-1)),
      nElectrons_(cfg.getUntrackedParameter<int> ("nElectrons",-1))
{
}

void FilterTT::beginJob() {
      nall=0;
      nsel=0;	

}

void FilterTT::endJob() {
     cout<<"********************************************************************"<<endl;
     cout<<"GEN LEVEL FILTERING"<<endl<<endl;
     cout<<"Total Analyzed =   "<<nall<<endl;
     cout<<"GEN Selection  =   "<<nsel<<endl;
     cout<<"********************************************************************"<<endl;
}

bool FilterTT::filter (Event & ev, const EventSetup &) {
  nall++;

  bool found=true;
  int NMuons=0;
  int NElectrons=0;

  std::auto_ptr<GenParticleCollection> cleanPart(new GenParticleCollection);


  edm::Handle< vector<reco::GenParticle> >pGenPart;
  if(ev.getByLabel("genParticles", pGenPart)){


  for( size_t i = 0; i < pGenPart->size(); ++ i ) {
        const reco::GenParticle& genpart = (*pGenPart)[i];
        if ( genpart.status()!=3) continue;
        std::cout<<genpart.status()<<"   "<<genpart.pdgId()<<std::endl;
        if ( fabs(genpart.pdgId())==13) NMuons++;
        if ( fabs(genpart.pdgId())==11) NElectrons++;
  }

  }
  if(nMuons_!=-1  && nMuons_!=NMuons) found=false;
  if(nElectrons_!=-1  && nElectrons_!=NElectrons) found=false;

  if (found) nsel++;



  return found;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FilterTT);
