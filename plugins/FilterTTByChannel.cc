#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <map>
#include <memory>

using namespace edm;
using namespace std;
using namespace reco;



class FilterTTByChannel : public edm::EDFilter {
 public:
  FilterTTByChannel(const edm::ParameterSet& pset);
  virtual ~FilterTTByChannel(){}
  bool filter(edm::Event& evt, const edm::EventSetup& es);
  virtual void beginJob();
  virtual void endJob();
 private:
  int nMuons_;
  int nElectrons_;
  int nTaus_;

  double nall;
  double nsel;

  double events2muon;
  double events2electron;
  double events2tau;
  double eventsmuonelectron;
  double eventsmuontau;
  double eventselectrontau;
  double events1muon;
  double events1electron;
  double events1tau;
  double eventsnolepton;


};

FilterTTByChannel::FilterTTByChannel(const edm::ParameterSet& pset) :
 nMuons_(pset.getUntrackedParameter<int> ("nMuons",-1)),
 nElectrons_(pset.getUntrackedParameter<int> ("nElectrons",-1)),
 nTaus_(pset.getUntrackedParameter<int> ("nTaus",-1))
{
}

void FilterTTByChannel::endJob() {
 cout<<"********************************************************************"<<endl;
 cout<<"GEN LEVEL FILTERING"<<endl<<endl;
 cout<<"Total Analyzed =   "<<nall<<endl;
 cout<<"GEN Selection  =   "<<nsel<<endl;
 cout<<"--------------------------------------------------------------------"<<endl;
 cout<<"    2Muons            : "<<events2muon<<"   "<<events2muon/nall<<endl;
 cout<<"    2Electrons        : "<<events2electron<<"   "<<events2electron/nall<<endl;
 cout<<"    2Taus             : "<<events2tau<<"   "<<events2tau/nall<<endl;
 cout<<"    1Muon, 1 Electron : "<<eventsmuonelectron<<"   "<<eventsmuonelectron/nall<<endl;
 cout<<"    1Muon, 1 Tau      : "<<eventsmuontau<<"   "<<eventsmuontau/nall<<endl;
 cout<<"    1Electron, 1 Tau  : "<<eventselectrontau<<"   "<<eventselectrontau/nall<<endl;
 cout<<"    1Muons            : "<<events1muon<<"   "<<events1muon/nall<<endl;
 cout<<"    1Electrons        : "<<events1electron<<"   "<<events1electron/nall<<endl;
 cout<<"    1Taus             : "<<events1tau<<"   "<<events1tau/nall<<endl;
 cout<<"    NoLeptons         : "<<eventsnolepton<<"   "<<eventsnolepton/nall<<endl;    
 cout<<"********************************************************************"<<endl;
}

void FilterTTByChannel::beginJob() {
 nall=0;
 nsel=0;
 events2muon=0;
 events2electron=0;
 events2tau=0;
 eventsmuonelectron=0;
 eventsmuontau=0;
 eventselectrontau=0;
 events1muon=0;
 events1electron=0;
 events1tau=0;
 eventsnolepton=0;  
}

bool FilterTTByChannel::filter(edm::Event& evt, const edm::EventSetup& es) {
 nall++;

 bool found=true;
 int NMuons=0;
 int NElectrons=0;
 int NTaus=0;

 double gen_met_px=0;
 double gen_met_py=0;

 std::auto_ptr<GenParticleCollection> cleanPart(new GenParticleCollection);


 edm::Handle< vector<reco::GenParticle> >pGenPart;
 if(evt.getByLabel("genParticles", pGenPart)){

//  std::cout<<nall<<std::endl;


  for( size_t i = 0; i < pGenPart->size(); ++ i ) {
   const reco::GenParticle& genpart = (*pGenPart)[i];

   if(genpart.status()==1 && (fabs(genpart.pdgId())==12 || fabs(genpart.pdgId())==14 || fabs(genpart.pdgId())== 16) ){
      gen_met_px+=genpart.px();
      gen_met_py+=genpart.py();
   }

   if(fabs(genpart.pdgId())==6){
    //        bool filled=false; 
  //  std::cout<<i<<"   "<<genpart.numberOfDaughters()<<"   "<<genpart.pdgId()<<"   "<<genpart.pt()<<std::endl;
    for(unsigned int j=0; j<genpart.numberOfDaughters(); j++){
     const reco::Candidate* topdaughter=genpart.daughter(j);
    // std::cout<<"\t"<<j<<"   "<<genpart.numberOfDaughters()<<"   "<<fabs(topdaughter->pdgId())<<std::endl;
     if(fabs(topdaughter->pdgId()) != 24) continue;
     for (unsigned int k=0; k<topdaughter->numberOfDaughters(); k++){
      //filled=true;
      const reco::Candidate* grandaughter=topdaughter->daughter(k);
    //  std::cout<<"\t\t"<<k<<"   "<<fabs(grandaughter->pdgId())<<std::endl; 
      if ( fabs(grandaughter->pdgId())==13) NMuons++;
      if ( fabs(grandaughter->pdgId())==11) NElectrons++;
      if ( fabs(grandaughter->pdgId())==15) NTaus++;
      if (fabs(grandaughter->pdgId()) == 24) {

       for (unsigned int l=0; l<grandaughter->numberOfDaughters(); l++){ 
        const reco::Candidate* ggrandaughter=grandaughter->daughter(l); 
    //    std::cout<<"\t\t\t"<<l<<"   "<<fabs(ggrandaughter->pdgId())<<std::endl;  
        if ( fabs(ggrandaughter->pdgId())==13) NMuons++; 
        if ( fabs(ggrandaughter->pdgId())==11) NElectrons++; 
        if ( fabs(ggrandaughter->pdgId())==15) NTaus++; 
       }

      }


     }
    }

   }
  }
 }

// std::cout<<NMuons<<"   "<<NElectrons<<"   "<<NTaus<<std::endl;


 if(nMuons_!=-1  && nMuons_!=NMuons) found=false;
 if(nElectrons_!=-1  && nElectrons_!=NElectrons) found=false;
 if(nTaus_!=-1  && nTaus_!=NTaus) found=false;

 if (found) nsel++;

 if(NMuons==2) events2muon++;
 else if(NElectrons==2) events2electron++;
 else if(NTaus==2) events2tau++;
 else if(NMuons==1 && NElectrons==1) eventsmuonelectron++;
 else if(NMuons==1 && NTaus==1) eventsmuontau++;
 else if(NTaus==1 && NElectrons==1) eventselectrontau++;
 else if(NMuons==1) events1muon++;
 else if(NElectrons==1) events1electron++;
 else if(NTaus==1) events1tau++;
 else eventsnolepton++;


 return found;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FilterTTByChannel);
