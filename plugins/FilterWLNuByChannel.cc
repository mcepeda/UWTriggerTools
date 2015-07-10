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



class FilterWLNuByChannel : public edm::EDFilter {
 public:
  FilterWLNuByChannel(const edm::ParameterSet& pset);
  virtual ~FilterWLNuByChannel(){}
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

FilterWLNuByChannel::FilterWLNuByChannel(const edm::ParameterSet& pset) :
 nMuons_(pset.getUntrackedParameter<int> ("nMuons",-1)),
 nElectrons_(pset.getUntrackedParameter<int> ("nElectrons",-1)),
 nTaus_(pset.getUntrackedParameter<int> ("nTaus",-1))
{
}

void FilterWLNuByChannel::endJob() {
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

void FilterWLNuByChannel::beginJob() {
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

bool FilterWLNuByChannel::filter(edm::Event& evt, const edm::EventSetup& es) {
 nall++;

 bool found=true;
 int NMuons=0;
 int NElectrons=0;
 int NTaus=0;
 int NNeutrinos=0;



 edm::Handle< vector<reco::GenParticle> >pGenPart;
 if(evt.getByLabel("genParticles", pGenPart)){

  for( size_t i = 0; i < pGenPart->size(); ++ i ) {
   const reco::GenParticle& genpart = (*pGenPart)[i];
   if(fabs(genpart.pdgId())!=24 && genpart.pdgId()!=23)continue;
   bool filled=false;
   if(genpart.numberOfDaughters()==1) continue;
   if(abs(genpart.mother(0)->pdgId()/100)==4 || abs(genpart.mother(0)->pdgId()/1000)==4) continue;  // random extra Ws coming from the bs?
   if(abs(genpart.mother(0)->pdgId()/100)==5 || abs(genpart.mother(0)->pdgId()/1000)== 5) continue;  // random extra Ws coming from the bs?

 //  std::cout<<nall<<"  "<<i<<"  "<<genpart.pdgId()<<"   "<<genpart.numberOfDaughters()<<"  "<<genpart.mother(0)->pdgId()<<std::endl; 
   for(unsigned int j=0; j<genpart.numberOfDaughters()&&!filled; j++){
    const reco::Candidate* Wdaughter=genpart.daughter(j);
 //  std::cout<<" \t "<<j<<"  "<<Wdaughter->pdgId()<<"   "<<Wdaughter->numberOfDaughters()<<std::endl;
    if ( fabs(Wdaughter->pdgId())==13) NMuons++;
    else if ( fabs(Wdaughter->pdgId())==11) NElectrons++;
    else if ( fabs(Wdaughter->pdgId())==15) NTaus++;
    else if ( fabs(Wdaughter->pdgId())==12 || fabs(Wdaughter->pdgId())==14 || fabs(Wdaughter->pdgId())==16) NNeutrinos++;
    else if (fabs(Wdaughter->pdgId()) == 24 || genpart.pdgId()==23) {
     for (unsigned int l=0; l<Wdaughter->numberOfDaughters(); l++){ 
      const reco::Candidate* gWdaughter=Wdaughter->daughter(l); 
//   std::cout<<" \t\t "<<l<<"   "<<gWdaughter->pdgId()<<"   "<<gWdaughter->numberOfDaughters()<<std::endl;
      if ( fabs(gWdaughter->pdgId())==13) NMuons++; 
      else if ( fabs(gWdaughter->pdgId())==11) NElectrons++; 
      else if ( fabs(gWdaughter->pdgId())==15) NTaus++; 
      else if ( fabs(gWdaughter->pdgId())==12 || fabs(gWdaughter->pdgId())==14 || fabs(gWdaughter->pdgId())==16) NNeutrinos++;
     }
    }
   if(NElectrons>0 || NMuons>0 || NTaus>0) {filled=true;}
   }

  }
 }

// std::cout<<NMuons<<"   "<<NElectrons<<"   "<<NTaus<<"  "<<NNeutrinos<<std::endl;


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
DEFINE_FWK_MODULE(FilterWLNuByChannel);
