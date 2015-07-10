/*
 * =====================================================================================
 *
 *       Filename:  MyGenMet.cc
 *
 *    Description:  Compute total SumET and SumHt from PF objects
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
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class MyGenMet : public edm::EDProducer {
  public:
    typedef std::vector<reco::LeafCandidate> LeafCandidateCollection;
    MyGenMet(const edm::ParameterSet& pset);
    virtual ~MyGenMet(){}
    void produce(edm::Event& evt, const edm::EventSetup& es);
  private:
    edm::InputTag src_;
};

MyGenMet::MyGenMet(const edm::ParameterSet& pset) {
  produces<LeafCandidateCollection>("mygenmet");
}
void MyGenMet::produce(edm::Event& evt, const edm::EventSetup& es) {

 double gen_met_px=0;
 double gen_met_py=0;

 edm::Handle< std::vector<reco::GenParticle> >pGenPart;
 if(!evt.getByLabel("genParticles", pGenPart)) {std::cout<<"NO GenParticles!!!!"<<std::endl; return;}

    for( size_t i = 0; i < pGenPart->size(); ++ i ) {
       const reco::GenParticle& genpart = (*pGenPart)[i];

          if(genpart.status()==1 && (fabs(genpart.pdgId())==12 || fabs(genpart.pdgId())==14 || fabs(genpart.pdgId())== 16) ){
                gen_met_px+=genpart.px();
                      gen_met_py+=genpart.py();
                         }
    }

  double MET=sqrt(gen_met_px*gen_met_px+gen_met_py*gen_met_py);
  double PHI = atan2(gen_met_py, gen_met_px) + 3.1415927;

  std::auto_ptr<LeafCandidateCollection> mygenmet(new LeafCandidateCollection);
  mygenmet->push_back(reco::LeafCandidate(0, reco::LeafCandidate::PolarLorentzVector(MET,0,PHI,0)));
  evt.put(mygenmet, "mygenmet");


}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MyGenMet);
