// system include files
#include <memory>
#include <math.h>
#include <vector>
#include <list>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

using namespace std;
using namespace edm;



class FilterSaturated : public edm::EDFilter {
 public:
  FilterSaturated(const edm::ParameterSet&);
  virtual void beginJob();
  virtual void endJob();
  bool filter(edm::Event&, const edm::EventSetup&);
 private:

  // ----------member data ---------------------------

  bool debug_;
  int  filterType_;

  InputTag ecalSrc_;
  InputTag hcalSrc_;

  Long_t nall, npass;

};


FilterSaturated::FilterSaturated(const edm::ParameterSet& iConfig) :
 debug_(iConfig.getUntrackedParameter<bool>("debug",false)),
 filterType_(iConfig.getUntrackedParameter<int>("filterType",1)),
 ecalSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ecalTag", edm::InputTag("ecalDigis:EcalTriggerPrimitives"))),
 hcalSrc_(iConfig.getUntrackedParameter<edm::InputTag>("hcalTag", edm::InputTag("hackHCALMIPs")))
{
}

void FilterSaturated::beginJob() {
 nall=0;
 npass=0;
}

void FilterSaturated::endJob() {
 cout<<"********************************************************************"<<endl;
 cout<<"Total Analyzed = "<<nall<<endl;
 cout<<"Passed = "<<npass<<endl;
 cout<<"********************************************************************"<<endl;
}

bool
FilterSaturated::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 nall++;
 edm::Handle<EcalTrigPrimDigiCollection> ecal;
 edm::Handle<HcalTrigPrimDigiCollection> hcal;
 iEvent.getByLabel(ecalSrc_, ecal);
 iEvent.getByLabel(hcalSrc_, hcal);
 edm::ESHandle<L1CaloHcalScale> hcalScale;
 iSetup.get<L1CaloHcalScaleRcd>().get(hcalScale);

   double satECALTPGs=0, satHCALTPGs=0;

    for (size_t i = 0; i < hcal->size(); ++i) {

     int ieta = (*hcal)[i].id().ieta();
     int iphi = (*hcal)[i].id().iphi();

     short absieta = std::abs((*hcal)[i].id().ieta());
     short zside = (*hcal)[i].id().zside();

     double energy = hcalScale->et(
       (*hcal)[i].SOI_compressedEt(), absieta, zside);

     if((*hcal)[i].SOI_compressedEt()==255) satHCALTPGs++;

     if (debug_ && (*hcal)[i].SOI_compressedEt() > 0 ) {
      std::cout << "hcal eta/phi=" << ieta << "/" << iphi
       << " et=" << (*hcal)[i].SOI_compressedEt()
       << " energy=" << energy<<std::endl;
     }
    }


    for (size_t i = 0; i < ecal->size(); ++i) {
     int ieta = (*ecal)[i].id().ieta();
     int iphi = (*ecal)[i].id().iphi();

     if((*ecal)[i].compressedEt()==255) satECALTPGs++;

     if(debug_  && (*ecal)[i].compressedEt()>0){
      std::cout << "ecal eta/phi=" << ieta << "/" << iphi
       << " et="<< (*ecal)[i].compressedEt() << " fg=" << (*ecal)[i].fineGrain()
       << std::endl;
     }
    }

 if(debug_) {
       std::cout<<" ECAL: "<<satECALTPGs<<"   HCAL: "<<satHCALTPGs<<std::endl;
       cout<<"==========================================="<<endl;
 }

 if(filterType_==0)  return true;
 if (filterType_==1) return (satECALTPGs==0); 
 if (filterType_==2) return (satHCALTPGs==0);
 if (filterType_==3) return (satECALTPGs==0 && satHCALTPGs==0);
 if (filterType_==4) return (satECALTPGs>0 || satHCALTPGs>0);

 return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FilterSaturated);
