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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloMipQuietRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


using namespace std;
using namespace edm;

class DebugRegions : public edm::EDAnalyzer {
 public:
  DebugRegions(const edm::ParameterSet&);
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
 private:
  //Note the physical definitions are here but not used in calculation
  double egPhysicalEt(const L1CaloEmCand& cand) const {
   return egLSB_*cand.rank();
  }

  //Note the physical definitions are here but not used in calculation
  double regionPhysicalEt(const L1CaloRegion& cand) const {
   return regionLSB_*cand.et();
  }

  //These are the definitions used in calculation below
  double egEt(const L1CaloEmCand& cand) const {
   return cand.rank();
  }

  //These are the definitions used in calculation below
  double regionEt(const L1CaloRegion& cand) const {
   return cand.et();
  }


  // Helper methods

  // ----------member data ---------------------------

  bool debug_;

  unsigned int puMult;
  bool puMultCorrect_;
  bool applyCalibration_;

  InputTag uctDigis_;
  InputTag ecalSrc_;
  InputTag hcalSrc_;

  double egLSB_;
  double regionLSB_;
  bool isMC_;

  TTree* tree;
  std::vector<Float_t>* pts2x1_;
  std::vector<Float_t>* pts4x4_;
  std::vector<Int_t>* rank2x1_;
  std::vector<Int_t>* rank4x4_;
  std::vector<Int_t>* etas4x4_;
  std::vector<Int_t>* phis4x4_;
  std::vector<Int_t>* tpgsECAL_;
  std::vector<Int_t>* tpgsSATECAL_;
  std::vector<Int_t>* sumTPECAL_;
  std::vector<Int_t>* tpgsHCAL_;
  std::vector<Int_t>* tpgsSATHCAL_;
  std::vector<Int_t>* sumTPHCAL_;

  std::vector<Int_t>* ptSATECALTP_;
  std::vector<Int_t>* etaSATECALTP_;
  std::vector<Int_t>* phiSATECALTP_;
  std::vector<Int_t>* ietaSATECALTP_;
  std::vector<Int_t>* iphiSATECALTP_;



  ULong64_t truepu_, pu_;
  ULong64_t event_;


  ULong64_t nall;
};


DebugRegions::DebugRegions(const edm::ParameterSet& iConfig) :
 debug_(iConfig.getUntrackedParameter<bool>("debug",false)),
 uctDigis_(iConfig.getUntrackedParameter<edm::InputTag>("uctDigisTag", edm::InputTag("uctDigis"))),
 ecalSrc_(iConfig.getUntrackedParameter<edm::InputTag>("ecalTag", edm::InputTag("ecalDigis:EcalTriggerPrimitives"))),
 hcalSrc_(iConfig.getUntrackedParameter<edm::InputTag>("hcalTag", edm::InputTag("hackHCALMIPs"))),
 egLSB_(iConfig.getParameter<double>("egammaLSB")),
 regionLSB_(iConfig.getParameter<double>("regionLSB")),
 isMC_(iConfig.getUntrackedParameter<bool>("isMC",true))
{
}

void DebugRegions::beginJob() {

 edm::Service<TFileService> fs;
 tree = fs->make<TTree>("Ntuple", "Ntuple");
 pts4x4_ = new std::vector<Float_t>();
 pts2x1_ = new std::vector<Float_t>();
 rank4x4_ = new std::vector<Int_t>();
 rank2x1_ = new std::vector<Int_t>();
 phis4x4_ = new std::vector<Int_t>();
 etas4x4_ = new std::vector<Int_t>();
 tpgsECAL_ = new std::vector<Int_t>();
 tpgsSATECAL_ = new std::vector<Int_t>();
 sumTPECAL_ = new std::vector<Int_t>();
 tpgsHCAL_ = new std::vector<Int_t>();
 tpgsSATHCAL_ = new std::vector<Int_t>();
 sumTPHCAL_ = new std::vector<Int_t>();

 ptSATECALTP_ = new std::vector<Int_t>();
 etaSATECALTP_ = new std::vector<Int_t>();
 phiSATECALTP_ = new std::vector<Int_t>();
 ietaSATECALTP_ = new std::vector<Int_t>();
 iphiSATECALTP_ = new std::vector<Int_t>();


 tree->Branch("pts2x1", "std::vector<float>", &pts2x1_);
 tree->Branch("pts4x4", "std::vector<float>", &pts4x4_);
 tree->Branch("rank2x1", "std::vector<int>", &rank2x1_);
 tree->Branch("rank4x4", "std::vector<int>", &rank4x4_);
 tree->Branch("etas4x4", "std::vector<int>", &etas4x4_);
 tree->Branch("phis4x4", "std::vector<int>", &phis4x4_);
 tree->Branch("tpgsSATECAL", "std::vector<int>", &tpgsSATECAL_);
 tree->Branch("tpgsECAL", "std::vector<int>", &tpgsECAL_);
 tree->Branch("sumTPECAL", "std::vector<int>", &sumTPECAL_);
 tree->Branch("tpgsSATHCAL", "std::vector<int>", &tpgsSATHCAL_);
 tree->Branch("tpgsHCAL", "std::vector<int>", &tpgsHCAL_);
 tree->Branch("sumTPHCAL", "std::vector<int>", &sumTPHCAL_);
 tree->Branch("ptSATECALTP", "std::vector<int>", &ptSATECALTP_);
 tree->Branch("etaSATECALTP", "std::vector<int>", &etaSATECALTP_);
 tree->Branch("phiSATECALTP", "std::vector<int>", &phiSATECALTP_);
 tree->Branch("ietaSATECALTP", "std::vector<int>", &ietaSATECALTP_);
 tree->Branch("iphiSATECALTP", "std::vector<int>", &iphiSATECALTP_);


 tree->Branch("truepu", &truepu_, "truepu/l");
 tree->Branch("pu", &pu_, "pu/l");

 tree->Branch("evt", &event_, "evt/l");
 tree->Branch("noEvt", &nall, "noEvt/l");

 nall=0;

}

void DebugRegions::endJob() {
 cout<<"********************************************************************"<<endl;
 cout<<"Total Analyzed = "<<nall<<endl;
 cout<<"********************************************************************"<<endl;




}

unsigned int  convertEta(int ieta){
 int testEta=-100;
 if(ieta>=1&&ieta<5){testEta=11;}
 if(ieta>=5&&ieta<9){testEta=12;}
 if(ieta>=9&&ieta<13){testEta=13;}
 if(ieta>=13&&ieta<17){testEta=14;}
 if(ieta>=17&&ieta<21){testEta=15;}
 if(ieta>=21&&ieta<25){testEta=16;}
 if(ieta>=25&&ieta<29){testEta=17;}
 if(ieta>=29&&ieta<30){testEta=18;}
 if(ieta>=30&&ieta<31){testEta=19;}
 if(ieta>=31&&ieta<32){testEta=20;}
 if(ieta>=32){testEta=21;}

 if(ieta<=-1&&ieta>-5){testEta=10;}
 if(ieta<=-5&&ieta>-9){testEta=9;}
 if(ieta<=-9&&ieta>-13){testEta=8;}
 if(ieta<=-13&&ieta>-17){testEta=7;}
 if(ieta<=-17&&ieta>-21){testEta=6;}
 if(ieta<=-21&&ieta>-25){testEta=5;}
 if(ieta<=-25&&ieta>-29){testEta=4;}
 if(ieta<=-29&&ieta>-30){testEta=3;}
 if(ieta<=-30&&ieta>-31){testEta=2;}
 if(ieta<=-31&&ieta>-32){testEta=1;}
 if(ieta<=-32){testEta=0;}

 return testEta;
}

unsigned int  convertPhi(int iphi){
 int testPhi=100;
 if(iphi==72||iphi==1||iphi==2||iphi==71) testPhi=0;
 else if(iphi<7) testPhi=1;
 else if(iphi<11) testPhi=2;
 else if(iphi<15) testPhi=3;
 else if(iphi<19) testPhi=4;
 else if(iphi<23) testPhi=5;
 else if(iphi<27) testPhi=6;
 else if(iphi<31) testPhi=7;
 else if(iphi<35) testPhi=8;
 else if(iphi<39) testPhi=9;
 else if(iphi<43) testPhi=10;
 else if(iphi<47) testPhi=11;
 else if(iphi<51) testPhi=12;
 else if(iphi<55) testPhi=13;
 else if(iphi<59) testPhi=14;
 else if(iphi<63) testPhi=15;
 else if(iphi<67) testPhi=16;
 else if(iphi<71) testPhi=17;

 return testPhi;
}


 void
DebugRegions::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 nall++;
 Handle<L1CaloRegionCollection> regions;
 Handle<L1CaloEmCollection> emCands;

 iEvent.getByLabel(uctDigis_, regions);
 iEvent.getByLabel(uctDigis_, emCands);

 edm::Handle<EcalTrigPrimDigiCollection> ecal;
 edm::Handle<HcalTrigPrimDigiCollection> hcal;

 iEvent.getByLabel(ecalSrc_, ecal);
 iEvent.getByLabel(hcalSrc_, hcal);

 edm::ESHandle<L1CaloHcalScale> hcalScale;
 iSetup.get<L1CaloHcalScaleRcd>().get(hcalScale);

 if(isMC_){
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
   int BX = PVI->getBunchCrossing();
   if(BX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
    truepu_ = PVI->getTrueNumInteractions();
    pu_ = PVI->getPU_NumInteractions();
   }
  }                
 }

 pts4x4_ ->clear();
 pts2x1_ ->clear(); 
 rank4x4_->clear(); 
 rank2x1_->clear();
 etas4x4_ ->clear();
 phis4x4_ ->clear();
 tpgsECAL_ ->clear();
 tpgsSATECAL_ ->clear();
 sumTPECAL_ ->clear();
 tpgsHCAL_ ->clear();
 tpgsSATHCAL_ ->clear();
 sumTPHCAL_ ->clear();
 ptSATECALTP_   ->clear();
 etaSATECALTP_  ->clear();
 phiSATECALTP_  ->clear();
 ietaSATECALTP_  ->clear();
 iphiSATECALTP_  ->clear();



 for(L1CaloRegionCollection::const_iterator region =
   regions->begin();
   region != regions->end(); region++){
  double regionET =  regionPhysicalEt(*region);
  //unsigned int regionEta = region->gctEta();
  // Only non-empty regions are corrected
  double energyECAL2x1=0;
  // Find associated 2x1 ECAL energy (EG are calibrated, we should not scale them up, it affects the isolation routines)
  // 2x1 regions have the MAX tower contained in the 4x4 region that its position points to.
  // This is to not break isolation.
  for(L1CaloEmCollection::const_iterator emCand =emCands->begin(); emCand != emCands->end(); emCand++){
   double et = egEt(*emCand);
   int countTPGs=0;
   int satTPGs=0;
   double sumTP=0, sumTP2=0;
   double countHCALTPGs=0, satHCALTPGs=0;
   double sumTPHCAL=0, sumTP2HCAL=0;

   if(et<10) continue;        

   if(emCand->regionId().iphi() == region->gctPhi() &&  emCand->regionId().ieta() == region->gctEta()) {
    energyECAL2x1=et;
    //                    if (energyECAL2x1>24) energyECAL2x1=24+(energyECAL2x1-24)*2;

    if(debug_){
     cout<<"--------------------------------"<<endl;
     cout<<"Position:  "<<emCand->regionId().ieta()<<"   "<<region->gctPhi()<<"    "<<emCand->regionId().rctEta()<<"  "<<emCand->regionId().rctPhi()<<endl;
     cout<<"2x1 --> RANK (output of RCT):"<<emCand->rank()<<"  --> PhysicalPt: "<<energyECAL2x1<<endl;
     cout<<"4x4 --> RANK (output of RCT):"<<region->et()<<"  --> PhysicalPt:"<<regionET<<endl;
    }

    for (size_t i = 0; i < hcal->size(); ++i) {

     int ieta = (*hcal)[i].id().ieta();
     int iphi = (*hcal)[i].id().iphi();

     unsigned int testEta=convertEta(ieta);
     unsigned int testPhi=convertPhi(iphi);

     if(region->gctPhi()!=testPhi || testEta!=region->gctEta()) continue;

     short absieta = std::abs((*hcal)[i].id().ieta());
     short zside = (*hcal)[i].id().zside();

     double energy = hcalScale->et(
       (*hcal)[i].SOI_compressedEt(), absieta, zside);

     countHCALTPGs++;
     if(energy>=127) satHCALTPGs++;

     sumTPHCAL+=energy;
     if(sumTPHCAL<127)  sumTP2+=energy;

     if (debug_) {
      std::cout << "hcal eta/phi=" << ieta << "/" << iphi
       << " et=" << (*hcal)[i].SOI_compressedEt()
       << " energy=" << energy
       << " fg=" << (*hcal)[i].SOI_fineGrain() << std::endl;
     }
    }


    for (size_t i = 0; i < ecal->size(); ++i) {
     int ieta = (*ecal)[i].id().ieta();
     int iphi = (*ecal)[i].id().iphi();


     unsigned int testEta=convertEta(ieta);
     unsigned int testPhi=convertPhi(iphi);

     if(region->gctPhi()!=testPhi || testEta!=region->gctEta()) continue;

     countTPGs++;
     if((*ecal)[i].compressedEt()==255) {
        satTPGs++;

        ptSATECALTP_->push_back((*ecal)[i].compressedEt());       
        ietaSATECALTP_->push_back(ieta);
        iphiSATECALTP_->push_back(iphi);
        etaSATECALTP_->push_back(testEta);
        phiSATECALTP_->push_back(testPhi);
     }

     sumTP+=(*ecal)[i].compressedEt();        
     if(debug_  && (*ecal)[i].compressedEt()>0){
      std::cout << "ecal eta/phi=" << ieta << "/" << iphi
       << " = (" << testEta <<"/" << testPhi << ") "
       << " et="<< (*ecal)[i].compressedEt() << " fg=" << (*ecal)[i].fineGrain()
       << std::endl;
     }
    }

    pts4x4_ ->push_back(regionET);
    pts2x1_ ->push_back(energyECAL2x1);
    phis4x4_ ->push_back(region->gctPhi());
    etas4x4_ ->push_back(region->gctEta());
    rank4x4_->push_back(region->et());
    rank2x1_->push_back(emCand->rank());
    tpgsECAL_ ->push_back(countTPGs);
    tpgsSATECAL_ ->push_back(satTPGs);
    tpgsHCAL_ ->push_back(countHCALTPGs);
    tpgsSATHCAL_ ->push_back(satHCALTPGs);
    sumTPECAL_ ->push_back(sumTP);
    sumTPHCAL_ ->push_back(sumTPHCAL);

    break;
   }
  }

 }


 tree->Fill();                 

 if(debug_)        cout<<"==========================================="<<endl;

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DebugRegions);
