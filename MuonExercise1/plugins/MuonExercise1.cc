// -*- C++ -*-
//
// Package:    CMSDASExercises/MuonExercise1
// Class:      MuonExercise1
// 
/**\class MuonExercise1 MuonExercise1.cc CMSDASExercises/MuonExercise1/plugins/MuonExercise1.cc

 Description: Short Muon exercise for CMSDAS

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Norbert Neumeister
//         Created:  Thu, 14 Dec 2017 09:31:13 GMT
//
//

// system include files
#include <memory>
#include <iomanip>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "Math/GenVector/VectorUtil.h"
#include "TRandom3.h"
#include "RoccoR.cc"

//
// class declaration
//

class MuonExercise1 : public edm::one::EDAnalyzer<edm::one::SharedResources> {

  public:

    explicit MuonExercise1(const edm::ParameterSet&);
    ~MuonExercise1();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:

    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
     
    // TFileService
    edm::Service<TFileService> fs;
  
    // miniAOD Collections
    edm::EDGetTokenT<std::vector<pat::Muon> > muonCollToken;
    //edm::EDGetTokenT<reco::GenParticleCollection> genCollToken;
    edm::EDGetTokenT<pat::PackedGenParticleCollection> genCollToken;

    // Histograms
    TH1D* h_ngen;     // number of generated muons
    TH1D* h_nrec;     // number of reconstructed muons
    TH1D* h_genpt;    // pt of generated muons 
    TH1D* h_recpt;    // pt of reconstructed muons
    TH1D* h_iso;      // isolation score for reconstructed muons     
    TH1D* h_geneta;   // eta value for generated muons
    TH1D* h_receta;   // eta value for reconstructed muons
    TH1D* h_genphi;    // phi value for generated muons
    TH1D* h_recphi;    // phi value for reconstructed muons
    TH1F* h_genmass;   // gen dimuon Z mass
    TH1F* h_recmass;   // reco dimuon Z mass
    TH1F* h_recmasscorr; // reco dimuon Z mass with rochester corrections

};

    RoccoR rc("/root://cmseos.fnal.gov///store/user/hats/2018/Muon/Exercises/MuonExercise2/plugin/rcdata.2016.v3");
    TRandom3 rnd(1234);


//
// constructors and destructor

 MuonExercise1::MuonExercise1(const edm::ParameterSet& iConfig) {

  usesResource("TFileService");
 
  edm::InputTag muonTag("slimmedMuons");
  edm::InputTag genPartTag("packedGenParticles");
  //edm::InputTag genPartTag("prunedGenParticles");
 
  muonCollToken = consumes<pat::MuonCollection>(muonTag);
  genCollToken = consumes<pat::PackedGenParticleCollection>(genPartTag);
  //genCollToken = consumes<reco::GenParticleCollection>(genPartTag);

  h_ngen  = fs->make<TH1D>("ngen", "Number of GEN muons", 10, 0.0, 10.0);
  h_nrec  = fs->make<TH1D>("nrec", "Number of RECO muons", 10, 0.0, 10.0);

  h_recpt = fs->make<TH1D>("recpt", "RECO pt", 100, 0.0, 200.0);
  h_genpt = fs->make<TH1D>("genpt", "GEN pt", 100, 0.0, 200.0);

  h_geneta = fs->make<TH1D>("geneta", "GEN eta", 20, -2.4, 2.4);
  h_receta = fs->make<TH1D>("receta", "RECO eta", 20, -2.4, 2.4);

  h_genphi = fs->make<TH1D>("genphi", "GEN phi", 20, -M_PI, M_PI);
  h_recphi = fs->make<TH1D>("recphi", "RECO phi", 20, -M_PI, M_PI); 

  h_genmass = fs->make<TH1F>("genmass", ";m_{#mu^{+}#mu^{-}};", 80, 70, 110);
  h_recmass = fs->make<TH1F>("recmass", ";m_{#mu^{+}#mu^{-}};", 80, 70, 110); 
  h_recmasscorr = fs->make<TH1F>("recmasscorr", ";m_{#mu^{+}#mu^{-}};", 80, 70, 110);
  
  h_iso = fs->make<TH1D>("iso", "Isolation Score", 100, 0.0, 2);

};

MuonExercise1::~MuonExercise1() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void MuonExercise1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  using namespace l1extra;


  double ngen(0.);
  double nrec(0.);

  ////////////////////////////////////////////
  //
  // RECO Muons
  //  
  ////////////////////////////////////////////
  
  edm::Handle <vector<pat::Muon>> muonColl;
  iEvent.getByToken(muonCollToken, muonColl);
  //nrec = muonColl->size();

  for (auto mu = muonColl->cbegin(); mu != muonColl->cend(); ++mu) {
    if ( fabs((*mu).eta()) < 2.4 && (*mu).pt() > 25 && (*mu).isMediumMuon()) {

	nrec++;	
       
	// PF Isolation score --------
	double chargedHadronIso = (*mu).pfIsolationR04().sumChargedHadronPt;
        double neutralHadronIso  = (*mu).pfIsolationR04().sumNeutralHadronEt;
        double photonIso  = (*mu).pfIsolationR04().sumPhotonEt;
        double chargedHadronIsoPU = (*mu).pfIsolationR04().sumPUPt;
        float relativeIsolationDBetaCorr = ( chargedHadronIso + std::max(0., neutralHadronIso + photonIso - 0.5*chargedHadronIsoPU) )/(*mu).pt();

	// Histograms
	//h_iso->Fill(relativeIsolationDBetaCorr);	

	if (relativeIsolationDBetaCorr < 0.15) {
		h_nrec->Fill(nrec);
		h_recpt->Fill((*mu).pt());
		h_receta->Fill((*mu).eta());
		h_recphi->Fill((*mu).phi());	
		h_iso->Fill(relativeIsolationDBetaCorr);
        }
     }
  }	
 
  
  cout << "Number of RECO muons: " << nrec << endl;

  ////////////////////////////////////////////
  //
  // GEN Muons
  //  
  ////////////////////////////////////////////

  edm::Handle <pat::PackedGenParticleCollection> genColl;
  //edm::Handle <reco::GenParticleCollection> genColl;
  iEvent.getByToken(genCollToken, genColl);
  //ngen = genColl->size();
  for (auto mu = genColl->cbegin(); mu != genColl->cend(); ++mu) {
    if ( abs((*mu).pdgId()) == 13 && fabs((*mu).eta()) < 2.4 && (*mu).pt() > 25 ) {
	
	ngen++;
	
	// Histograms
	h_ngen->Fill(ngen);
	h_genpt->Fill((*mu).pt());
	h_geneta->Fill((*mu).eta());
	h_genphi->Fill((*mu).phi());
    }
  }

  cout << "Number of GEN muons: " << ngen << endl;

  // for (auto it = genColl->cbegin(); it != genColl->cend(); ++it) {

    // put your code here
    
  //}




// Z Mass
// // find mu+
  for (auto mup = muonColl->cbegin(); mup != muonColl->cend(); ++mup) {
     if ( not (mup->charge() > 0 ) ) continue;
  for (auto mup = muonColl->cbegin(); mup != muonColl->cend(); ++mup) {
     if ( not (mup->charge() > 0 ) ) continue;
     if ( not (mup->isGlobalMuon()) ) continue;
     if ( not (mup->pt() > 20.0 ) ) continue;
     if ( fabs(mup->eta()) > 2.4 ) continue;
     if ( not (mup->chargedHadronIso() < 0.15 ) ) continue;

     std::cout << "Stops after line " << __LINE__ << std::endl;

     double mupSF = rc.kScaleFromGenMC(mup->charge(),
				       mup->pt(),
				       mup->eta(),
				       mup->phi(),
				       mup->innerTrack()->hitPattern().trackerLayersWithMeasurement(),
				       mup->genParticle()->pt(),
				       rnd.Rndm());

     std::cout << "Stops after line " << __LINE__ << std::endl;

     // find mu-
     for (auto mum = muonColl->cbegin(); mum != muonColl->cend(); ++mum) {
        if ( not (mum->charge() < 0 ) ) continue;
        if ( not (mum->isGlobalMuon()) ) continue;
        if ( not (mum->pt() > 20.0 ) ) continue;
        if ( fabs(mum->eta()) > 2.4 ) continue;
        if ( not (mum->chargedHadronIso() < 0.15 ) ) continue;

	std::cout << "Stops after line " << __LINE__ << std::endl;

        double mumSF = rc.kScaleFromGenMC(mum->charge(),
				          mum->pt(),
				          mum->eta(),
				          mum->phi(),
				          mum->innerTrack()->hitPattern().trackerLayersWithMeasurement(),
				          mum->genParticle()->pt(),
				          rnd.Rndm());

	std::cout << "Stops after line " << __LINE__ << std::endl;

        double diMuonRecMass = (mup->p4() + mum->p4()).M();
        if ( diMuonRecMass < 70 || diMuonRecMass > 110) continue; // only look around the Z peak
        h_recmass->Fill(diMuonRecMass);

	std::cout << "Stops after line " << __LINE__ << std::endl;

	auto mupp4 = mup->p4() * mupSF; 
	auto mump4 = mum->p4() * mumSF;

	std::cout << "Stops after line " << __LINE__ << std::endl;

	double diMuonRecMassCorr = (mupp4 + mump4).M();
	if ( diMuonRecMassCorr < 70 || diMuonRecMassCorr > 110) continue;
	h_recmasscorr->Fill(diMuonRecMassCorr);

	std::cout << "Stops after line " << __LINE__ << std::endl;

        int idxmup_Gen = -1;
        int idxmum_Gen = -1;

        // Gen matching
           for (auto genParticle = genColl->cbegin(); genParticle != genColl->cend(); ++genParticle) {
              const pat::PackedGenParticle& mcMuon = (*genParticle);
              if ( not (abs(mcMuon.pdgId()) == 13 ) ) continue; // make sure it is a muon
              if ( fabs(mcMuon.eta()) > 2.4 ) continue;
              if ( not (mcMuon.pt() > 1.5 ) ) continue;
              if ( deltaR(mcMuon, *(mup->innerTrack())) < 0.1 && mcMuon.charge() > 0 ) idxmup_Gen = std::distance(genColl->cbegin(), genParticle);
              if ( deltaR(mcMuon, *(mum->innerTrack())) < 0.1 && mcMuon.charge() < 0 ) idxmum_Gen = std::distance(genColl->cbegin(), genParticle);
           }
           if ( idxmup_Gen > -1 && idxmum_Gen > -1) {
              double diMuonRecMassGen = (genColl->at(idxmup_Gen).p4() + genColl->at(idxmum_Gen).p4()).M();
              h_genmass->Fill(diMuonRecMassGen);
           }
     }
   }
}
}
// ------------ method called once each job just before starting event loop  ------------
void MuonExercise1::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonExercise1::endJob() {
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonExercise1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise1);
