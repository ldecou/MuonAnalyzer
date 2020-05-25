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
#include <string>
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
#include "FWCore/Utilities/interface/Exception.h"
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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TString.h"
#include "TGraphAsymmErrors.h"

#include "TRandom3.h"
#include "RoccoR.h"
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

    std::string yearToProcess_;
     
    // TFileService
    edm::Service<TFileService> fs;
    
    std::vector<std::string> triggerList;
    std::map<TString, size_t> triggerIdxList; 
    
    // miniAOD Collections
    edm::EDGetTokenT<std::vector<pat::Muon> > muonCollToken;
    //edm::EDGetTokenT<reco::GenParticleCollection> genCollToken;
    edm::EDGetTokenT<pat::PackedGenParticleCollection> genCollToken;

    // Trigger Tokens
    edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjCollToken;
  
    edm::EDGetTokenT<std::vector<reco::Vertex> > vertexCollToken;

    // Trigger matching 
    bool matchTriggerObject(const pat::Muon& mu,
                            const std::map<TString, size_t> & trigList,
                            const edm::TriggerNames& trigNames,
                            edm::Handle<pat::TriggerObjectStandAloneCollection> trigObjs,
                            edm::Handle<edm::TriggerResults> trigBits);


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

    TH1D *h_nvtx;
    TH1D *h_rec_nvtx;

    TH1D *h_rec_loose_pt;
    TH1D *h_rec_loose_eta;
    TH1D *h_rec_loose_nvtx;
  
    TH1D *h_rec_medium_pt;
    TH1D *h_rec_medium_eta;
    TH1D *h_rec_medium_nvtx;
 
    TH1D *h_rec_tight_pt;
    TH1D *h_rec_tight_eta;
    TH1D *h_rec_tight_nvtx;
  
    TH1D *h_rec_tight_iso_pt;
    TH1D *h_rec_tight_iso_eta;
    TH1D *h_rec_tight_iso_nvtx;
  
    TH1D *h_rec_tight_iso_hlt_pt;
    TH1D *h_rec_tight_iso_hlt_eta;
    TH1D *h_rec_tight_iso_hlt_nvtx;

    // Rochester Corrections
    //RoccoR rc("/root://cmseos.fnal.gov///store/user/hats/2018/Muon/Exercises/MuonExercise2/plugin/rcdata.2016.v3");
    //TRandom3 rnd(1234);
    RoccoR rc;
    TRandom3 randomNum;

    // Efficiencies
    TGraphAsymmErrors *gae_rec_pt;
    TGraphAsymmErrors *gae_rec_eta;
    TGraphAsymmErrors *gae_rec_nvtx;
     
    TGraphAsymmErrors *gae_rec_loose_pt;
    TGraphAsymmErrors *gae_rec_loose_eta;
    TGraphAsymmErrors *gae_rec_loose_nvtx;
    
    TGraphAsymmErrors *gae_rec_medium_pt;
    TGraphAsymmErrors *gae_rec_medium_eta;
    TGraphAsymmErrors *gae_rec_medium_nvtx;
    
    TGraphAsymmErrors *gae_rec_tight_pt;
    TGraphAsymmErrors *gae_rec_tight_eta;
    TGraphAsymmErrors *gae_rec_tight_nvtx;
     
    TGraphAsymmErrors *gae_rec_tight_iso_pt;
    TGraphAsymmErrors *gae_rec_tight_iso_eta;
    TGraphAsymmErrors *gae_rec_tight_iso_nvtx;
     
    TGraphAsymmErrors *gae_rec_tight_iso_hlt_pt;
    TGraphAsymmErrors *gae_rec_tight_iso_hlt_eta;
    TGraphAsymmErrors *gae_rec_tight_iso_hlt_nvtx;


};
//
// constructors and destructor

 MuonExercise1::MuonExercise1(const edm::ParameterSet& iConfig) {

  usesResource("TFileService");
 
  edm::InputTag muonTag("slimmedMuons");
  edm::InputTag genPartTag("packedGenParticles");
  //edm::InputTag genPartTag("prunedGenParticles");
 
  edm::InputTag vertexTag("offlineSlimmedPrimaryVertices");
  edm::InputTag triggerTag("TriggerResults", "", "HLT");
  edm::InputTag trigObjTag("selectedPatTrigger");

  muonCollToken = consumes<pat::MuonCollection>(muonTag);
  genCollToken = consumes<pat::PackedGenParticleCollection>(genPartTag);
  //genCollToken = consumes<reco::GenParticleCollection>(genPartTag);

  vertexCollToken = consumes<std::vector<reco::Vertex> >(vertexTag);
  trigResultsToken = consumes<edm::TriggerResults>(triggerTag);
  trigObjCollToken = consumes<pat::TriggerObjectStandAloneCollection>(trigObjTag);

  // Histograms

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

  h_nvtx = fs->make<TH1D>("nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5);
  h_rec_nvtx = fs->make<TH1D>("rec_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5);

  h_rec_loose_pt   = fs->make<TH1D>("rec_loose_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0);
  h_rec_loose_eta  = fs->make<TH1D>("rec_loose_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4);
  h_rec_loose_nvtx = fs->make<TH1D>("rec_loose_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5);
 
  h_rec_medium_pt   = fs->make<TH1D>("rec_medium_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0);
  h_rec_medium_eta  = fs->make<TH1D>("rec_medium_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4);
  h_rec_medium_nvtx = fs->make<TH1D>("rec_medium_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5);
 
  h_rec_tight_pt   = fs->make<TH1D>("rec_tight_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0);
  h_rec_tight_eta  = fs->make<TH1D>("rec_tight_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4);
  h_rec_tight_nvtx = fs->make<TH1D>("rec_tight_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5);

  h_rec_tight_iso_pt   = fs->make<TH1D>("rec_tight_iso_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0);
  h_rec_tight_iso_eta  = fs->make<TH1D>("rec_tight_iso_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4);
  h_rec_tight_iso_nvtx = fs->make<TH1D>("rec_tight_iso_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5);
 
  h_rec_tight_iso_hlt_pt   = fs->make<TH1D>("rec_tight_iso_hlt_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0);
  h_rec_tight_iso_hlt_eta  = fs->make<TH1D>("rec_tight_iso_hlt_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4);
  h_rec_tight_iso_hlt_nvtx = fs->make<TH1D>("rec_tight_iso_hlt_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5);

  // Rochester corrections
  //if (yearToProcess_ == "2016") rc.init(edm::FileInPath("data/RoccoR2016.txt").fullPath()); 
  rc.init(edm::FileInPath("data/RoccoR2017.txt").fullPath()); 
  //else rc.init(edm::FileInPath("data/RoccoR2018.txt").fullPath()); 
  randomNum.SetSeed(1234);

  // Efficiencies
  unsigned int nbinpt   = h_recpt->GetNbinsX();
  unsigned int nbineta  = h_receta->GetNbinsX();
  unsigned int nbinnvtx = h_rec_nvtx->GetNbinsX();
 
  gae_rec_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_pt->SetName("eff_rec_pt");
  gae_rec_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_eta->SetName("eff_rec_eta");
  gae_rec_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_nvtx->SetName("eff_rec_nvtx");

  gae_rec_loose_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_loose_pt->SetName("eff_rec_loose_pt");
  gae_rec_loose_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_loose_eta->SetName("eff_rec_loose_eta");
  gae_rec_loose_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_loose_nvtx->SetName("eff_rec_loose_nvtx");

  gae_rec_medium_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_medium_pt->SetName("eff_rec_medium_pt");
  gae_rec_medium_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_medium_eta->SetName("eff_rec_medium_eta");
  gae_rec_medium_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_medium_nvtx->SetName("eff_rec_medium_nvtx");

  gae_rec_tight_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_tight_pt->SetName("eff_rec_tight_pt");
  gae_rec_tight_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_tight_eta->SetName("eff_rec_tight_eta");
  gae_rec_tight_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_tight_nvtx->SetName("eff_rec_tight_nvtx");

  gae_rec_tight_iso_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_tight_iso_pt->SetName("eff_rec_tight_iso_pt");
  gae_rec_tight_iso_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_tight_iso_eta->SetName("eff_rec_tight_iso_eta");
  gae_rec_tight_iso_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_tight_iso_nvtx->SetName("eff_rec_tight_iso_nvtx");

  gae_rec_tight_iso_hlt_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_tight_iso_hlt_pt->SetName("eff_rec_tight_iso_hlt_pt");
  gae_rec_tight_iso_hlt_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_tight_iso_hlt_eta->SetName("eff_rec_tight_iso_hlt_eta");
  gae_rec_tight_iso_hlt_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_tight_iso_hlt_nvtx->SetName("eff_rec_tight_iso_hlt_nvtx");
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



///////////////////
// Z Mass
////////////////////
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
     
     double mupSF(0.);

     if (mup->genParticle()) {
         
     mupSF = rc.kSpreadMC(mup->charge(),
                          mup->pt(),
                          mup->eta(),
                          mup->phi(),
                          mup->genParticle()->pt());
     }
     
     else  {
     mupSF = rc.kSmearMC(mup->charge(),
			 mup->pt(),
			 mup->eta(),
			 mup->phi(),
			 mup->innerTrack()->hitPattern().trackerLayersWithMeasurement(),
			 //mup->genParticle()->pt(),
			 randomNum.Rndm());
     }

     std::cout << "Stops after line " << __LINE__ << std::endl;

     // find mu-
     for (auto mum = muonColl->cbegin(); mum != muonColl->cend(); ++mum) {
        if ( not (mum->charge() < 0 ) ) continue;
        if ( not (mum->isGlobalMuon()) ) continue;
        if ( not (mum->pt() > 20.0 ) ) continue;
        if ( fabs(mum->eta()) > 2.4 ) continue;
        if ( not (mum->chargedHadronIso() < 0.15 ) ) continue;

	std::cout << "Stops after line " << __LINE__ << std::endl;

	double mumSF(0.);

	if (mum->genParticle()) {

	mumSF = rc.kSpreadMC(mum->charge(),
			     mum->pt(),
			     mum->eta(),
			     mum->phi(),
			     mum->genParticle()->pt());
	}

	else {
        mumSF = rc.kSmearMC(mum->charge(),
		            mum->pt(),
		            mum->eta(),
		            mum->phi(),
		            mum->innerTrack()->hitPattern().trackerLayersWithMeasurement(),
		            //mum->genParticle()->pt(),
		            randomNum.Rndm());
	}

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
              //const reco::GenParticle& mcMuon = (*genParticle);
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

////////////////////
// Efficiencies
////////////////////
  //edm::Handle <pat::PackedGenParticleCollection> packedgenColl;
  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(vertexCollToken, vertices);
  if(!vertices.isValid()) {
    throw cms::Exception("Vertex collection not valid!");
  }
  
  int nGoodVtx = 0;
  const reco::Vertex *goodVtx = nullptr; 
  for(std::vector<reco::Vertex>::const_iterator it=vertices->begin(), endVtx = vertices->end(); it!=endVtx; ++it) {
    if(!it->isFake() && it->ndof()>4 && it->position().Rho()<2. && std::abs(it->position().Z())<24.) {
      nGoodVtx++;
      if(goodVtx==nullptr) goodVtx = &(*it);
    }
  }
  // Require a good vertex
  if(nGoodVtx==0) return;
  
  // Retrieve the GenParticle collection and loop over it 
  //edm::Handle<reco::GenParticleCollection> genColl;
  //iEvent.getByToken(genCollToken, genColl);
  if(!genColl.isValid()) {
    throw cms::Exception("GenParticle collection not valid!");
  }
  std::vector<const reco::GenParticle*> genmus;	


	for(auto&& genPart : *(genColl.product())) {
  		// Check if it's a muon from Drell-Yan process
		if(genPart.isPromptFinalState() && std::abs(genPart.pdgId()) == 13) {
			// Only muons within acceptance and pt>20
  			if(genPart.pt() > 20. && std::abs(genPart.eta())<2.4) {
  				h_genpt->Fill(genPart.pt());
 				h_nvtx->Fill(nGoodVtx);
                        }
                }
        }
  if(genmus.size()==0) return;
  std::vector<const reco::GenParticle*>::const_iterator gmbeg = genmus.begin(), gmend = genmus.end();

  // --- Now prepare all the trigger objects --- 
  // Retrieve TriggerResults 
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(trigResultsToken, triggerBits);
  if(!triggerBits.isValid()) {
    throw cms::Exception("TriggerResults collection not valid!");
  }

  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerBits);

  // Find trigger indexes (only once)
  if(triggerIdxList.size()==0) {
    for(auto&& t : triggerList) {
      for(size_t i=0, n=triggerBits->size(); i<n; ++i) {
        if(TString(triggerNames.triggerName(i)).Contains((t+"_v").c_str())) {
          triggerIdxList[(t+"_v").c_str()] = i;
          break;
        }
      }
    }
  }

  // Retrieve TriggerObjects 
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(trigObjCollToken, triggerObjects);
  if(!triggerObjects.isValid()) {
    throw cms::Exception("TriggerObjectStandAloneCollection collection not valid!");
  }

  // Retrieve the pat::Muon collection and loop over it 
  //edm::Handle<pat::MuonCollection> muonColl;
  if(!muonColl.isValid()) {
    throw cms::Exception("Muon collection not valid!");
  }

  for(auto mu = muonColl->cbegin(); mu != muonColl->cend(); ++mu) {
    cout << "Muon PT: " << (*mu).pt() << endl;        
    // Let's skip muons that are standalone only
    // Check if it is matched to a GenParticle
    const reco::GenParticle *gmu = (*mu).genParticle();
    if(!(*mu).genParticle()) continue;

    // Check if the GenParticle is among the ones we selected 
    //if(std::find(gmbeg, gmend, gmu)==gmend) continue;

    // Fill the plots here 
    // (Let's not cut on pat::Muon pt and |eta| (already cut at GEN level))

    // -- ID --
    //  Fill ID plots 
    //  - Loose 
    if((*mu).isLooseMuon()) {
    cout << "Loose Muon PT: " << (*mu).pt() << endl;
    h_rec_loose_pt->Fill((*mu).pt());
    h_rec_loose_eta->Fill((*mu).eta());
    h_rec_loose_nvtx->Fill(nGoodVtx);    

    }

    //  - Medium 
    if((*mu).isMediumMuon()) {
    cout << "Medium Muon PT: " << (*mu).pt() << endl;
    h_rec_medium_pt->Fill((*mu).pt());
    h_rec_medium_eta->Fill((*mu).eta());
    h_rec_medium_nvtx->Fill(nGoodVtx);
    }

    //  - Tight
    if((*mu).isTightMuon(*goodVtx)==false) continue; // if it's not tight, nothing else to do 
    cout << "Tight Muon PT: " << (*mu).pt() << endl;
    h_rec_tight_pt->Fill((*mu).pt());
    h_rec_tight_eta->Fill((*mu).eta());
    h_rec_tight_nvtx->Fill(nGoodVtx);

    // Now isolation, only for tight muons
    if((*mu).isIsolationValid()==false) continue;
    const reco::MuonPFIsolation &pfR04 = (*mu).pfIsolationR04();

    // Calculate PF combined relative isolation with Delta-beta correction 
    double chargedHadronIso = (*mu).pfIsolationR04().sumChargedHadronPt;
    double neutralHadronIso  = (*mu).pfIsolationR04().sumNeutralHadronEt;
    double photonIso  = (*mu).pfIsolationR04().sumPhotonEt;
    double chargedHadronIsoPU = (*mu).pfIsolationR04().sumPUPt;

    double corriso = ( chargedHadronIso + std::max(0., neutralHadronIso + photonIso - 0.5*chargedHadronIsoPU) )/(*mu).pt();

    if(corriso>0.15) continue; // not isolated, nothing else to do
    h_rec_tight_iso_pt->Fill((*mu).pt());
    h_rec_tight_iso_eta->Fill((*mu).eta());
    h_rec_tight_iso_nvtx->Fill(nGoodVtx);
    
    // Finally, let's see if the isolated tight muon fired a trigger
    bool passTrigger = matchTriggerObject((*mu), triggerIdxList, triggerNames, triggerObjects, triggerBits);
    if(passTrigger==false) continue;
    h_rec_tight_iso_hlt_pt->Fill((*mu).pt());
    h_rec_tight_iso_hlt_eta->Fill((*mu).eta());
    h_rec_tight_iso_hlt_nvtx->Fill(nGoodVtx);

  } // end for(auto mu = muonColl->cbegin(); mu != muonColl->cend(); ++mu)



}
// ------------ method called once each job just before starting event loop  ------------
void MuonExercise1::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise1::endJob() {

   // Compute efficiencies
   try{gae_rec_pt  ->Divide(h_recpt  , h_genpt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_eta ->Divide(h_receta , h_geneta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_nvtx->Divide(h_rec_nvtx, h_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
 
   try{gae_rec_loose_pt  ->Divide(h_rec_loose_pt  , h_recpt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_loose_eta ->Divide(h_rec_loose_eta , h_receta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_loose_nvtx->Divide(h_rec_loose_nvtx, h_rec_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
 
   try{gae_rec_medium_pt  ->Divide(h_rec_medium_pt  , h_recpt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_medium_eta ->Divide(h_rec_medium_eta , h_receta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_medium_nvtx->Divide(h_rec_medium_nvtx, h_rec_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
 
   try{gae_rec_tight_pt  ->Divide(h_rec_tight_pt  , h_recpt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_tight_eta ->Divide(h_rec_tight_eta , h_receta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_tight_nvtx->Divide(h_rec_tight_nvtx, h_rec_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
 
   try{gae_rec_tight_iso_pt  ->Divide(h_rec_tight_iso_pt  , h_rec_tight_pt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_tight_iso_eta ->Divide(h_rec_tight_iso_eta , h_rec_tight_eta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_tight_iso_nvtx->Divide(h_rec_tight_iso_nvtx, h_rec_tight_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
 
   try{gae_rec_tight_iso_hlt_pt  ->Divide(h_rec_tight_iso_hlt_pt  , h_rec_tight_iso_pt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_tight_iso_hlt_eta ->Divide(h_rec_tight_iso_hlt_eta , h_rec_tight_iso_eta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}
   try{gae_rec_tight_iso_hlt_nvtx->Divide(h_rec_tight_iso_hlt_nvtx, h_rec_tight_iso_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}


   // Labels
   gae_rec_pt  ->GetXaxis()->SetTitle(h_recpt  ->GetXaxis()->GetTitle());
   gae_rec_eta ->GetXaxis()->SetTitle(h_receta ->GetXaxis()->GetTitle());
   gae_rec_nvtx->GetXaxis()->SetTitle(h_rec_nvtx->GetXaxis()->GetTitle());
 
   gae_rec_loose_pt  ->GetXaxis()->SetTitle(h_rec_loose_pt  ->GetXaxis()->GetTitle());
   gae_rec_loose_eta ->GetXaxis()->SetTitle(h_rec_loose_eta ->GetXaxis()->GetTitle());
   gae_rec_loose_nvtx->GetXaxis()->SetTitle(h_rec_loose_nvtx->GetXaxis()->GetTitle());
 
   gae_rec_medium_pt  ->GetXaxis()->SetTitle(h_rec_medium_pt  ->GetXaxis()->GetTitle());
   gae_rec_medium_eta ->GetXaxis()->SetTitle(h_rec_medium_eta ->GetXaxis()->GetTitle());
   gae_rec_medium_nvtx->GetXaxis()->SetTitle(h_rec_medium_nvtx->GetXaxis()->GetTitle());
 
   gae_rec_tight_pt  ->GetXaxis()->SetTitle(h_rec_tight_pt  ->GetXaxis()->GetTitle());
   gae_rec_tight_eta ->GetXaxis()->SetTitle(h_rec_tight_eta ->GetXaxis()->GetTitle());
   gae_rec_tight_nvtx->GetXaxis()->SetTitle(h_rec_tight_nvtx->GetXaxis()->GetTitle());
 
   gae_rec_tight_iso_pt  ->GetXaxis()->SetTitle(h_rec_tight_iso_pt  ->GetXaxis()->GetTitle());
   gae_rec_tight_iso_eta ->GetXaxis()->SetTitle(h_rec_tight_iso_eta ->GetXaxis()->GetTitle());
   gae_rec_tight_iso_nvtx->GetXaxis()->SetTitle(h_rec_tight_iso_nvtx->GetXaxis()->GetTitle());
 
   gae_rec_tight_iso_hlt_pt  ->GetXaxis()->SetTitle(h_rec_tight_iso_hlt_pt  ->GetXaxis()->GetTitle());
   gae_rec_tight_iso_hlt_eta ->GetXaxis()->SetTitle(h_rec_tight_iso_hlt_eta ->GetXaxis()->GetTitle());
   gae_rec_tight_iso_hlt_nvtx->GetXaxis()->SetTitle(h_rec_tight_iso_hlt_nvtx->GetXaxis()->GetTitle());
 
   gae_rec_pt  ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_eta ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_nvtx->GetYaxis()->SetTitle("Efficiency");
 
   gae_rec_loose_pt  ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_loose_eta ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_loose_nvtx->GetYaxis()->SetTitle("Efficiency");
 
   gae_rec_medium_pt  ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_medium_eta ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_medium_nvtx->GetYaxis()->SetTitle("Efficiency");
 
   gae_rec_tight_pt  ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_tight_eta ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_tight_nvtx->GetYaxis()->SetTitle("Efficiency");
 
   gae_rec_tight_iso_pt  ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_tight_iso_eta ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_tight_iso_nvtx->GetYaxis()->SetTitle("Efficiency");
 
   gae_rec_tight_iso_hlt_pt  ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_tight_iso_hlt_eta ->GetYaxis()->SetTitle("Efficiency");
   gae_rec_tight_iso_hlt_nvtx->GetYaxis()->SetTitle("Efficiency");

}

bool MuonExercise1::matchTriggerObject(const pat::Muon& mu,
                                              const std::map<TString, size_t>& trigList,
                                              const edm::TriggerNames& trigNames,
                                              edm::Handle<pat::TriggerObjectStandAloneCollection> trigObjs,
                                              edm::Handle<edm::TriggerResults> trigBits) {
   double dR = 0.1;
 
   for(auto&& t : trigList) {
     if(trigBits->accept(t.second)==false) continue;
     for(pat::TriggerObjectStandAlone obj : *trigObjs) {
       obj.unpackPathNames(trigNames);
       if(obj.hasPathName((t.first+"*").Data(), true, true)) {
         if(reco::deltaR(obj, *(mu.innerTrack())) < dR) {
           return true;
         }
       }
     } // end for(pat::TriggerObjectStandAlone obj : *trigObjs)
   } // end for(auto&& t : trigList)
 
   return false;
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
