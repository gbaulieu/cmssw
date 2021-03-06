// -*- C++ -*-
//
// Package:    PrintL1TkObjects
// Class:      PrintL1TkObjects
// 
/**\class PrintL1TkObjects PrintL1TkObjects.cc SLHCUpgradeSimulations/PrintL1TkObjects/src/PrintL1TkObjects.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Thu Nov 14 11:22:13 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"

#include "TFile.h"
#include "TH1F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1extra;



//
// class declaration
//

class PrintL1TkObjects : public edm::EDAnalyzer {
   public:

   typedef TTTrack< Ref_PixelDigi_ >  L1TkTrackType;
   typedef std::vector< L1TkTrackType >  L1TkTrackCollectionType;

      explicit PrintL1TkObjects(const edm::ParameterSet&);
      ~PrintL1TkObjects();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

	// to test the L1TkPrimaryVertex :
	edm::InputTag L1VtxInputTag;

	// for L1TrackEtmiss:
	edm::InputTag L1TkEtMissInputTag;

	// for L1TkEmParticles
        edm::InputTag L1TkPhotonsInputTag;

	// for L1TkElectrons
	edm::InputTag L1TkElectronsInputTag;

	// for L1TkJetParticles
        edm::InputTag L1TkJetsInputTag;

	// for L1TkHTMParticle
	edm::InputTag L1TkHTMInputTag;

        // for L1TkMuonParticle
        edm::InputTag L1TkMuonsInputTag;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PrintL1TkObjects::PrintL1TkObjects(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  L1VtxInputTag = iConfig.getParameter<edm::InputTag>("L1VtxInputTag") ;
  L1TkEtMissInputTag = iConfig.getParameter<edm::InputTag>("L1TkEtMissInputTag");
  L1TkElectronsInputTag = iConfig.getParameter<edm::InputTag>("L1TkElectronsInputTag");
  L1TkPhotonsInputTag = iConfig.getParameter<edm::InputTag>("L1TkPhotonsInputTag");
  L1TkJetsInputTag = iConfig.getParameter<edm::InputTag>("L1TkJetsInputTag");
  L1TkHTMInputTag = iConfig.getParameter<edm::InputTag>("L1TkHTMInputTag");
  L1TkMuonsInputTag = iConfig.getParameter<edm::InputTag>("L1TkMuonsInputTag");

}


PrintL1TkObjects::~PrintL1TkObjects()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PrintL1TkObjects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   std::cout << " ----  a new event ----- " << std::endl;

	// First, retrieve the generated primary vertex


  edm::Handle<edm::HepMCProduct> HepMCEvt;
  iEvent.getByLabel("generator",HepMCEvt);

     float zvtx_gen = -999;

 if (HepMCEvt.isValid() ) {

     const HepMC::GenEvent* MCEvt = HepMCEvt->GetEvent();
     const double mm=0.1;

     for ( HepMC::GenEvent::vertex_const_iterator ivertex = MCEvt->vertices_begin(); ivertex != MCEvt->vertices_end(); ++ivertex )
         {
             bool hasParentVertex = false;
 
             // Loop over the parents looking to see if they are coming from a production vertex
             for (
                 HepMC::GenVertex::particle_iterator iparent = (*ivertex)->particles_begin(HepMC::parents);
                 iparent != (*ivertex)->particles_end(HepMC::parents);
                 ++iparent
             )
                 if ( (*iparent)->production_vertex() )
                 {
                     hasParentVertex = true;
                     break;
                }
 
             // Reject those vertices with parent vertices
             if (hasParentVertex) continue;
 
             // Get the position of the vertex
             HepMC::FourVector pos = (*ivertex)->position();
	     zvtx_gen = pos.z()*mm; 

	     break;  // there should be one single primary vertex

          }  // end loop over gen vertices

     std::cout << " Generated zvertex : " << zvtx_gen << std::endl;
  }


	// ----------------------------------------------------------------------
	//
        // retrieve the L1 vertices
        
 edm::Handle<L1TkPrimaryVertexCollection> L1VertexHandle;
 iEvent.getByLabel(L1VtxInputTag,L1VertexHandle);
 std::vector<L1TkPrimaryVertex>::const_iterator vtxIter;
 
 if ( L1VertexHandle.isValid() ) {
     std::cout << " -----  L1TkPrimaryVertex objects   ----- " << std::endl;
     vtxIter = L1VertexHandle->begin();     // only one algorithm is run in the L1TkPrimaryVertexProducer
					    // (in contrast to earlier, under-dev, versions of the code)
        float z = vtxIter -> getZvertex();
        float sum = vtxIter -> getSum();
        std::cout << " a vertex with  zvtx " << z << " (cm) and SumPT " << sum << " (GeV) " << std::endl;
 }

	//
        // ----------------------------------------------------------------------
	// retrieve the EtMiss objects
	//

 edm::Handle<L1TkEtMissParticleCollection> L1TkEtMissHandle;
 iEvent.getByLabel(L1TkEtMissInputTag, L1TkEtMissHandle);
 std::vector<L1TkEtMissParticle>::const_iterator etmIter;

 if (L1TkEtMissHandle.isValid() ) {
    std::cout << " -----  L1TkEtMiss objects  -----  " << std::endl; 
    etmIter = L1TkEtMissHandle -> begin();	// idem: only one TrkMET now.
	float etmis = etmIter -> et();
	const edm::Ref< L1TkPrimaryVertexCollection > vtxRef = etmIter -> getVtxRef();
	float zvtx = vtxRef -> getZvertex();
        float etMissPU = etmIter -> etMissPU();
	std::cout << " ETmiss = " << etmis << " for zvtx = " << zvtx << " and ETmiss from PU = " << etMissPU << std::endl;
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the L1TkJetParticle objects
        //
        
 edm::Handle<L1TkJetParticleCollection> L1TkJetsHandle;
 iEvent.getByLabel(L1TkJetsInputTag, L1TkJetsHandle);
 std::vector<L1TkJetParticle>::const_iterator jetIter ;

 if ( L1TkJetsHandle.isValid() ) {
    std::cout << " -----   L1TkJetParticle  objects -----  " << std::endl;
    for (jetIter = L1TkJetsHandle -> begin(); jetIter != L1TkJetsHandle->end(); ++jetIter) {
        float et = jetIter -> pt();
        float phi = jetIter -> phi();
        float eta = jetIter -> eta();
        int bx = jetIter -> bx() ;
	float jetvtx = jetIter -> getJetVtx();
        const edm::Ref< L1JetParticleCollection > Jetref = jetIter -> getJetRef();
        float et_L1Jet = Jetref -> et();
	const std::vector< edm::Ptr< L1TkTrackType > > trkPtrs = jetIter -> getTrkPtrs() ;
	L1JetParticle::JetType type = Jetref -> type();

        std::cout << " a Jet candidate ET eta phi zvertex " << et << " " << eta << " " << phi << " " << jetvtx  << std::endl;
        std::cout << "                Calo  ET, typ " << et_L1Jet << " " << type << std::endl;
        std::cout << "                bx = " << bx << std::endl;
	std::cout << "      Tracks associated with the jet : " << std::endl;
        for (int it=0; it < (int)trkPtrs.size(); it++) {
            edm::Ptr< L1TkTrackType > aTrack = trkPtrs.at( it );
            std::cout << "             a track PT " << aTrack -> getMomentum().perp() << " eta " << 
		   aTrack -> getMomentum().eta() << " phi " << aTrack -> getMomentum().phi() << std::endl;
        }

    }
 }

        //
        // ----------------------------------------------------------------------
        // retrieve HT and HTM
	//

 edm::Handle<L1TkHTMissParticleCollection> L1TkHTMHandle;
 iEvent.getByLabel(L1TkHTMInputTag, L1TkHTMHandle);

 if ( L1TkHTMHandle.isValid() ) {
	std::cout << " -----  L1TkHTMissParticle: size (should be 1) = " << L1TkHTMHandle -> size() << std::endl;
	std::vector<L1TkHTMissParticle>::const_iterator HTMIter = L1TkHTMHandle -> begin();
	float HTT = HTMIter -> EtTotal();
	float HTM = HTMIter -> EtMiss();
	//float HTM_the_same = HTMIter -> et();

	// phi of the HTM vector :
	float phi = HTMIter -> phi();
	std::cout << " HTT = " << HTT << " HTM = " << HTM << " " << "phi(HTM) = " << phi << std::endl;
	 
	// access the L1TkJets used to build HT and HTM :
	const edm::RefProd< L1TkJetParticleCollection > jetCollRef = HTMIter -> getjetCollectionRef();
 	std::vector<L1TkJetParticle>::const_iterator jet = jetCollRef -> begin();
	std::cout << " ET of the first L1TkJet = " << jet -> et() << std::endl;
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the L1TkEmParticle objects
	//

 edm::Handle<L1TkEmParticleCollection> L1TkPhotonsHandle;
 iEvent.getByLabel(L1TkPhotonsInputTag, L1TkPhotonsHandle);
 std::vector<L1TkEmParticle>::const_iterator phoIter ;

 if ( L1TkPhotonsHandle.isValid() ) {
    std::cout << " -----   L1TkEmParticle  objects -----  " << std::endl;
    for (phoIter = L1TkPhotonsHandle -> begin(); phoIter != L1TkPhotonsHandle->end(); ++phoIter) {
	float et = phoIter -> pt();
	float phi = phoIter -> phi();
	float eta = phoIter -> eta();
	int bx = phoIter -> bx() ;
        float trkisol = phoIter -> getTrkIsol() ;
	const edm::Ref< L1EmParticleCollection > EGref = phoIter -> getEGRef();
	float et_L1Calo = EGref -> et();
	float eta_calo = EGref -> eta();
	float phi_calo = EGref -> phi();

	std::cout << " a photon candidate ET eta phi trkisol " << et << " " << eta << " " << phi << " " << trkisol << std::endl;
	std::cout << "                Calo  ET eta phi " << et_L1Calo << " " << eta_calo << " " << phi_calo << std::endl; 
	std::cout << "                bx = " << bx << std::endl;
    }
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the L1TkElectronParticle objects
        //

 edm::Handle<L1TkElectronParticleCollection> L1TkElectronsHandle;
 iEvent.getByLabel(L1TkElectronsInputTag, L1TkElectronsHandle);
 std::vector<L1TkElectronParticle>::const_iterator eleIter ;

 if ( L1TkElectronsHandle.isValid() ) {
    std::cout << " -----   L1TkElectronParticle  objects -----  " << std::endl;
    for (eleIter = L1TkElectronsHandle -> begin(); eleIter != L1TkElectronsHandle->end(); ++eleIter) {
        float et = eleIter -> pt();
        float phi = eleIter -> phi();
        float eta = eleIter -> eta();
	int bx = eleIter -> bx();
    	float trkisol = eleIter -> getTrkIsol() ;
	float ztr = eleIter -> getTrkzVtx() ;
        std::cout << "an electron candidate ET eta phi bx trkisol ztr " << et << " " << eta << " " << phi << " " << bx << " " << trkisol << " " << ztr << std::endl;

        const edm::Ref< L1EmParticleCollection > EGref = eleIter -> getEGRef();
        if ( EGref.isNonnull() ) {
           float et_L1Calo = EGref -> et();
           float eta_calo = EGref -> eta();
           float phi_calo = EGref -> phi();
           std::cout << "                Calo  ET eta phi " << et_L1Calo << " " << eta_calo << " " << phi_calo << std::endl;
	}
	else {
	    std::cout << " .... edm::Ref to EGamma is unvalid !! ?  " << std::endl;
	}

        const edm::Ptr< L1TkTrackType > TrkRef = eleIter -> getTrkPtr();
	if ( TrkRef.isNonnull() ) {
            float pt_track = TrkRef -> getMomentum().perp();
            float phi_track = TrkRef -> getMomentum().phi();
            float eta_track = TrkRef -> getMomentum().eta();
            float ztrack = TrkRef -> getPOCA().z();
            std::cout << "                Track PT eta phi ztr " << pt_track << " " << eta_track << " " << phi_track << " " << ztrack << std::endl;
	}
	else {
	    std::cout << " ... edm::Ptr to L1Tracks is unvalid (e.g. electron was matched to stubs) " << std::endl;
	}
    }
 }

        //  
        // ----------------------------------------------------------------------
        // retrieve the L1TkMuons
        //

 edm::Handle<L1TkMuonParticleCollection> L1TkMuonsHandle;
 iEvent.getByLabel(L1TkMuonsInputTag, L1TkMuonsHandle);
 std::vector<L1TkMuonParticle>::const_iterator muIter;

 if ( L1TkMuonsHandle.isValid() ) {
    std::cout << " -----   L1TkMuonPaticle objects ---- " << std::endl;
    for (muIter = L1TkMuonsHandle -> begin(); muIter != L1TkMuonsHandle->end(); ++muIter) {
        float pt = muIter -> pt();
        float eta = muIter -> eta();
        float phi = muIter -> phi();
        float zvtx = muIter -> getTrkzVtx();
        unsigned int quality = muIter -> quality();
        // access the quality via the reference :
        const edm::Ref< L1MuonParticleCollection >      MuRef = muIter -> getMuRef();
        unsigned int qualityBis = MuRef -> gmtMuonCand().quality();
        std::cout << " a muon candidate pt eta phi " << pt << " " << eta << " " << phi << " zvertex = " << zvtx << std::endl;
        std::cout << "   quality (the two qual flags are the same by definition) = " << quality << " " << qualityBis << std::endl;

    }
 }




}


// ------------ method called once each job just before starting event loop  ------------
void 
PrintL1TkObjects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PrintL1TkObjects::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PrintL1TkObjects::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PrintL1TkObjects::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PrintL1TkObjects::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PrintL1TkObjects::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PrintL1TkObjects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrintL1TkObjects);
