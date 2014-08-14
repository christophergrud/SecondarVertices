// -*- C++ -*-
//
// Package:    SecondaryVertices
// Class:      SecondaryVertices
// 
/**\class SecondaryVertices SecondaryVertices.cc MyAnalysis/SecondaryVertices/src/SecondaryVertices.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Grud
//         Created:  Mon Aug  4 12:39:44 CDT 2014
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
//PSimHitContainer includes PSimHit and Vector already
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "MyAnalysis/Vertex/interface/VertFitter.h"
#include "MyAnalysis/Vertex/interface/VertProducer.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TH1F.h"
#include "TH2F.h"

#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <TLorentzVector.h>
//
// class declaration
//

class SecondaryVertices : public edm::EDAnalyzer {
   public:
      explicit SecondaryVertices(const edm::ParameterSet&);
      ~SecondaryVertices();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

      edm::Service<TFileService> fs;

      TH1F *h1_nMatched;
      TH1F *h1_perEvent_eff;
      TH1F *h1_perEvent_numberOfReco;

      TH1F *h1_LL_pt;
      TH1F *h1_LL_eta;
      TH1F *h1_LL_rho;

      TH2F *h2_LLpT_vs_deltaR;

      TH1F *h1_matched_pt;
      TH1F *h1_matched_eta;
      TH1F *h1_matched_rho;
      TH1F *h1_matched_mass;

      TH1F *h1_eff_pt;
      TH1F *h1_eff_eta;
      TH1F *h1_eff_rho;

      TH1F *h1_chi2_1;
      TH1F *h1_chi2_1_5;
      TH1F *h1_chi2_5;

      TH1F *h1_deltaR;
      TH1F *h1_deltaPhi;
      TH1F *h1_deltaEta;
      
      TH1F *h1_deltaR_highMass;
      TH1F *h1_howMany_highMass;
      TH1F *h1_howMany_lowChi2_highMass;
      TH1F *h1_rho_highMass;

      TH1F *h1_lowChi2_mass;
      TH1F *h1_vertCand_mass;
      TH1F *h1_vertCand_rho;

      TH2F *h2_normal_pt_vs_chi2;
      TH2F *h2_mass_vs_dR;
      TH2F *h2_mass_vs_chi2;
      TH2F *h2_mass_vs_ndof;
      TH2F *h2_mass_vs_normalizedChi2;

      float howManyReco;
      float howManyLL;
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
SecondaryVertices::SecondaryVertices(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


SecondaryVertices::~SecondaryVertices()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SecondaryVertices::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  edm::Handle<std::vector <reco::VertexCompositeCandidate> > myVertexHandle;
  iEvent.getByLabel(edm::InputTag("generalVertCandidates:MyVertices"), myVertexHandle);

  edm::Handle<std::vector<reco::GenParticle>> genPartHandle;
  iEvent.getByLabel(edm::InputTag("genParticles"), genPartHandle);

  std::vector< TLorentzVector >   vLL;
  std::vector< TLorentzVector >   vVS;
  std::vector< TLorentzVector >   vMatchedLL;
  std::vector< TLorentzVector >   vMatchedVS;
  std::vector< double >           vRho;
  std::vector< double >           vVSRho;
  std::vector< double >           vMatchedRho;
  std::vector< double >           vRecoVertChi2;
  std::vector< double >           vMatchedVertChi2;

  int     numberOfMatched = 0;
  float   numberOfLL = 0;
  float   numberOfReco = 0;
  double  vChi2;
  double  vndof;
  double  vNormChi2;
  double  vMass;

  for (std::vector<reco::GenParticle>::const_iterator iParticle = genPartHandle->begin(); iParticle != genPartHandle->end(); ++iParticle){
    
    TLorentzVector tmpV;
    TLorentzVector posMuon;
    TLorentzVector negMuon;


    if (iParticle->status()!= 3) continue;

    if (iParticle->pdgId()==6001113){ 
      numberOfLL++;
      howManyLL++;
      tmpV.SetPtEtaPhiM(iParticle->pt(), iParticle->eta(), iParticle->phi(), iParticle->mass());
      h1_LL_pt->Fill(iParticle->pt());
      h1_LL_eta->Fill(iParticle->eta());
      for (size_t i = 0; i<iParticle->numberOfDaughters(); ++i){
        if (iParticle->daughter(i)->pdgId()==13){
          for (size_t j = 0; j<iParticle->daughter(i)->numberOfDaughters(); ++j){
            if (iParticle->daughter(i)->daughter(j)->pdgId()==13){
              double LL_x = iParticle->daughter(i)->daughter(j)->vx();
              double LL_y = iParticle->daughter(i)->daughter(j)->vy();
              double LL_rho = sqrt(LL_x*LL_x + LL_y*LL_y);
              h1_LL_rho->Fill(LL_rho);
              vRho.push_back(LL_rho);
              vLL.push_back(tmpV);
              posMuon.SetPtEtaPhiM(iParticle->daughter(i)->daughter(j)->pt(),iParticle->daughter(i)->daughter(j)->eta(),iParticle->daughter(i)->daughter(j)->phi(),iParticle->daughter(i)->daughter(j)->mass() );
            }
          }
        }
        if (iParticle->daughter(i)->pdgId()==-13){
          for (size_t j = 0; j<iParticle->daughter(i)->numberOfDaughters(); ++j){
            if (iParticle->daughter(i)->daughter(j)->pdgId()==-13){
              negMuon.SetPtEtaPhiM(iParticle->daughter(i)->daughter(j)->pt(),iParticle->daughter(i)->daughter(j)->eta(),iParticle->daughter(i)->daughter(j)->phi(),iParticle->daughter(i)->daughter(j)->mass() );
            }
          }
        }
      }      
      if (posMuon.Pt() > 0 && negMuon.Pt() > 0) {
        double muonDeltaR = reco::deltaR2(negMuon.Eta(), negMuon.Phi(), posMuon.Eta(), posMuon.Phi());
        h2_LLpT_vs_deltaR->Fill(iParticle->pt(), muonDeltaR);
      }
    }
  }

  for (std::vector<reco::VertexCompositeCandidate>::const_iterator it = myVertexHandle->begin(); it != myVertexHandle->end(); ++it) {
    double vert_x = it->vx();
    double vert_y = it->vy();
    double vert_rho = sqrt(vert_x*vert_x + vert_y*vert_y);

    vChi2       = it->vertexChi2();
    vndof       = it->vertexNdof();
    vNormChi2   = vChi2/vndof;
    vMass       = it->mass();
    
    howManyReco++;
    numberOfReco++;

    TLorentzVector tempV2;
    tempV2.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(), vMass);
    vVS.push_back(tempV2);
    vVSRho.push_back(vert_rho);
    vRecoVertChi2.push_back(vChi2);
    h1_vertCand_mass->Fill(tempV2.M());
    h1_vertCand_rho->Fill(vert_rho);
    h2_mass_vs_chi2->Fill(vMass, vChi2);
    h2_mass_vs_ndof->Fill(vMass, vndof);
    h2_mass_vs_normalizedChi2->Fill(vMass, vNormChi2);
  }

  std::cout<<"Rho = " <<vRho.size() << std::endl;
  std::cout<<"VS = " <<vVS.size() << std::endl;
  std::cout<<"LL = " <<vLL.size() << std::endl;
  std::cout<<"Chi2 = " <<vRecoVertChi2.size() << std::endl;
  std::cout<<"*******************************************"<<std::endl;

  for (size_t j = 0; j<vLL.size(); ++j){
    for (size_t k = 0; k<vVS.size(); ++k){
      if (reco::deltaR2(vLL.at(j).Eta(), vLL.at(j).Phi(), vVS.at(k).Eta(), vVS.at(k).Phi()) <= 7e-5){
        numberOfMatched++;
        vMatchedVS.push_back(vVS.at(k));
        vMatchedRho.push_back(vRho.at(j));
        vMatchedLL.push_back(vLL.at(j));
        vMatchedVertChi2.push_back(vRecoVertChi2.at(k));
      }
    }
  }

  double ptDiff;
  double normPtDiff;
  int nLowChiHighMass = 0;

  std::cout<<"Matched Rho = " <<vMatchedRho.size() << std::endl;
  std::cout<<"Matched VS = " <<vMatchedVS.size() << std::endl;
  std::cout<<"Matched LL = " <<vMatchedLL.size() << std::endl;
  std::cout<<"Matched Chi2 = " <<vMatchedVertChi2.size() << std::endl;
  std::cout<<"Matched Rho = " <<vMatchedRho.size() << std::endl;

  for (size_t i = 0; i<vMatchedVS.size(); ++i){
    h1_matched_rho->Fill(vMatchedRho.at(i));
    h1_eff_rho->Fill(vMatchedRho.at(i));
    h1_matched_pt->Fill(vMatchedLL.at(i).Pt());
    h1_eff_pt->Fill(vMatchedLL.at(i).Pt());
    h1_matched_eta->Fill(vMatchedLL.at(i).Eta());
    h1_eff_eta->Fill(vMatchedLL.at(i).Eta());
    h1_matched_mass->Fill(vMatchedVS.at(i).M());

    h1_deltaR->Fill(reco::deltaR2(vMatchedLL.at(i).Eta(), vMatchedLL.at(i).Phi(), vMatchedVS.at(i).Eta(), vMatchedVS.at(i).Phi()));
    h1_deltaEta->Fill(fabs(vMatchedLL.at(i).Eta() - vMatchedVS.at(i).Eta()));
    h1_deltaPhi->Fill(fabs(vMatchedLL.at(i).Phi() - vMatchedVS.at(i).Phi()));

    if (fabs(vMatchedVS.at(i).M() - 2) > 0.13) {
      h1_deltaR_highMass->Fill(reco::deltaR2(vMatchedLL.at(i).Eta(), vMatchedLL.at(i).Phi(), vMatchedVS.at(i).Eta(), vMatchedVS.at(i).Phi()));
      h1_howMany_highMass->Fill(numberOfMatched);
      h1_rho_highMass->Fill(vVSRho.at(i));
    }

    h2_mass_vs_dR->Fill(vMatchedVS.at(i).M(), reco::deltaR2(vMatchedLL.at(i).Eta(), vMatchedLL.at(i).Phi(), vMatchedVS.at(i).Eta(), vMatchedVS.at(i).Phi()));

    ptDiff      = fabs(vMatchedLL.at(i).Pt() - vMatchedVS.at(i).Pt());
    normPtDiff  = ptDiff/vMatchedLL.at(i).Pt();
    h2_normal_pt_vs_chi2->Fill(normPtDiff, vMatchedVertChi2.at(i));  

    if (vMatchedVertChi2.at(i)<2) {
      h1_lowChi2_mass->Fill(vMatchedVS.at(i).M());
      if (fabs(vMatchedVS.at(i).M() - 2) > 0.13) nLowChiHighMass++;
    }

    if (normPtDiff <= 0.01) h1_chi2_1->Fill(vMatchedVertChi2.at(i));
    else if (normPtDiff <= 0.05 && normPtDiff>0.01) h1_chi2_1_5->Fill(vMatchedVertChi2.at(i));
    else h1_chi2_5->Fill(vMatchedVertChi2.at(i));
  }

  h1_nMatched->Fill(numberOfMatched);
  h1_howMany_lowChi2_highMass->Fill(nLowChiHighMass);
  h1_perEvent_eff->Fill(numberOfReco/numberOfLL);
  h1_perEvent_numberOfReco->Fill(numberOfReco);
}


// ------------ method called once each job just before starting event loop  ------------
void 
SecondaryVertices::beginJob()
{
  howManyReco = 0;
  howManyLL   = 0;

  h1_LL_rho             = fs->make<TH1F>("h1_LL_rho",             "h1_LL_rho",            100, 0, 100);
  h1_LL_pt              = fs->make<TH1F>("h1_LL_pt",              "h1_LL_pt",             100, 0, 300);
  h1_LL_eta             = fs->make<TH1F>("h1_LL_eta",             "h1_LL_eta",            100, -8, 8);  

  h2_LLpT_vs_deltaR     = fs->make<TH2F>("h2_LLpT_vs_deltaR",     "h2_LLpT_vs_deltaR",    100, 0, 300, 100, 0, 5); 

  h1_perEvent_numberOfReco      = fs->make<TH1F>("h1_perEvent_numberOfReco",      "h1_perEvent_numberOfReco",         6, 0, 6);
  h1_perEvent_eff               = fs->make<TH1F>("h1_perEvent_eff",               "h1_perEvent_eff",                  100, 0, 4);

  h1_matched_rho        = fs->make<TH1F>("h1_matched_rho",        "h1_matched_rho",       100, 0, 100);
  h1_matched_pt         = fs->make<TH1F>("h1_matched_pt",         "h1_matched_pt",        100, 0, 300);
  h1_matched_eta        = fs->make<TH1F>("h1_matched_eta",        "h1_matched_eta",       100, -8, 8);
  h1_matched_mass       = fs->make<TH1F>("h1_matched_mass",       "h1_matched_mass",      100, 1, 3);

  h1_deltaR             = fs->make<TH1F>("h1_deltaR",             "h1_deltaR",            100, 0, 0.0001);
  h1_deltaEta           = fs->make<TH1F>("h1_deltaEta",           "h1_deltaEta",          100, 0, 0.02);
  h1_deltaPhi           = fs->make<TH1F>("h1_deltaPhi",           "h1_deltaPhi",          100, 0, 0.02);

  h1_eff_rho            = fs->make<TH1F>("h1_eff_rho",            "h1_eff_tho",           100, 0, 100);
  h1_eff_pt             = fs->make<TH1F>("h1_eff_pt",             "h1_eff_pt",            100, 0, 300);
  h1_eff_eta            = fs->make<TH1F>("h1_eff_eta",            "h1_eff_eta",           100, -8, 8);

  h1_chi2_5             = fs->make<TH1F>("h1_chi2_5",             "h1_chi2_5",            100, 0, 8);
  h1_chi2_1_5           = fs->make<TH1F>("h1_chi2_1_5",           "h1_chi2_1_5",          100, 0, 8);
  h1_chi2_1             = fs->make<TH1F>("h1_chi2_1",             "h1_chi2_1",            100, 0, 8);

  h1_nMatched                   = fs->make<TH1F>("h1_nMatched",                   "h1_nMatched",                      6, 0, 6);
  h1_vertCand_mass              = fs->make<TH1F>("h1_vertCand_mass",              "h1_vertCand_mass",                 110, 0, 110);
  h1_vertCand_rho               = fs->make<TH1F>("h1_vertCand_rho",               "h1_vertCand_rho",                  100, 0, 100);
  h1_lowChi2_mass               = fs->make<TH1F>("h1_lowChi2_mass",               "h1_lowChi2_mass",                  110, 0, 110);

  h1_deltaR_highMass            = fs->make<TH1F>("h1_deltaR_highMass",            "h1_deltaR_highMass",               100, 0, 0.1);
  h1_howMany_highMass           = fs->make<TH1F>("h1_howMany_highMass",           "h1_howMany_highMass",              6, 0, 6);
  h1_howMany_lowChi2_highMass   = fs->make<TH1F>("h1_howMany_lowChi2_highMass",   "h1_howMany_lowChi2_highMass",      6, 0, 6);
  h1_rho_highMass               = fs->make<TH1F>("h1_rho_highMass",               "h1_rho_highMass",                  100, 0, 100);
  
  h2_normal_pt_vs_chi2          = fs->make<TH2F>("h2_normal_pt_vs_chi2",          "h2_normal_pt_vs_chi2",             100, 0, 1, 50, 0, 10);
  h2_mass_vs_dR                 = fs->make<TH2F>("h2_mass_vs_dR",                 "h2_mass_vs_dR",                    110, 0, 110, 50, 0, 0.1); 
  h2_mass_vs_chi2               = fs->make<TH2F>("h2_mass_vs_chi2",               "h2_mass_vs_chi2",                  110, 0, 110, 50, 0, 7);
  h2_mass_vs_ndof               = fs->make<TH2F>("h2_mass_vs_ndof",               "h2_mass_vs_ndof",                  110, 0, 110, 10, 0, 2);
  h2_mass_vs_normalizedChi2     = fs->make<TH2F>("h2_mass_vs_normalizedChi2",     "h2_mass_vs_normalizedChi2",        110, 0, 110, 50, 0, 7);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SecondaryVertices::endJob() 
{
  h1_eff_eta->Divide(h1_LL_eta);
  h1_eff_pt->Divide(h1_LL_pt);
  h1_eff_rho->Divide(h1_LL_rho);

  std::cout <<"Number of reconstructed = "<< howManyReco <<std::endl;
  std::cout <<"Number of generated LLs = "<< howManyLL <<std::endl;
  std::cout <<"Percent of reconstruced = "<< howManyReco/howManyLL <<std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
SecondaryVertices::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SecondaryVertices::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SecondaryVertices::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SecondaryVertices::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SecondaryVertices::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SecondaryVertices);
