// include files
#include "PPSKinFilter.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

#include "FWCore/Framework/interface/MakerMacros.h"

// fill descriptions
//
void PPSKinFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("jetInputTag", edm::InputTag("hltAK4PFJetsCorrected"))
    ->setComment("input tag of the jet track collection");
  desc.add<edm::InputTag>("forwardProtonInputTag", edm::InputTag("ctppsProtons", "singleRP"))
    ->setComment("input tag of the forward proton collection");

  desc.add<std::string>("lhcInfoLabel", std::string(""))->setComment("label used for LHCInfo");
  
  desc.add<double>("maxDiffxi", 1.)->setComment("maximum relative deviation of RP xi from dijet xi");
  desc.add<double>("maxDiffm", 1.)->setComment("maximum relative deviation of RP m from dijet m");
  desc.add<double>("maxDiffy", 1.)->setComment("maximum absolute deviation of RP y from dijet y");

  desc.add<unsigned int>("nJets", 2)->setComment("number of jets to be used. Only 2 jets can be used for the xi matching option");

  desc.add<bool>("do_xi", true)->setComment("toggle to trigger on xi deviation");
  desc.add<bool>("do_my", false)->setComment("toggle to trigger on m,y deviation");

  descriptions.add("hltCTPPSKinematicFilter", desc);
  return;
}

// destructor and constructor
//
PPSKinFilter::~PPSKinFilter() = default;

PPSKinFilter::PPSKinFilter(const edm::ParameterSet& iConfig)
  :
  jetInputTag_(iConfig.getParameter<edm::InputTag>("jetInputTag")),
  jet_token_(consumes<reco::PFJetCollection>(jetInputTag_)),

  forwardProtonInputTag_(iConfig.getParameter<edm::InputTag>("forwardProtonInputTag")), 
  recoProtonSingleRPToken_(consumes<std::vector<reco::ForwardProton>>(forwardProtonInputTag_)),

  lhcInfoLabel_(iConfig.getParameter<std::string>("lhcInfoLabel")),

  maxDiffxi_(iConfig.getParameter<double>("maxDiffxi")),
  maxDiffm_(iConfig.getParameter<double>("maxDiffm")),
  maxDiffy_(iConfig.getParameter<double>("maxDiffy")),

  n_jets_(iConfig.getParameter<unsigned int>("nJets")),
  
  do_xi_(iConfig.getParameter<bool>("do_xi")),
  do_my_(iConfig.getParameter<bool>("do_my"))
{
}

// member functions
//
bool PPSKinFilter::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  
  float min45 = 1000.;
  float min56 = 1000.;
  float xi45 = 0.;
  float xi56 = 0.;

  bool kin_pass = false; 

  edm::ESHandle<LHCInfo> hLHCInfo;
  iSetup.get<LHCInfoRcd>().get(lhcInfoLabel_, hLHCInfo);
  float sqs = 2.*hLHCInfo->energy();
  
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jet_token_, jets); // get jet collection

  edm::Handle<std::vector<reco::ForwardProton>> recoSingleRPProtons;
  iEvent.getByToken(recoProtonSingleRPToken_, recoSingleRPProtons); // get RP proton collection

  if( jets->size() < n_jets_) return false; // cond for nr jets
  
  if(do_xi_) {

    float sum45 = 0, sum56 = 0, xi;
    
    for(unsigned int i = 0; i < n_jets_; i++){
      sum45 += (*jets)[i].energy() + (*jets)[i].pz();
      sum56 += (*jets)[i].energy() - (*jets)[i].pz();
    }

    xi45 = sum45/sqs; // get arm45 xi for n leading-pT jets
    xi56 = sum56/sqs; // get arm56 xi for n leading-pT jets

    for (const auto & proton : *recoSingleRPProtons) // cycle over proton tracks
      {
	if (proton.validFit()) // Check that the track fit is valid
	  {
	    xi = proton.xi(); // Get the proton xi
	       
	    CTPPSDetId rpId((*proton.contributingLocalTracks().begin())->getRPId()); // get RP ID (rpId.arm() is 0 for 45 and 1 for 56)

	    if(rpId.arm() == 0 && abs(xi-xi45) < min45) 
	      min45 = abs(xi-xi45);
	    if(rpId.arm() == 1 && abs(xi-xi56) < min56)
	      min56 = abs(xi-xi56);
	  }
      }
  }

  if(do_my_) {

    float m, y;
  
    // get the mass and rap of the n jets
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > j_sum;
    for(unsigned int i = 0; i < n_jets_; i++)
      j_sum = j_sum+(*jets)[i].p4();
    
    float mjet = j_sum.M();
    float yjet = j_sum.Rapidity();

    for (const auto & proton1 : *recoSingleRPProtons) // cycle over first RP (only arm45)
      {
	if(kin_pass) break; // leave cycle if a combination that meets cuts is found
	
	if (proton1.validFit()) {	       
	  CTPPSDetId rpId1((*proton1.contributingLocalTracks().begin())->getRPId()); // get RP ID (rpId.arm() is 0 for 45 and 1 for 56)	    
	  if(rpId1.arm() == 0) {
	    xi45 = proton1.xi();
	      
	    for (const auto & proton2 : *recoSingleRPProtons) // cycle over second RP (only arm56)
	      {
		if (proton2.validFit())  {	       
		  CTPPSDetId rpId2((*proton2.contributingLocalTracks().begin())->getRPId()); 
		  if(rpId2.arm() == 1) {
		    xi56 = proton2.xi();
		      
		    m = sqs * sqrt( xi45 * xi56);
		    y = 0.5 * log( xi45 / xi56);
		    if((abs(m-mjet)/mjet < maxDiffm_ || maxDiffm_ <= 0) && (abs(y-yjet) < maxDiffy_ || maxDiffy_ <= 0)) {
		      kin_pass = true;
		      break;
		    }
		  }
		}
	      }
	  }
	}
      }
  }

  if( ( min56 / xi56 > maxDiffxi_ || min45 / xi45 > maxDiffxi_ )  && maxDiffxi_ > 0 && do_xi_) return false; // cond for xi matching

  if(!kin_pass && do_my_ ) return false; // cond for m,y matching
  
  return true; // if none of the conds are met, event has passed the trigger
}

DEFINE_FWK_MODULE(PPSKinFilter);
