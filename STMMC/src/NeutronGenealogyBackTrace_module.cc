// Adapted from VirtualDetectorTree_module.cc
// For StepPointMCs in virtualdetectors, find all neutrons and trace back to their origin.
// Including two trees:
// 1. Tree with neutrons in given virtualdetector, "ttree_neutron"
//    with index in the file, hit time, PDG ID, virtualdetector ID, creation code, kinetic energy, positions x, y, and z
//    in branches "index", "time", "virtualdetectorId", "pdgId", "creationCode", "E", "x", "y", and "z" respectively.
// 2. Tree with the neutron genealogy, "ttree_neutron_genealogy"
//    with corresponding neutron index in the file, then hit time, PDG ID, creation code,kinetic energy, positions x, y, and z of each parent particle
//    in branches "index", "time", "pdgId", "creationCode", "E", "x", "y", and "z" respectively.
// Adapted by: Yongyi Wu

// stdlib includes
#include <cmath>
#include <iostream>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"


typedef cet::map_vector_key key_type;
typedef unsigned long VolumeId_type;

namespace mu2e {
  class NeutronGenealogyBackTrace : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> SimParticlemvTag{Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit NeutronGenealogyBackTrace(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      GlobalConstantsHandle<ParticleDataList> pdt;

      double mass = 0.0;
      // Data members for the ttree_neutron TTree
      int ttree_neutron_index = -1;
      double ttree_neutron_time = 0.0;
      VolumeId_type ttree_neutron_virtualdetectorId = 0;
      int ttree_neutron_creationCode = 0;
      int ttree_neutron_pdgId = 0;
      double ttree_neutron_E = 0.0;
      double ttree_neutron_x = 0.0;
      double ttree_neutron_y = 0.0;
      double ttree_neutron_z = 0.0;
      double ttree_neutron_startx = 0.0;
      double ttree_neutron_starty = 0.0;
      double ttree_neutron_startz = 0.0;
      double ttree_neutron_starttime = 0.0;
      // Data members for the ttree_neutron_genealogy TTree
      int ttree_neutron_genealogy_index = -1;
      double ttree_neutron_genealogy_startGlobalTime = 0.0;
      double ttree_neutron_genealogy_endGlobalTime = 0.0;
      int ttree_neutron_genealogy_pdgId = 0;
      int ttree_neutron_genealogy_creationCode = 0;
      double ttree_neutron_genealogy_E = 0.0;
      double ttree_neutron_genealogy_x = 0.0;
      double ttree_neutron_genealogy_y = 0.0;
      double ttree_neutron_genealogy_z = 0.0;

      TTree* ttree_neutron;
      TTree* ttree_neutron_genealogy;

      std::map<int, int> pdgIds; // <id, count>
  };

  NeutronGenealogyBackTrace::NeutronGenealogyBackTrace(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())) {
      art::ServiceHandle<art::TFileService> tfs;
      ttree_neutron = tfs->make<TTree>("ttree_neutron", "Neutrons in given virtual detector");
      ttree_neutron->Branch("index", &ttree_neutron_index, "index/I");
      ttree_neutron->Branch("time", &ttree_neutron_time, "time/D"); // ns
      ttree_neutron->Branch("virtualdetectorId", &ttree_neutron_virtualdetectorId, "virtualdetectorId/l");
      ttree_neutron->Branch("creationCode", &ttree_neutron_creationCode, "creationCode/I");
      ttree_neutron->Branch("pdgId", &ttree_neutron_pdgId, "pdgId/I");
      ttree_neutron->Branch("E", &ttree_neutron_E, "E/D"); // MeV
      ttree_neutron->Branch("x", &ttree_neutron_x, "x/D"); // mm
      ttree_neutron->Branch("y", &ttree_neutron_y, "y/D"); // mm
      ttree_neutron->Branch("z", &ttree_neutron_z, "z/D"); // mm
      ttree_neutron->Branch("startx", &ttree_neutron_startx, "startx/D"); // mm
      ttree_neutron->Branch("starty", &ttree_neutron_starty, "starty/D"); // mm
      ttree_neutron->Branch("startz", &ttree_neutron_startz, "startz/D"); // mm
      ttree_neutron->Branch("starttime", &ttree_neutron_starttime, "starttime/D"); // ns

      ttree_neutron_genealogy = tfs->make<TTree>("ttree_neutron_genealogy", "Neutron genealogy");
      ttree_neutron_genealogy->Branch("index", &ttree_neutron_genealogy_index, "index/I");
      ttree_neutron_genealogy->Branch("startGlobalTime", &ttree_neutron_genealogy_startGlobalTime, "startGlobalTime/D"); // ns
      ttree_neutron_genealogy->Branch("endGlobalTime", &ttree_neutron_genealogy_endGlobalTime, "endGlobalTime/D"); // ns
      ttree_neutron_genealogy->Branch("pdgId", &ttree_neutron_genealogy_pdgId, "pdgId/I");
      ttree_neutron_genealogy->Branch("creationCode", &ttree_neutron_genealogy_creationCode, "creationCode/I");
      ttree_neutron_genealogy->Branch("E", &ttree_neutron_genealogy_E, "E/D"); // MeV
      ttree_neutron_genealogy->Branch("x", &ttree_neutron_genealogy_x, "x/D"); // mm
      ttree_neutron_genealogy->Branch("y", &ttree_neutron_genealogy_y, "y/D"); // mm
      ttree_neutron_genealogy->Branch("z", &ttree_neutron_genealogy_z, "z/D"); // mm
    };

  void NeutronGenealogyBackTrace::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    if (StepPointMCs.empty())
      return;
    auto const& SimParticles = event.getProduct(SimParticlemvToken);
    if (SimParticles.empty())
      return;

    // Loop over all VD hits
    for (const StepPointMC& step : StepPointMCs) {
      // Get the associated particle
      auto particle = step.simParticle();

      // Extract the parameters
      if (particle->pdgId() != PDGCode::n0)
        continue;
      ttree_neutron_index ++;
      ttree_neutron_time = step.time();
      ttree_neutron_virtualdetectorId = step.virtualDetectorId();
      ttree_neutron_creationCode = particle->creationCode();
      ttree_neutron_pdgId = particle->pdgId();
      ttree_neutron_x = step.position().x();
      ttree_neutron_y = step.position().y();
      ttree_neutron_z = step.position().z();
      ttree_neutron_startx = particle->startPosition().x();
      ttree_neutron_starty = particle->startPosition().y();
      ttree_neutron_startz = particle->startPosition().z();
      ttree_neutron_starttime = particle->startGlobalTime();
      mass = pdt->particle(ttree_neutron_pdgId).mass();
      ttree_neutron_E = std::sqrt(step.momentum().mag2() + mass * mass) - mass; // Subtract the rest mass
      if (ttree_neutron_E < 0)
        throw cet::exception("LogicError", "Energy is negative");
      ttree_neutron->Fill();

      // Fill the genealogy tree
      while (particle->hasParent()) {
        particle = particle->parent();
        ttree_neutron_genealogy_index = ttree_neutron_index;
        ttree_neutron_genealogy_startGlobalTime = particle->startGlobalTime();
        ttree_neutron_genealogy_endGlobalTime = particle->endGlobalTime();
        ttree_neutron_genealogy_pdgId = particle->pdgId();
        ttree_neutron_genealogy_creationCode = particle->creationCode();
        ttree_neutron_genealogy_x = particle->startPosition().x();
        ttree_neutron_genealogy_y = particle->startPosition().y();
        ttree_neutron_genealogy_z = particle->startPosition().z();
        mass = pdt->particle(ttree_neutron_genealogy_pdgId).mass();
        ttree_neutron_genealogy_E = std::sqrt(particle->startMomentum().vect().mag2() + mass * mass) - mass; // Subtract the rest mass
        if (ttree_neutron_genealogy_E < 0)
          throw cet::exception("LogicError", "Energy is negative");
        ttree_neutron_genealogy->Fill();
      }

      // Generate the data summary
      if (pdgIds.find(ttree_neutron_pdgId) != pdgIds.end())
        pdgIds[ttree_neutron_pdgId] += 1;
      else
        pdgIds.emplace(std::make_pair(ttree_neutron_pdgId, 1));
    };

    return;
  };

  void NeutronGenealogyBackTrace::endJob() {
    mf::LogInfo log("Virtual detector tree summary");
    log << "========= Data summary =========\n";
    for (auto part : pdgIds)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "================================\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::NeutronGenealogyBackTrace)
