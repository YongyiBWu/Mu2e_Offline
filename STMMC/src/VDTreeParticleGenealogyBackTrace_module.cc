// Adapted from VirtualDetectorTree_module.cc
// For StepPointMCs in virtualdetectors, find all particles and trace back to their origin.
// Including two trees:
// 1. Tree with particles in given virtualdetector, "ttree"
//    with index in the file, hit time, PDG ID, virtualdetector ID, creation code, kinetic energy, positions x, y, and z
//    in branches "index", "time", "virtualdetectorId", "pdgId", "creationCode", "E", "x", "y", and "z" respectively.
// 2. Tree with the particle genealogy, "ttree_genealogy"
//    with corresponding particle index in the file, then hit time, PDG ID, creation code,kinetic energy, positions x, y, and z of each parent particle
//    in branches "index", "startGlobalTime", "endGlobalTime", "pdgId", "creationCode", "E", "x", "y", and "z" respectively.
// 3. Tree with the particle generation information, "ttree_source"
//    with the info of the first particle in the genealogy chain that has the same pdgId as the particle in the virtualdetector
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
  class VDTreeParticleGenealogyBackTrace : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> SimParticlemvTag{Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VDTreeParticleGenealogyBackTrace(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      GlobalConstantsHandle<ParticleDataList> pdt;

      double mass = 0.0;
      bool ttree_source_filled = false;
      // Data members for the ttree TTree
      int ttree_index = -1;
      double ttree_time = 0.0;
      VolumeId_type ttree_virtualdetectorId = 0;
      int ttree_creationCode = 0;
      int ttree_pdgId = 0;
      double ttree_E = 0.0;
      double ttree_x = 0.0;
      double ttree_y = 0.0;
      double ttree_z = 0.0;
      double ttree_px = 0.0;
      double ttree_py = 0.0;
      double ttree_pz = 0.0;
      double ttree_startx = 0.0;
      double ttree_starty = 0.0;
      double ttree_startz = 0.0;
      double ttree_startpx = 0.0;
      double ttree_startpy = 0.0;
      double ttree_startpz = 0.0;
      double ttree_starttime = 0.0;
      // Data members for the ttree_genealogy TTree
      int ttree_genealogy_index = -1;
      double ttree_genealogy_startGlobalTime = 0.0;
      double ttree_genealogy_endGlobalTime = 0.0;
      int ttree_genealogy_pdgId = 0;
      int ttree_genealogy_creationCode = 0;
      double ttree_genealogy_E = 0.0;
      double ttree_genealogy_startx = 0.0;
      double ttree_genealogy_starty = 0.0;
      double ttree_genealogy_startz = 0.0;
      double ttree_genealogy_startpx = 0.0;
      double ttree_genealogy_startpy = 0.0;
      double ttree_genealogy_startpz = 0.0;
      double ttree_genealogy_endx = 0.0;
      double ttree_genealogy_endy = 0.0;
      double ttree_genealogy_endz = 0.0;
      double ttree_genealogy_endpx = 0.0;
      double ttree_genealogy_endpy = 0.0;
      double ttree_genealogy_endpz = 0.0;
      // Data members for the ttree_source TTree
      int ttree_source_index = -1;
      double ttree_source_starttime = 0.0;
      double ttree_source_parentendtime = 0.0;
      int ttree_source_creationCode = 0;
      int ttree_source_parentpdgId = 0;
      double ttree_source_x = 0.0;
      double ttree_source_y = 0.0;
      double ttree_source_z = 0.0;

      TTree* ttree;
      TTree* ttree_genealogy;
      TTree* ttree_source;

      std::map<int, int> pdgIds; // <id, count>
  };

  VDTreeParticleGenealogyBackTrace::VDTreeParticleGenealogyBackTrace(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())) {
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>("ttree", "Particles in given virtual detector");
      ttree->Branch("index", &ttree_index, "index/I");
      ttree->Branch("time", &ttree_time, "time/D"); // ns
      ttree->Branch("virtualdetectorId", &ttree_virtualdetectorId, "virtualdetectorId/l");
      ttree->Branch("creationCode", &ttree_creationCode, "creationCode/I");
      ttree->Branch("pdgId", &ttree_pdgId, "pdgId/I");
      ttree->Branch("E", &ttree_E, "E/D"); // MeV
      ttree->Branch("x", &ttree_x, "x/D"); // mm
      ttree->Branch("y", &ttree_y, "y/D"); // mm
      ttree->Branch("z", &ttree_z, "z/D"); // mm
      ttree->Branch("px", &ttree_px, "px/D"); // MeV/c
      ttree->Branch("py", &ttree_py, "py/D"); // MeV/c
      ttree->Branch("pz", &ttree_pz, "pz/D"); // MeV/c
      ttree->Branch("startx", &ttree_startx, "startx/D"); // mm
      ttree->Branch("starty", &ttree_starty, "starty/D"); // mm
      ttree->Branch("startz", &ttree_startz, "startz/D"); // mm
      ttree->Branch("startpx", &ttree_startpx, "startpx/D"); // MeV/c
      ttree->Branch("startpy", &ttree_startpy, "startpy/D"); // MeV/c
      ttree->Branch("startpz", &ttree_startpz, "startpz/D"); // MeV/c
      ttree->Branch("starttime", &ttree_starttime, "starttime/D"); // ns

      ttree_genealogy = tfs->make<TTree>("ttree_genealogy", "Particle genealogy");
      ttree_genealogy->Branch("index", &ttree_genealogy_index, "index/I");
      ttree_genealogy->Branch("startGlobalTime", &ttree_genealogy_startGlobalTime, "startGlobalTime/D"); // ns
      ttree_genealogy->Branch("endGlobalTime", &ttree_genealogy_endGlobalTime, "endGlobalTime/D"); // ns
      ttree_genealogy->Branch("pdgId", &ttree_genealogy_pdgId, "pdgId/I");
      ttree_genealogy->Branch("creationCode", &ttree_genealogy_creationCode, "creationCode/I");
      ttree_genealogy->Branch("E", &ttree_genealogy_E, "E/D"); // MeV
      ttree_genealogy->Branch("startx", &ttree_genealogy_startx, "startx/D"); // mm
      ttree_genealogy->Branch("starty", &ttree_genealogy_starty, "starty/D"); // mm
      ttree_genealogy->Branch("startz", &ttree_genealogy_startz, "startz/D"); // mm
      ttree_genealogy->Branch("endx", &ttree_genealogy_endx, "endx/D"); // mm
      ttree_genealogy->Branch("endy", &ttree_genealogy_endy, "endy/D"); // mm
      ttree_genealogy->Branch("endz", &ttree_genealogy_endz, "endz/D"); // mm
      ttree_genealogy->Branch("startpx", &ttree_genealogy_startpx, "startpx/D"); // MeV/c
      ttree_genealogy->Branch("startpy", &ttree_genealogy_startpy, "startpy/D"); // MeV/c
      ttree_genealogy->Branch("startpz", &ttree_genealogy_startpz, "startpz/D"); // MeV/c
      ttree_genealogy->Branch("endpx", &ttree_genealogy_endpx, "endpx/D"); // MeV/c
      ttree_genealogy->Branch("endpy", &ttree_genealogy_endpy, "endpy/D"); // MeV/c
      ttree_genealogy->Branch("endpz", &ttree_genealogy_endpz, "endpz/D"); // MeV/c

      ttree_source = tfs->make<TTree>("ttree_source", "Particle source information");
      ttree_source->Branch("index", &ttree_source_index, "index/I");
      ttree_source->Branch("starttime", &ttree_source_starttime, "starttime/D"); // ns
      ttree_source->Branch("parentendtime", &ttree_source_parentendtime, "parentendtime/D"); // ns
      ttree_source->Branch("creationCode", &ttree_source_creationCode, "creationCode/I");
      ttree_source->Branch("parentpdgId", &ttree_source_parentpdgId, "parentpdgId/I");
      ttree_source->Branch("x", &ttree_source_x, "x/D"); // mm
      ttree_source->Branch("y", &ttree_source_y, "y/D"); // mm
      ttree_source->Branch("z", &ttree_source_z, "z/D"); // mm
    };

  void VDTreeParticleGenealogyBackTrace::analyze(const art::Event& event) {
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
      ttree_index ++;
      ttree_time = step.time();
      ttree_virtualdetectorId = step.virtualDetectorId();
      ttree_creationCode = particle->creationCode();
      ttree_pdgId = particle->pdgId();
      ttree_x = step.position().x();
      ttree_y = step.position().y();
      ttree_z = step.position().z();
      ttree_px = step.momentum().x();
      ttree_py = step.momentum().y();
      ttree_pz = step.momentum().z();
      ttree_startx = particle->startPosition().x();
      ttree_starty = particle->startPosition().y();
      ttree_startz = particle->startPosition().z();
      ttree_startpx = particle->startMomentum().x();
      ttree_startpy = particle->startMomentum().y();
      ttree_startpz = particle->startMomentum().z();
      ttree_starttime = particle->startGlobalTime();
      mass = pdt->particle(ttree_pdgId).mass();
      ttree_E = std::sqrt(step.momentum().mag2() + mass * mass) - mass; // Subtract the rest mass
      if (ttree_E < 0)
        throw cet::exception("LogicError", "Energy is negative");
      ttree->Fill();

      ttree_source_filled = false;
      ttree_source_parentpdgId = particle->pdgId();
      ttree_source_x = particle->startPosition().x();
      ttree_source_y = particle->startPosition().y();
      ttree_source_z = particle->startPosition().z();
      ttree_source_starttime = particle->startGlobalTime();
      ttree_source_creationCode = particle->creationCode();

      // Fill the genealogy tree
      while (particle->hasParent()) {
        particle = particle->parent();
        ttree_genealogy_index = ttree_index;
        ttree_genealogy_startGlobalTime = particle->startGlobalTime();
        ttree_genealogy_endGlobalTime = particle->endGlobalTime();
        ttree_genealogy_pdgId = particle->pdgId();
        ttree_genealogy_creationCode = particle->creationCode();
        ttree_genealogy_startx = particle->startPosition().x();
        ttree_genealogy_starty = particle->startPosition().y();
        ttree_genealogy_startz = particle->startPosition().z();
        ttree_genealogy_startpx = particle->startMomentum().x();
        ttree_genealogy_startpy = particle->startMomentum().y();
        ttree_genealogy_startpz = particle->startMomentum().z();
        ttree_genealogy_endx = particle->endPosition().x();
        ttree_genealogy_endy = particle->endPosition().y();
        ttree_genealogy_endz = particle->endPosition().z();
        ttree_genealogy_endpx = particle->endMomentum().x();
        ttree_genealogy_endpy = particle->endMomentum().y();
        ttree_genealogy_endpz = particle->endMomentum().z();
        mass = pdt->particle(ttree_genealogy_pdgId).mass();
        ttree_genealogy_E = std::sqrt(particle->startMomentum().vect().mag2() + mass * mass) - mass; // Subtract the rest mass
        if (ttree_genealogy_E < 0)
          throw cet::exception("LogicError", "Energy is negative");
        ttree_genealogy->Fill();

        if (!ttree_source_filled) {
          if (ttree_source_parentpdgId == particle->pdgId()) {
            ttree_source_x = particle->startPosition().x();
            ttree_source_y = particle->startPosition().y();
            ttree_source_z = particle->startPosition().z();
            ttree_source_starttime = particle->startGlobalTime();
            ttree_source_creationCode = particle->creationCode();
          }
          else {
            ttree_source_parentpdgId = particle->pdgId();
            ttree_source_parentendtime = particle->endGlobalTime();
            ttree_source_filled = true;
            ttree_source_index = ttree_index;
            ttree_source->Fill();
          }
        }
      }

      // Generate the data summary
      if (pdgIds.find(ttree_pdgId) != pdgIds.end())
        pdgIds[ttree_pdgId] += 1;
      else
        pdgIds.emplace(std::make_pair(ttree_pdgId, 1));
    };

    return;
  };

  void VDTreeParticleGenealogyBackTrace::endJob() {
    mf::LogInfo log("Virtual detector tree summary");
    log << "========= Data summary =========\n";
    for (auto part : pdgIds)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "================================\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VDTreeParticleGenealogyBackTrace)
