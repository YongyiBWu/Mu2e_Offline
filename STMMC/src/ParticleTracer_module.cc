// ParticleTracer_module.cc
// Traces particle trajectories that reach a set of StepPointMCs (e.g. a virtual detector),
// following their MC genealogy back to the primary, and writes the trajectories to a ROOT
// TTree.
//
// Workflow:
//   1) Filter StepPointMCs by a configurable set of PDG codes (the same filtering pattern
//      as VirtualDetectorTree_module.cc). This is the seed: "particles that reached here".
//   2) For each matching StepPointMC, take its SimParticle and walk parent() up to the
//      primary, collecting the matched particle AND every ancestor.
//   3) For each such particle, look up its MCTrajectory (keyed by art::Ptr<SimParticle>) and
//      write it as one TTree entry -- keeping mother trajectories even when the mother's PDG
//      does not match the filter. Ancestors whose trajectory was not stored (did not pass
//      the G4 trajectory cuts) are skipped.
//
// The filter therefore applies to StepPointMCs, not to trajectories. Each entry stores the
// per-point coordinates as std::vector<double> branches plus genealogy (immediate parent and
// the full ancestor chain up to the primary). A single MCTrajectory spans one SimParticle
// (creation vertex to stop/exit); ancestors' points live in their own (separate) entries.
//
// The SimParticleCollection (SimParticlemvTag) is consumed so it is loaded; the working
// art::Ptr<SimParticle> used for the genealogy walk and trajectory lookup comes from
// StepPointMC::simParticle(). This mirrors the STMMC package convention
// (VirtualDetectorTree/HPGeTree both take a SimParticlemvTag).
//
// Modelled on VirtualDetectorTree_module.cc (same package). For a comparison of how the
// MCTrajectoryCollection is read, see Offline/Print/src/MCTrajectoryPrinter.cc and its
// driver Offline/Print/src/PrintModule_module.cc, and Analyses/src/ReadMCTrajectories_module.cc.
//
// Yongyi Wu, Jul. 2026

// stdlib includes
#include <iostream>
#include <map>
#include <set>
#include <vector>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {
  class ParticleTracer : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{  Name("StepPointMCsTag"),   Comment("Tag identifying the StepPointMCs to seed from (the PDG filter applies here)")};
        fhicl::Atom<art::InputTag> SimParticlemvTag{ Name("SimParticlemvTag"),  Comment("Tag identifying the SimParticleCollection (loaded so the genealogy Ptrs resolve)")};
        fhicl::Atom<art::InputTag> MCTrajectoriesTag{Name("MCTrajectoriesTag"), Comment("Tag identifying the MCTrajectoryCollection"), art::InputTag("g4run")};
        fhicl::Sequence<int> pdgIDs{                 Name("pdgIDs"),             Comment("PDG codes of the StepPointMC particles to trace"), std::vector<int>{2112}};
        fhicl::Atom<bool> keepAllPdgIds{             Name("keepAllPdgIds"),      Comment("If true (or pdgIDs is empty), trace every StepPointMC particle"), false};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit ParticleTracer(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      // Fill the point/genealogy branches for one SimParticle's trajectory and Fill() the
      // tree. Returns true if the particle had a stored trajectory (and an entry was written).
      bool writeTrajectory(const art::Ptr<SimParticle>& simPtr,
                           const MCTrajectoryCollection& coll,
                           bool matchedStep);

      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      art::ProductToken<MCTrajectoryCollection> MCTrajectoriesToken;
      std::set<int> pdgIDs;
      bool keepAllPdgIds = false;
      GlobalConstantsHandle<ParticleDataList> pdt;

      // TTree + branch buffers (one entry per trajectory)
      TTree* ttree = nullptr;
      unsigned int run = 0, subRun = 0, event = 0, simId = 0;
      int pdgId = 0, nPoints = 0;
      bool matched = false; // true if this particle is a filtered StepPointMC particle (vs. an ancestor)
      std::vector<double> x, y, z, t, kE, E;
      // Genealogy
      int parentSimId = -1, parentPdgId = 0, creationCode = 0;
      bool isPrimary = false;
      std::vector<int> ancestorSimIds; // parent-first, ultimate primary last
      std::vector<int> ancestorPdgIds;

      std::map<int, int> vdHitCounts;          // <pdgId, number of matching StepPointMC (VD) hits>
      std::set<unsigned int> writtenThisEvent; // simIds already written in the current event
  };

  ParticleTracer::ParticleTracer(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())),
    MCTrajectoriesToken(consumes<MCTrajectoryCollection>(conf().MCTrajectoriesTag())),
    keepAllPdgIds(conf().keepAllPdgIds()) {
      const std::vector<int>& ids = conf().pdgIDs();
      pdgIDs.insert(ids.begin(), ids.end());

      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>("ttree", "MCTrajectory ttree");
      ttree->Branch("run", &run, "run/i");
      ttree->Branch("subRun", &subRun, "subRun/i");
      ttree->Branch("event", &event, "event/i");
      ttree->Branch("simId", &simId, "simId/i");
      ttree->Branch("pdgId", &pdgId, "pdgId/I");
      ttree->Branch("matched", &matched, "matched/O"); // true: seed StepPointMC particle; false: ancestor
      ttree->Branch("nPoints", &nPoints, "nPoints/I");
      ttree->Branch("x", &x);   // mm
      ttree->Branch("y", &y);   // mm
      ttree->Branch("z", &z);   // mm
      ttree->Branch("t", &t);   // ns
      ttree->Branch("kE", &kE); // MeV
      ttree->Branch("E", &E);   // MeV, total energy = kE + mass
      // Genealogy
      ttree->Branch("parentSimId", &parentSimId, "parentSimId/I"); // -1 if primary/unavailable
      ttree->Branch("parentPdgId", &parentPdgId, "parentPdgId/I");
      ttree->Branch("creationCode", &creationCode, "creationCode/I");
      ttree->Branch("isPrimary", &isPrimary, "isPrimary/O");
      ttree->Branch("ancestorSimIds", &ancestorSimIds); // parent-first, primary last
      ttree->Branch("ancestorPdgIds", &ancestorPdgIds);
    };

  bool ParticleTracer::writeTrajectory(const art::Ptr<SimParticle>& simPtr,
                                       const MCTrajectoryCollection& coll,
                                       bool matchedStep) {
    if (!simPtr.isAvailable()) {
      mf::LogWarning("ParticleTracer")
        << "SimParticle Ptr (key " << simPtr.key() << ") is unavailable"
        << (matchedStep ? " (matched StepPointMC particle)" : " (ancestor)")
        << "; skipping its trajectory. Genealogy may be truncated in this file.";
      return false;
    }

    // Only write each particle once per event.
    if (!writtenThisEvent.insert(simPtr.key()).second)
      return false;

    // Look up this particle's trajectory; absent if it didn't pass the G4 trajectory cuts.
    auto traj_it = coll.find(simPtr);
    if (traj_it == coll.end()) {
      mf::LogWarning("ParticleTracer")
        << "No MCTrajectory found for SimParticle (key " << simPtr.key()
        << ", pdgId " << simPtr->pdgId() << ")"
        << (matchedStep ? " (matched StepPointMC particle)" : " (ancestor)")
        << "; skipping. It likely did not pass the G4 trajectory cuts.";
      return false;
    }
    const MCTrajectory& traj = traj_it->second;

    const SimParticle& sim = *simPtr;
    pdgId = sim.pdgId();
    matched = matchedStep;
    const double mass = pdt->particle(pdgId).mass();

    simId = simPtr.key();
    nPoints = traj.size();
    x.clear(); y.clear(); z.clear(); t.clear(); kE.clear(); E.clear();
    for (auto const& p : traj.points()) {
      x.push_back(p.x());
      y.push_back(p.y());
      z.push_back(p.z());
      t.push_back(p.t());
      kE.push_back(p.kineticEnergy());
      E.push_back(p.kineticEnergy() + mass);
    }

    // Genealogy: immediate parent + full ancestor chain walked to the primary.
    isPrimary = sim.isPrimary();
    creationCode = sim.creationCode().id();
    parentSimId = -1;
    parentPdgId = 0;
    ancestorSimIds.clear();
    ancestorPdgIds.clear();
    art::Ptr<SimParticle> anc = sim.parent();
    bool first = true;
    while (anc.isAvailable()) {
      if (first) {
        parentSimId = anc.key();
        parentPdgId = anc->pdgId();
        first = false;
      }
      ancestorSimIds.push_back(anc.key());
      ancestorPdgIds.push_back(anc->pdgId());
      anc = anc->parent();
    }

    ttree->Fill();
    return true;
  };

  void ParticleTracer::analyze(const art::Event& evt) {
    auto const& steps = evt.getProduct(StepPointMCsToken);
    if (steps.empty())
      return;
    auto const& SimParticles = evt.getProduct(SimParticlemvToken);
    if (SimParticles.empty())
      return;
    auto const& trajectories = evt.getProduct(MCTrajectoriesToken);
    if (trajectories.empty())
      return;

    run = evt.id().run();
    subRun = evt.id().subRun();
    event = evt.id().event();
    writtenThisEvent.clear();

    for (const StepPointMC& step : steps) {
      const art::Ptr<SimParticle>& simPtr = step.simParticle();
      if (!simPtr.isAvailable()) {
        mf::LogWarning("ParticleTracer")
          << "StepPointMC SimParticle Ptr (key " << simPtr.key()
          << ") is unavailable; skipping this hit.";
        continue;
      }

      // The PDG filter applies to the StepPointMC particle (the seed).
      const int stepPdgId = simPtr->pdgId();
      if (!keepAllPdgIds && !pdgIDs.empty() &&
          pdgIDs.find(stepPdgId) == pdgIDs.end())
        continue;

      // Count this VD hit (one per matching StepPointMC, keyed by the hit particle's PDG;
      // not deduplicated -- this is the raw hit multiplicity on the StepPointMC volume).
      ++vdHitCounts[stepPdgId];

      // Write the matched particle's trajectory, then every ancestor up to the primary
      // (regardless of the ancestor's PDG).
      writeTrajectory(simPtr, trajectories, true);
      art::Ptr<SimParticle> anc = simPtr->parent();
      while (anc.isAvailable()) {
        writeTrajectory(anc, trajectories, false);
        anc = anc->parent();
      }
    }
    return;
  };

  void ParticleTracer::endJob() {
    mf::LogInfo log("ParticleTracer summary");
    log << "========= StepPointMC (VD) hits per PDG ID =========\n";
    for (auto const& part : vdHitCounts)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "===================================================\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ParticleTracer)
