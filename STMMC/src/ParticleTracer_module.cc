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

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"

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
        fhicl::Atom<bool> keepParticlesWithPositivePz{ Name("keepParticlesWithPositivePz"), Comment("If true, only seed from StepPointMCs whose momentum pz > 0"), false};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit ParticleTracer(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      typedef cet::map_vector_key key_type;
      typedef std::map<key_type, const MCTrajectory*> TrajIndex;

      // Fill the point/genealogy branches for one SimParticle's trajectory and Fill() the
      // tree. The trajectory is matched by SimParticle key (map_vector_key) rather than by
      // art::Ptr, since the StepPointMC/SimParticle collection and the MCTrajectory
      // collection are generally different products (e.g. compressed vs. g4run) whose Ptrs
      // do not compare equal. Genealogy is walked via the SimParticle parent art::Ptrs.
      // Outcome of an attempt to write one particle's trajectory.
      enum class WriteResult {
        kWritten,       // trajectory found and a TTree entry was filled
        kNoTrajectory,  // no MCTrajectory stored for this particle (failed G4 cuts / key mismatch)
        kDuplicate      // this particle was already written earlier in the event
      };
      WriteResult writeTrajectory(const art::Ptr<SimParticle>& simPtr,
                                  const TrajIndex& trajIndex,
                                  bool matchedStep);

      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      art::ProductToken<MCTrajectoryCollection> MCTrajectoriesToken;
      std::set<int> pdgIDs;
      bool keepAllPdgIds = false;
      bool keepParticlesWithPositivePz = false;
      GlobalConstantsHandle<ParticleDataList> pdt;

      // TTree + branch buffers (one entry per trajectory)
      TTree* ttree = nullptr;
      unsigned int run = 0, subRun = 0, event = 0, simId = 0;
      int pdgId = 0, nPoints = 0;
      bool matched = false; // true if this particle is a filtered StepPointMC particle (vs. an ancestor)
      std::vector<double> x, y, z, t, kE, E;
      // Per-trajectory momentum is not stored in MCTrajectoryPoint (only kE), so we record
      // the SimParticle start/end momentum (MeV/c) as a best-available substitute.
      double startPx = 0, startPy = 0, startPz = 0, startP = 0;
      double endPx = 0, endPy = 0, endPz = 0, endP = 0;
      // Genealogy
      int parentSimId = -1, parentPdgId = 0, creationCode = 0;
      bool isPrimary = false;
      std::vector<int> ancestorSimIds; // parent-first, ultimate primary last
      std::vector<int> ancestorPdgIds;

      std::map<int, int> vdHitCounts;          // <pdgId, number of matching StepPointMC (VD) hits>
      std::map<int, int> tracedBackCounts;     // <pdgId, matching hits whose particle had a stored trajectory (kWritten)>
      // Breakdown of matched hits that did NOT trace, so traced < hits can be explained:
      std::map<int, int> noTrajectoryCounts;   // <pdgId, matched hits with no stored MCTrajectory>
      std::map<int, int> duplicateCounts;      // <pdgId, matched hits whose particle was already written this event>
      std::set<key_type> writtenThisEvent;     // simIds already written in the current event
  };

  ParticleTracer::ParticleTracer(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())),
    MCTrajectoriesToken(consumes<MCTrajectoryCollection>(conf().MCTrajectoriesTag())),
    keepAllPdgIds(conf().keepAllPdgIds()),
    keepParticlesWithPositivePz(conf().keepParticlesWithPositivePz()) {
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
      // SimParticle start/end momentum (MeV/c); trajectory has no per-point momentum
      ttree->Branch("startPx", &startPx, "startPx/D");
      ttree->Branch("startPy", &startPy, "startPy/D");
      ttree->Branch("startPz", &startPz, "startPz/D");
      ttree->Branch("startP",  &startP,  "startP/D");  // magnitude
      ttree->Branch("endPx",   &endPx,   "endPx/D");
      ttree->Branch("endPy",   &endPy,   "endPy/D");
      ttree->Branch("endPz",   &endPz,   "endPz/D");
      ttree->Branch("endP",    &endP,    "endP/D");    // magnitude
      // Genealogy
      ttree->Branch("parentSimId", &parentSimId, "parentSimId/I"); // -1 if primary/unavailable
      ttree->Branch("parentPdgId", &parentPdgId, "parentPdgId/I");
      ttree->Branch("creationCode", &creationCode, "creationCode/I");
      ttree->Branch("isPrimary", &isPrimary, "isPrimary/O");
      ttree->Branch("ancestorSimIds", &ancestorSimIds); // parent-first, primary last
      ttree->Branch("ancestorPdgIds", &ancestorPdgIds);
    };

  ParticleTracer::WriteResult
  ParticleTracer::writeTrajectory(const art::Ptr<SimParticle>& simPtr,
                                  const TrajIndex& trajIndex,
                                  bool matchedStep) {
    const key_type key(simPtr.key());
    const SimParticle& sim = *simPtr;

    // Only write each particle once per event.
    if (!writtenThisEvent.insert(key).second) {
      mf::LogWarning("ParticleTracer")
        << "Duplicate: SimParticle (key " << key
        << ", pdgId " << sim.pdgId() << ")"
        << (matchedStep ? " (matched StepPointMC particle)" : " (ancestor)")
        << " already written this event; skipping (repeated VD hit / shared ancestor).";
      return WriteResult::kDuplicate;
    }

    // Look up this particle's trajectory by key; absent if it didn't pass the G4
    // trajectory cuts (or wasn't produced for this particle).
    auto traj_it = trajIndex.find(key);
    if (traj_it == trajIndex.end()) {
      mf::LogWarning("ParticleTracer")
        << "No MCTrajectory found for SimParticle (key " << key
        << ", pdgId " << sim.pdgId() << ")"
        << (matchedStep ? " (matched StepPointMC particle)" : " (ancestor)")
        << "; skipping. It likely did not pass the G4 trajectory cuts.";
      return WriteResult::kNoTrajectory;
    }
    const MCTrajectory& traj = *traj_it->second;

    pdgId = sim.pdgId();
    matched = matchedStep;
    const double mass = pdt->particle(pdgId).mass();

    simId = key.asUint();
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

    // SimParticle start/end momentum (MeV/c). HepLorentzVector: vect() = 3-momentum.
    const CLHEP::HepLorentzVector startMom = sim.startMomentum();
    const CLHEP::HepLorentzVector endMom   = sim.endMomentum();
    startPx = startMom.x(); startPy = startMom.y(); startPz = startMom.z();
    startP  = startMom.vect().mag();
    endPx   = endMom.x();   endPy   = endMom.y();   endPz   = endMom.z();
    endP    = endMom.vect().mag();

    // Genealogy: immediate parent + full ancestor chain walked to the primary via the
    // SimParticle parent art::Ptrs (key and pdgId read out correctly from these).
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
    return WriteResult::kWritten;
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

    // Index the trajectories by SimParticle key (map_vector_key). We match on the key
    // rather than the art::Ptr because the StepPointMC/SimParticle collection and the
    // MCTrajectory collection are generally different products (e.g. compressed vs. g4run)
    // whose Ptrs do not compare equal, even though the keys are preserved.
    TrajIndex trajIndex;
    for (auto const& i : trajectories)
      trajIndex.emplace(key_type(i.first.key()), &i.second); // i.first is the SimParticle Ptr (its key == simid)

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

      // Optional pz filter on the StepPointMC momentum at this hit.
      if (keepParticlesWithPositivePz && step.momentum().z() <= 0)
        continue;

      // Count this VD hit (one per matching StepPointMC, keyed by the hit particle's PDG;
      // not deduplicated -- this is the raw hit multiplicity on the StepPointMC volume).
      ++vdHitCounts[stepPdgId];

      // Write the matched particle's trajectory, then every ancestor up to the primary
      // (regardless of the ancestor's PDG), walking genealogy via the parent art::Ptrs.
      // Record why a matched hit did (or did not) trace, so traced < hits is explained.
      switch (writeTrajectory(simPtr, trajIndex, true)) {
        case WriteResult::kWritten:      ++tracedBackCounts[stepPdgId];  break;
        case WriteResult::kNoTrajectory: ++noTrajectoryCounts[stepPdgId]; break;
        case WriteResult::kDuplicate:    ++duplicateCounts[stepPdgId];   break;
      }

      art::Ptr<SimParticle> anc = simPtr->parent();
      while (anc.isAvailable()) {
        writeTrajectory(anc, trajIndex, false);
        anc = anc->parent();
      }
    }
    return;
  };

  void ParticleTracer::endJob() {
    mf::LogInfo log("ParticleTracer summary");
    log << "========= StepPointMC (VD) hits per PDG ID =========\n";
    log << "hits    = matching StepPointMCs (raw VD hit multiplicity)\n";
    log << "traced  = hits whose particle had a stored MCTrajectory (a TTree entry)\n";
    log << "Why traced < hits, per PDG:\n";
    log << "  noTraj = no MCTrajectory stored (failed G4 trajectory cuts or key mismatch)\n";
    log << "  dup    = same particle already written this event (repeated VD hit)\n";
    log << "  (hits should equal traced + noTraj + dup)\n";
    for (auto const& part : vdHitCounts) {
      const int pdg = part.first;
      auto get = [pdg](const std::map<int,int>& m) {
        auto it = m.find(pdg);
        return (it != m.end()) ? it->second : 0;
      };
      const int hits   = part.second;
      const int traced = get(tracedBackCounts);
      const int noTraj = get(noTrajectoryCounts);
      const int dup    = get(duplicateCounts);
      log << "PDGID " << pdg << ": hits " << hits
          << ", traced " << traced
          << ", noTraj " << noTraj
          << ", dup " << dup;
      if (traced + noTraj + dup != hits)
        log << "  [WARNING: components do not sum to hits]";
      log << "\n";
    }
    log << "===================================================\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ParticleTracer)
