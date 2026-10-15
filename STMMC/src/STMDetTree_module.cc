// STMDetTree_module.cc
//
// One detector-aware module covering BOTH STM crystals, carrying everything the energy
// deposition analyses need: per-count energy (ionizing AND non-ionizing), the genealogy
// root that identifies what caused the deposit, and the full ancestor chain for back-trace
// plots.
//
// This is the merge the FIXME at the top of LaBrTree_module.cc asks for. HPGeTree and
// LaBrTree are byte-for-byte the same algorithm with the x cut flipped, and HPGeTree's
// "Detector" parameter was meant to select between them but throws for "LaBr". Neither
// module records enough to answer "which source / which species caused this deposit", and
// both discard the non-ionizing half of every step. They are left in place - existing
// analyses depend on their exact output - and this module is additive.
//
// What it writes
// --------------
//   ttree            one row per COUNT: one time-clustered deposit in one crystal, i.e.
//                    one pulse the detector would digitize, aggregating however many
//                    StepPointMC steps fell inside the coincidence window. A count is NOT
//                    a step - a shower's tens of steps normally collapse into a single
//                    count, and nSteps records how many went into it. "Count" is used in
//                    the detector sense, as in the counts-per-second a rate plot shows.
//                    This is the analysis-facing tree, the counterpart of HPGeTree's.
//   ttree_genealogy  one row per (count, ancestor): the back-trace chain.
//                    Written when WriteGenealogy is true (the default).
//
// A count is identified by (run, subRun, event, detector, groupId, cluster):
//   detector  which crystal - 0 = HPGe, 1 = LaBr.
//   groupId   the SimParticle KEY of the particle the deposit is attributed to, i.e. the
//             one the grouping walk below arrived at. Under the default grouping that is
//             the genealogy root - the generator's VD116 primary - so an event normally
//             has a single groupId; under GroupByTopParentInDetector it is HPGeTree's
//             top-parent-within-STMDet and an event has many. It is a real SimParticle
//             key either way, so it also identifies the particle for any further lookup.
//   cluster   which time cluster of that particle's deposits (see below), 0-based.
// The same tuple is the join key into ttree_genealogy.
//
// Grouping: by (detector, genealogy root, TIME CLUSTER)
// ----------------------------------------------------
// Two corrections to HPGeTree's grouping, in order.
//
// First, HPGeTree walks up only while the parent ITSELF deposited in a crystal, so a
// secondary born in the shielding or the hall air - whose parent never touched a crystal -
// becomes its own entry. STMDetProvenance measured the effect: for the Neutrals source that
// splits 377 incident particles into 456 entries (+21%), ~1% for EleBeam. Walking to the
// genealogy root instead merges those splits, because every particle in an event descends
// from the one VD116 primary the generator emits.
//
// Second - and this is why the root alone is NOT the right grouping - a single VD116 primary
// can deposit energy at widely separated times: a prompt photon shower, then a neutron
// thermalising and capturing much later. A real detector integrates over its own recovery
// time and produces a SEPARATE WAVEFORM, hence a separate count, for each such arrival.
// Grouping by root alone would merge them into one row and undercount the rate.
//
// So each root's steps are sorted by time and split wherever the gap to the previous step
// exceeds that detector's coincidence window (gap-based clustering: a long shower stays one
// cluster as long as energy keeps arriving, while a late neutron starts a new one). Each
// cluster is one row - one count.
//
// The two windows differ by three orders of magnitude (HPGe ~50 us recovery, LaBr ~50 ns
// scintillation), so the same event can yield one HPGe count and several LaBr counts. How
// often that actually happens is reported per detector at endJob rather than assumed here.
//
// GroupByTopParentInDetector restores HPGeTree's original grouping (and disables clustering)
// for a direct comparison against the old modules.
//
// Energy: both halves of the Geant4 split
// ---------------------------------------
// StepPointMC::ionizingEdep() is totalEDep() - nonIonizingEDep(), and HPGeTree/LaBrTree sum
// only the first. The non-ionizing part is mostly nuclear-recoil displacement damage (NIEL),
// which matters for neutrons: a keV-scale Ge recoil puts well under half its energy into
// ionization, so a neutron-induced deposit is systematically smaller in E than the energy
// actually transferred. Both are stored here, plus the total, so the identity
// total = ionizing + nonIonizing is checkable rather than assumed.
//
// Detector assignment
// -------------------
// Every STMDet step is in one of the two crystals by construction: Mu2eWorld attaches the
// sensitive detector by exact name match on fLaBrPV and fHPGePV only, so nothing else can
// produce one. The only question is WHICH crystal, and the two methods differ in how they
// answer it:
//   * by VOLUME index (when the volume table is available) - exact.
//   * by the x cut HPGeTree uses - a proxy. The HPGe is rotated 45 degrees and offset by
//     offset_HPGe, so a crystal can extend across a plane fixed at the beam centre and have
//     steps attributed to the wrong detector.
// Which was used is recorded per row in "detectorByVolume" so a sample can never silently
// mix the two conventions. A step matching neither configured crystal volume is therefore
// unexpected: it is counted and reported, since it means the configured volume names do not
// match the geometry that was actually run.
//
// The genealogy tree
// ------------------
// For each emitted count, walk particle->parent() to the root and write one row per link,
// recording each ancestor's species, creation process, birth volume, and BOTH its start and
// end point in space and time. Both endpoints matter: a back-trace draws a segment per
// ancestor, and with start points alone the deepest ancestor has nothing to draw to.
//
// The chain is written once per EMITTED ROW rather than once per step, so a shower's tens of
// steps do not repeat it; every cluster of one root carries the same chain, which is correct
// since they share an ancestry, and lets a row be plotted without a second lookup.
//
// Yongyi Wu, Oct. 2026

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "TTree.h"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  class STMDetTree : public art::EDAnalyzer {
    public:
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      // Values written to the "detector" branch, plus the array bound that sizes the
      // per-detector tallies below. An int branch rather than a string keeps the cut
      // detector==0 cheap; the order must match DETECTOR_NAMES.
      enum DetectorId { kHPGe = 0, kLaBr = 1, kNDetectors = 2 };

      struct Config {
        fhicl::Atom<art::InputTag> stepPointMCsTag{
          Name("StepPointMCsTag"),
          Comment("STMDet StepPointMCs, e.g. \"compressSTMDet:STMDet\"")};
        fhicl::Atom<art::InputTag> simParticlemvTag{
          Name("SimParticlemvTag"),
          Comment("The matching SimParticle collection, e.g. \"compressSTMDet\"")};
        fhicl::OptionalAtom<art::InputTag> physicalVolumesTag{
          Name("PhysicalVolumesTag"),
          Comment("PhysicalVolumeInfoMultiCollection (SubRun), produced by g4run. When "
                  "given, steps are assigned to a crystal by VOLUME NAME, which is exact. "
                  "Omit it and the module falls back to the x cut HPGeTree uses, which can "
                  "mis-assign between the two crystals - see the header comment.")};
        fhicl::Sequence<std::string> hpgeCrystalVolumes{
          Name("HPGeCrystalVolumes"),
          Comment("Logical volume name(s) of the HPGe crystal. These are the volumes "
                  "Mu2eWorld makes sensitive for STMDet."),
          std::vector<std::string>({"fHPGePV"})};
        fhicl::Sequence<std::string> labrCrystalVolumes{
          Name("LaBrCrystalVolumes"),
          Comment("Logical volume name(s) of the LaBr crystal."),
          std::vector<std::string>({"fLaBrPV"})};
        fhicl::Atom<double> xBeamCentre{
          Name("xBeamCentre"),
          Comment("Fallback x split when no volume table is available: HPGe is x < this, "
                  "LaBr is x > this. Matches HPGeTree/LaBrTree."),
          -3904.0};

        // Coincidence windows. One count per time cluster, so these set how finely a single
        // incident particle's deposits are resolved into separate waveforms.
        fhicl::Atom<double> hpgeCoincidenceWindow{
          Name("HPGeCoincidenceWindow"),
          Comment("HPGe gap [ns] above which a new count starts. Default 50 us is the "
                  "detector's recovery time."),
          50000.0};
        fhicl::Atom<double> labrCoincidenceWindow{
          Name("LaBrCoincidenceWindow"),
          Comment("LaBr gap [ns] above which a new count starts. Default 50 ns is the "
                  "scintillation timescale, so LaBr resolves arrivals that HPGe merges."),
          50.0};

        fhicl::Atom<bool> groupByTopParentInDetector{
          Name("GroupByTopParentInDetector"),
          Comment("false (default): group by GENEALOGY ROOT, then split into time clusters. "
                  "true: reproduce HPGeTree's grouping exactly (walk up only while the "
                  "parent itself deposited in STMDet) and DISABLE time clustering, so the "
                  "output is directly comparable with the old modules."),
          false};
        fhicl::Atom<bool> writeGenealogy{
          Name("WriteGenealogy"),
          Comment("Write ttree_genealogy, one row per (count, ancestor) link, for back-trace "
                  "plots. Off makes the output considerably smaller when only the energy "
                  "tree is wanted."),
          true};
        fhicl::Atom<double> minEdep{
          Name("MinEdep"),
          Comment("Drop counts whose TOTAL (ionizing + non-ionizing) deposit is below this "
                  "[MeV]. 0 keeps everything; filter in the analysis instead unless the file "
                  "size demands otherwise."),
          0.0};
      };

      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit STMDetTree(const Parameters& conf);
      void beginSubRun(const art::SubRun& sr) override;
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      // One step reduced to what clustering and accumulation need, so the second pass does
      // not have to re-walk the genealogy.
      struct StepInfo {
        double time = 0.0;
        double ionizing = 0.0, nonIonizing = 0.0, total = 0.0;
        double x = 0.0, y = 0.0, z = 0.0;
      };

      int detectorForStep(const StepPointMC& step) const;
      // Kinetic energy from a momentum vector and the species' rest mass.
      double kineticEnergy(const CLHEP::Hep3Vector& p, int pdgId) const;
      void fillGenealogy(const art::Ptr<SimParticle>& particle);

      art::ProductToken<StepPointMCCollection> stepsToken_;
      art::ProductToken<SimParticleCollection> simsToken_;
      art::InputTag volumesTag_;
      std::vector<std::string> hpgeVolumes_, labrVolumes_;
      double xBeamCentre_ = -3904.0;
      double window_[kNDetectors] = {50000.0, 50.0};
      bool   groupByTopParent_ = false;
      bool   writeGenealogy_ = true;
      double minEdep_ = 0.0;

      std::vector<std::string> volumeNames_;
      // Volume indices resolved to each detector at beginSubRun; empty = fall back to x.
      // Keyed on StepPointMC::VolumeId_type so the lookup needs no width conversion.
      std::set<StepPointMC::VolumeId_type> hpgeIndices_, labrIndices_;
      bool haveVolumeAssignment_ = false;

      GlobalConstantsHandle<ParticleDataList> pdt_;

      TTree* ttree = nullptr;
      TTree* gtree = nullptr;

      // ---- ttree: one row per count ----
      int    runId_ = 0, subRunId_ = 0, eventId_ = 0;
      int    detector_ = 0;            // kHPGe / kLaBr
      bool   detectorByVolume_ = false;// true: assigned by volume; false: by the x cut
      int    groupId_ = 0;             // SimParticle key of the attributed particle
      int    cluster_ = 0;             // 0-based time cluster within (detector, groupId)
      int    nClusters_ = 0;           // how many this (detector, groupId) produced
      int    pdgId_ = 0;               // that particle's species
      int    creationCode_ = 0;
      bool   isPrimary_ = false;       // it is the generator's VD116 particle
      unsigned startVolumeIndex_ = 0;
      int    genDepth_ = 0;            // generations from the depositing particle to it

      double edep_ = 0.0;              // summed IONIZING deposit [MeV]
      double nonIonEdep_ = 0.0;        // summed NON-IONIZING deposit [MeV]
      double totalEdep_ = 0.0;         // their sum
      double time_ = 0.0;              // earliest step time in the cluster [ns]
      double lastTime_ = 0.0;          // latest; lastTime - time is the cluster's duration
      unsigned long nSteps_ = 0;
      double x_ = 0.0, y_ = 0.0, z_ = 0.0;   // energy-weighted mean deposit position [mm]

      // Where the attributed particle started, and with what - the VD116 kinematics the
      // back-trace plots annotate tracks with. Identical across a root's clusters.
      double startTime_ = 0.0;
      double startX_ = 0.0, startY_ = 0.0, startZ_ = 0.0;
      double startPx_ = 0.0, startPy_ = 0.0, startPz_ = 0.0;
      double startE_ = 0.0;

      // ---- ttree_genealogy: one row per (count, ancestor) ----
      int    gRunId_ = 0, gSubRunId_ = 0, gEventId_ = 0;
      int    gDetector_ = 0;
      int    gGroupId_ = 0;
      int    gCluster_ = 0;            // with the above, joins to exactly one ttree row
      int    gDepth_ = 0;              // 1 = parent, 2 = grandparent, ...
      int    gPdgId_ = 0;
      int    gCreationCode_ = 0;
      unsigned gStartVolumeIndex_ = 0;
      double gStartTime_ = 0.0, gEndTime_ = 0.0;
      double gStartX_ = 0.0, gStartY_ = 0.0, gStartZ_ = 0.0;
      double gEndX_ = 0.0, gEndY_ = 0.0, gEndZ_ = 0.0;
      double gStartPx_ = 0.0, gStartPy_ = 0.0, gStartPz_ = 0.0;
      double gStartE_ = 0.0;

      // ---- endJob tallies ----
      unsigned long nStepsSeen_ = 0, nStepsUnassigned_ = 0, nRows_ = 0, nDropped_ = 0;
      double totalIon_[kNDetectors] = {0.0, 0.0};
      double totalNonIon_[kNDetectors] = {0.0, 0.0};
      unsigned long nRowsPerDet_[kNDetectors] = {0, 0};
      // Attributed particles, and how many of them produced >1 cluster - i.e. where the
      // detector would see several waveforms from one of them. The reason clustering exists.
      unsigned long nGroupsPerDet_[kNDetectors] = {0, 0};
      unsigned long nSplitGroupsPerDet_[kNDetectors] = {0, 0};
      std::map<int, unsigned long> groupPdgRows_;
      // Volume indices of steps matching neither crystal, for the diagnostic at endJob.
      std::map<StepPointMC::VolumeId_type, unsigned long> unassignedVolumes_;
  };

  namespace {
    const char* const DETECTOR_NAMES[] = {"HPGe", "LaBr"};
  }

  STMDetTree::STMDetTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    stepsToken_(consumes<StepPointMCCollection>(conf().stepPointMCsTag())),
    simsToken_ (consumes<SimParticleCollection>(conf().simParticlemvTag())),
    hpgeVolumes_(conf().hpgeCrystalVolumes()),
    labrVolumes_(conf().labrCrystalVolumes()),
    xBeamCentre_(conf().xBeamCentre()),
    groupByTopParent_(conf().groupByTopParentInDetector()),
    writeGenealogy_(conf().writeGenealogy()),
    minEdep_(conf().minEdep())
  {
    window_[kHPGe] = conf().hpgeCoincidenceWindow();
    window_[kLaBr] = conf().labrCoincidenceWindow();
    for (int d = 0; d < kNDetectors; ++d)
      if (!(window_[d] > 0.0))
        throw cet::exception("STMDetTree")
          << DETECTOR_NAMES[d] << " coincidence window must be positive (got "
          << window_[d] << " ns). A zero or negative window would make every step its own "
          << "count.";

    art::InputTag vt;
    if (conf().physicalVolumesTag(vt)) volumesTag_ = vt;

    art::ServiceHandle<art::TFileService> tfs;

    ttree = tfs->make<TTree>("ttree", "STM crystal counts (time-clustered deposits)");
    ttree->Branch("run",    &runId_,    "run/I");
    ttree->Branch("subRun", &subRunId_, "subRun/I");
    ttree->Branch("event",  &eventId_,  "event/I");

    ttree->Branch("detector",         &detector_,         "detector/I");
    ttree->Branch("detectorByVolume", &detectorByVolume_, "detectorByVolume/O");

    ttree->Branch("groupId",          &groupId_,          "groupId/I");
    ttree->Branch("cluster",          &cluster_,          "cluster/I");
    ttree->Branch("nClusters",        &nClusters_,        "nClusters/I");
    ttree->Branch("pdgId",            &pdgId_,            "pdgId/I");
    ttree->Branch("creationCode",     &creationCode_,     "creationCode/I");
    ttree->Branch("isPrimary",        &isPrimary_,        "isPrimary/O");
    ttree->Branch("startVolumeIndex", &startVolumeIndex_, "startVolumeIndex/i");
    ttree->Branch("genDepth",         &genDepth_,         "genDepth/I");

    ttree->Branch("E",        &edep_,       "E/D");        // named E to match HPGeTree
    ttree->Branch("nonIonE",  &nonIonEdep_, "nonIonE/D");
    ttree->Branch("totalE",   &totalEdep_,  "totalE/D");
    ttree->Branch("time",     &time_,       "time/D");     // named time to match HPGeTree
    ttree->Branch("lastTime", &lastTime_,   "lastTime/D");
    ttree->Branch("nSteps",   &nSteps_,     "nSteps/l");
    ttree->Branch("x", &x_, "x/D");
    ttree->Branch("y", &y_, "y/D");
    ttree->Branch("z", &z_, "z/D");

    ttree->Branch("startTime", &startTime_, "startTime/D");
    ttree->Branch("startx",    &startX_,    "startx/D");
    ttree->Branch("starty",    &startY_,    "starty/D");
    ttree->Branch("startz",    &startZ_,    "startz/D");
    ttree->Branch("startpx",   &startPx_,   "startpx/D");
    ttree->Branch("startpy",   &startPy_,   "startpy/D");
    ttree->Branch("startpz",   &startPz_,   "startpz/D");
    ttree->Branch("startE",    &startE_,    "startE/D");

    if (writeGenealogy_) {
      gtree = tfs->make<TTree>("ttree_genealogy", "Ancestor chain per count");
      gtree->Branch("run",      &gRunId_,    "run/I");
      gtree->Branch("subRun",   &gSubRunId_, "subRun/I");
      gtree->Branch("event",    &gEventId_,  "event/I");
      gtree->Branch("detector", &gDetector_, "detector/I");
      gtree->Branch("groupId",  &gGroupId_,  "groupId/I");
      gtree->Branch("cluster",  &gCluster_,  "cluster/I");
      gtree->Branch("depth",    &gDepth_,    "depth/I");

      gtree->Branch("pdgId",            &gPdgId_,            "pdgId/I");
      gtree->Branch("creationCode",     &gCreationCode_,     "creationCode/I");
      gtree->Branch("startVolumeIndex", &gStartVolumeIndex_, "startVolumeIndex/i");

      gtree->Branch("startGlobalTime", &gStartTime_, "startGlobalTime/D");
      gtree->Branch("endGlobalTime",   &gEndTime_,   "endGlobalTime/D");
      gtree->Branch("startx", &gStartX_, "startx/D");
      gtree->Branch("starty", &gStartY_, "starty/D");
      gtree->Branch("startz", &gStartZ_, "startz/D");
      gtree->Branch("endx",   &gEndX_,   "endx/D");
      gtree->Branch("endy",   &gEndY_,   "endy/D");
      gtree->Branch("endz",   &gEndZ_,   "endz/D");
      gtree->Branch("startpx", &gStartPx_, "startpx/D");
      gtree->Branch("startpy", &gStartPy_, "startpy/D");
      gtree->Branch("startpz", &gStartPz_, "startpz/D");
      gtree->Branch("E",       &gStartE_,  "E/D");
    }
  }

  void STMDetTree::beginSubRun(const art::SubRun& sr) {
    volumeNames_.clear();
    hpgeIndices_.clear();
    labrIndices_.clear();
    haveVolumeAssignment_ = false;
    if (volumesTag_.label().empty()) return;

    art::Handle<PhysicalVolumeInfoMultiCollection> volsHandle;
    if (!sr.getByLabel(volumesTag_, volsHandle)) {
      mf::LogWarning("STMDetTree")
        << "PhysicalVolumesTag " << volumesTag_ << " not found in this SubRun; falling back "
        << "to the x cut for detector assignment. See detectorByVolume in the tree.";
      return;
    }

    for (const PhysicalVolumeInfoSingleStage& stage : *volsHandle) {
      for (const auto& entry : stage) {
        const unsigned idx = entry.first.asUint();
        if (volumeNames_.size() <= idx) volumeNames_.resize(idx + 1);
        volumeNames_[idx] = entry.second.name();
      }
    }

    for (StepPointMC::VolumeId_type i = 0; i < volumeNames_.size(); ++i) {
      for (const std::string& n : hpgeVolumes_)
        if (volumeNames_[i] == n) hpgeIndices_.insert(i);
      for (const std::string& n : labrVolumes_)
        if (volumeNames_[i] == n) labrIndices_.insert(i);
    }

    // Both crystals must resolve: assigning one by volume and the other by x would put the
    // two detectors on different definitions within a single file, which is worse than
    // using the proxy consistently.
    haveVolumeAssignment_ = !hpgeIndices_.empty() && !labrIndices_.empty();
    if (!haveVolumeAssignment_)
      mf::LogWarning("STMDetTree")
        << "Could not resolve both crystal volumes in the volume table (HPGe matched "
        << hpgeIndices_.size() << ", LaBr matched " << labrIndices_.size()
        << "); falling back to the x cut. Check HPGeCrystalVolumes / LaBrCrystalVolumes "
        << "against the geometry actually used.";
    else
      mf::LogInfo("STMDetTree")
        << "Detector assignment by volume: HPGe " << hpgeIndices_.size()
        << " volume(s), LaBr " << labrIndices_.size() << " volume(s).";
  }

  double STMDetTree::kineticEnergy(const CLHEP::Hep3Vector& p, int pdgId) const {
    // An unknown pdgId (an exotic ion, say) has no ParticleDataList entry; treat it as
    // massless rather than aborting the job, since this is a plotting convenience and the
    // momentum is recorded alongside it either way.
    double mass = 0.0;
    try {
      mass = pdt_->particle(pdgId).mass();
    } catch (...) {
      mass = 0.0;
    }
    return std::sqrt(p.mag2() + mass * mass) - mass;
  }

  // Which crystal a step is in. Returns -1 only when volume assignment is active and the
  // step matched neither configured crystal - which should not happen, since STMDet is
  // attached to those volumes alone, so it means the names are wrong for this geometry.
  // The x-cut fallback always returns one of the two: it cannot detect that case.
  int STMDetTree::detectorForStep(const StepPointMC& step) const {
    if (haveVolumeAssignment_) {
      const StepPointMC::VolumeId_type vid = step.volumeId();
      if (hpgeIndices_.count(vid)) return kHPGe;
      if (labrIndices_.count(vid)) return kLaBr;
      return -1;
    }
    return (step.position().x() < xBeamCentre_) ? kHPGe : kLaBr;
  }

  // Walk from `particle` up to the genealogy root, writing one row per ancestor. The caller
  // has already set the join-key branches (run/subRun/event/detector/groupId/cluster).
  //
  // Both the start and end point of every ancestor are recorded, not just the start: a
  // back-trace draws one segment per link, and the deepest ancestor would otherwise have no
  // endpoint to draw to.
  void STMDetTree::fillGenealogy(const art::Ptr<SimParticle>& particle) {
    if (!gtree) return;

    art::Ptr<SimParticle> cur = particle;
    int depth = 0;
    while (cur->hasParent() && cur->parent().isAvailable()) {
      cur = cur->parent();
      ++depth;
      if (depth > 10000)
        throw cet::exception("STMDetTree")
          << "Genealogy walk exceeded 10000 generations in event " << eventId_
          << "; the parent links are cyclic.";

      gDepth_            = depth;
      gPdgId_            = cur->pdgId();
      gCreationCode_     = static_cast<int>(cur->creationCode().id());
      gStartVolumeIndex_ = cur->startVolumeIndex();
      gStartTime_        = cur->startGlobalTime();
      gEndTime_          = cur->endGlobalTime();
      gStartX_ = cur->startPosition().x();
      gStartY_ = cur->startPosition().y();
      gStartZ_ = cur->startPosition().z();
      gEndX_   = cur->endPosition().x();
      gEndY_   = cur->endPosition().y();
      gEndZ_   = cur->endPosition().z();
      gStartPx_ = cur->startMomentum().x();
      gStartPy_ = cur->startMomentum().y();
      gStartPz_ = cur->startMomentum().z();
      gStartE_  = kineticEnergy(cur->startMomentum().vect(), gPdgId_);
      gtree->Fill();
    }
  }

  void STMDetTree::analyze(const art::Event& event) {
    auto const& steps = event.getProduct(stepsToken_);
    // A handle rather than the bare product: its ProductID is what the mis-pairing check
    // below compares each step's SimParticle Ptr against.
    auto const simsHandle = event.getValidHandle(simsToken_);
    auto const& sims = *simsHandle;
    if (steps.empty() || sims.empty()) return;

    runId_    = event.run();
    subRunId_ = event.subRun();
    eventId_  = event.event();

    // HPGeTree's set: particles with at least one STMDet step. Needed only for the
    // GroupByTopParentInDetector comparison mode.
    std::set<cet::map_vector_key> stmDetIds;
    if (groupByTopParent_)
      for (const StepPointMC& step : steps)
        stmDetIds.insert(sims.at(step.trackId()).id());

    // (detector, groupId) -> its steps, collected before clustering because the split
    // depends on the whole time-ordered set, not on any one step. An ordered map so rows
    // come out in a stable order regardless of step ordering.
    std::map<std::pair<int, int>, std::vector<StepInfo>> byGroup;
    std::map<int, art::Ptr<SimParticle>> groupPtrs;
    std::map<int, int> groupDepth;

    for (const StepPointMC& step : steps) {
      ++nStepsSeen_;
      const int det = detectorForStep(step);
      if (det < 0) {
        ++nStepsUnassigned_;
        ++unassignedVolumes_[step.volumeId()];
        continue;
      }

      art::Ptr<SimParticle> particle = step.simParticle();
      if (particle.isNull() || !particle.isAvailable()) {
        // Without the SimParticle there is nothing to group by; counted so the totals stay
        // honest rather than silently losing the energy.
        ++nStepsUnassigned_;
        continue;
      }

      // StepPointMCsTag and SimParticlemvTag must name the SAME compression chain. The STM
      // job runs several in parallel - compressSTMDet over the crystal (STMDet) and
      // virtualdetector steps, compressDetStepMCsSTM* over the VD101/89/90 extractions -
      // and each writes its OWN SimParticle collection, every one of them a full copy of
      // g4run's with the parent Ptrs remapped into itself. Pairing crystal steps with
      // another chain's SimParticles therefore resolves (or fails to resolve) against the
      // wrong collection, which is why this is checked rather than trusted: the symptom
      // would otherwise be silently wrong attribution, not a crash.
      if (particle.id() != simsHandle.id())
        throw cet::exception("STMDetTree")
          << "StepPointMCsTag steps carry SimParticle Ptrs into a different collection than "
          << "SimParticlemvTag names, so the two tags are from different compression "
          << "chains. Pair a chain's steps with its own SimParticles - e.g. "
          << "\"compressSTMDet:STMDet\" with \"compressSTMDet\", NOT with "
          << "\"compressDetStepMCsSTM\", which compresses the VD101 extraction and holds a "
          << "separate SimParticle collection.";

      // ---- find the particle this deposit is attributed to ----
      art::Ptr<SimParticle> group = particle;
      int depth = 0;
      if (groupByTopParent_) {
        // HPGeTree's walk: up only while the parent itself deposited in a crystal.
        if (group->hasParent() && group->parent().isAvailable() &&
            stmDetIds.count(group->parent()->id())) {
          while (group->parent()->hasParent() &&
                 group->parent()->parent().isAvailable() &&
                 stmDetIds.count(group->parent()->parent()->id())) {
            group = group->parent();
            ++depth;
          }
          group = group->parent();
          ++depth;
        }
      } else {
        // To the genealogy root, with no STMDet restriction.
        while (group->hasParent() && group->parent().isAvailable()) {
          group = group->parent();
          ++depth;
          if (depth > 10000)
            throw cet::exception("STMDetTree")
              << "Genealogy walk exceeded 10000 generations in event " << eventId_
              << "; the parent links are cyclic.";
        }
      }

      const int gid = static_cast<int>(group->id().asUint());
      if (!groupPtrs.count(gid)) {
        groupPtrs[gid]  = group;
        groupDepth[gid] = depth;
      }

      StepInfo s;
      s.time        = step.time();
      s.ionizing    = step.ionizingEdep();
      s.nonIonizing = step.nonIonizingEDep();
      s.total       = step.totalEDep();
      s.x = step.position().x();
      s.y = step.position().y();
      s.z = step.position().z();
      byGroup[{det, gid}].push_back(s);
    }

    // ---- cluster each group in time, and emit one row per cluster ----
    for (auto& [key, stepList] : byGroup) {
      const int det = key.first;
      const int gid = key.second;
      ++nGroupsPerDet_[det];

      // Steps arrive in StepPointMC order, which is not time order; the gap test is only
      // meaningful on a sorted list.
      std::sort(stepList.begin(), stepList.end(),
                [](const StepInfo& a, const StepInfo& b) { return a.time < b.time; });

      // Cluster boundaries: a new cluster starts wherever the gap to the previous step
      // exceeds this detector's window. In the HPGeTree-compatibility mode clustering is
      // disabled outright, so the old one-row-per-top-parent output is reproduced.
      const double window = window_[det];
      std::vector<std::pair<size_t, size_t>> clusters;   // [begin, end)
      size_t begin = 0;
      for (size_t i = 1; i < stepList.size(); ++i) {
        if (!groupByTopParent_ && (stepList[i].time - stepList[i - 1].time) > window) {
          clusters.emplace_back(begin, i);
          begin = i;
        }
      }
      clusters.emplace_back(begin, stepList.size());
      if (clusters.size() > 1) ++nSplitGroupsPerDet_[det];

      const art::Ptr<SimParticle>& group = groupPtrs[gid];
      const int gPdg = group->pdgId();

      for (size_t ci = 0; ci < clusters.size(); ++ci) {
        const size_t b = clusters[ci].first;
        const size_t e = clusters[ci].second;

        double ion = 0.0, nonIon = 0.0;
        double sumWx = 0.0, sumWy = 0.0, sumWz = 0.0, sumW = 0.0;
        for (size_t i = b; i < e; ++i) {
          ion    += stepList[i].ionizing;
          nonIon += stepList[i].nonIonizing;
          // Weight the mean position by total deposit, so a step that deposited nothing
          // does not pull the centroid.
          sumWx += stepList[i].total * stepList[i].x;
          sumWy += stepList[i].total * stepList[i].y;
          sumWz += stepList[i].total * stepList[i].z;
          sumW  += stepList[i].total;
        }

        totalEdep_ = ion + nonIon;
        if (minEdep_ > 0.0 && totalEdep_ < minEdep_) { ++nDropped_; continue; }

        detector_         = det;
        detectorByVolume_ = haveVolumeAssignment_;
        groupId_          = gid;
        cluster_          = static_cast<int>(ci);
        nClusters_        = static_cast<int>(clusters.size());
        pdgId_            = gPdg;
        creationCode_     = static_cast<int>(group->creationCode().id());
        isPrimary_        = group->isPrimary();
        startVolumeIndex_ = group->startVolumeIndex();
        genDepth_         = groupDepth[gid];

        edep_       = ion;
        nonIonEdep_ = nonIon;
        time_       = stepList[b].time;
        lastTime_   = stepList[e - 1].time;
        nSteps_     = e - b;

        if (sumW > 0.0) {
          x_ = sumWx / sumW;
          y_ = sumWy / sumW;
          z_ = sumWz / sumW;
        } else {
          x_ = y_ = z_ = 0.0;
        }

        startTime_ = group->startGlobalTime();
        startX_ = group->startPosition().x();
        startY_ = group->startPosition().y();
        startZ_ = group->startPosition().z();
        startPx_ = group->startMomentum().x();
        startPy_ = group->startMomentum().y();
        startPz_ = group->startMomentum().z();
        startE_  = kineticEnergy(group->startMomentum().vect(), gPdg);

        ttree->Fill();

        ++nRows_;
        ++nRowsPerDet_[det];
        totalIon_[det]    += ion;
        totalNonIon_[det] += nonIon;
        ++groupPdgRows_[gPdg];

        // One chain per emitted row. Written after the row so a dropped row writes no
        // chain, keeping the two trees joinable without a dangling-key check.
        if (writeGenealogy_) {
          gRunId_ = runId_; gSubRunId_ = subRunId_; gEventId_ = eventId_;
          gDetector_ = det; gGroupId_ = gid; gCluster_ = static_cast<int>(ci);
          fillGenealogy(group);
        }
      }
    }
  }

  void STMDetTree::endJob() {
    mf::LogInfo log("STMDetTree");
    log << "========== STMDet tree summary ==========\n"
        << "grouping                 : "
        << (groupByTopParent_
              ? "top parent within STMDet, NO time clustering (HPGeTree compatibility)"
              : "genealogy root, split into time clusters") << "\n"
        << "detector assignment      : "
        << (haveVolumeAssignment_ ? "by volume" : "by x cut (proxy)") << "\n";
    if (!groupByTopParent_)
      log << "coincidence windows      : HPGe " << window_[kHPGe] << " ns, LaBr "
          << window_[kLaBr] << " ns\n";
    log << "STMDet steps seen        : " << nStepsSeen_ << "\n"
        << "  not assigned to either : " << nStepsUnassigned_ << "\n"
        << "rows (counts) written    : " << nRows_ << "\n";
    if (nDropped_)
      log << "  dropped below MinEdep  : " << nDropped_ << "\n";
    for (int d = 0; d < kNDetectors; ++d) {
      log << "  " << DETECTOR_NAMES[d] << ": " << nRowsPerDet_[d] << " count(s) from "
          << nGroupsPerDet_[d] << " attributed particle(s), " << totalIon_[d]
          << " MeV ionizing, " << totalNonIon_[d] << " MeV non-ionizing\n";
      if (!groupByTopParent_)
        log << "        " << nSplitGroupsPerDet_[d]
            << " of them split into >1 count by the " << window_[d] << " ns window\n";
    }
    log << "counts by attributed-particle PDG:\n";
    for (const auto& [pdg, n] : groupPdgRows_)
      log << "  " << pdg << ": " << n << "\n";

    // STMDet is attached to the two crystal volumes alone, so every step should match one
    // of them. Any that did not means the configured names do not describe the geometry
    // that was run, and those steps' energy is missing from the trees.
    if (!unassignedVolumes_.empty()) {
      log << "\nWARNING: " << nStepsUnassigned_ << " step(s) matched neither configured "
          << "crystal volume, and their energy is NOT in the trees. STMDet is attached only "
          << "to the crystals, so this means HPGeCrystalVolumes / LaBrCrystalVolumes do not "
          << "match this geometry. Volume indices seen:\n";
      for (const auto& [idx, n] : unassignedVolumes_) {
        log << "  idx " << idx;
        if (idx < volumeNames_.size() && !volumeNames_[idx].empty())
          log << " (" << volumeNames_[idx] << ")";
        log << ": " << n << " step(s)\n";
      }
    }
    log << "=========================================\n";
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::STMDetTree)
