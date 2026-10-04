#ifndef STMGeom_InnerShielding_hh
#define STMGeom_InnerShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // ---- the updated, hand-stacked inner shielding ---------------------
  //
  // What lines the cavity the detectors sit in: six copper pieces and
  // eleven lead ones. Unlike the walls this is not a stack of layers --
  // every piece is placed on its own -- so there is no layer sequence,
  // just the three lists below.
  //
  // Some pieces are bored on the LaBr beam axis. That bore follows the
  // collimator rather than this section, so a piece carries only a flag
  // and the construction code reads the axis from the front shielding's
  // bore, as the right half's sheets and lining do.

  // A prism swept along one axis, given in the ExtShieldDownstream
  // form: a u/v cap polygon plus a sweep length.
  //
  // Each is placed by an ANCHOR on its cap origin rather than by a
  // centre, since the outlines are traced from the vertex nearest the
  // section reference. The sweep is centred on that point, as
  // G4ExtrudedSolid centres an extrusion on its placement point.
  struct InnerShieldingPrism {
    std::string         name;
    std::string         material;
    std::vector<double> uVerts;
    std::vector<double> vVerts;
    double              length;        // the sweep
    std::string         orientation;
    CLHEP::Hep3Vector   anchor;
    bool                bored;         // on the LaBr axis
  };

  // A plain box: the two copper plates, and the three lead pieces that
  // are turned 45 degrees about y and so cannot be stated as leadBrick
  // types.
  struct InnerShieldingBox {
    std::string       name;
    std::string       material;
    CLHEP::Hep3Vector halfDim;
    CLHEP::Hep3Vector center;
    std::string       orientation;
    bool              bored;           // on the LaBr axis
    // Set on the pieces whose nominal size would touch their
    // neighbours once turned. The construction code backs each
    // half-dimension off by its own nudge, so halfDim above stays the
    // nominal value the config states.
    bool              nudge;
  };

  // One of the standard lead bricks, placed individually.
  struct InnerShieldingBrick {
    int               type;            // index into LeadBrick's types
    std::string       orientation;
    CLHEP::Hep3Vector center;
    bool              bored;           // on the LaBr axis
  };

  class InnerShielding {
  public:

    // The earlier description. Its geometry is written out in
    // constructSTM.cc, so the class carries only the switch.
    InnerShielding(bool build):
      _build(build)
    {}

    // The updated, hand-stacked description. Taking a separate
    // constructor rather than defaulted arguments keeps each to the
    // parameters it actually has, so neither can be built with a
    // quantity that does not apply to it.
    //
    // Every position arrives resolved from STMMaker, fixed by one
    // reference -- the bottom wall's, moved in 0.5 in and up 6 in --
    // so nothing here was measured independently of the house.
    InnerShielding(bool build,
                   std::vector<InnerShieldingPrism> const & prisms,
                   std::vector<InnerShieldingBox> const & boxes,
                   std::vector<InnerShieldingBrick> const & bricks):
      _build(build),
      _prisms(prisms),
      _boxes(boxes),
      _bricks(bricks)
    {}

    bool   build()                               const { return _build;  }

    // ---- updated (hand-stacked) only; empty for the earlier one -----

    // The swept pieces: four copper, and on the lead side the angled
    // triangle, the bored trapezoid and the two clipped prisms.
    std::vector<InnerShieldingPrism> const & prisms() const {return _prisms;}

    // The boxes: two copper plates and the three angled lead pieces.
    std::vector<InnerShieldingBox> const & boxes() const {return _boxes;}

    // The standard lead bricks placed individually.
    std::vector<InnerShieldingBrick> const & bricks() const {return _bricks;}

    InnerShielding() {}
  private:

    bool               _build;

    std::vector<InnerShieldingPrism> _prisms;
    std::vector<InnerShieldingBox>   _boxes;
    std::vector<InnerShieldingBrick> _bricks;

  };

}

#endif/*STMGeom_InnerShielding_hh*/
