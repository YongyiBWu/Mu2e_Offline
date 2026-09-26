#ifndef STMGeom_FrontShieldingLeft_hh
#define STMGeom_FrontShieldingLeft_hh

// The left half of the STM front shielding: the section beside the
// spot-size collimator.
//
// Twenty-two pieces at five depths -- four poly sheets, three groups
// of lead bricks, and a lead prism. Nothing here is bored: the beam
// passes through the right half.
//
// The brick groups are BrickWalls, the same as the right half's lead
// layers, but they are not layers of one wall. Each is a separate
// stack at its own depth, so they are listed rather than sequenced,
// and the poly here is held as plain sheets because none of it
// interleaves with a lead layer the way the right half's does.
//
// Every position arrives resolved from STMMaker. The section keys off
// the right half -- x from FrontShieldingRight.brickEndX, y from the
// baseplate, z from the right half's last lead layer -- so nothing
// here is measured independently.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // A poly sheet. These carry no bores -- the beam is in the right
  // half -- and they need no orientation: the sheets do not all lie
  // the same way up, but halfDim is already in the Mu2e frame, so
  // which of the three is the 1 in thickness is read off it rather
  // than from a rotation applied afterwards.
  struct FrontShieldingLeftSheet {
    std::string       material;
    CLHEP::Hep3Vector halfDim;
    CLHEP::Hep3Vector center;
  };

  // The lead prism in the corner at the inner sheet's +x edge. Its
  // outline is given in the ExtShieldDownstream form, a u/v polygon
  // plus a sweep length, like the right half's copper lining.
  struct FrontShieldingLeftPrism {
    std::string         material;
    std::vector<double> uVerts;
    std::vector<double> vVerts;
    double              length;
    // Unlike a sheet, this one does need its orientation: the outline
    // is a polygon in its own u/v frame and the rotation is what puts
    // that frame into the Mu2e one.
    std::string         orientation;
    // Placed by an ANCHOR, not a centre. The outline is drawn so that
    // the local origin is the prism's RIGHT ANGLE, the corner its two
    // legs meet at, and that corner is what butts the inner sheet --
    // its +x face and its -z face. A triangle's centroid lines up
    // with nothing, so the corner is the only useful handle.
    //
    // The sweep is centred on the anchor, as G4ExtrudedSolid centres
    // an extrusion on its placement point, so the anchor's y is the
    // mid-plane of the 6 in and not its top edge.
    CLHEP::Hep3Vector   anchor;
  };

  class FrontShieldingLeft {
  public:

    FrontShieldingLeft(bool build,
                       std::vector<FrontShieldingLeftSheet> const & sheets,
                       std::vector<BrickWall> const & brickGroups,
                       FrontShieldingLeftPrism const & prism
                       ) :
      _build(build),
      _sheets(sheets),
      _brickGroups(brickGroups),
      _prism(prism)
    {
    }

    bool build() const {return _build;}

    // The four poly sheets: side, outer, inner, edge, in that order.
    std::vector<FrontShieldingLeftSheet> const & sheets() const {return _sheets;}

    // The three brick groups: the nine-brick block, the edge course,
    // and the grid behind the inner sheet, in that order. Each is a
    // BrickWall in its own plane at its own depth.
    std::vector<BrickWall> const & brickGroups() const {return _brickGroups;}

    FrontShieldingLeftPrism const & prism() const {return _prism;}

    // Genreflex can't do persistency of vector<FrontShieldingLeft>
    // without a default constructor
    FrontShieldingLeft() {}

  private:

    bool _build;

    std::vector<FrontShieldingLeftSheet> _sheets;
    std::vector<BrickWall>               _brickGroups;
    FrontShieldingLeftPrism              _prism;
  };

}

#endif/*STMGeom_FrontShieldingLeft_hh*/
