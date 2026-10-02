#ifndef STMGeom_LeftWall_hh
#define STMGeom_LeftWall_hh

// The left wall of the STM shield house: the wall along -x, running
// downstream from the front shielding.
//
// The mirror of the right wall, with the same five kinds of layer and
// the same sheet shapes, but no L-shaped edge piece: 26 pieces. Nothing
// here is bored.
//
// Layers are listed outside to inside and walk +x from the reference,
// which is the -x face of the outermost layer, so every piece is on the
// same side of it. The lead layers are BrickWalls with courses running
// along z, as in the right wall.
//
// Every position arrives resolved from STMMaker. The wall is fixed by
// that reference and by the baseplate, so each piece follows from
// butting against its neighbour rather than from a measurement of its
// own.
//
// This is deliberately its own class rather than a generalisation of
// RightWall: the two differ by the edge prism today, and whether a
// shared side-wall description is the right shape is better decided
// once the top and inner walls are in.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // A flat sheet: the two poly pairs and the copper sheet. No bores and
  // no orientation -- halfDim is already in the Mu2e frame, so which of
  // the three is the thickness is read off it rather than from a
  // rotation applied afterwards.
  //
  // These are not all of one material, so each carries the name its
  // volume takes, set by STMMaker where which-sheet-is-which is known.
  struct LeftWallSheet {
    std::string       name;
    std::string       material;
    CLHEP::Hep3Vector halfDim;
    CLHEP::Hep3Vector center;
  };

  class LeftWall {
  public:

    LeftWall(bool build,
             std::vector<BrickWall> const & leadLayers,
             std::vector<LeftWallSheet> const & sheets
             ) :
      _build(build),
      _leadLayers(leadLayers),
      _sheets(sheets)
    {
    }

    bool build() const {return _build;}

    // The two lead layers, outside to inside: four courses of three,
    // then three courses of three.
    std::vector<BrickWall> const & leadLayers() const {return _leadLayers;}

    // The flat sheets, in the order the layers are built: the outer
    // longwall and top edge, the inner pair, then the copper sheet.
    std::vector<LeftWallSheet> const & sheets() const {return _sheets;}

    // Genreflex can't do persistency of vector<LeftWall> without a
    // default constructor
    LeftWall() {}

  private:

    bool _build;

    std::vector<BrickWall>     _leadLayers;
    std::vector<LeftWallSheet> _sheets;
  };

}

#endif/*STMGeom_LeftWall_hh*/
