#ifndef STMGeom_SSCFrontShield_hh
#define STMGeom_SSCFrontShield_hh

// The shielding in front of the STM spot-size collimator: a wall of
// lead bricks, an aluminium shelf, and two polyethylene blocks, one
// of them bored for the beam.
//
// The bricks themselves are described by LeadBrick, which is shared
// with every other structure that stacks them; this class says only
// where they go.
//
// Positions arrive here already resolved. The geometry file writes
// every center as if the structure sat exactly on the SSC axis, and
// STMMaker adds the offsets that displace it -- the whole-structure
// offsetX/offsetY, the LaBr side brick offset, and the poly2 offset --
// so nothing downstream has to know which offset applies to which
// piece.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class SSCFrontShield {
  public:

    SSCFrontShield(bool build,
                   std::vector<CLHEP::Hep3Vector> const & brick2x4x8Center,
                   std::vector<std::string>       const & brick2x4x8Orientation,
                   std::vector<CLHEP::Hep3Vector> const & brick2x4x16Center,
                   std::vector<std::string>       const & brick2x4x16Orientation,
                   std::string const & shelfMaterial,
                   CLHEP::Hep3Vector const & shelfDim,
                   CLHEP::Hep3Vector const & shelfCenter,
                   std::string const & poly1Material,
                   CLHEP::Hep3Vector const & poly1Dim,
                   CLHEP::Hep3Vector const & poly1Center,
                   double poly1BoreR, double poly1BoreDX, double poly1BoreDY,
                   std::string const & poly2Material,
                   CLHEP::Hep3Vector const & poly2Dim,
                   CLHEP::Hep3Vector const & poly2Center
                   ) :
      _build(build),
      _brick2x4x8Center(brick2x4x8Center),
      _brick2x4x8Orientation(brick2x4x8Orientation),
      _brick2x4x16Center(brick2x4x16Center),
      _brick2x4x16Orientation(brick2x4x16Orientation),
      _shelfMaterial(shelfMaterial),
      _shelfDim(shelfDim),
      _shelfCenter(shelfCenter),
      _poly1Material(poly1Material),
      _poly1Dim(poly1Dim),
      _poly1Center(poly1Center),
      _poly1BoreR(poly1BoreR),
      _poly1BoreDX(poly1BoreDX),
      _poly1BoreDY(poly1BoreDY),
      _poly2Material(poly2Material),
      _poly2Dim(poly2Dim),
      _poly2Center(poly2Center)
    {
    }

    bool build() const {return _build;}

    // Where the bricks go. Dimensions and material come from LeadBrick.
    std::vector<CLHEP::Hep3Vector> const & brick2x4x8Center()  const {return _brick2x4x8Center;}
    std::vector<std::string>       const & brick2x4x8Orientation() const {return _brick2x4x8Orientation;}
    std::vector<CLHEP::Hep3Vector> const & brick2x4x16Center() const {return _brick2x4x16Center;}
    std::vector<std::string>       const & brick2x4x16Orientation() const {return _brick2x4x16Orientation;}

    // Aluminium shelf, carrying the bricks.
    std::string const & shelfMaterial() const {return _shelfMaterial;}
    CLHEP::Hep3Vector const & shelfDim()    const {return _shelfDim;}
    CLHEP::Hep3Vector const & shelfCenter() const {return _shelfCenter;}

    // The large bored poly block. The bore runs along z through the
    // full depth, offset from the block centre by (boreDX, boreDY) so
    // that it lands on the SSC axis.
    std::string const & poly1Material() const {return _poly1Material;}
    CLHEP::Hep3Vector const & poly1Dim()    const {return _poly1Dim;}
    CLHEP::Hep3Vector const & poly1Center() const {return _poly1Center;}
    double poly1BoreR()  const {return _poly1BoreR;}
    double poly1BoreDX() const {return _poly1BoreDX;}
    double poly1BoreDY() const {return _poly1BoreDY;}

    // The small poly block in the beam opening.
    std::string const & poly2Material() const {return _poly2Material;}
    CLHEP::Hep3Vector const & poly2Dim()    const {return _poly2Dim;}
    CLHEP::Hep3Vector const & poly2Center() const {return _poly2Center;}

    // Genreflex can't do persistency of vector<SSCFrontShield> without
    // a default constructor
    SSCFrontShield() {}

  private:

    bool _build;

    std::vector<CLHEP::Hep3Vector> _brick2x4x8Center;
    std::vector<std::string>       _brick2x4x8Orientation;
    std::vector<CLHEP::Hep3Vector> _brick2x4x16Center;
    std::vector<std::string>       _brick2x4x16Orientation;

    std::string _shelfMaterial;
    CLHEP::Hep3Vector _shelfDim;
    CLHEP::Hep3Vector _shelfCenter;

    std::string _poly1Material;
    CLHEP::Hep3Vector _poly1Dim;
    CLHEP::Hep3Vector _poly1Center;
    double      _poly1BoreR;
    double      _poly1BoreDX;
    double      _poly1BoreDY;

    std::string _poly2Material;
    CLHEP::Hep3Vector _poly2Dim;
    CLHEP::Hep3Vector _poly2Center;
  };

}

#endif/*STMGeom_SSCFrontShield_hh*/
