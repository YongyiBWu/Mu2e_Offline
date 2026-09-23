#ifndef STMGeom_LeadBrick_hh
#define STMGeom_LeadBrick_hh

// The standard lead bricks the STM shield house is stacked from.
//
// "wear" is taken off each face, so a brick is (dx - 2*wear) on a side
// while staying centred where it was placed. A stack therefore keeps
// its nominal pitch and the wear opens as gaps between bricks rather
// than displacing anything.
//
// Author: Yongyi Wu

#include <string>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LeadBrick {
  public:

    LeadBrick(CLHEP::Hep3Vector const & dim2x4x8,
              CLHEP::Hep3Vector const & dim2x4x16,
              double wear,
              std::string const & material
              ) :
      _dim2x4x8(dim2x4x8),
      _dim2x4x16(dim2x4x16),
      _wear(wear),
      _material(material)
    {
    }

    // Nominal outside dimensions, before wear.
    CLHEP::Hep3Vector const & dim2x4x8()  const {return _dim2x4x8;}
    CLHEP::Hep3Vector const & dim2x4x16() const {return _dim2x4x16;}

    // The same, with the wear taken off each face. This is what the
    // G4Box should be built from.
    CLHEP::Hep3Vector worn2x4x8()  const {return worn(_dim2x4x8);}
    CLHEP::Hep3Vector worn2x4x16() const {return worn(_dim2x4x16);}

    double              wear()     const {return _wear;}
    std::string const & material() const {return _material;}

    // Genreflex can't do persistency of vector<LeadBrick> without a
    // default constructor
    LeadBrick() {}

  private:

    CLHEP::Hep3Vector worn(CLHEP::Hep3Vector const & d) const {
      return CLHEP::Hep3Vector(d.x() - 2.*_wear,
                               d.y() - 2.*_wear,
                               d.z() - 2.*_wear);
    }

    CLHEP::Hep3Vector _dim2x4x8;
    CLHEP::Hep3Vector _dim2x4x16;
    double            _wear;
    std::string       _material;
  };

}

#endif/*STMGeom_LeadBrick_hh*/
