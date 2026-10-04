#ifndef STMGeom_LeadBrick_hh
#define STMGeom_LeadBrick_hh

// The standard lead bricks the STM shield house is stacked from.
//
// They are a shared primitive rather than a property of any one
// structure: the same sizes recur throughout the shielding, so their
// dimensions, material and wear live here and every structure that
// stacks them refers to this one object.
//
// The sizes are held as a numbered list rather than as named members,
// so that a course anywhere in the house can be written as a list of
// type numbers and a new size is a new entry rather than a change to
// every caller. A type is a SHAPE only -- it carries no orientation,
// because the same brick lies flat in one wall and stands on end in
// another, and no bores, because whether a hole passes through a given
// brick depends on where that brick sits.
//
// Wear is taken off each face, so a brick shrinks while staying
// centred where it was placed. A stack therefore keeps its nominal
// pitch and the wear opens as gaps between bricks rather than
// displacing anything. The as-delivered sizes are exact imperial.
//
// The wear is not isotropic: one value applies to whichever of the
// brick's own axes ends up vertical, another to the other two. That
// means it depends on how the brick is turned, so worn() takes the
// placement rotation and the construction code needs a solid per
// (type, orientation) pair rather than one per type.
//
// Author: Yongyi Wu

#include <cmath>
#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LeadBrick {
  public:

    LeadBrick(std::vector<std::string> const & names,
              std::vector<CLHEP::Hep3Vector> const & dims,
              double wearY,
              double wearXZ,
              std::string const & material
              ) :
      _names(names),
      _dims(dims),
      _wearY(wearY),
      _wearXZ(wearXZ),
      _material(material)
    {
    }

    // How many sizes are defined. Types are numbered from 1.
    size_t nTypes() const {return _dims.size();}

    // The name a type was given in the geometry file, e.g. "2x4x16".
    std::string const & name(int type) const {return _names.at(type-1);}

    // Nominal outside dimensions of a type, before wear.
    CLHEP::Hep3Vector const & dim(int type) const {return _dims.at(type-1);}

    // The same with the wear taken off each face, in the brick's OWN
    // frame. This is what the G4Box should be built from.
    //
    // Which of the brick's axes gets the vertical wear depends on how
    // it is turned, so the placement rotation is needed: the local
    // axis that rot maps closest to Mu2e y takes wearY and the other
    // two take wearXZ. For the 90 degree turns used here each local
    // axis maps exactly onto a Mu2e axis, so the choice is unambiguous.
    CLHEP::Hep3Vector worn(int type, CLHEP::HepRotation const & rot) const {
      CLHEP::Hep3Vector const & d = _dims.at(type-1);
      double w[3];
      for (int k = 0; k < 3; ++k) {
        const CLHEP::Hep3Vector axis =
          rot * CLHEP::Hep3Vector(k == 0, k == 1, k == 2);
        w[k] = (std::abs(axis.y()) > 0.5) ? _wearY : _wearXZ;
      }
      return CLHEP::Hep3Vector(d.x() - 2.*w[0],
                               d.y() - 2.*w[1],
                               d.z() - 2.*w[2]);
    }

    double              wearY()    const {return _wearY;}
    double              wearXZ()   const {return _wearXZ;}
    std::string const & material() const {return _material;}

    // Genreflex can't do persistency of vector<LeadBrick> without a
    // default constructor
    LeadBrick() {}

  private:

    std::vector<std::string>       _names;
    std::vector<CLHEP::Hep3Vector> _dims;
    double                         _wearY;
    double                         _wearXZ;
    std::string                    _material;
  };

}

#endif/*STMGeom_LeadBrick_hh*/
