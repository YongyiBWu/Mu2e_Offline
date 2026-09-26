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
// "wear" is taken off each face, so a brick is (dx - 2*wear) on a side
// while staying centred where it was placed. A stack therefore keeps
// its nominal pitch and the wear opens as gaps between bricks rather
// than displacing anything. The as-delivered sizes are exact imperial.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LeadBrick {
  public:

    LeadBrick(std::vector<std::string> const & names,
              std::vector<CLHEP::Hep3Vector> const & dims,
              double wear,
              std::string const & material
              ) :
      _names(names),
      _dims(dims),
      _wear(wear),
      _material(material)
    {
    }

    // How many sizes are defined. Types are numbered from 1.
    size_t nTypes() const {return _dims.size();}

    // The name a type was given in the geometry file, e.g. "2x4x16".
    std::string const & name(int type) const {return _names.at(type-1);}

    // Nominal outside dimensions of a type, before wear.
    CLHEP::Hep3Vector const & dim(int type) const {return _dims.at(type-1);}

    // The same with the wear taken off each face. This is what the
    // G4Box should be built from.
    CLHEP::Hep3Vector worn(int type) const {
      CLHEP::Hep3Vector const & d = _dims.at(type-1);
      return CLHEP::Hep3Vector(d.x() - 2.*_wear,
                               d.y() - 2.*_wear,
                               d.z() - 2.*_wear);
    }

    double              wear()     const {return _wear;}
    std::string const & material() const {return _material;}

    // Genreflex can't do persistency of vector<LeadBrick> without a
    // default constructor
    LeadBrick() {}

  private:

    std::vector<std::string>       _names;
    std::vector<CLHEP::Hep3Vector> _dims;
    double                         _wear;
    std::string                    _material;
  };

}

#endif/*STMGeom_LeadBrick_hh*/
