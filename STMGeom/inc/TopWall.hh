#ifndef STMGeom_TopWall_hh
#define STMGeom_TopWall_hh

// The roof of the STM shield house.
//
// Three layers worked downward in -y from the top of the stack: a pair
// of borated poly sheets, a layer of lead bricks laid flat, and a pair
// of aluminium plates. 29 pieces, none of them bored.
//
// The lead layer is a BrickWall like the side walls', but laid flat:
// its columns run side by side in x and each column runs along z, so
// the 2 in of a brick is the layer's thickness in y. Two of the five
// columns start 3 in along z from the rest, which is why the columns
// carry their own offsets rather than sharing one origin.
//
// Every position arrives resolved from STMMaker. The wall is fixed by
// one reference -- inboard of the right wall in x, at that wall's
// reference z, and at the top of the stack in y -- so each piece
// follows from butting against its neighbour.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // A flat sheet or plate: the two poly sheets and the two aluminium
  // plates. No bores and no orientation -- halfDim is already in the
  // Mu2e frame, so which of the three is the thickness is read off it
  // rather than from a rotation applied afterwards.
  //
  // These are not all of one material, so each carries the name its
  // volume takes, set by STMMaker where which-piece-is-which is known.
  struct TopWallSheet {
    std::string       name;
    std::string       material;
    CLHEP::Hep3Vector halfDim;
    CLHEP::Hep3Vector center;
  };

  class TopWall {
  public:

    TopWall(bool build,
            std::vector<BrickWall> const & leadLayers,
            std::vector<TopWallSheet> const & sheets
            ) :
      _build(build),
      _leadLayers(leadLayers),
      _sheets(sheets)
    {
    }

    bool build() const {return _build;}

    // The lead. One layer, held as a vector so the construction code
    // walks it the same way it walks the side walls'.
    std::vector<BrickWall> const & leadLayers() const {return _leadLayers;}

    // The flat pieces, in the order the layers are built: the two poly
    // sheets, then the two aluminium plates.
    std::vector<TopWallSheet> const & sheets() const {return _sheets;}

    // Genreflex can't do persistency of vector<TopWall> without a
    // default constructor
    TopWall() {}

  private:

    bool _build;

    std::vector<BrickWall>    _leadLayers;
    std::vector<TopWallSheet> _sheets;
  };

}

#endif/*STMGeom_TopWall_hh*/
