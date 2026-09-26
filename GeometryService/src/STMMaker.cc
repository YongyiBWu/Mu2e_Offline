//
// Construct and return Stopping Target Monitor (STM)
//
// Author: Anthony Palladino
// Update: Haichuan Cao Sept 2023
//
// Notes
// See mu2e-doc-XXXX for naming conventions etc.

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <sstream>
#include <algorithm>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "Offline/GeneralUtilities/inc/OrientationResolver.hh"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Offline/GeometryService/inc/STMMaker.hh"
#include "Offline/STMGeom/inc/STM.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"

using namespace std;

namespace mu2e {

  namespace {

    // How far a box reaches along one direction once it has been
    // turned.
    //
    // The hand-stacked sections describe a piece by its dimensions in
    // its own frame plus an orientation code, so those dimensions do
    // not say how much of any Mu2e axis the piece occupies until the
    // rotation is applied. Reaching for dx because a piece "looks
    // long in x" is the recurring mistake here: a 2x4x16 brick laid on
    // its side spans 2 in along a course whose dx says 16, and the
    // copper lining's 12.7 mm depth is its outline's v range under one
    // orientation and its sweep under another.
    //
    // dim is in the piece's own frame; rot is the placement rotation;
    // dir is the Mu2e direction asked about, and should be a unit
    // vector.
    double spanAlong(CLHEP::Hep3Vector const & dim,
                     CLHEP::HepRotation const & rot,
                     CLHEP::Hep3Vector const & dir) {
      double span = 0.;
      for (int k = 0; k < 3; ++k) {
        const CLHEP::Hep3Vector axis =
          rot * CLHEP::Hep3Vector(k == 0, k == 1, k == 2);
        span += std::abs(axis.dot(dir)) * dim[k];
      }
      return span;
    }

  } // anonymous namespace

  // Constructor that gets information from the config file instead of
  // from arguments.
  STMMaker::STMMaker(SimpleConfig const & _config,
                     double solenoidOffset)
  {
    // if( ! _config.getBool("hasSTM",false) ) return;

    // create an empty STM
    _stm = unique_ptr<STM>(new STM());

    // access its object through a reference

    STM & stm = *_stm.get();

    parseConfig(_config);

    // One switch for the whole downstream geometry rather than a flag per
    // component; see STM::handstacked().
    stm._handstacked = _handstacked;

    // The shared standard lead bricks. 
    if (_handstacked) {
      stm._pLeadBrickParams = std::unique_ptr<LeadBrick>
              (new LeadBrick(_leadBrickNames,
                             _leadBrickDims,
                             _leadBrickWear,
                             _leadBrickMaterial));
    }

    // now create the specific components

    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;
    const CLHEP::Hep3Vector &dsP( ds->position() );

    //Create a reference position (most things in the STM geometry will be defined w.r.t. this position)
    //Our reference z is the downstream edge of the CRV-D as it was in crv_counters_v09.txt.
    //It is now set as a user variable in STM_v08.txt without changing its value.
    const CLHEP::Hep3Vector _STMMOffsetInMu2e(dsP.x(), 0.0, _stmReferenceZ );
    const CLHEP::HepRotation _magnetRotation = CLHEP::HepRotation::IDENTITY;
    double magnetZOffset = _magnetUpStrSpace+_magnetHalfLength;
    //calculate the magnet position assuming the shield pipe is flush to the wall
    if(_config.getBool("stm.magnet.usePipeAsOrigin", false)) {
      magnetZOffset = _magnetHalfLength + _shieldDnStrWallGap + 2.*_shieldPipeHalfLength + _shieldUpStrWallGap;
      if(!_shieldMatchPipeBlock) magnetZOffset += 2.*_shieldDnStrWallHalfLength;
    }
    const CLHEP::Hep3Vector _magnetOffsetInMu2e  = _STMMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,magnetZOffset);
    const CLHEP::Hep3Vector _magnetHoleOffset  = CLHEP::Hep3Vector(_magnetHoleXOffset,_magnetHoleYOffset, 0.);

    //if (_magnetBuild){
      stm._pSTMMagnetParams = std::unique_ptr<PermanentMagnet>
        (new PermanentMagnet(_magnetBuild,
                             _magnetHalfWidth,
                             _magnetHalfHeight,
                             _magnetHalfLength,
                             _magnetHoleHalfWidth,
                             _magnetHoleHalfHeight,
                             _magnetOffsetInMu2e,
                             _magnetRotation,
                             _magnetHoleOffset,
                             _magnetMaterial,
                             _magnetHasLiner,
                             _magnetField,
                             _magnetFieldVisible
                            ));
    //}


    const CLHEP::HepRotation _FOVCollimatorRotation     = CLHEP::HepRotation::IDENTITY;
    CLHEP::Hep3Vector  _FOVCollimatorOffsetInMu2e = _magnetOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_magnetHalfLength + _FOVCollimatorUpStrSpace + _FOVCollimatorHalfLength + _shieldPipeUpStrAirGap);
    // If we actually don't want to build the magnet, subtract off the offsets related to the magnet.
    // (We can't just set _magnetHalfLength = 0 in config because it is needed in various parts of constructSTM.cc
    // including to make a G4Box, which cannot have length 0...)
    if(_config.getBool("stm.magnet.build") == false) {
      _FOVCollimatorOffsetInMu2e -= CLHEP::Hep3Vector(0.0, 0.0, 2*_magnetHalfLength);
    }

    //if (_FOVCollimatorBuild){
      stm._pSTMFOVCollimatorParams = std::unique_ptr<STMCollimator>
        (new STMCollimator(_FOVCollimatorBuild,
                           _FOVCollimatorHalfWidth,
                           _FOVCollimatorHalfHeight,
                           _FOVCollimatorHalfLength,
                           _FOVCollimatorLinerBuild,
                           _FOVCollimatorLinerHalfWidth,
                           _FOVCollimatorLinerHalfHeight,
                           _FOVCollimatorLinerHalfLength,
                           _FOVCollimatorLinerCutOutHalfLength,
                           _FOVCollimatorHole1xOffset,
                           _FOVCollimatorHole1RadiusUpStr,
                           _FOVCollimatorHole1RadiusDnStr,
                           _FOVCollimatorHole1LinerBuild,
                           _FOVCollimatorHole1LinerThickness,
                           _FOVCollimatorHole2Build,
                           _FOVCollimatorHole2xOffset,
                           _FOVCollimatorHole2RadiusUpStr,
                           _FOVCollimatorHole2RadiusDnStr,
                           _FOVCollimatorHole2LinerBuild,
                           _FOVCollimatorHole2LinerThickness,
                           _FOVCollimatorOffsetInMu2e,
                           _FOVCollimatorRotation,
                           _FOVCollimatorMaterial,
                           _FOVCollimatorLinerMaterial,
                           _FOVCollimatorHoleLinerMaterial
                          ));
    //}


    const CLHEP::HepRotation _pipeRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _pipeOffsetInMu2e = _STMMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_magnetUpStrSpace+_magnetHalfLength);
    //if (_pipeBuild){
      stm._pSTMTransportPipeParams = std::unique_ptr<TransportPipe>
        (new TransportPipe(_pipeBuild,
                           _pipeRadiusIn,
                           _pipeRadiusOut,
                           _pipeMaterial,
                           _pipeGasMaterial,
                           _pipeUpStrSpace,
                           _pipeDnStrHalfLength,
                           _pipeUpStrWindowMaterial,
                           _pipeUpStrWindowHalfLength,
                           _pipeDnStrWindowMaterial,
                           _pipeDnStrWindowHalfLength,
                           _pipeFlangeHalfLength,
                           _pipeFlangeOverhangR,
                           _pipeOffsetInMu2e,
                           _pipeRotation
                          ));
    //}

    double _magnetTableTopHalfWidth = 0.0;
    if ( _magnetBuild && !_FOVCollimatorBuild) _magnetTableTopHalfWidth = _magnetHalfWidth;
    if (!_magnetBuild &&  _FOVCollimatorBuild) _magnetTableTopHalfWidth = _FOVCollimatorHalfWidth;
    if ( _magnetBuild &&  _FOVCollimatorBuild) _magnetTableTopHalfWidth = std::max(_magnetHalfWidth,_FOVCollimatorHalfWidth);
    _magnetTableTopHalfWidth += _magnetTableTopExtraWidth;

    double _magnetTableTopHalfLength = 0.0;
    if ( _magnetBuild )        _magnetTableTopHalfLength += _magnetHalfLength;
    if ( _pipeBuild )          _magnetTableTopHalfLength += _pipeDnStrHalfLength;
    if ( _FOVCollimatorBuild ) _magnetTableTopHalfLength += 0.5*_FOVCollimatorUpStrSpace+_FOVCollimatorHalfLength;
    if ( _shieldBuild ) {
      _magnetTableTopHalfLength += 0.5*_shieldDnStrSpace+_shieldDnStrWallHalfLength;
      if (!_magnetBuild) {
        // Previous versions (<STM_v09) didn't need to include the shield pipe length in here.
        // Now we will include it, otherwise the table is tiny
        _magnetTableTopHalfLength += _shieldPipeHalfLength;
      }
    }
    _magnetTableTopHalfLength += _magnetTableTopExtraLength;

    const CLHEP::HepRotation _magnetTableRotation     = CLHEP::HepRotation::IDENTITY;
    CLHEP::Hep3Vector  _magnetTableOffsetInMu2e = _STMMOffsetInMu2e - CLHEP::Hep3Vector(0.0,_magnetHalfHeight+_magnetTableTopHalfHeight,0.0);
    if ( _magnetBuild )        _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,_magnetUpStrSpace+_magnetHalfLength);
    if ( _pipeBuild )          _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,_pipeDnStrHalfLength);
    if ( _FOVCollimatorBuild ) _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,0.5*_FOVCollimatorUpStrSpace+_FOVCollimatorHalfLength);
    if ( _shieldBuild ) {
      _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,-0.5*_shieldDnStrSpace-_shieldDnStrWallHalfLength + _shieldPipeUpStrAirGap);
      if (!_magnetBuild) {
        // Previous versions (<STM_v09) didn't need to include the shield pipe length in here.
        // Now we will include it, otherwise the table is tiny
        _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,_shieldPipeHalfLength+0.5*_shieldDnStrSpace+_shieldDnStrWallHalfLength+_shieldPipeUpStrAirGap);
      }
    }


    //if (_magnetTableBuild && (_magnetBuild||_FOVCollimatorBuild) ){
      stm._pSTMMagnetSupportTableParams = std::unique_ptr<SupportTable>
        (new SupportTable( _magnetTableBuild,
                           _magnetTableTopHalfWidth,
                           _magnetTableTopHalfHeight,
                           _magnetTableTopHalfLength,
                           _magnetTableLegRadius,
                           _magnetTableOffsetInMu2e,
                           _magnetTableRotation,
                           _magnetTableMaterial
                          ));
    //}


      ////////////////////////////////////
      // STM Downstream Area
      //
    //The STM geometry must fit inside the detector hall, so find the z of the East hall wall
    GeomHandle<Mu2eHall> hall;
    const double z_hall_inside_max = hall->getWallExtentz("dsArea",1)/CLHEP::mm;//the integer allows you to specify which side of which wall you want the z for: 1 = west side of east wall (i.e. the z of the inside surface of the east wall)
    const CLHEP::Hep3Vector BeamAxisAtEastWallInMu2e(dsP.x(), 0.0, z_hall_inside_max );
    const double yExtentLow = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
    const CLHEP::Hep3Vector FloorAtEastWallInMu2e = BeamAxisAtEastWallInMu2e - CLHEP::Hep3Vector(0.0, yExtentLow, 0.0);

    // Define the envelope w.r.t the floor at the east wall
    const CLHEP::HepRotation _stmDnStrEnvRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector _stmDnStrEnvPositionInMu2e = FloorAtEastWallInMu2e + CLHEP::Hep3Vector(0.0, +_stmDnStrEnvHalfHeight, -_stmDnStrEnvHalfLength);
    stm._pSTMDnStrEnvParams = std::unique_ptr<STMDownstreamEnvelope>
      (new STMDownstreamEnvelope(_stmDnStrEnvBuild,
                                 _stmDnStrEnvHalfWidth,
                                 _stmDnStrEnvHalfHeight,
                                 _stmDnStrEnvHalfLength,
                                 _stmDnStrEnvPositionInMu2e,
                                 _stmDnStrEnvRotation,
                                 _stmDnStrEnvMaterial
                                 ));

    const CLHEP::HepRotation _SSCollimatorRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _SSCollimatorOffsetInMu2e = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0.0,0.,-_stmZAllowed+_SSCollimatorHalfLength);
    //if(_SSCollimatorBuild){
      stm._pSTMSSCollimatorParams = std::unique_ptr<STMCollimator>
        (new STMCollimator(_SSCollimatorBuild,
                           _SSCollimatorHalfWidth,
                           _SSCollimatorHalfHeight,
                           _SSCollimatorHalfLength,
                           _SSCollimatorLinerBuild,
                           _SSCollimatorLinerHalfWidth,
                           _SSCollimatorLinerHalfHeight,
                           _SSCollimatorLinerHalfLength,
                           _SSCollimatorLinerCutOutHalfLength,
                           _SSCollimatorHole1xOffset,
                           _SSCollimatorHole1RadiusUpStr,
                           _SSCollimatorHole1RadiusDnStr,
                           _SSCollimatorHole1LinerBuild,
                           _SSCollimatorHole1LinerThickness,
                           _SSCollimatorHole2Build,
                           _SSCollimatorHole2xOffset,
                           _SSCollimatorHole2RadiusUpStr,
                           _SSCollimatorHole2RadiusDnStr,
                           _SSCollimatorHole2LinerBuild,
                           _SSCollimatorHole2LinerThickness,
                           _SSCollimatorOffsetInMu2e,
                           _SSCollimatorRotation,
                           _SSCollimatorMaterial,
                           _SSCollimatorLinerMaterial,
                           _SSCollimatorHoleLinerMaterial
                          ));
    //}


    double _detectorTableTopHalfWidth = _SSCollimatorHalfWidth + _detectorTableTopExtraWidth;
    double _detectorTableTopHalfLength = 0.5*_stmZAllowed - 1.0;
    const CLHEP::HepRotation _detectorTableRotation = CLHEP::HepRotation::IDENTITY;
    CLHEP::Hep3Vector  _detectorTableOffsetInMu2e = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0.0,-_SSCollimatorHalfHeight-_detectorTableTopHalfHeight, -_stmZAllowed+_detectorTableTopHalfLength);

    //if (_detectorTableBuild && _SSCollimatorBuild ){
      stm._pSTMDetectorSupportTableParams = std::unique_ptr<SupportTable>
        (new SupportTable( _detectorTableBuild,
                           _detectorTableTopHalfWidth,
                           _detectorTableTopHalfHeight,
                           _detectorTableTopHalfLength,
                           _detectorTableLegRadius,
                           _detectorTableOffsetInMu2e,
                           _detectorTableRotation,
                           _detectorTableMaterial
                          ));
    //}


    const CLHEP::HepRotation _detector1Rotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _detector1OffsetInMu2e  = _SSCollimatorOffsetInMu2e + CLHEP::Hep3Vector(_detector1xOffset, 0.0, _SSCollimatorHalfLength+_detector1CanUpStrSpace+_detector1CanHalfLength);
    //if (_detector1Build){
      stm._pSTMDetector1Params = std::unique_ptr<GeDetector>
        (new GeDetector(_detector1Build,
                        _detector1CrystalMaterial,
                        _detector1CrystalRadiusIn,
                        _detector1CrystalRadiusOut,
                        _detector1CrystalHalfLength,
                        _detector1CanMaterial,
                        _detector1CanRadiusIn,
                        _detector1CanRadiusOut,
                        _detector1CanHalfLength,
                        _detector1CanUpStrWindowMaterial,
                        _detector1CanUpStrWindowHalfLength,
                        _detector1CanGasMaterial,
                        _detector1OffsetInMu2e,
                        _detector1Rotation
                       ));
      //}

    const CLHEP::HepRotation _detector2Rotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _detector2OffsetInMu2e  = _SSCollimatorOffsetInMu2e + CLHEP::Hep3Vector(_detector2xOffset, 0.0, _SSCollimatorHalfLength+_detector2CanUpStrSpace+_detector2CanHalfLength);
    //if (_detector2Build){
      stm._pSTMDetector2Params = std::unique_ptr<GeDetector>
        (new GeDetector(_detector2Build,
                        _detector2CrystalMaterial,
                        _detector2CrystalRadiusIn,
                        _detector2CrystalRadiusOut,
                        _detector2CrystalHalfLength,
                        _detector2CanMaterial,
                        _detector2CanRadiusIn,
                        _detector2CanRadiusOut,
                        _detector2CanHalfLength,
                        _detector2CanUpStrWindowMaterial,
                        _detector2CanUpStrWindowHalfLength,
                        _detector2CanGasMaterial,
                        _detector2OffsetInMu2e,
                        _detector2Rotation
                       ));
      //}


    const CLHEP::HepRotation _shieldRotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _shieldOffsetInMu2e  = _FOVCollimatorOffsetInMu2e + CLHEP::Hep3Vector(0.0, 0.0, -_FOVCollimatorHalfLength);
    //if (_shieldBuild){
      stm._pSTMShieldPipeParams = std::unique_ptr<ShieldPipe>
        (new ShieldPipe(_shieldBuild,
                        _shieldRadiusIn,
                        _shieldHasLiner,
                        _shieldLinerWidth,
                        _shieldRadiusOut,
                        _shieldPipeHalfLength,
                        _shieldMaterialLiner,
                        _shieldMaterial,
                        _shieldMatchPipeBlock,
                        _shieldUpStrSpace,
                        _shieldDnStrSpace,
                        _shieldDnStrWallHalfLength,
                        _shieldDnStrWallHoleRadius,
                        _shieldDnStrWallHalfHeight,
                        _shieldDnStrWallHalfWidth,
                        _shieldDnStrWallGap,
                        _shieldDnStrWallMaterial,
                        _shieldBuildMatingBlock,
                        _shieldPipeUpStrAirGap,
                        _shieldOffsetInMu2e, //This is upstream edge of FOV collimator for now.
                        _shieldRotation
                       ));
    //}

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /// The geometries below were updated by Haichuan Cao in Sept. 2023

    const CLHEP::Hep3Vector  _STMShieldingRef = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0., 0., -_STM_SSCFrontToWall);

   ////////////////////////////////////////////////////////////////
   //Tungsten Spot-Size Collimator

    const CLHEP::HepRotation _STM_SSCRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _STM_SSCOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(0, 0, _STM_SSCWdepth_f/2);

    if (_handstacked) {
      // The updated SSC: a single symmetric block, bores stepped in z.
          stm._pSTM_SSCParams = std::unique_ptr<STM_SSC>
                  (new STM_SSC(_STM_SSCBuild,
                               _STM_SSCVDBuild,
                               _STM_SSCW_width,
                               _STM_SSCW_height,
                               _STM_SSCWdepth_f,
                               _STM_SSCWdepth_b,
                               _STM_SSCr_LaBr_f,
                               _STM_SSCr_HPGe_f,
                               _STM_SSCr_LaBr_b,
                               _STM_SSCr_HPGe_b,
                               _STM_SSCoffset_Spot,
                               _STM_SSCleak,
                               _STM_SSCFrontToWall,
                               _STM_SSCZGap,
                               _STM_SSCZGapBack,
                               _STM_SSCboreToBase,
                               _STM_SSCOffsetInMu2e,
                               _STM_SSCRotation,
                               _STM_SSCMaterial));
    } else {
          stm._pSTM_SSCParams = std::unique_ptr<STM_SSC>
                  (new STM_SSC(_STM_SSCBuild,
                               _STM_SSCVDBuild,
                               _STM_SSCdelta_WlR,
                               _STM_SSCdelta_WlL,
                               _STM_SSCW_middle,
                               _STM_SSCW_height,
                               _STM_SSCWdepth_f,
                               _STM_SSCWdepth_b,
                               _STM_SSCAperture_HPGe1,
                               _STM_SSCAperture_HPGe2,
                               _STM_SSCAperture_LaBr1,
                               _STM_SSCAperture_LaBr2,
                               _STM_SSCoffset_Spot,
                               _STM_SSCleak,
                               _STM_SSCFrontToWall,
                               _STM_SSCZGap,
                               _STM_SSCZGapBack,
                               _STM_SSCOffsetInMu2e,
                               _STM_SSCRotation,
                               _STM_SSCMaterial));
    }

   ////////////////////////////////////////////////////////////////
   //Spot-Size Collimator Support

    const CLHEP::HepRotation _SSCSupportRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _SSCSupportOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(0, 0, _STM_SSCZGap/2);

    if (_handstacked) {
      // The updated support: a five-plate steel cradle around the collimator.
          stm._pSSCSupportParams = std::unique_ptr<SSCSupport>
                  (new SSCSupport(_SSCSupportBuild,
                               _SSCSupportdepth,
                               _SSCSupportside_T,
                               _SSCSupportside_H,
                               _SSCSupportplate_base_T,
                               _SSCSupportbottom_T,
                               _SSCSupporttop_T,
                               _SSCSupportMaterial,
                               _SSCSupportOffsetInMu2e,
                               _SSCSupportRotation));
    } else {
          stm._pSSCSupportParams = std::unique_ptr<SSCSupport>
                  (new SSCSupport(_SSCSupportBuild,
                               _SSCSupporttable_L,
                               _SSCSupporttable_H,
                               _SSCSupporttable_T,
                               _SSCSupportleg_L,
                               _SSCSupportleg_H,
                               _SSCSupportleg_T,
                               _SSCSupportbase_L,
                               _SSCSupportbase_H,
                               _SSCSupportbase_T,
                               _SSCSupportwall_L,
                               _SSCSupportwall_H,
                               _SSCSupportwall_T,
                               _SSCSupporthole_H,
                               _SSCSupporthole_T,
                               _SSCSupportFLeadStand_L,
                               _SSCSupportFLeadStand_H,
                               _SSCSupportFLeadStand_T,
                               _SSCSupportFLeadShim_H,
                               _SSCSupportFLeadShim_T,
                               _SSCSupportFAluminumShim_T,
                               _SSCSupportFAluminumExtra_L,
                               _SSCSupportFAluminumExtra_H,
                               _SSCSupportOffsetInMu2e,
                               _SSCSupportRotation));
    }

   ////////////////////////////////////////////////////////////////
   //STM Front Shielding
   //
   // The updated geometry splits this into a left and a right half.
   // The right one -- the section behind the collimator -- is built
   // here in place of the earlier single description; its lead courses
   // and poly sheets were expanded from the config's type lists in
   // parseConfig, so what is passed on are finished placements.

    // The front shielding's depth and width, used further down by the
    // Bottom, Left, Right and Top sections. Those are still the earlier
    // description and are due to be replaced, so these stay zero under
    // the updated geometry and are set only by the branch that needs
    // them.
    double _FrontS_Thickness = 0.;
    double _FrontS_Length    = 0.;

    if (_handstacked) {
      stm._pFrontShieldingRightParams = std::unique_ptr<FrontShieldingRight>
              (new FrontShieldingRight(_FrontShieldingRightBuild,
                                       _FrontShieldingRightLeadLayers,
                                       _FrontShieldingRightSheets,
                                       _FrontShieldingRightBlockMaterial,
                                       _FrontShieldingRightBlockHalfDim,
                                       _FrontShieldingRightBlockCenter,
                                       _FrontShieldingRightPipes,
                                       _FrontShieldingRightPlate,
                                       _FrontShieldingRightBackZ));

      stm._pFrontShieldingLeftParams = std::unique_ptr<FrontShieldingLeft>
              (new FrontShieldingLeft(_FrontShieldingLeftBuild,
                                      _FrontShieldingLeftSheets,
                                      _FrontShieldingLeftBrickGroups,
                                      _FrontShieldingLeftPrism));
    } else {

      const CLHEP::HepRotation _FrontSRotation   = CLHEP::HepRotation::IDENTITY;
      _FrontS_Thickness = _FrontStungstendepth + _FrontSLeakForSSC + _FrontSleaddepth2*3 + _FrontSBPdepth*2 + _FrontScopperdepth;
      _FrontS_Length    = _FrontStungstenlength + 2*_FrontSLeakForSSC + _FrontSfPb_lengthL + _FrontSfPb_lengthR;

      const CLHEP::Hep3Vector  _FrontSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(0, 0, _FrontS_Thickness/2);

          stm._pSTMFrontShieldingParams = std::unique_ptr<FrontShielding>
          (new FrontShielding(_FrontShieldingBuild,
                              _FrontSHeightofRoom,
                              _FrontStungstenlength,
                              _FrontStungstendepth,
                              _FrontSleaddepth1,
                              _FrontSleaddepth2,
                              _FrontSaluminumdepth,
                              _FrontScopperdepth,
                              _FrontSBPdepth,
                              _FrontSfPb_lengthL,
                              _FrontSfPb_lengthR,
                              _FrontSGapForTop,
                              _FrontSLeakForSSC,
                              _FrontSCopperL,
                              _FrontS_H,
                              _FrontSHole_r,
                              _FrontSOffsetInMu2e,
                              _FrontSRotation));
    }

   ////////////////////////////////////////////////////////////////
   //STM HPGe Detector

    CLHEP::HepRotation rotHPGe = CLHEP::HepRotation::IDENTITY;
    rotHPGe.rotateY(45*CLHEP::degree);

    const CLHEP::HepRotation _HPGeRotation     =  rotHPGe;
    const CLHEP::Hep3Vector  _HPGeOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(-_STM_SSCoffset_Spot + _HPGeoffset_HPGe, 0., _HPGeZ_HPGe);
          stm._pSTMHPGeDetectorParams = std::unique_ptr<HPGeDetector>
          (new HPGeDetector(_HPGeBuild,
                            _HPGecrystalMaterial,
                            _HPGeholeMaterial,
                            _HPGewindowMaterial,
                            _HPGewallMaterial,
                            _HPGecapsuleMaterial,
                            _HPGeEndcapR,
                            _HPGeEndcapL,
                            _HPGeCrystalR,
                            _HPGeCrystalL,
                            _HPGeZ_HPGe,
                            _HPGeHoleR,
                            _HPGeHoleL,
                            _HPGeCapsule_Wallthick,
                            _HPGeCapsule_Windowthick,
                            _HPGeCapsule_Endthick,
                            _HPGeCapsule_Walllength,
                            _HPGeWindowD,
                            _HPGeEndcapD,
                            _HPGeAirD,
                            _HPGeoffset_HPGe,
                            _HPGeOffsetInMu2e,
                            _HPGeRotation));


   ////////////////////////////////////////////////////////////////
   //STM LaBr Detector

    const CLHEP::HepRotation _LaBrRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _LaBrOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(_STM_SSCoffset_Spot + _LaBroffset_LaBr, 0., _LaBrZ_LaBr);
          stm._pSTMLaBrDetectorParams = std::unique_ptr<LaBrDetector>
          (new LaBrDetector(_LaBrBuild,
                            _LaBrcrystalMaterial,
                            _LaBrwindowMaterial,
                            _LaBrwallMaterial,
                            _LaBrEndcapR,
                            _LaBrEndcapL,
                            _LaBrCrystalR,
                            _LaBrCrystalL,
                            _LaBrZ_LaBr,
                            _LaBrWindowD,
                            _LaBrEndcapD,
                            _LaBrAirD,
                            _LaBroffset_LaBr,
                            _LaBrOffsetInMu2e,
                            _LaBrRotation));


   ////////////////////////////////////////////////////////////////
   //STM Bottom Shielding

      const double _BottomS_Thickness = _BottomSleaddepth*2 + _BottomSBPdepth*2 + _BottomScopperdepth;

      const double B_dX = - _FrontS_Length - _BottomSfloor_Zlength/4 + 264.12;
      const double B_dY = -_FrontSHeightofRoom/2 - _BottomS_Thickness/2;
      const double B_dZ = _FrontS_Thickness - 12.7 + _BottomSfloor_Zlength/2;

      const CLHEP::HepRotation _BottomSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _BottomSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(B_dX, B_dY, B_dZ);


          stm._pSTMBottomShieldingParams = std::unique_ptr<BottomShielding>
          (new BottomShielding(_BottomShieldingBuild,
                               _BottomSfloor_Zlength,
                               _BottomSFront_LB,
                               _BottomSFront_LB_inner,
                               _BottomSleaddepth,
                               _BottomScopperdepth,
                               _BottomSBPdepth,
                               _BottomSOffsetInMu2e,
                               _BottomSRotation));


   ////////////////////////////////////////////////////////////////
   //STM Left Shielding

      const double _LeftS_Thickness = _LeftSleaddepth*2 + _LeftSBPdepth*2 + _LeftScopperdepth;

      const double L_dX = _LeftSXmin + _LeftS_Thickness/2;
      const double L_dY = 0;
      const double L_dZ = _FrontS_Thickness + _LeftS_Length/2;

      const CLHEP::HepRotation _LeftSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _LeftSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(L_dX, L_dY, L_dZ);

          stm._pSTMLeftShieldingParams = std::unique_ptr<LeftShielding>
          (new LeftShielding(_LeftShieldingBuild,
                             _LeftS_Length,
                             _LeftSleaddepth,
                             _LeftScopperdepth,
                             _LeftSBPdepth,
                             _LeftSXmin,
                             _LeftSOffsetInMu2e,
                             _LeftSRotation));


   ////////////////////////////////////////////////////////////////
   //STM Right Shielding

      const double _RightS_Thickness = (_RightSleaddepth*2 + _RightSBPdepth*2 + _RightScopperdepth)*sqrt(2);

      const double R_dX = _RightSXmax - _RightS_Thickness/2;
      const double R_dY = 0;
      const double R_dZ = _FrontS_Thickness + _RightS_Length/2;

      const CLHEP::HepRotation _RightSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _RightSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(R_dX, R_dY, R_dZ);

          stm._pSTMRightShieldingParams = std::unique_ptr<RightShielding>
          (new RightShielding(_RightShieldingBuild,
                              _RightS_Length,
                              _RightSleaddepth,
                              _RightScopperdepth,
                              _RightSBPdepth,
                              _RightSXmax,
                              _RightSOffsetInMu2e,
                              _RightSRotation));


   ////////////////////////////////////////////////////////////////
   //STM Top Shielding

      const double _TopS_Thickness = _TopSleaddepth*2 + _TopSBPdepth*2 + _TopScopperdepth;

      const double _T_Xlength1 = _TopSFront_LT + _TopSBarLeft + _TopSBarLeft + _TopSGapRight + _TopSGapRight;

      const double T_dX = -(_T_Xlength1 + _TopSXlength)/4 + _TopSBarLeft + _TopSGapLeft + _TopSFront_LT + (_STM_SSCdelta_WlL - _STM_SSCdelta_WlR)/2;
      const double T_dY = _FrontSHeightofRoom/2 +  _TopS_Thickness/2;
      const double T_dZ = _FrontS_Thickness + _TopSZlength/2 - _TopSZHole;

      const CLHEP::HepRotation _TopSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _TopSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(T_dX, T_dY, T_dZ);


          stm._pSTMTopShieldingParams = std::unique_ptr<TopShielding>
          (new TopShielding(_TopShieldingBuild,
                            _TopShieldingSkirtBuild,
                            _TopLiftBeam_L,
                            _TopLiftBeam_H,
                            _TopLiftBeam_T,
                            _TopLiftBeam_Xmove,
                            _TopSZlength,
                            _TopSXlength,
                            _TopSFront_LT,
                            _TopTFZlength,
                            _TopTFXlength,
                            _TopTBZlength,
                            _TopScontainerdepth,
                            _TopSleaddepth,
                            _TopScopperdepth,
                            _TopSBPdepth,
                            _TopSZHole,
                            _TopSBarLeft,
                            _TopSBarRight,
                            _TopSGapLeft,
                            _TopSGapRight,
                            _TopSLeak,
                            _TopSOffsetInMu2e,
                            _TopSRotation));


   ////////////////////////////////////////////////////////////////
   //STM Inner Shielding
          stm._pSTMInnerShieldingParams = std::unique_ptr<InnerShielding>
          (new InnerShielding(_InnerShieldingBuild));

   ////////////////////////////////////////////////////////////////
   //STM Back Shielding


      const double Back_dZ = _BottomSfloor_Zlength + _BackSBPThick/2;

      const CLHEP::HepRotation _BackSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _BackSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(_BackS_dX, _BackS_dY, Back_dZ);

          stm._pSTMBackShieldingParams = std::unique_ptr<BackShielding>
          (new BackShielding(_BackShieldingBuild,
                              _BackSBPThick,
                              _BackSBPLength,
                              _BackSBPHeight,
                              _BackS_dX,
                              _BackS_dY,
                              _BackSPipeGap,
                              _BackSOffsetInMu2e,
                              _BackSRotation));

   ////////////////////////////////////////////////////////////////
   //STM Electronic Shielding
          stm._pSTMElectronicShieldingParams = std::unique_ptr<ElectronicShielding>
          (new ElectronicShielding(_ElectronicShieldingBuild,
                                   _ElectronicSSiGridX,
                                   _ElectronicSSiGridY,
                                   _ElectronicSSiGridZ,
                                   _ElectronicSSiXcenter,
                                   _ElectronicSSiYcenter,
                                   _ElectronicSSiZcenter,
                                   _ElectronicSConcreteT,
                                   _ElectronicSGapToSi));

   ////////////////////////////////////////////////////////////////
   //STM Absorber Shielding, or the SSC front shield that replaces it

    if (_handstacked) {
      // The updated geometry puts a stacked lead brick wall, a shelf and
      // two poly blocks where the earlier description had a single
      // absorber block.
      stm._pSSCFrontShieldParams = std::unique_ptr<SSCFrontShield>
              (new SSCFrontShield(_SSCFrontShieldBuild,
                                  _SSCFrontShieldBrick2x4x8Center,
                                  _SSCFrontShieldBrick2x4x8Orientation,
                                  _SSCFrontShieldBrick2x4x16Center,
                                  _SSCFrontShieldBrick2x4x16Orientation,
                                  _SSCFrontShieldShelfMaterial,
                                  _SSCFrontShieldShelfDim,
                                  _SSCFrontShieldShelfCenter,
                                  _SSCFrontShieldPoly1Material,
                                  _SSCFrontShieldPoly1Dim,
                                  _SSCFrontShieldPoly1Center,
                                  _SSCFrontShieldPoly1BoreR,
                                  _SSCFrontShieldPoly1BoreDX,
                                  _SSCFrontShieldPoly1BoreDY,
                                  _SSCFrontShieldPoly2Material,
                                  _SSCFrontShieldPoly2Dim,
                                  _SSCFrontShieldPoly2Center));
    } else {
          stm._pSTMSTM_AbsorberParams = std::unique_ptr<STM_Absorber>
          (new STM_Absorber(_STM_AbsorberBuild,
                            _STM_Absorber_hW,
                            _STM_Absorber_hH,
                            _STM_Absorber_hT,
                            _STM_Absorber_GaptoSSC));
    }

  }

  void STMMaker::parseConfig( SimpleConfig const & _config ){

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "stmMagnetField", "stm.magnet.field");


    _verbosityLevel            = _config.getInt("stm.verbosityLevel",0);
    _stmZAllowed               = _config.getDouble("stm.z.allowed");

    // The updated, hand-stacked downstream shielding, taken from the NX STEP
    // export of the shield house. One switch for the whole downstream
    // geometry: those pieces are dimensioned against a shared datum and a
    // shared stack chain, so they ship together rather than one at a time.
    _handstacked               = _config.getBool("stm.downstream.handstacked",false);

    // The standard lead bricks. Read here rather than with any one
    // structure: most of the shield house is stacked from them, so they
    // are a shared primitive that many components refer to.
    if (_handstacked) {
      // Read by number, so adding a size to the geometry file is
      // enough: nothing here or downstream names a particular one.
      const int nBrickTypes = _config.getInt("stm.leadBrick.typeN");
      for (int t = 1; t <= nBrickTypes; ++t) {
        const std::string name =
          _config.getString("stm.leadBrick.type" + std::to_string(t));
        _leadBrickNames.push_back(name);
        _leadBrickDims.push_back(
          CLHEP::Hep3Vector(_config.getDouble("stm.leadBrick." + name + ".dx"),
                            _config.getDouble("stm.leadBrick." + name + ".dy"),
                            _config.getDouble("stm.leadBrick." + name + ".dz")));
      }
      _leadBrickWear      = _config.getDouble("stm.leadBrick.wear");
      _leadBrickMaterial  = _config.getString("stm.leadBrick.material");
    }

    _stmReferenceZ             = _config.getDouble("stm.referenceZ");  //was previously calculated automatically based on the location of the CRV-D

    _magnetBuild               = _config.getBool(  "stm.magnet.build",false);
    _magnetUpStrSpace          = _config.getDouble("stm.magnet.UpStrSpace");
    _magnetHalfLength          = _config.getDouble("stm.magnet.halfLength");
    _magnetHalfWidth           = _config.getDouble("stm.magnet.halfWidth");
    _magnetHalfHeight          = _config.getDouble("stm.magnet.halfHeight");
    _magnetHoleHalfWidth       = _config.getDouble("stm.magnet.holeHalfWidth");
    _magnetHoleHalfHeight      = _config.getDouble("stm.magnet.holeHalfHeight");
    _magnetHoleXOffset         = _config.getDouble("stm.magnet.holeXOffset", 0.);
    _magnetHoleYOffset         = _config.getDouble("stm.magnet.holeYOffset", 0.);
    _magnetMaterial            = _config.getString("stm.magnet.material");
    _magnetHasLiner            = _config.getBool("stm.magnet.hasLiner", true);
    _magnetField               = _config.getDouble("stm.magnet.field");
    //_magnetFieldVisible        = _config.getBool(  "stm.magnet.fieldVisible",false);
    _magnetFieldVisible        = geomOptions->isVisible("stmMagnetField");

    _FOVCollimatorBuild            = _config.getBool(  "stm.FOVcollimator.build");
    _FOVCollimatorMaterial         = _config.getString("stm.FOVcollimator.material");
    _FOVCollimatorUpStrSpace       = _config.getDouble("stm.FOVcollimator.UpStrSpace");
    _FOVCollimatorHalfWidth        = _config.getDouble("stm.FOVcollimator.halfWidth");
    _FOVCollimatorHalfHeight       = _config.getDouble("stm.FOVcollimator.halfHeight");
    _FOVCollimatorHalfLength       = _config.getDouble("stm.FOVcollimator.halfLength");
    _FOVCollimatorLinerBuild       = _config.getBool(  "stm.FOVcollimator.liner.build");
    _FOVCollimatorLinerMaterial    = _config.getString("stm.FOVcollimator.liner.material");
    _FOVCollimatorLinerHalfWidth   = _config.getDouble("stm.FOVcollimator.liner.halfWidth");
    _FOVCollimatorLinerHalfHeight  = _config.getDouble("stm.FOVcollimator.liner.halfHeight");
    _FOVCollimatorLinerHalfLength  = _config.getDouble("stm.FOVcollimator.liner.halfLength");
    _FOVCollimatorLinerCutOutHalfLength  = _config.getDouble("stm.FOVcollimator.liner.cutOutHalfLength");
    _FOVCollimatorHole1xOffset     = _config.getDouble("stm.FOVcollimator.hole1.xoffset");
    _FOVCollimatorHole1RadiusUpStr = _config.getDouble("stm.FOVcollimator.hole1.radiusUpStr");
    _FOVCollimatorHole1RadiusDnStr = _config.getDouble("stm.FOVcollimator.hole1.radiusDnStr");
    _FOVCollimatorHole1LinerBuild     = _config.getBool(  "stm.FOVcollimator.hole1.liner.build");
    _FOVCollimatorHole1LinerThickness = _config.getDouble("stm.FOVcollimator.hole1.liner.thickness");
    _FOVCollimatorHole2Build       = _config.getBool(  "stm.FOVcollimator.hole2.build");
    _FOVCollimatorHole2xOffset     = _config.getDouble("stm.FOVcollimator.hole2.xoffset");
    _FOVCollimatorHole2RadiusUpStr = _config.getDouble("stm.FOVcollimator.hole2.radiusUpStr");
    _FOVCollimatorHole2RadiusDnStr = _config.getDouble("stm.FOVcollimator.hole2.radiusDnStr");
    _FOVCollimatorHole2LinerBuild     = _config.getBool(  "stm.FOVcollimator.hole2.liner.build");
    _FOVCollimatorHole2LinerThickness = _config.getDouble("stm.FOVcollimator.hole2.liner.thickness");
    _FOVCollimatorHoleLinerMaterial= _config.getString("stm.FOVcollimator.hole.liner.material");

    _pipeBuild                 = _config.getBool(  "stm.pipe.build");
    _pipeRadiusIn              = _config.getDouble("stm.pipe.rIn");
    _pipeRadiusOut             = _config.getDouble("stm.pipe.rOut");
    _pipeMaterial              = _config.getString("stm.pipe.material");
    _pipeGasMaterial           = _config.getString("stm.pipe.gas.material");
    _pipeUpStrSpace            = _config.getDouble("stm.pipe.UpStrSpace");
    _pipeDnStrHalfLength       = _config.getDouble("stm.pipe.DnStrHalfLength");
    _pipeUpStrWindowMaterial   = _config.getString("stm.pipe.UpStrWindow.material");
    _pipeUpStrWindowHalfLength = _config.getDouble("stm.pipe.UpStrWindow.halfLength");
    _pipeDnStrWindowMaterial   = _config.getString("stm.pipe.DnStrWindow.material");
    _pipeDnStrWindowHalfLength = _config.getDouble("stm.pipe.DnStrWindow.halfLength");
    _pipeFlangeHalfLength      = _config.getDouble("stm.pipe.flange.halfLength");
    _pipeFlangeOverhangR       = _config.getDouble("stm.pipe.flange.overhangR");

    _magnetTableBuild          = _config.getBool(  "stm.magnet.stand.build",false);
    _magnetTableMaterial       = _config.getString("stm.magnet.stand.material");
    _magnetTableTopExtraWidth  = _config.getDouble("stm.magnet.stand.topExtraWidth");
    _magnetTableTopExtraLength = _config.getDouble("stm.magnet.stand.topExtraLength");
    _magnetTableTopHalfHeight  = _config.getDouble("stm.magnet.stand.topHalfHeight");
    _magnetTableLegRadius      = _config.getDouble("stm.magnet.stand.legRadius");

    _SSCollimatorBuild            = _config.getBool(  "stm.SScollimator.build");
    _SSCollimatorMaterial         = _config.getString("stm.SScollimator.material");
    _SSCollimatorUpStrSpace       = _config.getDouble("stm.SScollimator.UpStrSpace");
    _SSCollimatorHalfWidth        = _config.getDouble("stm.SScollimator.halfWidth");
    _SSCollimatorHalfHeight       = _config.getDouble("stm.SScollimator.halfHeight");
    _SSCollimatorHalfLength       = _config.getDouble("stm.SScollimator.halfLength");
    _SSCollimatorLinerBuild       = _config.getBool(  "stm.SScollimator.liner.build");
    _SSCollimatorLinerMaterial    = _config.getString("stm.SScollimator.liner.material");
    _SSCollimatorLinerHalfWidth   = _config.getDouble("stm.SScollimator.liner.halfWidth");
    _SSCollimatorLinerHalfHeight  = _config.getDouble("stm.SScollimator.liner.halfHeight");
    _SSCollimatorLinerHalfLength  = _config.getDouble("stm.SScollimator.liner.halfLength");
    _SSCollimatorLinerCutOutHalfLength  = _config.getDouble("stm.SScollimator.liner.cutOutHalfLength");
    _SSCollimatorHole1xOffset     = _config.getDouble("stm.SScollimator.hole1.xoffset");
    _SSCollimatorHole1RadiusUpStr = _config.getDouble("stm.SScollimator.hole1.radiusUpStr");
    _SSCollimatorHole1RadiusDnStr = _config.getDouble("stm.SScollimator.hole1.radiusDnStr");
    _SSCollimatorHole1LinerBuild     = _config.getBool(  "stm.SScollimator.hole1.liner.build");
    _SSCollimatorHole1LinerThickness = _config.getDouble("stm.SScollimator.hole1.liner.thickness");
    _SSCollimatorHole2Build       = _config.getBool(  "stm.SScollimator.hole2.build");
    _SSCollimatorHole2xOffset     = _config.getDouble("stm.SScollimator.hole2.xoffset");
    _SSCollimatorHole2RadiusUpStr = _config.getDouble("stm.SScollimator.hole2.radiusUpStr");
    _SSCollimatorHole2RadiusDnStr = _config.getDouble("stm.SScollimator.hole2.radiusDnStr");
    _SSCollimatorHole2LinerBuild     = _config.getBool(  "stm.SScollimator.hole2.liner.build");
    _SSCollimatorHole2LinerThickness = _config.getDouble("stm.SScollimator.hole2.liner.thickness");
    _SSCollimatorHoleLinerMaterial = _config.getString("stm.SScollimator.hole.liner.material");

    _detectorTableBuild          = _config.getBool(  "stm.detector.stand.build",false);
    _detectorTableMaterial       = _config.getString("stm.detector.stand.material");
    _detectorTableTopExtraWidth  = _config.getDouble("stm.detector.stand.topExtraWidth");
    _detectorTableTopExtraLength = _config.getDouble("stm.detector.stand.topExtraLength");
    _detectorTableTopHalfHeight  = _config.getDouble("stm.detector.stand.topHalfHeight");
    _detectorTableLegRadius      = _config.getDouble("stm.detector.stand.legRadius");

    _detector1Build                    = _config.getBool(  "stm.det1.build",false);
    _detector1CrystalMaterial          = _config.getString("stm.det1.material");
    _detector1CrystalRadiusIn          = _config.getDouble("stm.det1.rIn");
    _detector1CrystalRadiusOut         = _config.getDouble("stm.det1.rOut");
    _detector1CrystalHalfLength        = _config.getDouble("stm.det1.halfLength");
    _detector1xOffset                  = _config.getDouble("stm.det1.xoffset");
    _detector1CanMaterial              = _config.getString("stm.det1.can.material");
    _detector1CanRadiusIn              = _config.getDouble("stm.det1.can.rIn");
    _detector1CanRadiusOut             = _config.getDouble("stm.det1.can.rOut");
    _detector1CanHalfLength            = _config.getDouble("stm.det1.can.halfLength");
    _detector1CanUpStrSpace            = _config.getDouble("stm.det1.can.UpStrSpace");
    _detector1CanUpStrWindowMaterial   = _config.getString("stm.det1.can.UpStrWindowMaterial");
    _detector1CanUpStrWindowHalfLength = _config.getDouble("stm.det1.can.UpStrWindowHalfLength");
    _detector1CanGasMaterial           = _config.getString("stm.det1.can.gas");

    _detector2Build                    = _config.getBool(  "stm.det2.build",false);
    _detector2CrystalMaterial          = _config.getString("stm.det2.material");
    _detector2CrystalRadiusIn          = _config.getDouble("stm.det2.rIn");
    _detector2CrystalRadiusOut         = _config.getDouble("stm.det2.rOut");
    _detector2CrystalHalfLength        = _config.getDouble("stm.det2.halfLength");
    _detector2xOffset                  = _config.getDouble("stm.det2.xoffset");
    _detector2CanMaterial              = _config.getString("stm.det2.can.material");
    _detector2CanRadiusIn              = _config.getDouble("stm.det2.can.rIn");
    _detector2CanRadiusOut             = _config.getDouble("stm.det2.can.rOut");
    _detector2CanHalfLength            = _config.getDouble("stm.det2.can.halfLength");
    _detector2CanUpStrSpace            = _config.getDouble("stm.det2.can.UpStrSpace");
    _detector2CanUpStrWindowMaterial   = _config.getString("stm.det2.can.UpStrWindowMaterial");
    _detector2CanUpStrWindowHalfLength = _config.getDouble("stm.det2.can.UpStrWindowHalfLength");
    _detector2CanGasMaterial           = _config.getString("stm.det2.can.gas");

    _shieldBuild                = _config.getBool(  "stm.shield.build",false);
    _shieldRadiusIn             = _config.getDouble("stm.shield.rIn");
    _shieldHasLiner             = _config.getBool(  "stm.shield.hasLiner", true /*true default for backwards compatibility*/);
    _shieldLinerWidth           = _config.getDouble("stm.shield.widthLiner");
    _shieldRadiusOut            = _config.getDouble("stm.shield.rOut");
    _shieldPipeHalfLength       = _config.getDouble("stm.shield.pipe.halfLength");
    _shieldMaterialLiner        = _config.getString("stm.shield.materialLiner");
    _shieldMaterial             = _config.getString("stm.shield.material");
    _shieldMatchPipeBlock       = _config.getBool  ("stm.shield.matchPipeBlock", false);
    _shieldUpStrSpace           = _config.getDouble("stm.shield.UpStrSpace");
    _shieldDnStrSpace           = _config.getDouble("stm.shield.DnStrSpace");
    _shieldDnStrWallHalfLength  = _config.getDouble("stm.shield.DnStrWall.halfLength");
    _shieldDnStrWallHoleRadius  = _config.getDouble("stm.shield.DnStrWall.holeRadius", -1.);
    _shieldDnStrWallHalfHeight  = _config.getDouble("stm.shield.DnStrWall.halfHeight", -1.);
    _shieldDnStrWallHalfWidth   = _config.getDouble("stm.shield.DnStrWall.halfWidth", -1.);
    _shieldDnStrWallGap         = _config.getDouble("stm.shield.DnStrWall.gap", 0.);
    _shieldUpStrWallGap         = _config.getDouble("stm.shield.UpStrWall.gap", 0.); //only if using pipe as origin
    _shieldDnStrWallMaterial    = _config.getString("stm.shield.DnStrWall.material", _shieldMaterial);
    _shieldBuildMatingBlock     = _config.getBool("stm.shield.matingBlock.build", true); // default to true because that is what older versions did
    _shieldPipeUpStrAirGap      = _config.getDouble("stm.shield.pipe.upStrAirGap", 0); // default to 0 for backwards compatibility

    _stmDnStrEnvBuild       = _config.getBool("stm.downstream.build");
    _stmDnStrEnvHalfLength  = _config.getDouble("stm.downstream.halfLength");
    _stmDnStrEnvHalfWidth   = _config.getDouble("stm.downstream.halfWidth");
    _stmDnStrEnvHalfHeight  = - _config.getDouble("yOfFloorSurface.below.mu2eOrigin");
    _stmDnStrEnvMaterial    = _config.getString("stm.downstream.material");

    _STM_SSCBuild          = _config.getBool(  "stm.STM_SSC.build");
    _STM_SSCVDBuild        = _config.getBool(  "stm.STM_SSC.VDbuild");
    _STM_SSCW_height       = _config.getDouble("stm.STM_SSC.W_height");
    _STM_SSCWdepth_f       = _config.getDouble("stm.STM_SSC.Wdepth_f");
    _STM_SSCWdepth_b       = _config.getDouble("stm.STM_SSC.Wdepth_b");

    // Each description reads only the keys it has: the earlier SSC is a
    // middle block plus two wings with the bores given as aperture areas,
    // the updated one a single block with the bores given as radii. Reading
    // them apart means a geometry file need not carry the other's parameters.
    if (_handstacked) {
      _STM_SSCW_width      = _config.getDouble("stm.STM_SSC.W_width");
      _STM_SSCr_LaBr_f     = _config.getDouble("stm.STM_SSC.r_LaBr_f");
      _STM_SSCr_HPGe_f     = _config.getDouble("stm.STM_SSC.r_HPGe_f");
      _STM_SSCr_LaBr_b     = _config.getDouble("stm.STM_SSC.r_LaBr_b");
      _STM_SSCr_HPGe_b     = _config.getDouble("stm.STM_SSC.r_HPGe_b");
      // The height of the bore above the baseplate: the datum every
      // hand-stacked piece downstream measures its own height from.
      _STM_SSCboreToBase   = _config.getDouble("stm.STM_SSC.boreToBase");
    } else {
      _STM_SSCdelta_WlR      = _config.getDouble("stm.STM_SSC.delta_WlR");
      _STM_SSCdelta_WlL      = _config.getDouble("stm.STM_SSC.delta_WlL");
      _STM_SSCW_middle       = _config.getDouble("stm.STM_SSC.W_middle");
      _STM_SSCAperture_HPGe1 = _config.getDouble("stm.STM_SSC.Aperture_HPGe1");
      _STM_SSCAperture_HPGe2 = _config.getDouble("stm.STM_SSC.Aperture_HPGe2");
      _STM_SSCAperture_LaBr1 = _config.getDouble("stm.STM_SSC.Aperture_LaBr1");
      _STM_SSCAperture_LaBr2 = _config.getDouble("stm.STM_SSC.Aperture_LaBr2");
    }
    _STM_SSCoffset_Spot    = _config.getDouble("stm.STM_SSC.offset_Spot");
    _STM_SSCleak           = _config.getDouble("stm.STM_SSC.leak");
    _STM_SSCFrontToWall    = _config.getDouble("stm.STM_SSC.FrontToWall");
    _STM_SSCZGap           = _config.getDouble("stm.STM_SSC.ZGap");
    _STM_SSCZGapBack       = _config.getDouble("stm.STM_SSC.ZGapBack");
    _STM_SSCMaterial       = _config.getString("stm.STM_SSC.material");

    _SSCSupportBuild          = _config.getBool(  "stm.SSCSupport.build");

    // The updated support is a five-plate steel cradle, not the earlier
    // table/legs/base/walls plus shims, so the two read different keys.
    if (_handstacked) {
      _SSCSupportdepth        = _config.getDouble("stm.SSCSupport.depth");
      _SSCSupportside_T       = _config.getDouble("stm.SSCSupport.side_T");
      _SSCSupportside_H       = _config.getDouble("stm.SSCSupport.side_H");
      _SSCSupportplate_base_T = _config.getDouble("stm.SSCSupport.base_T");
      _SSCSupportbottom_T     = _config.getDouble("stm.SSCSupport.bottom_T");
      _SSCSupporttop_T        = _config.getDouble("stm.SSCSupport.top_T");
      _SSCSupportMaterial     = _config.getString("stm.SSCSupport.material");
    } else {
      _SSCSupporttable_L      = _config.getDouble("stm.SSCSupport.table_L");
      _SSCSupporttable_H      = _config.getDouble("stm.SSCSupport.table_H");
      _SSCSupporttable_T      = _config.getDouble("stm.SSCSupport.table_T");
      _SSCSupportleg_L      = _config.getDouble("stm.SSCSupport.leg_L");
      _SSCSupportleg_H      = _config.getDouble("stm.SSCSupport.leg_H");
      _SSCSupportleg_T      = _config.getDouble("stm.SSCSupport.leg_T");
      _SSCSupportbase_L      = _config.getDouble("stm.SSCSupport.base_L");
      _SSCSupportbase_H      = _config.getDouble("stm.SSCSupport.base_H");
      _SSCSupportbase_T      = _config.getDouble("stm.SSCSupport.base_T");
      _SSCSupportwall_L      = _config.getDouble("stm.SSCSupport.wall_L");
      _SSCSupportwall_H      = _config.getDouble("stm.SSCSupport.wall_H");
      _SSCSupportwall_T      = _config.getDouble("stm.SSCSupport.wall_T");
      _SSCSupporthole_H      = _config.getDouble("stm.SSCSupport.hole_H");
      _SSCSupporthole_T      = _config.getDouble("stm.SSCSupport.hole_T");
      _SSCSupportFLeadStand_L      = _config.getDouble("stm.SSCSupport.FLeadStand_L");
      _SSCSupportFLeadStand_H      = _config.getDouble("stm.SSCSupport.FLeadStand_H");
      _SSCSupportFLeadStand_T      = _config.getDouble("stm.SSCSupport.FLeadStand_T");
      _SSCSupportFLeadShim_H       = _config.getDouble("stm.SSCSupport.FLeadShim_H");
      _SSCSupportFLeadShim_T       = _config.getDouble("stm.SSCSupport.FLeadShim_T");
      _SSCSupportFAluminumShim_T   = _config.getDouble("stm.SSCSupport.FAluminumShim_T");
      _SSCSupportFAluminumExtra_L  = _config.getDouble("stm.SSCSupport.FAluminumExtra_L");
      _SSCSupportFAluminumExtra_H  = _config.getDouble("stm.SSCSupport.FAluminumExtra_H");
    }

    // ---- SSC front shield ---------------------------------------------
    //
    // The geometry file writes every centre as if the structure sat
    // exactly on the SSC axis. Three offsets displace it, and they are
    // resolved here so that nothing downstream has to know which one
    // applies to which piece:
    //
    //   offsetX/offsetY        the whole structure
    //   LaBrSideBrickOffsetX   the bricks, except 2x4x8 centres 1, 2 and
    //                          4, which are constrained by poly2 instead
    //   poly2.offsetX          poly2, and those same three bricks
    if (_handstacked) {

      _SSCFrontShieldBuild = _config.getBool("stm.SSCFrontShield.build");

      const double fsOffX      = _config.getDouble("stm.SSCFrontShield.offsetX");
      const double fsOffY      = _config.getDouble("stm.SSCFrontShield.offsetY");
      const double fsLaBrOffX  = _config.getDouble("stm.SSCFrontShield.LaBrSideBrickOffsetX");
      const double fsPoly2OffX = _config.getDouble("stm.SSCFrontShield.poly2.offsetX");

      const CLHEP::Hep3Vector fsShift(fsOffX, fsOffY, 0.);

      // Read a numbered list of centres, shifting each by the structure
      // offset plus whichever x offset that piece follows.
      auto readCenters = [&](std::string const & base, int n,
                             std::vector<double> const & extraX) {
        std::vector<CLHEP::Hep3Vector> out;
        out.reserve(n);
        for (int i = 1; i <= n; ++i) {
          std::ostringstream key;
          key << base << ".center" << i;
          std::vector<double> c;
          _config.getVectorDouble(key.str(), c, 3);
          out.push_back(CLHEP::Hep3Vector(c[0] + extraX.at(i-1), c[1], c[2]) + fsShift);
        }
        return out;
      };
      // A structure may give one orientation for all of its bricks, as
      // the 2x4x16 wall does, or one per brick.
      auto readOrientations = [&](std::string const & base, int n) {
        std::vector<std::string> out;
        out.reserve(n);
        const std::string shared = _config.getString(base + ".orientation", "");
        for (int i = 1; i <= n; ++i) {
          if (!shared.empty()) { out.push_back(shared); continue; }
          std::ostringstream key;
          key << base << ".orientation" << i;
          out.push_back(_config.getString(key.str(), "000"));
        }
        return out;
      };

      const int n16 = _config.getInt("stm.SSCFrontShield.brick2x4x16.n");
      _SSCFrontShieldBrick2x4x16Center =
        readCenters("stm.SSCFrontShield.brick2x4x16", n16,
                    std::vector<double>(n16, fsLaBrOffX));
      _SSCFrontShieldBrick2x4x16Orientation =
        readOrientations("stm.SSCFrontShield.brick2x4x16", n16);

      // Bricks 1, 2 and 4 sit against poly2 rather than the LaBr side
      // stack, so they follow poly2's offset instead.
      const int n8 = _config.getInt("stm.SSCFrontShield.brick2x4x8.n");
      std::vector<double> extraX8(n8, fsLaBrOffX);
      for (int i : {1, 2, 4}) {
        if (i <= n8) extraX8[i-1] = fsPoly2OffX;
      }
      _SSCFrontShieldBrick2x4x8Center =
        readCenters("stm.SSCFrontShield.brick2x4x8", n8, extraX8);
      _SSCFrontShieldBrick2x4x8Orientation =
        readOrientations("stm.SSCFrontShield.brick2x4x8", n8);

      // The shelf and poly1 take only the whole-structure offset.
      _SSCFrontShieldShelfMaterial = _config.getString("stm.SSCFrontShield.shelf.material");
      _SSCFrontShieldShelfDim = CLHEP::Hep3Vector(_config.getDouble("stm.SSCFrontShield.shelf.dx"),
                                                  _config.getDouble("stm.SSCFrontShield.shelf.dy"),
                                                  _config.getDouble("stm.SSCFrontShield.shelf.dz"));
      {
        std::vector<double> c;
        _config.getVectorDouble("stm.SSCFrontShield.shelf.center", c, 3);
        _SSCFrontShieldShelfCenter = CLHEP::Hep3Vector(c[0], c[1], c[2]) + fsShift;
      }

      _SSCFrontShieldPoly1Material = _config.getString("stm.SSCFrontShield.poly1.material");
      _SSCFrontShieldPoly1Dim = CLHEP::Hep3Vector(_config.getDouble("stm.SSCFrontShield.poly1.dx"),
                                                  _config.getDouble("stm.SSCFrontShield.poly1.dy"),
                                                  _config.getDouble("stm.SSCFrontShield.poly1.dz"));
      {
        std::vector<double> c;
        _config.getVectorDouble("stm.SSCFrontShield.poly1.center", c, 3);
        _SSCFrontShieldPoly1Center = CLHEP::Hep3Vector(c[0], c[1], c[2]) + fsShift;
      }
      _SSCFrontShieldPoly1BoreR  = _config.getDouble("stm.SSCFrontShield.poly1.boreR");
      _SSCFrontShieldPoly1BoreDX = _config.getDouble("stm.SSCFrontShield.poly1.boreDX");
      _SSCFrontShieldPoly1BoreDY = _config.getDouble("stm.SSCFrontShield.poly1.boreDY");

      // poly2 carries its own x offset on top of the structure offset.
      _SSCFrontShieldPoly2Material = _config.getString("stm.SSCFrontShield.poly2.material");
      _SSCFrontShieldPoly2Dim = CLHEP::Hep3Vector(_config.getDouble("stm.SSCFrontShield.poly2.dx"),
                                                  _config.getDouble("stm.SSCFrontShield.poly2.dy"),
                                                  _config.getDouble("stm.SSCFrontShield.poly2.dz"));
      {
        std::vector<double> c;
        _config.getVectorDouble("stm.SSCFrontShield.poly2.center", c, 3);
        _SSCFrontShieldPoly2Center = CLHEP::Hep3Vector(c[0] + fsPoly2OffX, c[1], c[2]) + fsShift;
      }
    }


    if (_handstacked) {
      // ---- Front shielding, right ---------------------------------------
      //
      // The config gives a pattern, not positions: courses as lists of
      // brick types, a layer sequence in depth, and a map saying which
      // bores pass through which brick. Everything is expanded here, so
      // the construction code receives finished placements.
      //
      // Three constraints fix the section. brickEndX pins the +x end of
      // every course; the bore sits boreToBase above the baseplate with
      // courses on a 4 in pitch from it; and wallToCradleGap sets the
      // wall's front face behind the cradle. Pieces then butt together.

      _FrontShieldingRightBuild = _config.getBool("stm.FrontShieldingRight.build");

      const double fsrBrickEndX = _config.getDouble("stm.FrontShieldingRight.brickEndX");
      const double fsrGap       = _config.getDouble("stm.FrontShieldingRight.wallToCradleGap");
      const double fsrPitch     = _config.getDouble("stm.FrontShieldingRight.coursePitch");
      const std::string fsrBrickOrient =
        _config.getString("stm.FrontShieldingRight.brickOrientation");

      // The wall starts behind the cradle and grows downstream.
      const double fsrFrontZ = _SSCSupportdepth + fsrGap;

      // The bores follow the collimator rather than the wall, so they
      // are placed from offset_Spot with only a correction of their
      // own. Their y is the beam plane.
      std::vector<BrickWallBore> fsrBores;
      const int fsrBoreN = _config.getInt("stm.FrontShieldingRight.boreN");
      for (int i = 1; i <= fsrBoreN; ++i) {
        std::ostringstream base;
        base << "stm.FrontShieldingRight.bore" << i;
        BrickWallBore b;
        b.axis   = _config.getString(base.str() + "Axis");
        b.radius = _config.getDouble(base.str() + "R");
        b.offset = CLHEP::Hep3Vector(_config.getDouble(base.str() + "OffsetX"),
                                     _config.getDouble(base.str() + "OffsetY"), 0.);
        const double sign = (b.axis == "LaBr") ? +1. : -1.;
        b.center = CLHEP::Hep3Vector(sign*_STM_SSCoffset_Spot + b.offset.x(),
                                     b.offset.y(), 0.);
        fsrBores.push_back(b);
      }

      // Which bores go through which brick: {layer, course, position,
      // bore}, one entry per hole.
      struct FSRBoreMapEntry { int layer, course, position, bore; };
      std::vector<FSRBoreMapEntry> fsrBoreMap;
      const int fsrBoreMapN = _config.getInt("stm.FrontShieldingRight.leadBoreMapN");
      for (int i = 1; i <= fsrBoreMapN; ++i) {
        std::ostringstream key;
        key << "stm.FrontShieldingRight.leadBoreMap" << i;
        std::vector<int> e;
        _config.getVectorInt(key.str(), e, 4);
        fsrBoreMap.push_back({e[0], e[1], e[2], e[3]});
      }

      // Courses, numbered so that n = 0 is the one the beam passes
      // through.
      std::vector<int> fsrCourse;
      const int fsrCourseN = _config.getInt("stm.FrontShieldingRight.courseN");
      for (int i = 1; i <= fsrCourseN; ++i) {
        std::ostringstream key;
        key << "stm.FrontShieldingRight.course" << i;
        fsrCourse.push_back(_config.getInt(key.str()));
      }

      // Where the courses actually sit.
      //
      // They stack up from the baseplate -- block, then course after
      // course -- so their height comes from that chain, not from the
      // beam. The beam's own height is boreToBase above the same
      // baseplate, and the two coincide only if the chain happens to
      // leave half a course on the beam plane. That is true of this
      // section as built (236.000 - 83.600 - 101.600 = 50.800, half a
      // course) but it is arithmetic, not a constraint, so the offset
      // is derived here and checked against the bores below.
      const int fsrCoursesBelow =
        std::count_if(fsrCourse.begin(), fsrCourse.end(),
                      [](int n) { return n < 0; });
      const double fsrBlockDyForCourses =
        _config.getDouble("stm.FrontShieldingRight.block.dy");
      const double fsrCourseOffsetY =
        -_STM_SSCboreToBase + fsrBlockDyForCourses
        + fsrCoursesBelow*fsrPitch + fsrPitch/2;

      // How much of the course a brick of this type takes up. Rotation
      // is also taken into consideration
      const CLHEP::Hep3Vector fsrCourseDir(-1., 0., 0.);  // courses run toward -x
      const CLHEP::Hep3Vector fsrPitchDir  ( 0., 1., 0.);  // courses stack upward
      const CLHEP::Hep3Vector fsrDepthDir  ( 0., 0., 1.);  // layers face downstream
      CLHEP::HepRotation fsrBrickRot(CLHEP::HepRotation::IDENTITY);
      {
        OrientationResolver OR;
        OR.getRotationFromOrientation(fsrBrickRot, fsrBrickOrient);
      }
      auto fsrTypeSpan = [&](int type, CLHEP::Hep3Vector const & dir) {
        if (type < 1 || type > int(_leadBrickDims.size())) {
          throw cet::exception("GEOM")
            << "STMMaker: FrontShieldingRight names brick type " << type
            << ", but stm.leadBrick.typeN defines only "
            << _leadBrickDims.size() << ".\n";
        }
        return spanAlong(_leadBrickDims[type-1], fsrBrickRot, dir);
      };
      auto fsrTypeWidth = [&](int type) { return fsrTypeSpan(type, fsrCourseDir); };

      // Walk the layer sequence in depth. A lead layer expands into a
      // BrickWall; a sheet layer into plates held by the section.
      const int fsrLayerN = _config.getInt("stm.FrontShieldingRight.layerN");
      double fsrZ = fsrFrontZ;
      int fsrLeadLayer = 0;

      const std::string fsrSheetMat =
        _config.getString("stm.FrontShieldingRight.BPlayer.material");
      const double fsrSheetT  = _config.getDouble("stm.FrontShieldingRight.BPlayer.thickness");
      const double fsrSheetLo = _config.getDouble("stm.FrontShieldingRight.BPlayer.dyLower");
      const double fsrSheetUp = _config.getDouble("stm.FrontShieldingRight.BPlayer.dyUpper");

      // The courses span this, and the sheets share their footprint.
      double fsrCourseWidth = 0.;

      // The last lead layer's front face and thickness. The left half
      // shares this plane, so it is recorded as the loop runs rather
      // than worked back out of fsrZ afterwards, which by then has
      // walked past the copper lining as well.
      double fsrLastLeadZ = 0.;
      double fsrLastLeadT = 0.;

      for (int L = 1; L <= fsrLayerN; ++L) {
        std::ostringstream key;
        key << "stm.FrontShieldingRight.layer" << L;
        const std::string kind = _config.getString(key.str());

        if (kind == "Pb") {
          ++fsrLeadLayer;

          // How deep the layer is: the bricks' own thickness, rotated. Read
          // from the first brick of the first course, since a layer is
          // one brick deep and they all lie the same way up.
          double thick = 0.;

          std::vector<BrickWallBrick> bricks;
          for (size_t c = 0; c < fsrCourse.size(); ++c) {
            std::ostringstream ckey;
            ckey << "stm.FrontShieldingRight.leadLayer" << fsrLeadLayer
                 << "Course" << c+1;
            std::vector<int> types;
            _config.getVectorInt(ckey.str(), types);
            if (types.empty()) {
              throw cet::exception("GEOM")
                << "STMMaker: " << ckey.str() << " is empty.\n";
            }

            // The layer is one brick deep, so its thickness is that of
            // any of its bricks projected onto the depth direction.
            // Only assigned once, which is fine as all layers have the 
            // same thickness.
            if (thick == 0.) thick = fsrTypeSpan(types[0], fsrDepthDir);

            // Butt the types end to end, flush at brickEndX.
            double total = 0.;
            for (int t : types) total += fsrTypeWidth(t);
            fsrCourseWidth = total;
            // x start at the -x end, increase over the following iteration
            double x = fsrBrickEndX - total;
            for (size_t p = 0; p < types.size(); ++p) {
              const double w = fsrTypeWidth(types[p]);
              BrickWallBrick b;
              b.type        = types[p];
              b.orientation = fsrBrickOrient;
              // Final position: no offset is added later. Unlike
              // SSCFrontShield, this section carries no whole-structure
              // shift -- it is pinned by brickEndX, boreToBase and
              // wallToCradleGap, so it is moved by tuning one of those
              // rather than by displacing finished coordinates. The
              // only offsets here are per-bore, and they are already in
              // the bore centres above.
              b.center      = CLHEP::Hep3Vector(x + w/2,
                                                fsrCourseOffsetY + fsrCourse[c]*fsrPitch,
                                                fsrZ + thick/2);
              for (auto const & e : fsrBoreMap) {
                if (e.layer == fsrLeadLayer && e.course == int(c)+1
                    && e.position == int(p)+1) {
                  if (e.bore < 1 || e.bore > fsrBoreN) {
                    throw cet::exception("GEOM")
                      << "STMMaker: FrontShieldingRight bore map names bore "
                      << e.bore << ", but only " << fsrBoreN << " are defined.\n";
                  }
                  // The holes are pinned to the beam axes while
                  // brickEndX positions the courses, so the two can be
                  // moved apart. Check the hole really lands on the
                  // brick the map claims: otherwise the subtraction
                  // would clip nothing and the beam would see solid
                  // lead, with nothing to say so.
                  const BrickWallBore & bore = fsrBores[e.bore-1];
                  const double xlo = x, xhi = x + w;
                  const double h    = fsrTypeSpan(types[p], fsrPitchDir);
                  const double ylo  = b.center.y() - h/2;
                  const double yhi  = b.center.y() + h/2;
                  // Both axes matter. In x the courses are pinned by
                  // brickEndX while the holes follow the beam, and in y
                  // the courses stack from the baseplate while the beam
                  // sits boreToBase above it -- so either can drift
                  // away from the other.
                  if (bore.center.x() - bore.radius < xlo ||
                      bore.center.x() + bore.radius > xhi ||
                      bore.center.y() - bore.radius < ylo ||
                      bore.center.y() + bore.radius > yhi) {
                    throw cet::exception("GEOM")
                      << "STMMaker: FrontShieldingRight bore " << e.bore
                      << " (" << bore.axis << ") does not fit the brick named by"
                      << " leadBoreMap {layer " << e.layer << ", course "
                      << e.course << ", position " << e.position << "}.\n"
                      << "The hole spans x = [" << bore.center.x() - bore.radius
                      << ", " << bore.center.x() + bore.radius
                      << "], y = [" << bore.center.y() - bore.radius
                      << ", " << bore.center.y() + bore.radius << "] mm;\n"
                      << "that brick spans x = [" << xlo << ", " << xhi
                      << "], y = [" << ylo << ", " << yhi << "] mm.\n"
                      << "Check stm.FrontShieldingRight.brickEndX ("
                      << fsrBrickEndX << " mm), the course stack "
                      << "(boreToBase " << _STM_SSCboreToBase
                      << ", block " << fsrBlockDyForCourses
                      << ", pitch " << fsrPitch << " mm), and this bore's own "
                      << "offset (" << bore.offset.x() << ", " << bore.offset.y()
                      << " mm) against the bore map.\n";
                  }
                  b.bores.push_back(e.bore);
                }
              }
              bricks.push_back(b);
              x += w;
            }
          }

          _FrontShieldingRightLeadLayers.push_back(
            BrickWall(_FrontShieldingRightBuild,
                      CLHEP::Hep3Vector(fsrBrickEndX, 0., fsrZ),
                      fsrCourseDir, fsrPitchDir, fsrDepthDir,
                      bricks, fsrBores));
          fsrLastLeadZ = fsrZ;
          fsrLastLeadT = thick;
          fsrZ += thick;

        } else {
          // A sheet layer: two plates split at the beam, sharing the
          // courses' footprint. The lower one carries the bores.
          const double top = (fsrCourse.back() + 1)*fsrPitch - fsrPitch/2;
          const double bot = top - (fsrSheetLo + fsrSheetUp);

          FrontShieldingRightSheet lower;
          lower.material    = fsrSheetMat;
          lower.orientation = "000";
          lower.halfDim     = CLHEP::Hep3Vector(fsrCourseWidth/2, fsrSheetLo/2, fsrSheetT/2);
          lower.center      = CLHEP::Hep3Vector(fsrBrickEndX - fsrCourseWidth/2,
                                                bot + fsrSheetLo/2,
                                                fsrZ + fsrSheetT/2);
          for (int i = 0; i < fsrBoreN; ++i) {
            lower.bores.push_back(i+1);
            lower.boreRadius.push_back(fsrBores[i].radius);
          }
          _FrontShieldingRightSheets.push_back(lower);

          FrontShieldingRightSheet upper;
          upper.material    = fsrSheetMat;
          upper.orientation = "000";
          upper.halfDim     = CLHEP::Hep3Vector(fsrCourseWidth/2, fsrSheetUp/2, fsrSheetT/2);
          upper.center      = CLHEP::Hep3Vector(fsrBrickEndX - fsrCourseWidth/2,
                                                bot + fsrSheetLo + fsrSheetUp/2,
                                                fsrZ + fsrSheetT/2);
          _FrontShieldingRightSheets.push_back(upper);

          fsrZ += fsrSheetT;
        }
      }

      // The blocks the wall stands on: they span the courses' width and
      // their undersides land on the baseplate.
      _FrontShieldingRightBlockMaterial =
        _config.getString("stm.FrontShieldingRight.block.material");
      const double fsrBlockDx = _config.getDouble("stm.FrontShieldingRight.block.dx");
      const double fsrBlockDy = _config.getDouble("stm.FrontShieldingRight.block.dy");
      const double fsrBlockDz = _config.getDouble("stm.FrontShieldingRight.block.dz");
      _FrontShieldingRightBlockHalfDim =
        CLHEP::Hep3Vector(fsrBlockDx/2, fsrBlockDy/2, fsrBlockDz/2);

      const int fsrBlockN = _config.getInt("stm.FrontShieldingRight.block.n");
      const double fsrBlockY = -_STM_SSCboreToBase + fsrBlockDy/2;
      for (int i = 0; i < fsrBlockN; ++i) {
        _FrontShieldingRightBlockCenter.push_back(
          CLHEP::Hep3Vector(fsrBrickEndX - fsrBlockDx/2 - i*fsrBlockDx,
                            fsrBlockY,
                            fsrFrontZ + fsrBlockDz/2));
      }

      // The two pipes, each on its bore, starting at the wall front.
      const std::string fsrPipeMat = _config.getString("stm.FrontShieldingRight.pipe.material");
      const double fsrPipeRInLaBr  = _config.getDouble("stm.FrontShieldingRight.pipe.rInLaBr");
      const double fsrPipeRInHPGe  = _config.getDouble("stm.FrontShieldingRight.pipe.rInHPGe");
      const double fsrPipeROutLaBr = _config.getDouble("stm.FrontShieldingRight.pipe.rOutLaBr");
      const double fsrPipeROutHPGe = _config.getDouble("stm.FrontShieldingRight.pipe.rOutHPGe");
      const double fsrPipeLenLaBr  = _config.getDouble("stm.FrontShieldingRight.pipe.lenLaBr");
      const double fsrPipeLenHPGe  = _config.getDouble("stm.FrontShieldingRight.pipe.lenHPGe");
      for (int i = 0; i < fsrBoreN; ++i) {
        const bool laBr = (fsrBores[i].axis == "LaBr");
        FrontShieldingRightPipe p;
        p.material   = fsrPipeMat;
        p.rIn        = laBr ? fsrPipeRInLaBr  : fsrPipeRInHPGe;
        p.rOut       = laBr ? fsrPipeROutLaBr : fsrPipeROutHPGe;
        p.halfLength = (laBr ? fsrPipeLenLaBr : fsrPipeLenHPGe)/2;
        p.center     = CLHEP::Hep3Vector(fsrBores[i].center.x(),
                                         fsrBores[i].center.y(),
                                         fsrFrontZ + p.halfLength);
        _FrontShieldingRightPipes.push_back(p);
      }

      // The copper plate, placed by the anchor corner of its outline
      // rather than by a centre.
      _FrontShieldingRightPlate.material =
        _config.getString("stm.FrontShieldingRight.copperLining.material");
      _config.getVectorDouble("stm.FrontShieldingRight.copperLining.UVerts",
                              _FrontShieldingRightPlate.uVerts);
      _config.getVectorDouble("stm.FrontShieldingRight.copperLining.VVerts",
                              _FrontShieldingRightPlate.vVerts);
      _FrontShieldingRightPlate.length =
        _config.getDouble("stm.FrontShieldingRight.copperLining.length");
      _FrontShieldingRightPlate.orientation =
        _config.getString("stm.FrontShieldingRight.copperLining.orientation");
      const double fsrCuLiningFromEnd =
        _config.getDouble("stm.FrontShieldingRight.copperLining.originFromBrickEnd");
      // The lining is not centred on the beam: its underside sits a
      // stated height above the baseplate, on the same chain that
      // fixes the bore height. The sweep runs along y and an extruded
      // solid is centred on its placement point, so the anchor's y is
      // the mid-plane -- bottom plus half the sweep -- not the edge.
      const double fsrCuLiningBottomToBase =
        _config.getDouble("stm.FrontShieldingRight.copperLining.bottomToBase");
      const double fsrCuLiningY =
        -_STM_SSCboreToBase + fsrCuLiningBottomToBase
        + _FrontShieldingRightPlate.length/2;
      _FrontShieldingRightPlate.anchor =
        CLHEP::Hep3Vector(fsrBrickEndX - fsrCuLiningFromEnd, fsrCuLiningY, fsrZ);
      // The plate's holes are the same bores as everything else on the
      // beamline, so it names them by id: whatever offset has moved a
      // bore off the collimator axis has already been applied there.
      const int fsrCuLiningHoles = _config.getInt("stm.FrontShieldingRight.copperLining.nHoles");
      if (fsrCuLiningHoles > fsrBoreN) {
        throw cet::exception("GEOM")
          << "STMMaker: stm.FrontShieldingRight.copperLining.nHoles is "
          << fsrCuLiningHoles << " but only " << fsrBoreN << " bores are defined.\n";
      }
      for (int i = 1; i <= fsrCuLiningHoles; ++i) {
        std::ostringstream key;
        key << "stm.FrontShieldingRight.copperLining.holeRadius" << i;
        _FrontShieldingRightPlate.bores.push_back(i);
        _FrontShieldingRightPlate.holeRadius.push_back(_config.getDouble(key.str()));
      }

      // The section's downstream face: the lining's front, plus however
      // much of the lining lies along the depth direction.
      {
        CLHEP::HepRotation liningRot(CLHEP::HepRotation::IDENTITY);
        OrientationResolver OR;
        OR.getRotationFromOrientation(liningRot, _FrontShieldingRightPlate.orientation);

        auto range = [](std::vector<double> const & v) {
          return v.empty() ? 0. : *std::max_element(v.begin(), v.end())
                                - *std::min_element(v.begin(), v.end());
        };
        // In the outline's own frame u and v span the polygon and the
        // sweep runs along its z.
        const CLHEP::Hep3Vector liningDim(range(_FrontShieldingRightPlate.uVerts),
                                          range(_FrontShieldingRightPlate.vVerts),
                                          _FrontShieldingRightPlate.length);
        _FrontShieldingRightBackZ =
          fsrZ + spanAlong(liningDim, liningRot, fsrDepthDir);
      }

      // ---- Front shielding, left ----------------------------------
      //
      // The half beside the collimator. Nothing here is bored and no
      // poly interleaves with a lead layer, so it is three brick
      // groups, four sheets and a prism, each placed against the
      // right half rather than measured on its own.
      //
      // Every x is a centre measured back from brickEndX and so is
      // negative; every y is a height above the baseplate. The
      // per-bore offsets do not enter: they move the beam, and these
      // pieces are anchored by brickEndX and the baseplate.

      _FrontShieldingLeftBuild = _config.getBool("stm.FrontShieldingLeft.build");

      const std::string fslPolyMat =
        _config.getString("stm.FrontShieldingLeft.BPmaterial");

      auto fslSheet = [&](std::string const & name,
                          double dx, double dy, double dz,
                          double z, double bottomToBase) {
        FrontShieldingLeftSheet s;
        s.material = fslPolyMat;
        s.halfDim  = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
        s.center   = CLHEP::Hep3Vector(
          fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft." + name + ".fromEndX"),
          -_STM_SSCboreToBase + bottomToBase + dy/2,
          z + dz/2);
        return s;
      };

      // The side slab stands on the baseplate, spanning the depth of
      // the right half's first four layers from its front face.
      const double fslSideDx = _config.getDouble("stm.FrontShieldingLeft.sidePoly.dx");
      const double fslSideDy = _config.getDouble("stm.FrontShieldingLeft.sidePoly.dy");
      const double fslSideDz = _config.getDouble("stm.FrontShieldingLeft.sidePoly.dz");
      _FrontShieldingLeftSheets.push_back(
        fslSheet("sidePoly", fslSideDx, fslSideDy, fslSideDz, fsrFrontZ, 0.));

      // The outer sheet is the layer in front of the brick block, and
      // the block shares the plane of the right half's last lead
      // layer, so the sheet sits one thickness ahead of that plane.
      const double fslOuterDx = _config.getDouble("stm.FrontShieldingLeft.outerPoly.dx");
      const double fslOuterDy = _config.getDouble("stm.FrontShieldingLeft.outerPoly.dy");
      const double fslOuterDz = _config.getDouble("stm.FrontShieldingLeft.outerPoly.dz");
      _FrontShieldingLeftSheets.push_back(
        fslSheet("outerPoly", fslOuterDx, fslOuterDy, fslOuterDz,
                 fsrLastLeadZ - fslOuterDz, 0.));

      // The inner sheet and the edge sheet share a depth -- the layer
      // behind the last lead layer -- but not a height.
      const double fslInnerDx = _config.getDouble("stm.FrontShieldingLeft.innerPoly.dx");
      const double fslInnerDy = _config.getDouble("stm.FrontShieldingLeft.innerPoly.dy");
      const double fslInnerDz = _config.getDouble("stm.FrontShieldingLeft.innerPoly.dz");
      const double fslInnerBase =
        _config.getDouble("stm.FrontShieldingLeft.innerPoly.bottomToBase");
      const double fslInnerZ = fsrLastLeadZ + fsrLastLeadT;
      _FrontShieldingLeftSheets.push_back(
        fslSheet("innerPoly", fslInnerDx, fslInnerDy, fslInnerDz,
                 fslInnerZ, fslInnerBase));

      const double fslEdgeDx = _config.getDouble("stm.FrontShieldingLeft.edgePoly.dx");
      const double fslEdgeDy = _config.getDouble("stm.FrontShieldingLeft.edgePoly.dy");
      const double fslEdgeDz = _config.getDouble("stm.FrontShieldingLeft.edgePoly.dz");
      const double fslEdgeBase =
        _config.getDouble("stm.FrontShieldingLeft.edgePoly.bottomToBase");
      _FrontShieldingLeftSheets.push_back(
        fslSheet("edgePoly", fslEdgeDx, fslEdgeDy, fslEdgeDz,
                 fslInnerZ, fslEdgeBase));

      // A brick group: one BrickWall, its courses butted end to end
      // from a corner. The same expansion serves all three groups
      // here, which differ only in where they start and which way
      // their courses run. Orientation is per course, since these
      // groups mix bricks lying flat with bricks stood on end.
      auto fslBrickGroup = [&](std::vector<std::vector<int> > const & courses,
                               std::vector<std::string> const & orientations,
                               CLHEP::Hep3Vector const & origin,
                               CLHEP::Hep3Vector const & courseDir,
                               CLHEP::Hep3Vector const & pitchDir,
                               CLHEP::Hep3Vector const & depthDir,
                               std::string const & what) {
        std::vector<BrickWallBrick> bricks;
        double pitchAt = 0.;
        for (size_t c = 0; c < courses.size(); ++c) {
          if (courses[c].empty()) {
            throw cet::exception("GEOM")
              << "STMMaker: FrontShieldingLeft " << what << " course "
              << c+1 << " is empty.\n";
          }
          CLHEP::HepRotation rot(CLHEP::HepRotation::IDENTITY);
          {
            OrientationResolver OR;
            OR.getRotationFromOrientation(rot, orientations[c]);
          }
          auto span = [&](int type, CLHEP::Hep3Vector const & dir) {
            if (type < 1 || type > int(_leadBrickDims.size())) {
              throw cet::exception("GEOM")
                << "STMMaker: FrontShieldingLeft " << what
                << " names brick type " << type
                << ", but stm.leadBrick.typeN defines only "
                << _leadBrickDims.size() << ".\n";
            }
            return spanAlong(_leadBrickDims[type-1], rot, dir);
          };

          // How far this course steps along the pitch. A course is one
          // brick deep, so any of its bricks gives the step -- read
          // per course rather than once, because these groups mix
          // brick sizes from course to course.
          const double step = span(courses[c][0], pitchDir);

          double along = 0.;
          for (int t : courses[c]) {
            const double w = span(t, courseDir);
            BrickWallBrick b;
            b.type        = t;
            b.orientation = orientations[c];
            b.center      = origin
                          + courseDir*(along + w/2)
                          + pitchDir *(pitchAt + step/2)
                          + depthDir *(span(t, depthDir)/2);
            bricks.push_back(b);
            along += w;
          }
          pitchAt += step;
        }
        return BrickWall(_FrontShieldingLeftBuild, origin,
                         courseDir, pitchDir, depthDir,
                         bricks, std::vector<BrickWallBore>());
      };

      // The nine-brick block. Its columns run up y, so a "course"
      // here is a column: the course direction is +y and the pitch
      // steps along -x, away from the right half. It stands on the
      // baseplate.
      //
      // Its +x face is flush with the side slab's +x face, which is the
      // right half's -x course edge, but going through the slab is
      // what keeps them together if the slab is ever moved or
      // retuned. The slab stands alongside the block in z, spanning
      // the right half's first four layers.
      const double fslBlockEndX =
        fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft.sidePoly.fromEndX")
        + fslSideDx/2;
      {
        std::vector<std::vector<int> > columns;
        std::vector<std::string>       orientations;
        const int n = _config.getInt("stm.FrontShieldingLeft.blockColumnN");
        for (int i = 1; i <= n; ++i) {
          std::ostringstream key, okey;
          key  << "stm.FrontShieldingLeft.blockColumn" << i;
          okey << "stm.FrontShieldingLeft.blockColumn" << i << "Orientation";
          std::vector<int> types;
          _config.getVectorInt(key.str(), types);
          columns.push_back(types);
          orientations.push_back(_config.getString(okey.str()));
        }
        _FrontShieldingLeftBrickGroups.push_back(
          fslBrickGroup(columns, orientations,
                        CLHEP::Hep3Vector(fslBlockEndX, -_STM_SSCboreToBase,
                                          fsrLastLeadZ),
                        CLHEP::Hep3Vector(0., 1., 0.),   // a column runs up
                        CLHEP::Hep3Vector(-1., 0., 0.),  // columns step -x
                        CLHEP::Hep3Vector(0., 0., 1.),
                        "outer block"));
      }

      // The two bricks beside the edge sheet. They are level with it
      // rather than stacked on it, so they take its height, and they
      // sit one layer along +z from it.
      {
        std::vector<int> types;
        _config.getVectorInt("stm.FrontShieldingLeft.edgeCourse", types);
        const std::string orient =
          _config.getString("stm.FrontShieldingLeft.edgeCourseOrientation");
        // Flush with the sheet's +x end, running -x, the same way the
        // block and the right half's courses run.
        const double xStart =
          fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft.edgePoly.fromEndX")
          + fslEdgeDx/2;
        _FrontShieldingLeftBrickGroups.push_back(
          fslBrickGroup(std::vector<std::vector<int> >(1, types),
                        std::vector<std::string>(1, orient),
                        CLHEP::Hep3Vector(xStart,
                                          -_STM_SSCboreToBase + fslEdgeBase,
                                          fslInnerZ + fslEdgeDz),
                        CLHEP::Hep3Vector(-1., 0., 0.),
                        CLHEP::Hep3Vector(0., 1., 0.),
                        CLHEP::Hep3Vector(0., 0., 1.),   // +z of the sheet
                        "edge course"));
      }

      // The grid behind the inner sheet: courses running -x from the
      // sheet's +x end, stacking up from the baseplate.
      {
        std::vector<std::vector<int> > courses;
        std::vector<std::string>       orientations;
        const int n = _config.getInt("stm.FrontShieldingLeft.innerGridCourseN");
        const std::string orient =
          _config.getString("stm.FrontShieldingLeft.innerGridOrientation");
        for (int i = 1; i <= n; ++i) {
          std::ostringstream key;
          key << "stm.FrontShieldingLeft.innerGridCourse" << i;
          std::vector<int> types;
          _config.getVectorInt(key.str(), types);
          courses.push_back(types);
          orientations.push_back(orient);
        }
        const double xStart =
          fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft.innerPoly.fromEndX")
          + fslInnerDx/2;
        _FrontShieldingLeftBrickGroups.push_back(
          fslBrickGroup(courses, orientations,
                        CLHEP::Hep3Vector(xStart,
                                          -_STM_SSCboreToBase + fslInnerBase,
                                          fslInnerZ + fslInnerDz),
                        CLHEP::Hep3Vector(-1., 0., 0.),
                        CLHEP::Hep3Vector(0., 1., 0.),
                        CLHEP::Hep3Vector(0., 0., 1.),   // +z of the sheet
                        "inner grid"));
      }

      // The prism, placed by its right angle rather than by a centre.
      _FrontShieldingLeftPrism.material =
        _config.getString("stm.FrontShieldingLeft.triangle.material");
      _config.getVectorDouble("stm.FrontShieldingLeft.triangle.UVerts",
                              _FrontShieldingLeftPrism.uVerts);
      _config.getVectorDouble("stm.FrontShieldingLeft.triangle.VVerts",
                              _FrontShieldingLeftPrism.vVerts);
      _FrontShieldingLeftPrism.length =
        _config.getDouble("stm.FrontShieldingLeft.triangle.length");
      _FrontShieldingLeftPrism.orientation =
        _config.getString("stm.FrontShieldingLeft.triangle.orientation");
      // The right angle butts the inner sheet's +x and -z faces, and
      // the prism's top is level with that sheet's top. An extruded
      // solid is centred on its placement point, so the anchor's y is
      // half a sweep below that top rather than at it.
      {
        const double innerX =
          fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft.innerPoly.fromEndX");
        const double innerTop =
          -_STM_SSCboreToBase + fslInnerBase + fslInnerDy;
        _FrontShieldingLeftPrism.anchor =
          CLHEP::Hep3Vector(innerX + fslInnerDx/2,
                            innerTop - _FrontShieldingLeftPrism.length/2,
                            fslInnerZ);
      }
    }
    else{
      _FrontShieldingBuild   = _config.getBool(  "stm.FrontShielding.build");
      _FrontSHeightofRoom    = _config.getDouble("stm.FrontShielding.HeightofRoom");
      _FrontStungstenlength  = _config.getDouble("stm.FrontShielding.tungstenlength");
      _FrontStungstendepth   = _config.getDouble("stm.FrontShielding.tungstendepth");
      _FrontSleaddepth1      = _config.getDouble("stm.FrontShielding.leaddepth1");
      _FrontSleaddepth2      = _config.getDouble("stm.FrontShielding.leaddepth2");
      _FrontSaluminumdepth   = _config.getDouble("stm.FrontShielding.aluminumdepth");
      _FrontScopperdepth     = _config.getDouble("stm.FrontShielding.copperdepth");
      _FrontSBPdepth         = _config.getDouble("stm.FrontShielding.BPdepth");
      _FrontSfPb_lengthL     = _config.getDouble("stm.FrontShielding.fPb_lengthL");
      _FrontSfPb_lengthR     = _config.getDouble("stm.FrontShielding.fPb_lengthR");
      _FrontSGapForTop       = _config.getDouble("stm.FrontShielding.GapForTop");
      _FrontSLeakForSSC      = _config.getDouble("stm.FrontShielding.LeakForSSC");
      _FrontSCopperL         = _config.getDouble("stm.FrontShielding.CopperL");
      _FrontS_H              = _config.getDouble("stm.FrontShielding.FrontS_H");
      _FrontSHole_r          = _config.getDouble("stm.FrontShielding.FrontSHole_r");
    }


    _HPGeBuild                = _config.getBool("stm.HPGe.build");
    _HPGecrystalMaterial      = _config.getString("stm.HPGe.crystalMaterial");
    _HPGeholeMaterial         = _config.getString("stm.HPGe.holeMaterial");
    _HPGewindowMaterial       = _config.getString("stm.HPGe.windowMaterial");
    _HPGewallMaterial         = _config.getString("stm.HPGe.wallMaterial");
    _HPGecapsuleMaterial      = _config.getString("stm.HPGe.capsuleMaterial");
    _HPGeEndcapR              = _config.getDouble("stm.HPGe.EndcapR");
    _HPGeEndcapL              = _config.getDouble("stm.HPGe.EndcapL");
    _HPGeCrystalR             = _config.getDouble("stm.HPGe.CrystalR");
    _HPGeCrystalL             = _config.getDouble("stm.HPGe.CrystalL");
    _HPGeZ_HPGe               = _config.getDouble("stm.HPGe.Z_HPGe");
    _HPGeHoleR                = _config.getDouble("stm.HPGe.HoleR");
    _HPGeHoleL                = _config.getDouble("stm.HPGe.HoleL");
    _HPGeCapsule_Wallthick    = _config.getDouble("stm.HPGe.Capsule_Wallthick");
    _HPGeCapsule_Windowthick  = _config.getDouble("stm.HPGe.Capsule_Windowthick");
    _HPGeCapsule_Endthick     = _config.getDouble("stm.HPGe.Capsule_Endthick");
    _HPGeCapsule_Walllength   = _config.getDouble("stm.HPGe.Capsule_Walllength");
    _HPGeWindowD              = _config.getDouble("stm.HPGe.WindowD");
    _HPGeEndcapD              = _config.getDouble("stm.HPGe.EndcapD");
    _HPGeAirD                 = _config.getDouble("stm.HPGe.AirD");
    _HPGeoffset_HPGe          = _config.getDouble("stm.HPGe.offset_HPGe");

    _LaBrBuild                = _config.getBool("stm.LaBr.build");
    _LaBrcrystalMaterial      = _config.getString("stm.LaBr.crystalMaterial");
    _LaBrwindowMaterial       = _config.getString("stm.LaBr.windowMaterial");
    _LaBrwallMaterial         = _config.getString("stm.LaBr.wallMaterial");
    _LaBrEndcapR              = _config.getDouble("stm.LaBr.EndcapR");
    _LaBrEndcapL              = _config.getDouble("stm.LaBr.EndcapL");
    _LaBrCrystalR             = _config.getDouble("stm.LaBr.CrystalR");
    _LaBrCrystalL             = _config.getDouble("stm.LaBr.CrystalL");
    _LaBrZ_LaBr               = _config.getDouble("stm.LaBr.Z_LaBr");
    _LaBrWindowD              = _config.getDouble("stm.LaBr.WindowD");
    _LaBrEndcapD              = _config.getDouble("stm.LaBr.EndcapD");
    _LaBrAirD                 = _config.getDouble("stm.LaBr.AirD");
    _LaBroffset_LaBr          = _config.getDouble("stm.LaBr.offset_LaBr");


    _BottomShieldingBuild  = _config.getBool("stm.BottomShielding.build");
    _BottomSfloor_Zlength  = _config.getDouble("stm.BottomShielding.floor_Zlength");
    _BottomSFront_LB       = _config.getDouble("stm.BottomShielding.Front_LB");
    _BottomSFront_LB_inner = _config.getDouble("stm.BottomShielding.Front_LB_inner");
    _BottomSleaddepth      = _config.getDouble("stm.BottomShielding.leaddepth");
    _BottomScopperdepth    = _config.getDouble("stm.BottomShielding.copperdepth");
    _BottomSBPdepth        = _config.getDouble("stm.BottomShielding.BPdepth");


    _LeftShieldingBuild   = _config.getBool("stm.LeftShielding.build");
    _LeftS_Length         = _config.getDouble("stm.LeftShielding.Length");
    _LeftSleaddepth       = _config.getDouble("stm.LeftShielding.leaddepth");
    _LeftScopperdepth     = _config.getDouble("stm.LeftShielding.copperdepth");
    _LeftSBPdepth         = _config.getDouble("stm.LeftShielding.BPdepth");
    _LeftSXmin            = _config.getDouble("stm.LeftShielding.Left_Xmin");

    _RightShieldingBuild  = _config.getBool("stm.RightShielding.build");
    _RightS_Length        = _config.getDouble("stm.RightShielding.Length");
    _RightSleaddepth      = _config.getDouble("stm.RightShielding.leaddepth");
    _RightScopperdepth    = _config.getDouble("stm.RightShielding.copperdepth");
    _RightSBPdepth        = _config.getDouble("stm.RightShielding.BPdepth");
    _RightSXmax           = _config.getDouble("stm.RightShielding.Right_Xmax");

    _TopShieldingBuild      = _config.getBool("stm.TopShielding.build");
    _TopShieldingSkirtBuild = _config.getBool("stm.TopShielding.Skirtbuild");
    _TopLiftBeam_L        = _config.getDouble("stm.TopShielding.LiftBeam_L");
    _TopLiftBeam_H        = _config.getDouble("stm.TopShielding.LiftBeam_H");
    _TopLiftBeam_T        = _config.getDouble("stm.TopShielding.LiftBeam_T");
    _TopLiftBeam_Xmove    = _config.getDouble("stm.TopShielding.LiftBeam_Xmove");
    _TopSZlength          = _config.getDouble("stm.TopShielding.Zlength");
    _TopSXlength          = _config.getDouble("stm.TopShielding.Xlength");
    _TopSFront_LT         = _config.getDouble("stm.TopShielding.Front_LT");
    _TopTFZlength         = _config.getDouble("stm.TopShielding.TFZlength");
    _TopTFXlength         = _config.getDouble("stm.TopShielding.TFXlength");
    _TopTBZlength         = _config.getDouble("stm.TopShielding.TBZlength");
    _TopScontainerdepth   = _config.getDouble("stm.TopShielding.containerdepth");
    _TopSleaddepth        = _config.getDouble("stm.TopShielding.leaddepth");
    _TopScopperdepth      = _config.getDouble("stm.TopShielding.copperdepth");
    _TopSBPdepth          = _config.getDouble("stm.TopShielding.BPdepth");
    _TopSZHole            = _config.getDouble("stm.TopShielding.Zhole");
    _TopSBarLeft          = _config.getDouble("stm.TopShielding.BarLeft");
    _TopSBarRight         = _config.getDouble("stm.TopShielding.BarRight");
    _TopSGapLeft          = _config.getDouble("stm.TopShielding.GapLeft");
    _TopSGapRight         = _config.getDouble("stm.TopShielding.GapRight");
    _TopSLeak             = _config.getDouble("stm.TopShielding.Leak");

    _BackShieldingBuild   = _config.getBool("stm.BackShielding.build");
    _BackSBPThick         = _config.getDouble("stm.BackShielding.BPThick");
    _BackSBPLength        = _config.getDouble("stm.BackShielding.BPLength");
    _BackSBPHeight        = _config.getDouble("stm.BackShielding.BPHeight");
    _BackS_dX             = _config.getDouble("stm.BackShielding.BackS_dX");
    _BackS_dY             = _config.getDouble("stm.BackShielding.BackS_dY");
    _BackSPipeGap         = _config.getDouble("stm.BackShielding.ShieldingPipeGap");


    _InnerShieldingBuild        = _config.getBool("stm.InnerShielding.build");



    _ElectronicShieldingBuild      = _config.getBool("stm.ElectronicShielding.build");
    _ElectronicSSiGridX    = _config.getDouble("stm.ElectronicShielding.SiGridX");
    _ElectronicSSiGridY    = _config.getDouble("stm.ElectronicShielding.SiGridY");
    _ElectronicSSiGridZ    = _config.getDouble("stm.ElectronicShielding.SiGridZ");
    _ElectronicSSiXcenter  = _config.getDouble("stm.ElectronicShielding.SiXcenter");
    _ElectronicSSiYcenter  = _config.getDouble("stm.ElectronicShielding.SiYcenter");
    _ElectronicSSiZcenter  = _config.getDouble("stm.ElectronicShielding.SiZcenter");
    _ElectronicSConcreteT  = _config.getDouble("stm.ElectronicShielding.ConcreteT");
    _ElectronicSGapToSi    = _config.getDouble("stm.ElectronicShielding.GapToSi");


    // The SSC front shield replaces the absorber in the updated geometry,
    // so only one description's keys are read.
    if (!_handstacked) {
      _STM_AbsorberBuild          = _config.getBool("stm.STM_Absorber.build");
      _STM_Absorber_hW            = _config.getDouble("stm.STM_Absorber.hW");
      _STM_Absorber_hH            = _config.getDouble("stm.STM_Absorber.hH");
      _STM_Absorber_hT            = _config.getDouble("stm.STM_Absorber.hT");
      _STM_Absorber_GaptoSSC      = _config.getDouble("stm.STM_Absorber.GaptoSSC");
    }

  }
} // namespace mu2e
