//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file optical/LXe/include/LXeMainVolume.hh
/// \brief Definition of the LXeMainVolume class
//
#ifndef LXeAccelerator_h
#define LXeAccelerator_h 1

#include "LXeDetectorConstruction.hh"
#include "globals.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "ElectricFieldSetup.hh"
#include "F03FieldSetup.hh"
#include "G4Cache.hh"

class G4Box;
class G4LogicalVolume;
class G4Sphere;
class G4Tubs;
class G4VPhysicalVolume;
class G4UniformMagField;
class ElectricFieldSetup;
class F03FieldSetup;

class LXeAccelerator : public G4PVPlacement
{
 public:
 
  LXeAccelerator(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo,
                LXeDetectorConstruction* c);
  
  G4LogicalVolume* getlogicShape6() const {return flogicShape6;}
  G4LogicalVolume* getlogicShape10() const {return flogicShape10;}
  G4LogicalVolume* getlogicShape1() const {return flogicShape1;}
  G4LogicalVolume* getlogicShape4() const {return flogicShape4;}
  G4LogicalVolume* getlogicShape8() const {return flogicShape8;}
  G4LogicalVolume* getlogicShape12() const {return flogicShape12;}
  
 private:

  LXeDetectorConstruction* fConstructor;
  
  G4LogicalVolume* flogicShape6;
  G4LogicalVolume* flogicShape10;
  G4LogicalVolume* flogicShape1;
  G4LogicalVolume* flogicShape4;
  G4LogicalVolume* flogicShape8;
  G4LogicalVolume* flogicShape12;
};

#endif
