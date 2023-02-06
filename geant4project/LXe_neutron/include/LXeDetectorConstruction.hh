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
/// \file optical/LXe/include/LXeDetectorConstruction.hh
/// \brief Definition of the LXeDetectorConstruction class
//
//
#ifndef LXeDetectorConstruction_h
#define LXeDetectorConstruction_h 1

#include "LXeDetectorMessenger.hh"

#include "G4Cache.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Types.hh"
#include "ElectricFieldSetup.hh"
#include "F03FieldSetup.hh"
#include "G4Timer.hh"

class LXeMainVolume;
class LXePMTSD;
class LXeScintSD;
class LXeAccelerator;
class LXeModule1;
class LXeModule2;
class LXeModule3;
class LXeModuleCCL;

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4Material;
class G4MaterialPropertiesTable;
class G4Sphere;
class G4Tubs;
class G4VPhysicalVolume;

class LXeDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  LXeDetectorConstruction();
  ~LXeDetectorConstruction();
  
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_0;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_0;

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;
  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  // Functions to modify the geometry
  void SetDimensions(G4ThreeVector);
  void SetHousingThickness(G4double);
  void SetNX(G4int);
  void SetNY(G4int);
  void SetNZ(G4int);
  void SetPMTRadius(G4double);
  void SetDefaults();
  void SetSaveThreshold(G4int);

  // Get values
  G4int GetNX() const { return fNx; };
  G4int GetNY() const { return fNy; };
  G4int GetNZ() const { return fNz; };
  G4int GetSaveThreshold() const { return fSaveThreshold; };
  G4double GetScintX() const { return fScint_x; }
  G4double GetScintY() const { return fScint_y; }
  G4double GetScintZ() const { return fScint_z; }
  G4double GetHousingThickness() const { return fD_mtl; }
  G4double GetPMTRadius() const { return fOuterRadius_pmt; }
  G4double GetSlabZ() const { return fSlab_z; }
  G4double GetModeratorThickness() const { return fD_moderator; }
  
 
  void SetSphereOn(G4bool);
  static G4bool GetSphereOn() { return fSphereOn; }

  void SetHousingReflectivity(G4double);
  G4double GetHousingReflectivity() const { return fRefl; }

  void SetWLSSlabOn(G4bool b);
  G4bool GetWLSSlabOn() const { return fWLSslab; }

  void SetMainVolumeOn(G4bool b);
  G4bool GetMainVolumeOn() const { return fMainVolumeOn; }

  void SetNFibers(G4int n);
  G4int GetNFibers() const { return fNfibers; }

  void SetMainScintYield(G4double);
  void SetWLSScintYield(G4double);

 private:
  void DefineMaterials();

  LXeDetectorMessenger* fDetectorMessenger;

  G4Box* fExperimentalHall_box;
  G4LogicalVolume* fExperimentalHall_log;
  G4VPhysicalVolume* fExperimentalHall_phys;
  
  G4Box* fWall_box;
  G4LogicalVolume* fWall_log;
  G4VPhysicalVolume* fWall_phys;
  
  G4Box* fTarget_box;
  G4LogicalVolume* fTarget_log;
  G4VPhysicalVolume* fTarget_phys;

  // Materials & Elements
  G4Material* fLXe;
  G4Material* fAl;
  G4Element* fN;
  G4Element* fO;
  G4Material* fAir;
  G4Material* fVacuum;
  G4Element* fC;
  G4Element* fH;
  G4Material* fGlass;
  G4Material* fPstyrene;
  G4Material* fPMMA;
  G4Material* fPethylene1;
  G4Material* fPethylene2;
  G4Element* fLi6;
  G4Element* fCe;
  G4Element* fNa;
  G4Element* fAll;
  G4Element* fMg;
  G4Element* fSi;
  G4Element* fK;
  G4Element* fCa;
  G4Element* fFe;
  G4Element* fP;
  G4Element* fS;
  G4Element* fTi;
  G4Element* fMn;
  G4Element* fZn;
  G4Element* fZr;
  G4Element* fBa;
  G4Element* fPb;
  G4Element* fSr;
  G4Material* fGlassSci;
  G4Material* fConcrete;
  G4Element* fBe;
  G4Material* fBerillium;

  // Geometry
  G4double fScint_x;
  G4double fScint_y;
  G4double fScint_z;
  
  G4double fD_mtl;
  G4int fNx;
  G4int fNy;
  G4int fNz;
  G4int fSaveThreshold;
  G4double fOuterRadius_pmt;
  G4int fNfibers;
  static G4bool fSphereOn;
  G4double fRefl;
  G4bool fWLSslab;
  G4bool fMainVolumeOn;
  G4double fSlab_z;
  G4double fD_moderator;
  
  G4LogicalVolume* fScoringVolume;

  LXeMainVolume* fMainVolume;

  G4MaterialPropertiesTable* fLXe_mt;
  G4MaterialPropertiesTable* fMPTPStyrene;
  G4MaterialPropertiesTable* fMPTGlassSci;

  // Sensitive Detectors
  G4Cache<LXeScintSD*> fScint_SD;
  G4Cache<LXePMTSD*> fPmt_SD;
  // LXeAccelerator Fields
  const LXeAccelerator* fLXeAccelerator;
  //static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_0;
  //static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_0;
  static G4ThreadLocal F03FieldSetup* femFieldSetup1_0;
  static G4ThreadLocal F03FieldSetup* femFieldSetup2_0;
  static G4ThreadLocal F03FieldSetup* femFieldSetup3_0;
  static G4ThreadLocal F03FieldSetup* femFieldSetup4_0;
  // LXeModule1 Fields
  const LXeModule1* fLXeModule1;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_1;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_1;
  static G4ThreadLocal F03FieldSetup* femFieldSetup1_1;
  static G4ThreadLocal F03FieldSetup* femFieldSetup2_1;
  // LXeModule2 Fields
  const LXeModule2* fLXeModule2;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_2;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_2;
  static G4ThreadLocal F03FieldSetup* femFieldSetup1_2;
  static G4ThreadLocal F03FieldSetup* femFieldSetup2_2;
  // LXeModule3 Fields
  const LXeModule3* fLXeModule3;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_3;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_3;
  static G4ThreadLocal F03FieldSetup* femFieldSetup1_3;
  static G4ThreadLocal F03FieldSetup* femFieldSetup2_3;
  // LXeModuleCCL1 Fields
  const LXeModuleCCL* fLXeModuleCCL1;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_4;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_4;
  static G4ThreadLocal F03FieldSetup* femFieldSetup1_4;
  static G4ThreadLocal F03FieldSetup* femFieldSetup2_4;
  // LXeModuleCCL2 Fields
  const LXeModuleCCL* fLXeModuleCCL2;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_5;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_5;
  static G4ThreadLocal F03FieldSetup* femFieldSetup1_5;
  static G4ThreadLocal F03FieldSetup* femFieldSetup2_5;
  // LXeModuleCCL3 Fields
  const LXeModuleCCL* fLXeModuleCCL3;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_6;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_6;
  static G4ThreadLocal F03FieldSetup* femFieldSetup1_6;
  static G4ThreadLocal F03FieldSetup* femFieldSetup2_6;
  // LXeModuleCCL4 Fields
  const LXeModuleCCL* fLXeModuleCCL4;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup1_7;
  static G4ThreadLocal ElectricFieldSetup* felFieldSetup2_7;
  static G4ThreadLocal F03FieldSetup* femFieldSetup1_7;
  static G4ThreadLocal F03FieldSetup* femFieldSetup2_7;
};

#endif
