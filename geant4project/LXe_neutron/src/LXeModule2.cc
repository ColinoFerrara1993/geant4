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
/// \file optical/LXe/src/LXeMainVolume.cc
/// \brief Implementation of the LXeMainVolume class
//
//
#include "LXeModule2.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "ElectricFieldSetup.hh"
#include "F03FieldSetup.hh"
#include "G4TransportationManager.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Trd.hh"
#include "G4PVPlacement.hh"
#include "G4AutoDelete.hh"
#include "G4AssemblyVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeModule2::LXeModule2(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                             G4LogicalVolume* pMotherLogical, G4bool pMany,
                             G4int pCopyNo, LXeDetectorConstruction* c)
  // Pass info to the G4PVPlacement constructor
  : G4PVPlacement(pRot, tlate,
                  // Temp logical volume must be created here
                  new G4LogicalVolume(new G4Box("temp5", 40.*cm, 40.*cm, 1.79*m),
                                      G4Material::GetMaterial("Vacuum"), "temp5",
                                      0, 0, 0),
                  "housing5", pMotherLogical, pMany, pCopyNo)
  , fConstructor(c)
{
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 20.*cm, env_sizeZ = 1.*m;
  G4Material* env_mat = G4Material::GetMaterial("Air");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  //Material - Air vacuum
  G4Material* vacuum = G4Material::GetMaterial("Vacuum");
  
  //
  // Shape 0 - Quadrupole 1
  //
  G4Material* shape0_mat = nist->FindOrBuildMaterial("G4_Al");

  //
  // Shape 5 - Module 0
  //
  G4Material* shape5_mat = nist->FindOrBuildMaterial("G4_Cu");
  G4double shape5_rmin =  2.*mm, shape5_rmax = 3.*cm;
  G4double shape5_hz = 3.*cm;
  G4double shape5_phimin = 0.*deg, shape5_phimax = 360.*deg;
  G4ThreeVector pos5 = G4ThreeVector(0*cm, 0*cm, -50.*cm + shape5_hz );
  G4Tubs* solidShape5 =    
    new G4Tubs("Shape5",                     //its name
      shape5_rmin, shape5_rmax, shape5_hz, shape5_phimin, shape5_phimax);//its size
                
  G4LogicalVolume* logicShape5 =                         
    new G4LogicalVolume(solidShape5,         //its solid
                        shape5_mat,          //its material
                        "Shape5");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos5,                    //at position
                    logicShape5,             //its logical volume
                    "Shape5",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  
  //
  // Shape 6 - Beam pipe 4
  //
  G4double shape6_rmin =  0.*mm, shape6_rmax = 2.*mm;
  G4double shape6_hz = 3.*cm;
  G4double shape6_phimin = 0.*deg, shape6_phimax = 360.*deg;
  G4ThreeVector pos6 = G4ThreeVector(0*cm, 0*cm, -50.*cm + shape5_hz);
  G4Tubs* solidShape6 =    
    new G4Tubs("Shape6",                     //its name
      shape6_rmin, shape6_rmax, shape6_hz, shape6_phimin, shape6_phimax);//its size
                
  flogicShape6 =                         
    new G4LogicalVolume(solidShape6,         //its solid
                        vacuum,              //its material
                        "Shape6");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos6,                    //at position
                    flogicShape6,             //its logical volume
                    "Shape6",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  

  
  //
  // Beam pipe between Module 0 and Quadrupole 1
  //
  //
  // Beam pipe envelope
  //
  G4Material* beampipe1_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4double beampipe1_rmin =  3.*mm, beampipe1_rmax = 3.2*mm;
  G4double beampipe1_hz = 0.5*cm;
  G4double beampipe1_phimin = 0.*deg, beampipe1_phimax = 360.*deg;
  
  //Position
  G4ThreeVector beampipe1_pos = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + beampipe1_hz);
  
  G4Tubs* beampipe1_solid =    
    new G4Tubs("beampipe1solid", 
      beampipe1_rmin, beampipe1_rmax, beampipe1_hz, beampipe1_phimin, beampipe1_phimax);
  
  G4LogicalVolume* beampipe1_logic =                         
    new G4LogicalVolume(beampipe1_solid,     //its solid
                        beampipe1_mat,       //its material
                        "beampipe1logic");   //its name

  new G4PVPlacement(0,                       //no rotation
                    beampipe1_pos,           //at position
                    beampipe1_logic,         //its logical volume
                    "beampipe1physical",     //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //     
  // Vacuum inside beampipe
  //
  G4double vacuum1_rmin =  0.*mm, vacuum1_rmax = 3.*mm;
  G4double vacuum1_hz = 0.5*cm;
  G4double vacuum1_phimin = 0.*deg, vacuum1_phimax = 360.*deg;
  
  //Position
  G4ThreeVector vacuum1_pos = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + vacuum1_hz);
  
  G4Tubs* vacuum1_solid =    
    new G4Tubs("vacuum1solid", 
      vacuum1_rmin, vacuum1_rmax, vacuum1_hz, vacuum1_phimin, vacuum1_phimax);
  
  G4LogicalVolume* vacuum1_logic =                         
    new G4LogicalVolume(vacuum1_solid,       //its solid
                        vacuum,              //its material
                        "vacuum1logic");     //its name

  new G4PVPlacement(0,                       //no rotation
                    vacuum1_pos,             //at position
                    vacuum1_logic,           //its logical volume
                    "vacuum1physical",       //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );
  
  //
  // Shape 7 - Quadrupole 3
  //
  G4double shape7_rmin =  3.*mm, shape7_rmax = 1.*cm;
  G4double shape7_hz = 1.*cm;
  G4double shape7_phimin = 0.*deg, shape7_phimax = 360.*deg;
  G4ThreeVector pos7 = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + shape7_hz);

  G4Tubs* solidShape7 =    
    new G4Tubs("Shape7",                     //its name
      shape7_rmin, shape7_rmax, shape7_hz, shape7_phimin, shape7_phimax);//its size
                
  G4LogicalVolume* logicShape7 =                         
    new G4LogicalVolume(solidShape7,         //its solid
                        shape0_mat,          //its material
                        "Shape7");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos7,                    //at position
                    logicShape7,             //its logical volume
                    "Shape7",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // Shape 8 - Beam pipe 5
  //
  G4double shape8_rmin =  0.*mm, shape8_rmax = 3.*mm;
  G4double shape8_hz = 1.*cm;
  G4double shape8_phimin = 0.*deg, shape8_phimax = 360.*deg;
  G4ThreeVector pos8 = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + shape7_hz);

  G4Tubs* solidShape8 =    
    new G4Tubs("Shape8",                     //its name
      shape8_rmin, shape8_rmax, shape8_hz, shape8_phimin, shape8_phimax);//its size
                
  flogicShape8 =                         
    new G4LogicalVolume(solidShape8,         //its solid
                        vacuum,              //its material
                        "Shape8");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos8,                    //at position
                    flogicShape8,             //its logical volume
                    "Shape8",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // Beam pipe between Quadrupole 3 and Module 1
  //
  //
  // Beam pipe envelope
  //
  G4Material* beampipe2_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4double beampipe2_rmin =  3.*mm, beampipe2_rmax = 3.2*mm;
  G4double beampipe2_hz = 0.5*cm;
  G4double beampipe2_phimin = 0.*deg, beampipe2_phimax = 360.*deg;
  
  //Position
  G4ThreeVector beampipe2_pos = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + beampipe2_hz);
  
  G4Tubs* beampipe2_solid =    
    new G4Tubs("beampipe2solid", 
      beampipe2_rmin, beampipe2_rmax, beampipe2_hz, beampipe2_phimin, beampipe2_phimax);
  
  G4LogicalVolume* beampipe2_logic =                         
    new G4LogicalVolume(beampipe2_solid,     //its solid
                        beampipe2_mat,       //its material
                        "beampipe2logic");   //its name

  new G4PVPlacement(0,                       //no rotation
                    beampipe2_pos,           //at position
                    beampipe2_logic,         //its logical volume
                    "beampipe2physical",     //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //     
  // Vacuum inside beampipe
  //
  G4double vacuum2_rmin =  0.*mm, vacuum2_rmax = 3.*mm;
  G4double vacuum2_hz = 0.5*cm;
  G4double vacuum2_phimin = 0.*deg, vacuum2_phimax = 360.*deg;
  
  //Position
  G4ThreeVector vacuum2_pos = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + vacuum2_hz);
  
  G4Tubs* vacuum2_solid =    
    new G4Tubs("vacuum2solid", 
      vacuum2_rmin, vacuum2_rmax, vacuum2_hz, vacuum2_phimin, vacuum2_phimax);
  
  G4LogicalVolume* vacuum2_logic =                         
    new G4LogicalVolume(vacuum2_solid,       //its solid
                        vacuum,              //its material
                        "vacuum2logic");     //its name

  new G4PVPlacement(0,                       //no rotation
                    vacuum2_pos,             //at position
                    vacuum2_logic,           //its logical volume
                    "vacuum2physical",       //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );
  
  //                              //
  // -- Modulo 1 --- Quadrupole 3 //
  //                              //
                    
  //
  // Shape 9 - Module 1
  //
  G4Material* shape9_mat = nist->FindOrBuildMaterial("G4_Cu");
  G4double shape9_rmin =  2.*mm, shape9_rmax = 3.*cm;
  G4double shape9_hz = 3.*cm;
  G4double shape9_phimin = 0.*deg, shape9_phimax = 360.*deg;
  G4ThreeVector pos9 = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*vacuum2_hz + shape9_hz);
  G4Tubs* solidShape9 =    
    new G4Tubs("Shape9",                     //its name
      shape9_rmin, shape9_rmax, shape9_hz, shape9_phimin, shape9_phimax);//its size
                
  G4LogicalVolume* logicShape9 =                         
    new G4LogicalVolume(solidShape9,         //its solid
                        shape9_mat,          //its material
                        "Shape9");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos9,                    //at position
                    logicShape9,             //its logical volume
                    "Shape9",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  
  //
  // Shape 10 - Beam pipe 6
  //
  G4double shape10_rmin =  0.*mm, shape10_rmax = 2.*mm;
  G4double shape10_hz = 3.*cm;
  G4double shape10_phimin = 0.*deg, shape10_phimax = 360.*deg;
  G4ThreeVector pos10 = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*vacuum2_hz + shape9_hz);
  G4Tubs* solidShape10 =    
    new G4Tubs("Shape10",                     //its name
      shape10_rmin, shape10_rmax, shape10_hz, shape10_phimin, shape10_phimax);//its size
                
  flogicShape10 =                         
    new G4LogicalVolume(solidShape10,        //its solid
                        vacuum,              //its material
                        "Shape10");          //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos10,                   //at position
                    flogicShape10,            //its logical volume
                    "Shape10",               //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  
  
  
  //
  // Beam pipe between Module 1 and Quadrupole 4
  //
  //
  // Beam pipe envelope
  //
  G4Material* beampipe3_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4double beampipe3_rmin =  3.*mm, beampipe3_rmax = 3.2*mm;
  G4double beampipe3_hz = 0.5*cm;
  G4double beampipe3_phimin = 0.*deg, beampipe3_phimax = 360.*deg;
  
  //Position
  G4ThreeVector beampipe3_pos = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*beampipe2_hz + 2*shape9_hz + beampipe3_hz);
  
  G4Tubs* beampipe3_solid =    
    new G4Tubs("beampipe3solid", 
      beampipe3_rmin, beampipe3_rmax, beampipe3_hz, beampipe3_phimin, beampipe3_phimax);
  
  G4LogicalVolume* beampipe3_logic =                         
    new G4LogicalVolume(beampipe3_solid,     //its solid
                        beampipe3_mat,       //its material
                        "beampipe3logic");   //its name

  new G4PVPlacement(0,                       //no rotation
                    beampipe3_pos,           //at position
                    beampipe3_logic,         //its logical volume
                    "beampipe3physical",     //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //     
  // Vacuum inside beampipe
  //
  G4double vacuum3_rmin =  0.*mm, vacuum3_rmax = 3.*mm;
  G4double vacuum3_hz = 0.5*cm;
  G4double vacuum3_phimin = 0.*deg, vacuum3_phimax = 360.*deg;
  
  //Position
  G4ThreeVector vacuum3_pos = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*vacuum2_hz + 2*shape9_hz + vacuum3_hz);
  
  G4Tubs* vacuum3_solid =    
    new G4Tubs("vacuum3solid", 
      vacuum3_rmin, vacuum3_rmax, vacuum3_hz, vacuum3_phimin, vacuum3_phimax);
  
  G4LogicalVolume* vacuum3_logic =                         
    new G4LogicalVolume(vacuum3_solid,       //its solid
                        vacuum,              //its material
                        "vacuum3logic");     //its name

  new G4PVPlacement(0,                       //no rotation
                    vacuum3_pos,             //at position
                    vacuum3_logic,           //its logical volume
                    "vacuum3physical",       //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //
  // Shape 11 - Quadrupole 4
  //
  G4double shape11_rmin =  3.*mm, shape11_rmax = 1.*cm;
  G4double shape11_hz = 1.*cm;
  G4double shape11_phimin = 0.*deg, shape11_phimax = 360.*deg;
  G4ThreeVector pos11 = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*vacuum2_hz + 2*shape9_hz + 2*vacuum3_hz + shape11_hz);

  G4Tubs* solidShape11 =    
    new G4Tubs("Shape11",                    //its name
      shape11_rmin, shape11_rmax, shape11_hz, shape11_phimin, shape11_phimax);//its size
                
  G4LogicalVolume* logicShape11 =                         
    new G4LogicalVolume(solidShape11,        //its solid
                        shape0_mat,          //its material
                        "Shape11");          //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos11,                   //at position
                    logicShape11,            //its logical volume
                    "Shape11",               //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // Shape 12 - Beam pipe 7
  //
  G4double shape12_rmin =  0.*mm, shape12_rmax = 3.*mm;
  G4double shape12_hz = 1.*cm;
  G4double shape12_phimin = 0.*deg, shape12_phimax = 360.*deg;
  G4ThreeVector pos12 = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*vacuum2_hz + 2*shape9_hz + 2*vacuum3_hz + shape11_hz);

  G4Tubs* solidShape12 =    
    new G4Tubs("Shape12",                     //its name
      shape12_rmin, shape12_rmax, shape12_hz, shape12_phimin, shape12_phimax);//its size
                
  flogicShape12 =                         
    new G4LogicalVolume(solidShape12,        //its solid
                        vacuum,              //its material
                        "Shape12");          //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos12,                   //at position
                    flogicShape12,            //its logical volume
                    "Shape12",               //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                    
  //
  // Beam pipe between Quadrupole 4 and Module 2
  //
  //
  // Beam pipe envelope
  //
  G4Material* beampipe4_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4double beampipe4_rmin =  3.*mm, beampipe4_rmax = 3.2*mm;
  G4double beampipe4_hz = 0.5*cm;
  G4double beampipe4_phimin = 0.*deg, beampipe4_phimax = 360.*deg;
  
  //Position
  G4ThreeVector beampipe4_pos = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*beampipe2_hz + 2*shape9_hz + 2*beampipe3_hz + 2*shape11_hz + beampipe4_hz);
  
  G4Tubs* beampipe4_solid =    
    new G4Tubs("beampipe4solid", 
      beampipe4_rmin, beampipe4_rmax, beampipe4_hz, beampipe4_phimin, beampipe4_phimax);
  
  G4LogicalVolume* beampipe4_logic =                         
    new G4LogicalVolume(beampipe4_solid,     //its solid
                        beampipe4_mat,       //its material
                        "beampipe4logic");   //its name

  new G4PVPlacement(0,                       //no rotation
                    beampipe4_pos,           //at position
                    beampipe4_logic,         //its logical volume
                    "beampipe4physical",     //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //     
  // Vacuum inside beampipe
  //
  G4double vacuum4_rmin =  0.*mm, vacuum4_rmax = 3.*mm;
  G4double vacuum4_hz = 0.5*cm;
  G4double vacuum4_phimin = 0.*deg, vacuum4_phimax = 360.*deg;
  
  //Position
  G4ThreeVector vacuum4_pos = G4ThreeVector(0*cm, 0*cm, -50.*cm + 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*vacuum2_hz + 2*shape9_hz + 2*vacuum3_hz + 2*shape11_hz + beampipe4_hz);
  
  G4Tubs* vacuum4_solid =    
    new G4Tubs("vacuum4solid", 
      vacuum4_rmin, vacuum4_rmax, vacuum4_hz, vacuum4_phimin, vacuum4_phimax);
  
  G4LogicalVolume* vacuum4_logic =                         
    new G4LogicalVolume(vacuum4_solid,       //its solid
                        vacuum,              //its material
                        "vacuum4logic");     //its name

  new G4PVPlacement(0,                       //no rotation
                    vacuum4_pos,             //at position
                    vacuum4_logic,           //its logical volume
                    "vacuum4physical",       //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  
  //
  // Construct Assembled Volume
  //
  G4AssemblyVolume* assembly = new G4AssemblyVolume();

  assembly->AddPlacedVolume(logicShape5,     pos5,          0);
  assembly->AddPlacedVolume(flogicShape6,     pos6,          0);
  assembly->AddPlacedVolume(beampipe1_logic, beampipe1_pos, 0);
  assembly->AddPlacedVolume(vacuum1_logic,   vacuum1_pos,   0);
  assembly->AddPlacedVolume(logicShape7,     pos7,          0);
  assembly->AddPlacedVolume(flogicShape8,     pos8,          0);
  assembly->AddPlacedVolume(beampipe2_logic, beampipe2_pos, 0);
  assembly->AddPlacedVolume(vacuum2_logic,   vacuum2_pos,   0);
  assembly->AddPlacedVolume(logicShape9,     pos9,          0);
  assembly->AddPlacedVolume(flogicShape10,    pos10,         0);
  assembly->AddPlacedVolume(beampipe3_logic, beampipe3_pos, 0);
  assembly->AddPlacedVolume(vacuum3_logic,   vacuum3_pos,   0);
  assembly->AddPlacedVolume(logicShape11,    pos11,         0);
  assembly->AddPlacedVolume(flogicShape12,    pos12,         0);
  assembly->AddPlacedVolume(beampipe4_logic, beampipe4_pos, 0);
  assembly->AddPlacedVolume(vacuum4_logic,   vacuum4_pos,   0);
  
  G4double z_pos = 2*shape5_hz + 2*vacuum1_hz + 2*shape7_hz + 2*vacuum2_hz + 2*shape9_hz + 2*vacuum3_hz + 2*shape11_hz + 2*beampipe4_hz;
  
  for (int i = 0; i < 4; i++)
  { 
     G4ThreeVector position = G4ThreeVector(0*cm, 0*cm, z_pos*(i+1));     
     assembly->MakeImprint(logicEnv, position, 0);  
  }
 
  SetLogicalVolume(logicEnv);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
