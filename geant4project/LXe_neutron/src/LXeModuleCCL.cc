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
#include "LXeModuleCCL.hh"

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
// CADMesh //
#include "CADMesh.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeModuleCCL::LXeModuleCCL(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                             G4LogicalVolume* pMotherLogical, G4bool pMany,
                             G4int pCopyNo, LXeDetectorConstruction* c)
  // Pass info to the G4PVPlacement constructor
  : G4PVPlacement(pRot, tlate,
                  // Temp logical volume must be created here
                  new G4LogicalVolume(new G4Box("temp6", 1.*m, 1.*m, 1.*m),
                                      G4Material::GetMaterial("Vacuum"), "temp6",
                                      0, 0, 0),
                  "housing6", pMotherLogical, pMany, pCopyNo)
  , fConstructor(c)
{
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 25.*cm;
  G4double world_sizeZ  = 65.5*cm;
  G4Material* world_mat = G4Material::GetMaterial("Air");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
  /*                     
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0., 0., 0.),//at (0,0,0)
                    logicWorld,                //its logical volume
                    "Envelope",              //its name
                    0,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  */

  //Material - Air vacuum
  G4Material* vacuum = G4Material::GetMaterial("Vacuum");


  
  //
  // Shape 1 - Beam pipe 1
  //
  G4double shape1_rmin =  0.*mm, shape1_rmax = 4.*mm;
  G4double shape1_hz = 12.5*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4ThreeVector pos1 = G4ThreeVector(0*cm, 0*cm, -20.25*cm);
  G4Tubs* solidShape1 =    
    new G4Tubs("Shape1",                     //its name
      shape1_rmin, shape1_rmax, shape1_hz, shape1_phimin, shape1_phimax);//its size
                
  flogicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        vacuum,              //its material
                        "Shape1");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    flogicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Shape 1
  //
  // CADMesh import function
  G4Material* shape_mat = nist->FindOrBuildMaterial("G4_Cu");
   
  auto cubo_mesh = CADMesh::TessellatedMesh::FromSTL("./PiastraSCL.STL");//import function
  //auto cubo_mesh = CADMesh::TessellatedMesh::FromSTL("./Modulo_0_SCDTL.STL");//import function              
  auto logicShape = 
    new G4LogicalVolume( 
                        cubo_mesh->GetSolid()//its solid
                       ,shape_mat           //its material
                       ,"Shape"             //its name
                       );
  
  
  G4ThreeVector pos0 = G4ThreeVector(-6.05*cm, -11.3*cm, -32.75*cm);  
          
  //                 
  //Assembly
  //
  G4AssemblyVolume* assembly = new G4AssemblyVolume();
  G4double theta = 180.*deg;
  assembly->AddPlacedVolume(logicShape, pos0, 0);
  
  for (int i = 0; i < 16; i++)
  { 
     G4ThreeVector pos = G4ThreeVector(0*cm, 0*cm, 1.6*i*cm);
     G4RotationMatrix* matrix = new G4RotationMatrix();
     matrix->rotateZ(theta*i);  
     assembly->MakeImprint(logicWorld, pos, matrix);  
  }
  
  //
  // Beam pipe between Shape 1 and Quadrupole 1
  //
  //
  // Beam pipe envelope
  //
  G4Material* beampipe1_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4double beampipe1_rmin =  4.*mm, beampipe1_rmax = 4.2*mm;
  G4double beampipe1_hz = 1.5*cm;
  G4double beampipe1_phimin = 0.*deg, beampipe1_phimax = 360.*deg;
  
  //Position
  G4ThreeVector beampipe1_pos = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + beampipe1_hz);
  
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
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //     
  // Vacuum inside beampipe
  //
  G4double vacuum1_rmin =  0.*mm, vacuum1_rmax = 4.*mm;
  G4double vacuum1_hz = 1.5*cm;
  G4double vacuum1_phimin = 0.*deg, vacuum1_phimax = 360.*deg;
  
  //Position
  G4ThreeVector vacuum1_pos = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + vacuum1_hz);
  
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
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );
                   
                   
  //
  // Shape 2 - Quadrupole 1
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_Al");
  G4double shape2_rmin =  3.*mm, shape2_rmax = 1.*cm;
  G4double shape2_hz = 1.*cm;
  G4double shape2_phimin = 0.*deg, shape2_phimax = 360.*deg;
  G4ThreeVector pos2 = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + shape2_hz);

  G4Tubs* solidShape2 =    
    new G4Tubs("Shape2",                     //its name
      shape2_rmin, shape2_rmax, shape2_hz, shape2_phimin, shape2_phimax);//its size
                
  G4LogicalVolume* logicShape2 =                         
    new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // Shape 3 - Quadrupole's beam pipe
  //
  G4double shape3_rmin =  0.*mm, shape3_rmax = 3.*mm;
  G4double shape3_hz = 1.*cm;
  G4double shape3_phimin = 0.*deg, shape3_phimax = 360.*deg;
  G4ThreeVector pos3 = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + shape3_hz);

  G4Tubs* solidShape3 =    
    new G4Tubs("Shape3",                     //its name
      shape3_rmin, shape3_rmax, shape3_hz, shape3_phimin, shape3_phimax);//its size
                
  flogicShape3 =                         
    new G4LogicalVolume(solidShape3,         //its solid
                        vacuum,              //its material
                        "Shape3");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos3,                    //at position
                    flogicShape3,             //its logical volume
                    "Shape3",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                    
  //
  // Beam pipe between Quadrupole 1 and Shape 4
  //
  //
  // Beam pipe envelope
  //
  G4Material* beampipe2_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4double beampipe2_rmin =  4.*mm, beampipe2_rmax = 4.2*mm;
  G4double beampipe2_hz = 1.5*cm;
  G4double beampipe2_phimin = 0.*deg, beampipe2_phimax = 360.*deg;
  
  //Position
  G4ThreeVector beampipe2_pos = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + beampipe2_hz);
  
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
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //     
  // Vacuum inside beampipe
  //
  G4double vacuum2_rmin =  0.*mm, vacuum2_rmax = 4.*mm;
  G4double vacuum2_hz = 1.5*cm;
  G4double vacuum2_phimin = 0.*deg, vacuum2_phimax = 360.*deg;
  
  //Position
  G4ThreeVector vacuum2_pos = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + vacuum2_hz);
  
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
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );
                   
                   
  //
  //Assembly 2
  //
  G4ThreeVector pos_assembly2 = G4ThreeVector(-6.05*cm, -11.3*cm, -20.3*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + 2*vacuum2_hz);              
  G4AssemblyVolume* assembly2 = new G4AssemblyVolume();
  assembly2->AddPlacedVolume(logicShape, pos_assembly2, 0);
  
  for (int i = 0; i < 16; i++)
  { 
     G4ThreeVector pos = G4ThreeVector(0*cm, 0*cm, 1.6*i*cm);
     G4RotationMatrix* matrix = new G4RotationMatrix();
     matrix->rotateZ(180.*deg*i);  
     assembly2->MakeImprint(logicWorld, pos, matrix);  
  }
  
  //
  // Beam pipe 2 - Module 2
  //
  G4double shape4_rmin =  0.*mm, shape4_rmax = 4.*mm;
  G4double shape4_hz = 12.5*cm;
  G4double shape4_phimin = 0.*deg, shape4_phimax = 360.*deg;
  G4ThreeVector pos4 = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + 2*vacuum2_hz + shape4_hz);
  G4Tubs* solidShape4 =    
    new G4Tubs("Shape4",                     //its name
      shape4_rmin, shape4_rmax, shape4_hz, shape4_phimin, shape4_phimax);//its size
                
  flogicShape4 =                         
    new G4LogicalVolume(solidShape4,         //its solid
                        vacuum,              //its material
                        "Shape4");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    flogicShape4,             //its logical volume
                    "Shape4",                //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // Beam pipe between Shape4 and Quadrupole 2
  //
  //
  // Beam pipe envelope
  //
  G4Material* beampipe3_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4double beampipe3_rmin =  4.*mm, beampipe3_rmax = 4.2*mm;
  G4double beampipe3_hz = 1.5*cm;
  G4double beampipe3_phimin = 0.*deg, beampipe3_phimax = 360.*deg;
  
  //Position
  G4ThreeVector beampipe3_pos = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + 2*vacuum2_hz + 2*shape4_hz + beampipe3_hz);
  
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
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //     
  // Vacuum inside beampipe
  //
  G4double vacuum3_rmin =  0.*mm, vacuum3_rmax = 4.*mm;
  G4double vacuum3_hz = 1.5*cm;
  G4double vacuum3_phimin = 0.*deg, vacuum3_phimax = 360.*deg;
  
  //Position
  G4ThreeVector vacuum3_pos = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + 2*vacuum2_hz + 2*shape4_hz + beampipe3_hz);
  
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
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );
                   
  //
  // Shape 5 - Quadrupole 2
  //
  G4Material* shape5_mat = nist->FindOrBuildMaterial("G4_Al");
  G4double shape5_rmin =  3.*mm, shape5_rmax = 1.*cm;
  G4double shape5_hz = 1.*cm;
  G4double shape5_phimin = 0.*deg, shape5_phimax = 360.*deg;
  G4ThreeVector pos5 = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + 2*vacuum2_hz + 2*shape4_hz + 2*beampipe3_hz + shape5_hz);

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
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // Shape 6 - Quadrupole's beam pipe
  //
  G4double shape6_rmin =  0.*mm, shape6_rmax = 3.*mm;
  G4double shape6_hz = 1.*cm;
  G4double shape6_phimin = 0.*deg, shape6_phimax = 360.*deg;
  G4ThreeVector pos6 = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + 2*vacuum2_hz + 2*shape4_hz + 2*beampipe3_hz + shape6_hz);

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
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                    
  //
  // Beam pipe over Quadrupole 2
  //
  //
  // Beam pipe envelope
  //
  G4Material* beampipe4_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4double beampipe4_rmin =  4.*mm, beampipe4_rmax = 4.2*mm;
  G4double beampipe4_hz = 1.25*cm;
  G4double beampipe4_phimin = 0.*deg, beampipe4_phimax = 360.*deg;
  
  //Position
  G4ThreeVector beampipe4_pos = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + 2*vacuum2_hz + 2*shape4_hz + 2*beampipe3_hz + 2*shape6_hz + beampipe4_hz);
  
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
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );

  //     
  // Vacuum inside beampipe
  //
  G4double vacuum4_rmin =  0.*mm, vacuum4_rmax = 4.*mm;
  G4double vacuum4_hz = 1.25*cm;
  G4double vacuum4_phimin = 0.*deg, vacuum4_phimax = 360.*deg;
  
  //Position
  G4ThreeVector vacuum4_pos = G4ThreeVector(0*cm, 0*cm, -20.25*cm + shape1_hz + 2*vacuum1_hz + 2*shape3_hz + 2*vacuum2_hz + 2*shape4_hz + 2*beampipe3_hz + 2*shape6_hz + beampipe4_hz);
  
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
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps            //overlaps checking
                   );
  
  SetLogicalVolume(logicWorld);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
