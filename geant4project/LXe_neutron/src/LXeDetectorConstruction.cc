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
/// \file optical/LXe/src/LXeDetectorConstruction.cc
/// \brief Implementation of the LXeDetectorConstruction class
//
//
#include "LXeDetectorConstruction.hh"

#include "LXeDetectorMessenger.hh"
#include "LXeMainVolume.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"
#include "LXeWLSSlab.hh"
#include "LXeAccelerator.hh"
#include "LXeModule1.hh"
#include "LXeModule2.hh"
#include "LXeModule3.hh"
#include "LXeModuleCCL.hh"

#include "G4Types.hh"
#include "ElectricFieldSetup.hh"
#include "F03FieldSetup.hh"

#include "LXeSteppingAction.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4VisAttributes.hh"
#include "G4AutoDelete.hh"
#include "G4Isotope.hh"

using namespace CLHEP;

G4bool LXeDetectorConstruction::fSphereOn = true;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::LXeDetectorConstruction()
  : fLXe_mt(nullptr)
  , fMPTPStyrene(nullptr)
  , fScoringVolume(nullptr)
  , fMPTGlassSci(nullptr)
{
  fExperimentalHall_box  = nullptr;
  fExperimentalHall_log  = nullptr;
  fExperimentalHall_phys = nullptr;

  fLXe = fAl = fAir = fVacuum = fGlass = fGlassSci = nullptr;
  fPstyrene = fPMMA = fPethylene1 = fPethylene2 = nullptr;
  fConcrete = nullptr;

  fN = fO = fC = fH = fLi6 = fCe = nullptr;
  fAll = nullptr;
  fNa = fMg = fSi = fK = fCa = fFe = fP = fS = fTi = fMn = fZn = fZr = fBa = fPb = fSr = nullptr;
  fBe = nullptr;
  fIsoBe = nullptr;

  fSaveThreshold = 0;
  SetDefaults();

  DefineMaterials();
  fDetectorMessenger = new LXeDetectorMessenger(this);

  felFieldSetup1_0 = 0;
  felFieldSetup2_0 = 0;
  femFieldSetup1_0 = 0;
  femFieldSetup2_0 = 0;
  femFieldSetup3_0 = 0;
  femFieldSetup4_0 = 0;
  
  felFieldSetup1_1 = 0;
  felFieldSetup2_1 = 0;
  femFieldSetup1_1 = 0;
  femFieldSetup2_1 = 0;
  
  felFieldSetup1_2 = 0;
  felFieldSetup2_2 = 0;
  femFieldSetup1_2 = 0;
  femFieldSetup2_2 = 0;
  
  felFieldSetup1_3 = 0;
  felFieldSetup2_3 = 0;
  femFieldSetup1_3 = 0;
  femFieldSetup2_3 = 0;
  
  felFieldSetup1_4 = 0;
  felFieldSetup2_4 = 0;
  femFieldSetup1_4 = 0;
  femFieldSetup2_4 = 0;
  
  felFieldSetup1_5 = 0;
  felFieldSetup2_5 = 0;
  femFieldSetup1_5 = 0;
  femFieldSetup2_5 = 0;
  
  felFieldSetup1_6 = 0;
  felFieldSetup2_6 = 0;
  femFieldSetup1_6 = 0;
  femFieldSetup2_6 = 0;
  
  felFieldSetup1_7 = 0;
  felFieldSetup2_7 = 0;
  femFieldSetup1_7 = 0;
  femFieldSetup2_7 = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::~LXeDetectorConstruction()
{
  if(fMainVolume)
  {
    delete fMainVolume;
  }
  delete fLXe_mt;
  delete fDetectorMessenger;
  delete fMPTPStyrene;
  delete fScoringVolume;
  delete fMPTGlassSci;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::DefineMaterials()
{
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  G4int polyPMMA = 1;
  G4int nC_PMMA  = 3 + 2 * polyPMMA;
  G4int nH_PMMA  = 6 + 2 * polyPMMA;

  G4int polyeth = 1;
  G4int nC_eth  = 2 * polyeth;
  G4int nH_eth  = 4 * polyeth;

  //***Elements
  fH = new G4Element("H", "H", z = 1., a = 1.01 * g / mole);
  fC = new G4Element("C", "C", z = 6., a = 12.01 * g / mole);
  fN = new G4Element("N", "N", z = 7., a = 14.01 * g / mole);
  fO = new G4Element("O", "O", z = 8., a = 16.00 * g / mole);
  fLi6 = new G4Element("Li6", "Li6", z = 3., a = 6.00 * g / mole);
  fCe = new G4Element("Ce", "Ce", z = 58., a = 140.12 * g / mole);
  
  fAll = new G4Element("Al", "Al", z = 13., a = 26.98 * g / mole);

  fNa = new G4Element("Na", "Na", z = 11., a = 22.99 * g / mole);
  fMg = new G4Element("Mg", "Mg", z = 12., a = 24.30 * g / mole);
  fSi = new G4Element("Si", "Si", z = 14., a = 28.08 * g / mole);
  fK = new G4Element("K", "K", z = 19., a = 39.10 * g / mole);
  fCa = new G4Element("Ca", "Ca", z = 20., a = 40.078 * g / mole);
  fFe = new G4Element("Fe", "Fe", z = 26., a = 55.84 * g / mole);
  fP = new G4Element("P", "P", z = 15., a = 30.97 * g / mole);
  fS = new G4Element("S", "S", z = 16., a = 32.065 * g / mole);
  fTi = new G4Element("Ti", "Ti", z = 22., a = 47.867 * g / mole);
  fMn = new G4Element("Mn", "Mn", z = 25., a = 54.94 * g / mole);
  fZn = new G4Element("Zn", "Zn", z = 30., a = 65.409 * g / mole);
  fZr = new G4Element("Zr", "Zr", z = 40., a = 91.224 * g / mole);
  fBa = new G4Element("Ba", "Ba", z = 56., a = 137.327 * g / mole);
  fPb = new G4Element("Pb", "Pb", z = 82., a = 207.2 * g / mole);
  fSr = new G4Element("Sr", "Sr", z = 38., a = 87.62 * g / mole);
  
  fIsoBe = new G4Isotope("_Be_", z = 4., a = 9., 9.01 * g / mole);
  fBe = new G4Element("Be", "Be", 1);
  fBe->AddIsotope(fIsoBe, 100. * perCent);

  //***Materials
  // Liquid Xenon
  fLXe = new G4Material("LXe", z = 54., a = 131.29 * g / mole,
                        density = 3.020 * g / cm3);
  // Aluminum
  fAl = new G4Material("Al", z = 13., a = 26.98 * g / mole,
                       density = 2.7 * g / cm3);
  // Vacuum
  fVacuum = new G4Material("Vacuum", density = universe_mean_density, 3, kStateGas,
                           0.1 * kelvin, 1.e-7 * pascal);
  fVacuum->AddElement(fH, 1. * perCent);
  fVacuum->AddElement(fN, 69. * perCent);
  fVacuum->AddElement(fO, 30. * perCent);
  // Air
  fAir = new G4Material("Air", density = 1.29 * mg / cm3, 3);
  fAir->AddElement(fH, 1. * perCent);
  fAir->AddElement(fN, 69. * perCent);
  fAir->AddElement(fO, 30. * perCent);
  // Glass
  fGlass = new G4Material("Glass", density = 1.032 * g / cm3, 2);
  fGlass->AddElement(fC, 91.533 * perCent);
  fGlass->AddElement(fH, 8.467 * perCent);
  // Glass scintillator
  fGlassSci = new G4Material("GlassSci", density = 2.40 * g / cm3, 4);
  fGlassSci->AddElement(fC, 84.18 * perCent);
  fGlassSci->AddElement(fH, 7.82 * perCent);
  fGlassSci->AddElement(fLi6, 7.5 * perCent);
  fGlassSci->AddElement(fCe, 0.5 * perCent);
  // Polystyrene
  fPstyrene = new G4Material("Polystyrene", density = 1.03 * g / cm3, 2);
  fPstyrene->AddElement(fC, 8);
  fPstyrene->AddElement(fH, 8);
  // Fiber(PMMA)
  fPMMA = new G4Material("PMMA", density = 1190. * kg / m3, 3);
  fPMMA->AddElement(fH, nH_PMMA);
  fPMMA->AddElement(fC, nC_PMMA);
  fPMMA->AddElement(fO, 2);
  // Cladding(polyethylene)
  fPethylene1 = new G4Material("Pethylene1", density = 1200. * kg / m3, 2);
  fPethylene1->AddElement(fH, nH_eth);
  fPethylene1->AddElement(fC, nC_eth);
  // Double cladding(flourinated polyethylene)
  fPethylene2 = new G4Material("Pethylene2", density = 1400. * kg / m3, 2);
  fPethylene2->AddElement(fH, nH_eth);
  fPethylene2->AddElement(fC, nC_eth);
  // Concrete composition
  fConcrete = new G4Material("Concrete", density = 2.42 * g / cm3, 19);
  fConcrete->AddElement(fO, 49.2875 * perCent);
  fConcrete->AddElement(fC, 5.62 * perCent);
  fConcrete->AddElement(fH, 0.6 * perCent);
  fConcrete->AddElement(fNa, 0.453 * perCent);
  fConcrete->AddElement(fMg, 0.663 * perCent);
  fConcrete->AddElement(fAll, 2.063 * perCent);
  fConcrete->AddElement(fSi, 18.867 * perCent);
  fConcrete->AddElement(fK, 0.656 * perCent);
  fConcrete->AddElement(fCa, 20.091 * perCent);
  fConcrete->AddElement(fFe, 1.118 * perCent);
  fConcrete->AddElement(fP, 0.048 * perCent);
  fConcrete->AddElement(fS, 0.012 * perCent);
  fConcrete->AddElement(fTi, 0.347 * perCent);
  fConcrete->AddElement(fMn, 0.0387 * perCent);
  fConcrete->AddElement(fZn, 0.0241 * perCent);
  fConcrete->AddElement(fZr, 0.0074 * perCent);
  fConcrete->AddElement(fBa, 0.0179 * perCent);
  fConcrete->AddElement(fPb, 0.0464 * perCent);
  fConcrete->AddElement(fSr, 0.04 * perCent);
  
  // Berillium/Lithium material
  fBerillium = new G4Material("Berillium", density = 1.848 * g / cm3, 1);
  fBerillium->AddElement(fBe, 100. * perCent);

  //***Material properties tables

  std::vector<G4double> lxe_Energy = { 7.0 * eV, 7.07 * eV, 7.14 * eV };

  std::vector<G4double> lxe_SCINT = { 0.1, 1.0, 0.1 };
  std::vector<G4double> lxe_RIND  = { 1.59, 1.57, 1.54 };
  std::vector<G4double> lxe_ABSL  = { 35. * cm, 35. * cm, 35. * cm };
  fLXe_mt = new G4MaterialPropertiesTable();
  fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT1", lxe_Energy, lxe_SCINT);
  fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT2", lxe_Energy, lxe_SCINT);
  fLXe_mt->AddProperty("RINDEX", lxe_Energy, lxe_RIND);
  fLXe_mt->AddProperty("ABSLENGTH", lxe_Energy, lxe_ABSL);
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
  fLXe_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
  fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20. * ns);
  fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
  fLXe->SetMaterialPropertiesTable(fLXe_mt);

  // Set the Birks Constant for the LXe scintillator
  fLXe->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  std::vector<G4double> glass_RIND      = { 1.49, 1.49, 1.49 };
  std::vector<G4double> glass_AbsLength = { 420. * cm, 420. * cm, 420. * cm };
  G4MaterialPropertiesTable* glass_mt   = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH", lxe_Energy, glass_AbsLength);
  glass_mt->AddProperty("RINDEX", lxe_Energy, glass_RIND);
  fGlass->SetMaterialPropertiesTable(glass_mt);

  std::vector<G4double> vacuum_Energy  = { 2.0 * eV, 7.0 * eV, 7.14 * eV };
  std::vector<G4double> vacuum_RIND    = { 1., 1., 1. };
  G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND);
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
  fAir->SetMaterialPropertiesTable(vacuum_mt);  // Give air the same rindex

  std::vector<G4double> wls_Energy = { 2.00 * eV, 2.87 * eV, 2.90 * eV,
                                       3.47 * eV };

  std::vector<G4double> rIndexPstyrene = { 1.58, 1.58, 1.58, 1.58 };
  std::vector<G4double> absorption1    = { 250. * cm, 250. * cm, 250. * cm, 250. * cm };
  std::vector<G4double> scintilFast    = { 0.0, 0.0, 1.0, 1.0 };
  fMPTPStyrene = new G4MaterialPropertiesTable();
  fMPTPStyrene->AddProperty("RINDEX", wls_Energy, rIndexPstyrene);
  fMPTPStyrene->AddProperty("ABSLENGTH", wls_Energy, absorption1);
  fMPTPStyrene->AddProperty("SCINTILLATIONCOMPONENT1", wls_Energy, scintilFast);
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
  fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
  fMPTPStyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.4 * ns);
  fPstyrene->SetMaterialPropertiesTable(fMPTPStyrene);

  // Set the Birks Constant for the Polystyrene scintillator
  fPstyrene->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  
  //New Material
  std::vector<G4double> rIndex_GlassSci = { 1.59, 1.59, 1.59, 1.59 };
  std::vector<G4double> AbsLength_GlassSci    = { 420. * cm, 420. * cm, 420. * cm, 420. * cm };
  std::vector<G4double> scintilFast_GlassSci    = { 0.0, 0.0, 1.0, 1.0 };
  fMPTGlassSci = new G4MaterialPropertiesTable();
  fMPTGlassSci->AddProperty("RINDEX", wls_Energy, rIndex_GlassSci);
  fMPTGlassSci->AddProperty("ABSLENGTH", wls_Energy, AbsLength_GlassSci);
  fMPTGlassSci->AddProperty("SCINTILLATIONCOMPONENT1", wls_Energy, scintilFast_GlassSci);
  fMPTGlassSci->AddProperty("SCINTILLATIONCOMPONENT2", wls_Energy, scintilFast_GlassSci);
  fMPTGlassSci->AddConstProperty("SCINTILLATIONYIELD", 10000. / MeV);
  fMPTGlassSci->AddConstProperty("RESOLUTIONSCALE", 1.0);
  fMPTGlassSci->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 18. * ns);
  fMPTGlassSci->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
  fMPTGlassSci->AddConstProperty("SCINTILLATIONYIELD1", 1.);
  fMPTGlassSci->AddConstProperty("SCINTILLATIONYIELD2", 0.);
  fGlassSci->SetMaterialPropertiesTable(fMPTGlassSci);

  // Set the Birks Constant for the GlassSci scintillator
  fGlassSci->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  

  std::vector<G4double> RefractiveIndexFiber = { 1.6, 1.6, 1.6, 1.6 };
  std::vector<G4double> AbsFiber    = { 9.0 * m, 9.0 * m, 0.1 * mm, 0.1 * mm };
  std::vector<G4double> EmissionFib = { 1.0, 1.0, 0.0, 0.0 };
  G4MaterialPropertiesTable* fiberProperty = new G4MaterialPropertiesTable();
  fiberProperty->AddProperty("RINDEX", wls_Energy, RefractiveIndexFiber);
  fiberProperty->AddProperty("WLSABSLENGTH", wls_Energy, AbsFiber);
  fiberProperty->AddProperty("WLSCOMPONENT", wls_Energy, EmissionFib);
  fiberProperty->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
  fPMMA->SetMaterialPropertiesTable(fiberProperty);

  std::vector<G4double> RefractiveIndexClad1 = { 1.49, 1.49, 1.49, 1.49 };
  G4MaterialPropertiesTable* clad1Property   = new G4MaterialPropertiesTable();
  clad1Property->AddProperty("RINDEX", wls_Energy, RefractiveIndexClad1);
  clad1Property->AddProperty("ABSLENGTH", wls_Energy, AbsFiber);
  fPethylene1->SetMaterialPropertiesTable(clad1Property);

  std::vector<G4double> RefractiveIndexClad2 = { 1.42, 1.42, 1.42, 1.42 };
  G4MaterialPropertiesTable* clad2Property   = new G4MaterialPropertiesTable();
  clad2Property->AddProperty("RINDEX", wls_Energy, RefractiveIndexClad2);
  clad2Property->AddProperty("ABSLENGTH", wls_Energy, AbsFiber);
  fPethylene2->SetMaterialPropertiesTable(clad2Property);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LXeDetectorConstruction::Construct()
{
  // The experimental hall walls are all 1m away from housing walls
  G4double expHall_x = fScint_x + fD_mtl + 1.1 * m;
  G4double expHall_y = fScint_y + fD_mtl + 1.1 * m;
  G4double expHall_z = ((1.79 + 6. + 4.*1.31 + 0.67 + 3.2 + 8.)/4.)*m;
  
  G4double wall_x = expHall_x + 1.20*m;
  G4double wall_y = expHall_y + 1.20*m;
  G4double wall_z = expHall_z + 1.20*m;
  
  //Create wall over experimental hall
  fWall_box = 
    new G4Box("wall_box", wall_x, wall_y, wall_z);
  fWall_log =
    new G4LogicalVolume(fWall_box, fConcrete, "wall_log", 0, 0, 0);
  fWall_phys = 
    new G4PVPlacement(0, G4ThreeVector(), fWall_log, "wall", 0, false, 0);
  
  //Create experimental hall
  fExperimentalHall_box =
    new G4Box("expHall_box", expHall_x, expHall_y, expHall_z);
  fExperimentalHall_log =
    new G4LogicalVolume(fExperimentalHall_box, fAir, "expHall_log", 0, 0, 0);
  fExperimentalHall_phys = 
    new G4PVPlacement(0, G4ThreeVector(), fExperimentalHall_log, "expHall", fWall_log, false, 0);
    
  fExperimentalHall_log->SetVisAttributes(G4Colour(0.8, 0.8, 0.8));
  fWall_log->SetVisAttributes(G4Colour(0.8, 0.8, 0.8));
  
  //Create experimental target of Berillium/Lithium
  G4double target_x = 4.*cm;
  G4double target_y = 4.*cm;
  G4double target_z = 0.04*cm;
  
  fTarget_box =
    new G4Box("target_box", target_x, target_y, target_z);
  fTarget_log = 
    new G4LogicalVolume(fTarget_box, fBerillium, "target_log", 0, 0, 0);
  fTarget_phys = 
    new G4PVPlacement(0, G4ThreeVector(0., 0., 4.1875*m - 4.6*m + 4.0*m + 0.72*m - 6.0*m), fTarget_log, "target", fExperimentalHall_log, false, 0);

  // Place the main volume
  if(fMainVolumeOn)
  {
    fMainVolume = new LXeMainVolume(0, G4ThreeVector(0. *m, 0. *m , 4.1875*m - 4.6*m + 4.0*m + 1.6*m - 6.0*m), fExperimentalHall_log, false, 0, this);
  }

  fLXeAccelerator = static_cast<const LXeAccelerator*>(new LXeAccelerator(0, G4ThreeVector( 0. , 0. , -6.60*m + 4.0*m), fExperimentalHall_log, false, 0, this));

  //fLXeModule1 = static_cast<const LXeModule1*>(new LXeModule1(0, G4ThreeVector( 0. , 0. , 0.395*m - 5.60*m + 4.0*m), fExperimentalHall_log, false, 0, this));
  
  //fLXeModule2 = static_cast<const LXeModule2*>(new LXeModule2(0, G4ThreeVector( 0. , 0. , 0.395*m - 4.6*m + 4.0*m), fExperimentalHall_log, false, 0, this));
  
  //fLXeModule3 = static_cast<const LXeModule3*>(new LXeModule3(0, G4ThreeVector( 0. , 0. , 1.395*m - 4.6*m + 4.0*m), fExperimentalHall_log, false, 0, this));
  
  //fLXeModuleCCL1 = static_cast<const LXeModuleCCL*>(new LXeModuleCCL(0, G4ThreeVector( 0. , 0. , 2.2225*m - 4.6*m + 4.0*m), fExperimentalHall_log, false, 0, this));
  
  //fLXeModuleCCL2 = static_cast<const LXeModuleCCL*>(new LXeModuleCCL(0, G4ThreeVector( 0. , 0. , 2.8775*m - 4.6*m + 4.0*m), fExperimentalHall_log, false, 0, this));
  
  //fLXeModuleCCL3 = static_cast<const LXeModuleCCL*>(new LXeModuleCCL(0, G4ThreeVector( 0. , 0. , 3.5325*m - 4.6*m + 4.0*m), fExperimentalHall_log, false, 0, this));
  
  //fLXeModuleCCL4 = static_cast<const LXeModuleCCL*>(new LXeModuleCCL(0, G4ThreeVector( 0. , 0. , 4.1875*m - 4.6*m + 4.0*m), fExperimentalHall_log, false, 0, this));

  // Place the WLS slab
  if(fWLSslab)
  {
    G4VPhysicalVolume* slab = new LXeWLSSlab(
      0, G4ThreeVector(0., 0., -fScint_z / 2. - fSlab_z - 1. * cm),
      fExperimentalHall_log, false, 0, this);

    // Surface properties for the WLS slab
    G4OpticalSurface* scintWrap = new G4OpticalSurface("ScintWrap");

    new G4LogicalBorderSurface("ScintWrap", slab, fExperimentalHall_phys,
                               scintWrap);

    scintWrap->SetType(dielectric_metal);
    scintWrap->SetFinish(polished);
    scintWrap->SetModel(glisur);

    std::vector<G4double> pp           = { 2.0 * eV, 3.5 * eV };
    std::vector<G4double> reflectivity = { 1.0, 1.0 };
    std::vector<G4double> efficiency   = { 0.0, 0.0 };

    G4MaterialPropertiesTable* scintWrapProperty =
      new G4MaterialPropertiesTable();

    scintWrapProperty->AddProperty("REFLECTIVITY", pp, reflectivity);
    scintWrapProperty->AddProperty("EFFICIENCY", pp, efficiency);
    scintWrap->SetMaterialPropertiesTable(scintWrapProperty);
  }
  
  fScoringVolume = fExperimentalHall_log;

  return fWall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::ConstructSDandField()
{
  if(!fMainVolume)
    return;

  // PMT SD

  LXePMTSD* pmt = fPmt_SD.Get();
  if(!pmt)
  {
    // Created here so it exists as pmts are being placed
    G4cout << "Construction /LXeDet/pmtSD" << G4endl;
    LXePMTSD* pmt_SD = new LXePMTSD("/LXeDet/pmtSD");
    fPmt_SD.Put(pmt_SD);

    pmt_SD->InitPMTs();
    pmt_SD->SetPmtPositions(fMainVolume->GetPmtPositions());
  }
  else
  {
    pmt->InitPMTs();
    pmt->SetPmtPositions(fMainVolume->GetPmtPositions());
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(fPmt_SD.Get());
  // sensitive detector is not actually on the photocathode.
  // processHits gets done manually by the stepping action.
  // It is used to detect when photons hit and get absorbed & detected at the
  // boundary to the photocathode (which doesn't get done by attaching it to a
  // logical volume.
  // It does however need to be attached to something or else it doesn't get
  // reset at the begining of events

  SetSensitiveDetector(fMainVolume->GetLogPhotoCath(), fPmt_SD.Get());

  // Scint SD

  if(!fScint_SD.Get())
  {
    G4cout << "Construction /LXeDet/scintSD" << G4endl;
    LXeScintSD* scint_SD = new LXeScintSD("/LXeDet/scintSD");
    fScint_SD.Put(scint_SD);
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(fScint_SD.Get());
  SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());

  // LXeAccelerator
  if(!felFieldSetup1_0){
  static G4LogicalVolume* flogicShape6 = fLXeAccelerator->getlogicShape6();
  felFieldSetup1_0 = new ElectricFieldSetup();                               
  felFieldSetup1_0->SetLocalFieldValue(G4ThreeVector(0.,0.,(4.8*megavolt/m)));
  flogicShape6->SetFieldManager(felFieldSetup1_0->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup1_0); 
  }
  if(!felFieldSetup2_0){
  static G4LogicalVolume* flogicShape10 = fLXeAccelerator->getlogicShape10();
  felFieldSetup2_0 = new ElectricFieldSetup();                               
  felFieldSetup2_0->SetLocalFieldValue(G4ThreeVector(0.,0.,(4.8*megavolt/m)));
  flogicShape10->SetFieldManager(felFieldSetup2_0->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup2_0); 
  }
  if(!femFieldSetup1_0){
  static G4LogicalVolume* flogicShape1 = fLXeAccelerator->getlogicShape1();
  femFieldSetup1_0 = new F03FieldSetup();                               
  femFieldSetup1_0->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape1->SetFieldManager(femFieldSetup1_0->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup1_0);
  }
  if(!femFieldSetup2_0){
  static G4LogicalVolume* flogicShape4 = fLXeAccelerator->getlogicShape4();
  femFieldSetup2_0 = new F03FieldSetup();                               
  femFieldSetup2_0->SetLocalFieldValue(196.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape4->SetFieldManager(femFieldSetup2_0->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup2_0);
  }
  if(!femFieldSetup3_0){
  static G4LogicalVolume* flogicShape8 = fLXeAccelerator->getlogicShape8();
  femFieldSetup3_0 = new F03FieldSetup();                               
  femFieldSetup3_0->SetLocalFieldValue(179.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape8->SetFieldManager(femFieldSetup3_0->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup3_0);
  }
  if(!femFieldSetup4_0){
  static G4LogicalVolume* flogicShape12 = fLXeAccelerator->getlogicShape12();
  femFieldSetup4_0 = new F03FieldSetup();                               
  femFieldSetup4_0->SetLocalFieldValue(186.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape12->SetFieldManager(femFieldSetup4_0->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup4_0);
  }
  /*
  // LXeModule1
  if(!felFieldSetup1_1){
  static G4LogicalVolume* flogicShape6 = fLXeModule1->getlogicShape6();
  felFieldSetup1_1 = new ElectricFieldSetup();                               
  //felFieldSetup1_1->SetLocalFieldValue(G4ThreeVector(0.,0.,8.*megavolt/m));
  flogicShape6->SetFieldManager(felFieldSetup1_1->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup1_1); 
  }
  if(!felFieldSetup2_1){
  static G4LogicalVolume* flogicShape10 = fLXeModule1->getlogicShape10();
  felFieldSetup2_1 = new ElectricFieldSetup();                               
  //felFieldSetup2_1->SetLocalFieldValue(G4ThreeVector(0.,0.,8.*megavolt/m));
  flogicShape10->SetFieldManager(felFieldSetup2_1->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup2_1); 
  }
  if(!femFieldSetup1_1){
  static G4LogicalVolume* flogicShape8 = fLXeModule1->getlogicShape8();
  femFieldSetup1_1 = new F03FieldSetup();                               
  femFieldSetup1_1->SetLocalFieldValue(179.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape8->SetFieldManager(femFieldSetup1_1->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup1_1);
  }
  if(!femFieldSetup2_1){
  static G4LogicalVolume* flogicShape12 = fLXeModule1->getlogicShape12();
  femFieldSetup2_1 = new F03FieldSetup();                               
  femFieldSetup2_1->SetLocalFieldValue(186.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape12->SetFieldManager(femFieldSetup2_1->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup2_1);
  }
  
  // LXeModule2
  if(!felFieldSetup1_2){
  static G4LogicalVolume* flogicShape6 = fLXeModule2->getlogicShape6();
  felFieldSetup1_2 = new ElectricFieldSetup();                               
  //felFieldSetup1_2->SetLocalFieldValue(G4ThreeVector(0.,0.,11.*megavolt/m));
  flogicShape6->SetFieldManager(felFieldSetup1_2->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup1_2); 
  }
  if(!felFieldSetup2_2){
  static G4LogicalVolume* flogicShape10 = fLXeModule2->getlogicShape10();
  felFieldSetup2_2 = new ElectricFieldSetup();                               
  //felFieldSetup2_2->SetLocalFieldValue(G4ThreeVector(0.,0.,11.*megavolt/m));
  flogicShape10->SetFieldManager(felFieldSetup2_2->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup2_2); 
  }
  if(!femFieldSetup1_2){
  static G4LogicalVolume* flogicShape8 = fLXeModule2->getlogicShape8();
  femFieldSetup1_2 = new F03FieldSetup();                               
  femFieldSetup1_2->SetLocalFieldValue(179.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape8->SetFieldManager(femFieldSetup1_2->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup1_2);
  }
  if(!femFieldSetup2_2){
  static G4LogicalVolume* flogicShape12 = fLXeModule2->getlogicShape12();
  femFieldSetup2_2 = new F03FieldSetup();                               
  femFieldSetup2_2->SetLocalFieldValue(186.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape12->SetFieldManager(femFieldSetup2_2->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup2_2);
  }
  
  // LXeModule3
  if(!felFieldSetup1_3){
  static G4LogicalVolume* flogicShape6 = fLXeModule3->getlogicShape6();
  felFieldSetup1_3 = new ElectricFieldSetup();                               
  //felFieldSetup1_3->SetLocalFieldValue(G4ThreeVector(0.,0.,13.4*megavolt/m));
  flogicShape6->SetFieldManager(felFieldSetup1_3->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup1_3); 
  }
  if(!felFieldSetup2_3){
  static G4LogicalVolume* flogicShape10 = fLXeModule3->getlogicShape10();
  felFieldSetup2_3 = new ElectricFieldSetup();                               
  //felFieldSetup2_3->SetLocalFieldValue(G4ThreeVector(0.,0.,13.4*megavolt/m));
  flogicShape10->SetFieldManager(felFieldSetup2_3->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup2_3); 
  }
  if(!femFieldSetup1_3){
  static G4LogicalVolume* flogicShape8 = fLXeModule3->getlogicShape8();
  femFieldSetup1_3 = new F03FieldSetup();                               
  femFieldSetup1_3->SetLocalFieldValue(179.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape8->SetFieldManager(femFieldSetup1_3->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup1_3);
  }
  if(!femFieldSetup2_3){
  static G4LogicalVolume* flogicShape12 = fLXeModule3->getlogicShape12();
  femFieldSetup2_3 = new F03FieldSetup();                               
  femFieldSetup2_3->SetLocalFieldValue(186.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape12->SetFieldManager(femFieldSetup2_3->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup2_3);
  }
  
  // LXeModuleCCL1
  if(!felFieldSetup1_4){
  static G4LogicalVolume* flogicShape1 = fLXeModuleCCL1->getlogicShape1();
  felFieldSetup1_4 = new ElectricFieldSetup();                               
  //felFieldSetup1_4->SetLocalFieldValue(G4ThreeVector(0.,0.,12.*megavolt/m));
  flogicShape1->SetFieldManager(felFieldSetup1_4->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup1_4); 
  }
  if(!felFieldSetup2_4){
  static G4LogicalVolume* flogicShape4 = fLXeModuleCCL1->getlogicShape4();
  felFieldSetup2_4 = new ElectricFieldSetup();                               
  //felFieldSetup2_4->SetLocalFieldValue(G4ThreeVector(0.,0.,12.*megavolt/m));
  flogicShape4->SetFieldManager(felFieldSetup2_4->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup2_4); 
  }
  if(!femFieldSetup1_4){
  static G4LogicalVolume* flogicShape3 = fLXeModuleCCL1->getlogicShape3();
  femFieldSetup1_4 = new F03FieldSetup();                               
  femFieldSetup1_4->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape3->SetFieldManager(femFieldSetup1_4->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup1_4);
  }
  if(!femFieldSetup2_4){
  static G4LogicalVolume* flogicShape6 = fLXeModuleCCL1->getlogicShape6();
  femFieldSetup2_4 = new F03FieldSetup();                               
  femFieldSetup2_4->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape6->SetFieldManager(femFieldSetup2_4->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup2_4);
  }
  
  // LXeModuleCCL2
  if(!felFieldSetup1_5){
  static G4LogicalVolume* flogicShape1 = fLXeModuleCCL2->getlogicShape1();
  felFieldSetup1_5 = new ElectricFieldSetup();                               
  //felFieldSetup1_5->SetLocalFieldValue(G4ThreeVector(0.,0.,12.*megavolt/m));
  flogicShape1->SetFieldManager(felFieldSetup1_5->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup1_5); 
  }
  if(!felFieldSetup2_5){
  static G4LogicalVolume* flogicShape4 = fLXeModuleCCL2->getlogicShape4();
  felFieldSetup2_5 = new ElectricFieldSetup();                               
  //felFieldSetup2_5->SetLocalFieldValue(G4ThreeVector(0.,0.,12.*megavolt/m));
  flogicShape4->SetFieldManager(felFieldSetup2_5->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup2_5); 
  }
  if(!femFieldSetup1_5){
  static G4LogicalVolume* flogicShape3 = fLXeModuleCCL2->getlogicShape3();
  femFieldSetup1_5 = new F03FieldSetup();                               
  femFieldSetup1_5->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape3->SetFieldManager(femFieldSetup1_5->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup1_5);
  }
  if(!femFieldSetup2_5){
  static G4LogicalVolume* flogicShape6 = fLXeModuleCCL2->getlogicShape6();
  femFieldSetup2_5 = new F03FieldSetup();                               
  femFieldSetup2_5->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape6->SetFieldManager(femFieldSetup2_5->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup2_5);
  }
  
  // LXeModuleCCL3
  if(!felFieldSetup1_6){
  static G4LogicalVolume* flogicShape1 = fLXeModuleCCL3->getlogicShape1();
  felFieldSetup1_6 = new ElectricFieldSetup();                               
  //felFieldSetup1_6->SetLocalFieldValue(G4ThreeVector(0.,0.,12.*megavolt/m));
  flogicShape1->SetFieldManager(felFieldSetup1_6->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup1_6); 
  }
  if(!felFieldSetup2_6){
  static G4LogicalVolume* flogicShape4 = fLXeModuleCCL3->getlogicShape4();
  felFieldSetup2_6 = new ElectricFieldSetup();                               
  //felFieldSetup2_6->SetLocalFieldValue(G4ThreeVector(0.,0.,12.*megavolt/m));
  flogicShape4->SetFieldManager(felFieldSetup2_6->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup2_6); 
  }
  if(!femFieldSetup1_6){
  static G4LogicalVolume* flogicShape3 = fLXeModuleCCL3->getlogicShape3();
  femFieldSetup1_6 = new F03FieldSetup();                               
  femFieldSetup1_6->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape3->SetFieldManager(femFieldSetup1_6->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup1_6);
  }
  if(!femFieldSetup2_6){
  static G4LogicalVolume* flogicShape6 = fLXeModuleCCL3->getlogicShape6();
  femFieldSetup2_6 = new F03FieldSetup();                               
  femFieldSetup2_6->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape6->SetFieldManager(femFieldSetup2_6->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup2_6);
  }
  
  // LXeModuleCCL4
  if(!felFieldSetup1_7){
  static G4LogicalVolume* flogicShape1 = fLXeModuleCCL4->getlogicShape1();
  felFieldSetup1_7 = new ElectricFieldSetup();                               
  //felFieldSetup1_7->SetLocalFieldValue(G4ThreeVector(0.,0.,12.*megavolt/m));
  flogicShape1->SetFieldManager(felFieldSetup1_7->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup1_7); 
  }
  if(!felFieldSetup2_7){
  static G4LogicalVolume* flogicShape4 = fLXeModuleCCL4->getlogicShape4();
  felFieldSetup2_7 = new ElectricFieldSetup();                               
  //felFieldSetup2_7->SetLocalFieldValue(G4ThreeVector(0.,0.,12.*megavolt/m));
  flogicShape4->SetFieldManager(felFieldSetup2_7->GetLocalFieldManager(), true );
  G4AutoDelete::Register(felFieldSetup2_7); 
  }
  if(!femFieldSetup1_7){
  static G4LogicalVolume* flogicShape3 = fLXeModuleCCL4->getlogicShape3();
  femFieldSetup1_7 = new F03FieldSetup();                               
  femFieldSetup1_7->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 0.*deg);
  flogicShape3->SetFieldManager(femFieldSetup1_7->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup1_7);
  }
  if(!femFieldSetup2_7){
  static G4LogicalVolume* flogicShape6 = fLXeModuleCCL4->getlogicShape6();
  femFieldSetup2_7 = new F03FieldSetup();                               
  femFieldSetup2_7->SetLocalFieldValue(160.*tesla/(1.*m), G4ThreeVector() , 90.*deg);
  flogicShape6->SetFieldManager(femFieldSetup2_7->GetLocalFieldManager(), true );
  G4AutoDelete::Register(femFieldSetup2_7);
  }
  */
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDimensions(G4ThreeVector dims)
{
  fScint_x = dims[0];
  fScint_y = dims[1];
  fScint_z = dims[2];
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingThickness(G4double d_mtl)
{
  fD_mtl = d_mtl;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNX(G4int nx)
{
  fNx = nx;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNY(G4int ny)
{
  fNy = ny;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNZ(G4int nz)
{
  fNz = nz;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetPMTRadius(G4double outerRadius_pmt)
{
  fOuterRadius_pmt = outerRadius_pmt;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDefaults()
{
  // Resets to default values
  fD_mtl = 0.0635 * cm;
  fD_moderator = 3. * cm;

  fScint_x = 6. * cm;
  fScint_y = 6. * cm;
  fScint_z = 20. * cm;

  fNx = 1;
  fNy = 1;
  fNz = 0;

  fOuterRadius_pmt = 2.3 * cm;

  fSphereOn = false;
  fRefl     = 1.0;

  fNfibers      = 1;
  fWLSslab      = false;
  fMainVolumeOn = true;
  fMainVolume   = nullptr;
  fSlab_z       = 2.5 * mm;

  G4UImanager::GetUIpointer()->ApplyCommand(
    "/LXe/detector/scintYieldFactor 1.");

  if(fLXe_mt)
    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
  if(fMPTPStyrene)
    fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetSphereOn(G4bool b)
{
  fSphereOn = b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingReflectivity(G4double r)
{
  fRefl = r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetWLSSlabOn(G4bool b)
{
  fWLSslab = b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainVolumeOn(G4bool b)
{
  fMainVolumeOn = b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNFibers(G4int n)
{
  fNfibers = n;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainScintYield(G4double y)
{
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", y / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetWLSScintYield(G4double y)
{
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", y / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetSaveThreshold(G4int save)
{
  // Sets the save threshold for the random number seed. If the number of
  // photons generated in an event is lower than this, then save the seed for
  // this event in a file called run###evt###.rndm

  fSaveThreshold = save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
