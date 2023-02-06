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
/// \file field/field03/src/F03FieldSetup.cc
/// \brief Implementation of the F03FieldSetup class
//
//
//
//
//   Field Setup class implementation.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ElectricFieldSetup.hh"
#include "ElectricFieldMessenger.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectricFieldSetup::ElectricFieldSetup()
 : fFieldManager(0),
   fLocalFieldManager(0),
   fChordFinder(0),
   fLocalChordFinder(0),
   fEquation(0),
   fLocalEquation(0),
   fElectricField(0),
   fLocalElectricField(0),
   fStepper(0),
   fLocalStepper(0),
   fIntgrDriver(0),
   fLocalIntgrDriver(0),
   fFieldMessenger(0)
{
  fElectricField = new G4UniformElectricField(G4ThreeVector(0.,
                                                            0.,
                                                            0.*megavolt/m));
  fLocalElectricField = new G4UniformElectricField(G4ThreeVector(0.,
                                                                 0.,
                                                                 0.*megavolt/m));

  fFieldMessenger = new ElectricFieldMessenger(this);
 
  fEquation = new G4EqMagElectricField(fElectricField);
  fLocalEquation = new G4EqMagElectricField(fLocalElectricField);
 
  fMinStep     = 0.0000001*mm ; // minimal step of 1 mm is default
  fStepperType = 5 ;       // ClassicalRK4 is default stepper
  
  fmissDistance       = 0.00001*mm;
  fdeltaIntersection  = 0.00001*mm;
  fdeltaOneStep       = 0.00001*mm;

  fFieldManager = GetGlobalFieldManager();
  fLocalFieldManager = new G4FieldManager();

  UpdateField();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectricFieldSetup::~ElectricFieldSetup()
{
  delete fChordFinder;
  delete fLocalChordFinder;
  delete fStepper;
  delete fLocalStepper;
  delete fEquation;
  delete fLocalEquation;
  delete fElectricField;
  delete fLocalElectricField;
  delete fFieldMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectricFieldSetup::UpdateField()
{
  // It must be possible to call 'again' - e.g. to choose an alternative stepper
  //   has been chosen, or in case other changes have been made.

  // 1. First clean up previous state.
  delete fChordFinder;
  fChordFinder = nullptr;
  delete fLocalChordFinder;
  fLocalChordFinder = nullptr;
  
  fIntgrDriver = nullptr;
  fLocalIntgrDriver = nullptr;

  G4cout<<"ElectricFieldSetup::UpdateField> The minimal step is equal to "
        << fMinStep/mm <<" mm"<<G4endl;
  G4cout<<"                            Stepper Type chosen = " << fStepperType
        << G4endl;

  // 2. Create the steppers ( Note: this also deletes the previous ones. )
  CreateSteppers();
  
  assert(fStepper != nullptr);
  
  // 3. Create the chord finder(s)
  if( fStepper ) {
     fIntgrDriver = new G4MagInt_Driver(fMinStep,
                                        fStepper,
                                        fStepper->GetNumberOfVariables());
     if( fIntgrDriver ){ 
        fChordFinder = new G4ChordFinder(fIntgrDriver);
     }
  }
  
  assert(fLocalStepper != nullptr);
  
  // 4. Create the local chord finder(s)
  if( fLocalStepper ) {
     fLocalIntgrDriver = new G4MagInt_Driver(fMinStep,
                                             fLocalStepper,
                                             fLocalStepper->GetNumberOfVariables());
     if( fLocalIntgrDriver ){ 
        fLocalChordFinder = new G4ChordFinder(fLocalIntgrDriver);
     }
  }
  
  fLocalChordFinder->SetDeltaChord( fmissDistance );
  fLocalFieldManager->SetDeltaIntersection( fdeltaIntersection );
  fLocalFieldManager->SetDeltaOneStep( fdeltaOneStep );

  // 5. Set Chords
  fFieldManager->SetChordFinder(fChordFinder);
  fLocalFieldManager->SetChordFinder(fLocalChordFinder);

  // 6. Ensure that the field is updated (in Field manager & equation)
  fFieldManager->SetDetectorField(fElectricField);
  fLocalFieldManager->SetDetectorField(fLocalElectricField);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectricFieldSetup::CreateSteppers()
{
  delete fStepper;
  fStepper= nullptr;

  delete fLocalStepper; 
  fLocalStepper= nullptr;
  
  const G4int nvar = 8;

  switch ( fStepperType )
  {
    case 0:
      fStepper = new G4ExplicitEuler( fEquation , nvar );
      fLocalStepper = new G4ExplicitEuler( fLocalEquation , nvar );
      G4cout<<"G4ExplicitEuler is called"<<G4endl;
      break;
    case 1:
      fStepper = new G4ImplicitEuler( fEquation , nvar );
      fLocalStepper = new G4ImplicitEuler( fLocalEquation , nvar );
      G4cout<<"G4ImplicitEuler is called"<<G4endl;
      break;
    case 2:
      fStepper = new G4SimpleRunge( fEquation , nvar );
      fLocalStepper = new G4SimpleRunge( fLocalEquation , nvar );
      G4cout<<"G4SimpleRunge is called"<<G4endl;
      break;
    case 3:
      fStepper = new G4SimpleHeum( fEquation , nvar );
      fLocalStepper = new G4SimpleHeum( fLocalEquation , nvar );
      G4cout<<"G4SimpleHeum is called"<<G4endl;
      break;
    case 4:
      fStepper = new G4ClassicalRK4( fEquation , nvar);
      fLocalStepper = new G4ClassicalRK4( fLocalEquation , nvar );
      G4cout<<"G4ClassicalRK4 is called"<<G4endl;
      break;    
    case 5:
      fStepper = new G4CashKarpRKF45( fEquation , nvar );
      fLocalStepper = new G4CashKarpRKF45( fLocalEquation , nvar );
      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
      break;
    default:
      fStepper = new G4ClassicalRK4( fEquation , nvar );
      fLocalStepper = new G4ClassicalRK4( fLocalEquation , nvar );
      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;
      break;
  }
  
  if( fIntgrDriver )
      fIntgrDriver->RenewStepperAndAdjust( fStepper );
      
  if( fLocalIntgrDriver )
      fLocalIntgrDriver->RenewStepperAndAdjust( fLocalStepper );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectricFieldSetup::SetFieldZValue(G4double fieldStrength)
{
  G4ThreeVector fieldSetVec(0.0, 0.0, fieldStrength);
  SetFieldValue( fieldSetVec );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectricFieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  if(fElectricField) delete fElectricField;

  if(fieldVector != G4ThreeVector(0.,0.,0.))
  {
    fElectricField = new  G4UniformElectricField(fieldVector);
  }
  else
  {
    // If the new field's value is Zero, then
    // setting the pointer to zero ensures
    // that it is not used for propagation.
    fElectricField = 0;
  }

  // Either
  //   - UpdateField() to reset all (ChordFinder, Equation);
  // UpdateField();
  //     or simply update the field manager & equation of motion
  //     with pointer to new field
  GetGlobalFieldManager()->SetDetectorField(fElectricField);
  fEquation->SetFieldObj( fElectricField );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectricFieldSetup::SetLocalFieldValue(G4ThreeVector fieldVector)
{
  if(fLocalElectricField) delete fLocalElectricField;

  if(fieldVector != G4ThreeVector(0.,0.,0.))
  {
    fLocalElectricField = new G4UniformElectricField(fieldVector);
  }
  else
  {
    // If the new field's value is Zero, then
    // setting the pointer to zero ensures
    // that it is not used for propagation.
    fLocalElectricField = 0;
  }

  // Either
  //   - UpdateField() to reset all (ChordFinder, Equation);
  //UpdateField();
  //     or simply update the field manager & equation of motion
  //     with pointer to new field
  GetLocalFieldManager()->SetDetectorField(fLocalElectricField);
  fEquation->SetFieldObj( fLocalElectricField );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FieldManager* ElectricFieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
                                  ->GetFieldManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ElectricFieldSetup::GetConstantFieldValue(G4ElectricField* ElectricField) const
{
  if ( ! ElectricField ) return G4ThreeVector();

  static G4double fieldValue[6],  position[4];
  position[0] = position[1] = position[2] = position[3] = 0.0;

  ElectricField->GetFieldValue( position, fieldValue);
  G4ThreeVector fieldVec(fieldValue[0], fieldValue[1], fieldValue[2]);

  return fieldVec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
