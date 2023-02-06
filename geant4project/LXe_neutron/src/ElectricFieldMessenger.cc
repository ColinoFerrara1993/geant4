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
/// \file field/field03/src/F03FieldMessenger.cc
/// \brief Implementation of the F03FieldMessenger class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ElectricFieldMessenger.hh"

#include "ElectricFieldSetup.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectricFieldMessenger::ElectricFieldMessenger(ElectricFieldSetup* fieldSetup)
 : G4UImessenger(),
   fEMfieldSetup(fieldSetup),
   fFieldDir(0),
   fStepperCmd(0),
   fEleFieldZCmd(0),
   fEleFieldCmd(0),
   fLocalEleFieldCmd(0),
   fMinStepCmd(0),
   fUpdateCmd(0)
{
  fFieldDir = new G4UIdirectory("/electricfield/");
  fFieldDir->SetGuidance("Electric field tracking control.");

  fStepperCmd = new G4UIcmdWithAnInteger("/electricfield/setStepperType",this);
  fStepperCmd->SetGuidance("Select stepper type for electric field");
  fStepperCmd->SetParameterName("choice",true);
  fStepperCmd->SetDefaultValue(4);
  fStepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/electricfield/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
 
  fEleFieldZCmd = new G4UIcmdWithADoubleAndUnit("/electricfield/setFieldZ",this);
  fEleFieldZCmd->SetGuidance("Define global electric field.");
  fEleFieldZCmd->SetGuidance("Global electric field will be in Z direction.");
  fEleFieldZCmd->SetParameterName("Ez",false,false);
  fEleFieldZCmd->SetDefaultUnit("megavolt/m");
  fEleFieldZCmd->AvailableForStates(G4State_Idle);
 
  fEleFieldCmd = new G4UIcmdWith3VectorAndUnit("/electricfield/setField",this);
  fEleFieldCmd->SetGuidance("Define global electric field.");
  fEleFieldCmd->SetParameterName("Ex","Ey","Ez",false,false);
  fEleFieldCmd->SetDefaultUnit("megavolt/m");
  fEleFieldCmd->AvailableForStates(G4State_Idle);
 
  fLocalEleFieldCmd = new G4UIcmdWith3VectorAndUnit("/electricfield/setLocalField",this);
  fLocalEleFieldCmd->SetGuidance("Define local electric field.");
  fLocalEleFieldCmd->SetParameterName("Elx","Ely","Elz",false,false);
  fLocalEleFieldCmd->SetDefaultUnit("megavolt/m");
  fLocalEleFieldCmd->AvailableForStates(G4State_Idle);
 
  fMinStepCmd = new G4UIcmdWithADoubleAndUnit("/electricfield/setMinStep",this);
  fMinStepCmd->SetGuidance("Define minimal step");
  fMinStepCmd->SetGuidance("Electric field will be in Z direction.");
  fMinStepCmd->SetParameterName("min step",false,false);
  fMinStepCmd->SetDefaultUnit("mm");
  fMinStepCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectricFieldMessenger::~ElectricFieldMessenger()
{
  delete fStepperCmd;
  delete fEleFieldZCmd;
  delete fEleFieldCmd;
  delete fLocalEleFieldCmd;
  delete fMinStepCmd;
  delete fFieldDir;
  delete fUpdateCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectricFieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{
  if( command == fStepperCmd )
    fEMfieldSetup->SetStepperType(fStepperCmd->GetNewIntValue(newValue));
  if( command == fUpdateCmd )
    fEMfieldSetup->UpdateField();
  if( command == fEleFieldZCmd )
  {
    fEMfieldSetup->SetFieldZValue(fEleFieldZCmd->GetNewDoubleValue(newValue));
    // Check the value
    G4cout << "Set global field value to " <<
      fEMfieldSetup->GetGlobalFieldValue() / (megavolt/m) << " MegaVolt/m " << G4endl;
  }
  if( command == fEleFieldCmd )
  {
    fEMfieldSetup->SetFieldValue(fEleFieldCmd->GetNew3VectorValue(newValue));
    // Check the value
    G4cout << "Set global field value to " <<
      fEMfieldSetup->GetGlobalFieldValue() / (megavolt/m) << " MegaVolt/m " << G4endl;
  }
  if( command == fLocalEleFieldCmd )
  {
    fEMfieldSetup->SetLocalFieldValue(fLocalEleFieldCmd->GetNew3VectorValue(newValue));
    fEMfieldSetup->UpdateField();
    // Check the value
    G4cout << "Set local field value to " <<
      fEMfieldSetup->GetLocalFieldValue() / (megavolt/m) << " MegaVolt/m " << G4endl;
  }
  if( command == fMinStepCmd )
    fEMfieldSetup->SetMinStep(fMinStepCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
