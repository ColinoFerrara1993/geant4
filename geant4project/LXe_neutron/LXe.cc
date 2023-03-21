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
/// \file optical/LXe/LXe.cc
/// \brief Main program of the optical/LXe example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LXeActionInitialization.hh"
#include "LXeDetectorConstruction.hh"
#include "LXeSteppingAction.hh"

#include "QGS_BIC.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC_AllHP.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4MTRunManager.hh"

#include "ElectricFieldSetup.hh"
#include "F03FieldSetup.hh"
#include "G4Timer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup1_0;
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup2_0;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup1_0;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup2_0;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup3_0;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup4_0;

G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup1_1;
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup2_1;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup1_1;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup2_1;

G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup1_2;
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup2_2;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup1_2;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup2_2;

G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup1_3;
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup2_3;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup1_3;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup2_3;

G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup1_4;
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup2_4;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup1_4;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup2_4;

G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup1_5;
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup2_5;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup1_5;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup2_5;

G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup1_6;
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup2_6;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup1_6;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup2_6;

G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup1_7;
G4ThreadLocal ElectricFieldSetup* LXeDetectorConstruction::felFieldSetup2_7;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup1_7;
G4ThreadLocal F03FieldSetup* LXeDetectorConstruction::femFieldSetup2_7;

G4ThreadLocal G4double LXeSteppingAction::ftime;

int main(int argc, char** argv)
{
  // detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if(argc == 1)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  //auto runManager = G4RunManagerFactory::CreateRunManager();
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(16);
  
  LXeDetectorConstruction* det = new LXeDetectorConstruction();
  runManager->SetUserInitialization(det);

  G4VModularPhysicsList* physicsList = new QGSP_BIC_AllHP;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  auto opticalParams               = G4OpticalParameters::Instance();

  opticalParams->SetWLSTimeProfile("delta");

  opticalParams->SetScintYieldFactor(1.0);
  opticalParams->SetScintExcitationRatio(0.0);
  opticalParams->SetScintTrackSecondariesFirst(true);
  opticalParams->SetScintEnhancedTimeConstants(true);

  opticalParams->SetCerenkovMaxPhotonsPerStep(100);
  opticalParams->SetCerenkovMaxBetaChange(10.0);
  opticalParams->SetCerenkovTrackSecondariesFirst(true);

  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new LXeActionInitialization(det));

  // initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(ui)
  {
    // interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    if(ui->IsGUI())
    {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }
  else
  {
    // batch mode
    G4String command  = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  // job termination
  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
