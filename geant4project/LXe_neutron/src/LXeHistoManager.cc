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
/// \file optical/LXe/src/LXeHistoManager.cc
/// \brief Implementation of the LXeHistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LXeHistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeHistoManager::LXeHistoManager()
  : fFileName("lxe")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeHistoManager::~LXeHistoManager() { delete G4AnalysisManager::Instance(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeHistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in LXeHistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);  // enable inactivation of histograms

  // Define histogram indices, titles

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins   = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // 0
  analysisManager->CreateH1("E kin neutrons in Module0", "Neutron spectrum in Module0", nbins, vmin, 10.);
  // 1
  analysisManager->CreateH1("hits per event", "hits per event", nbins, vmin, 10000.);
  // 2
  analysisManager->CreateH1("hits above threshold", "hits per event above threshold", nbins, vmin, 2.);
  // 3
  analysisManager->CreateH1("scintillation", "scintillation photons per event", nbins, vmin, 15000.);
  // 4
  analysisManager->CreateH1("Cerenkov", "Cerenkov photons per event", nbins, vmin, vmax);
  // 5
  analysisManager->CreateH1("absorbed", "absorbed photons per event", nbins, vmin, 15000.);
  // 6
  analysisManager->CreateH1("boundary absorbed", "photons absorbed at boundary per event", nbins, vmin, vmax);
  // 7
  analysisManager->CreateH1("E dep", "energy deposition in scintillator per event", nbins, vmin, 3.);  
  // 8
  analysisManager->CreateH1("E kin protons in cavities", "Proton spectrum in cavities - Kinetic Energy ", nbins, vmin, 10.);
  // 9
  analysisManager->CreateH1("E kin neutrons in moderator", "Neutron spectrum in moderator - Kinetic Energy ", nbins, vmin, 1.);
  // 10
  analysisManager->CreateH1("E kin neutrons in exp Hall", "Neutron spectrum in exp Hall - Kinetic Energy ", 200, vmin, 10.);
  // 11
  analysisManager->CreateH1("Z position of protons in target volume", "Proton position in target - Z Pos", 100, -170.30, -168.30);
  // 12
  analysisManager->CreateH1("Angle distribution of neutron from target", "Angle distribution of neutron from target", 60, -180., +180.);

  // Create all histograms as inactivated
  for(G4int i = 0; i < analysisManager->GetNofH1s(); ++i)
  {
    analysisManager->SetH1Activation(i, true);
  }
}












