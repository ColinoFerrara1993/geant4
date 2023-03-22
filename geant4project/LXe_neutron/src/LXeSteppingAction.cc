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
/// \file optical/LXe/src/LXeSteppingAction.cc
/// \brief Implementation of the LXeSteppingAction class
//
//
#include "LXeSteppingAction.hh"
#include "LXeDetectorConstruction.hh"

#include "LXeEventAction.hh"
#include "LXePMTSD.hh"
#include "LXeSteppingMessenger.hh"
#include "LXeTrajectory.hh"
#include "LXeUserTrackInformation.hh"

#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

#include "g4root.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSteppingAction::LXeSteppingAction(LXeEventAction* ea)
  : fOneStepPrimaries(false)
  , fEventAction(ea)
  , fScoringVolume(0)
  , fkinEnergy(0.)
{
  fSteppingMessenger = new LXeSteppingMessenger(this);

  fExpectedNextStatus = Undefined;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSteppingAction::~LXeSteppingAction() { delete fSteppingMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  //measure energy deposition in scoring volume
  if (!fScoringVolume) { 
    const LXeDetectorConstruction* detectorConstruction
      = static_cast<const LXeDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = theStep->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  G4StepPoint* prePoint    = theStep->GetPreStepPoint();
  G4VPhysicalVolume* prePV = prePoint->GetPhysicalVolume();
  
  G4StepPoint* postPoint  = theStep->GetPostStepPoint();
  G4VPhysicalVolume* postPV = postPoint->GetPhysicalVolume();
  
  // check if we are in scoring volume
  if (prePV->GetName() == "housing1") {

    // collect energy deposited in this step
    G4double edepStep = theStep->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edepStep);
  
    //Retrieve from step the track
    G4Track* track1 = theStep->GetTrack();
  
    //Pointer to dynamic particle
    const G4DynamicParticle* dynamicParticle1 = track1->GetDynamicParticle();
  
    //Get Kinetic Energy after the step
    G4double kinEnergy1 = dynamicParticle1->GetKineticEnergy();

    //Retrive particle definition
    G4ParticleDefinition* particle1 = dynamicParticle1->GetDefinition();
  
    //Get name of particle
    G4String particleName1 = particle1->GetParticleName();
  
    if (particleName1 == "proton") {
    
      if (track1->GetParentID()==0) {
  
        //Analysis
        auto analysisManager = G4AnalysisManager::Instance();
  
        //fill
        analysisManager->FillH1(8, kinEnergy1);
      
      }
      
    } else if (particleName1 == "neutron") {
    
             if (postPV->GetName() == "housing1") {
             
              //Analysis
              auto analysisManager = G4AnalysisManager::Instance();
  
              //fill
              analysisManager->FillH1(9, kinEnergy1);
              
             }
        
           } 
  
  } else if (prePV->GetName() == "housing2"){
      
     //Retrieve from step the track
     G4Track* track2 = theStep->GetTrack();
  
     //Pointer to dynamic particle
     const G4DynamicParticle* dynamicParticle2 = track2->GetDynamicParticle();
  
     //Get Kinetic Energy after the step
     fkinEnergy = dynamicParticle2->GetKineticEnergy();
     
     //Retrive particle definition
     G4ParticleDefinition* particle2 = dynamicParticle2->GetDefinition();
  
     //Get name of particle
     G4String particleName2 = particle2->GetParticleName();
  
     if (particleName2 == "neutron") {
     
       if (track2->GetParentID() > 0) {
  
         //Analysis
         auto analysisManager = G4AnalysisManager::Instance();
  
         //fill
         analysisManager->FillH1(0, fkinEnergy);
       
       }
        
     }
  
    } else if (prePV->GetName() == "target"){
    
        //Retrieve from step the track
        G4Track* track3 = theStep->GetTrack();
  
        //Pointer to dynamic particle
        const G4DynamicParticle* dynamicParticle3 = track3->GetDynamicParticle();
  
        //Get Kinetic Energy after the step
        G4double kinEnergy3 = dynamicParticle3->GetKineticEnergy();
     
        //Retrive particle definition
        G4ParticleDefinition* particle3 = dynamicParticle3->GetDefinition();
  
        //Get name of particle
        G4String particleName3 = particle3->GetParticleName();
  
        if (particleName3 == "neutron") {
     
              if (postPV->GetName() == "expHall") {
              
                //Get angle
                G4double _x = track3->GetMomentumDirection().x();
                G4double _y = track3->GetMomentumDirection().y();
                G4double _z = track3->GetMomentumDirection().z();
                G4double ang = std::atan(std::sqrt((_x)*(_x) + (_y)*(_y))/(_z));
  
                //Analysis
                auto analysisManager = G4AnalysisManager::Instance();
  
                //fill
                analysisManager->FillH1(10, kinEnergy3);
                
                analysisManager->FillH1(12, ang*180./(pi));
       
              }
        
        } else if (particleName3 == "proton") {

                 // Analysis
                 auto analysisManager = G4AnalysisManager::Instance();

                 // Position
                 G4double z_p = track3->GetPosition().z();

                 // Fill
                 analysisManager->FillH1(11, z_p);
                
               }
      } else {
    
            //Retrieve from step the track
            G4Track* track4 = theStep->GetTrack();
  
            //Pointer to dynamic particle
            const G4DynamicParticle* dynamicParticle4 = track4->GetDynamicParticle();
  
            //Get Kinetic Energy after the step
            G4double kinEnergy4 = dynamicParticle4->GetKineticEnergy();
     
            //Retrive particle definition
            G4ParticleDefinition* particle4 = dynamicParticle4->GetDefinition();
            
            // Get mass
            //G4double mass = particle4->GetPDGMass();
  
            //Get name of particle
            G4String particleName4 = particle4->GetParticleName();
  
            if (particleName4 == "proton") {
     
              if (track4->GetParentID()==0) {
              
                if (kinEnergy4 > 3.) {
                
                  //ftime = track4->GetGlobalTime();
                  //G4double fmomentum_z = track4->GetMomentum().z();
                  //G4double fvelocity_z = fmomentum_z/mass;
                  //LXeDetectorConstruction::felFieldSetup1_0->SetLocalFieldValue(G4ThreeVector(0.,0.,(+24.0*megavolt/m)*sin(2.*3.14*(fvelocity_z/(2*0.06*m))*ftime + 40.*deg)));
                  //LXeDetectorConstruction::felFieldSetup2_0->SetLocalFieldValue(G4ThreeVector(0.,0.,(+24.0*megavolt/m)*sin(2.*3.14*(fvelocity_z/(2*0.06*m))*ftime + 40.*deg)));
                  //LXeDetectorConstruction::felFieldSetup1_0->UpdateField();
                  //LXeDetectorConstruction::felFieldSetup2_0->UpdateField();
  
                  //if (kinEnergy4 > 3. && kinEnergy4 < 14.){
                  
                    //G4double frequency = fvelocity_z/(2*0.06*m);
                    
                    //Analysis
                    auto analysisManager = G4AnalysisManager::Instance();
  
                    //fill
                    //analysisManager->FillH1(11, kinEnergy4);
                    
                  //}
                
                }
       
              }
        
            }    
    
          }
          
  G4Track* theTrack = theStep->GetTrack();

  if(theTrack->GetCurrentStepNumber() == 1)
    fExpectedNextStatus = Undefined;

  LXeUserTrackInformation* trackInformation =
    (LXeUserTrackInformation*) theTrack->GetUserInformation();

  G4StepPoint* thePrePoint    = theStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

  G4StepPoint* thePostPoint    = theStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  G4OpBoundaryProcessStatus boundaryStatus           = Undefined;
  static G4ThreadLocal G4OpBoundaryProcess* boundary = nullptr;

  // find the boundary process only once
  if(!boundary)
  {
    G4ProcessManager* pm =
      theStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses    = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    for(G4int i = 0; i < nprocesses; ++i)
    {
      if((*pv)[i]->GetProcessName() == "OpBoundary")
      {
        boundary = (G4OpBoundaryProcess*) (*pv)[i];
        break;
      }
    }
  }

  if(theTrack->GetParentID() == 0)
  {
    // This is a primary track

    G4TrackVector* fSecondary = fpSteppingManager->GetfSecondary();
    G4int tN2ndariesTot       = fpSteppingManager->GetfN2ndariesAtRestDoIt() +
                          fpSteppingManager->GetfN2ndariesAlongStepDoIt() +
                          fpSteppingManager->GetfN2ndariesPostStepDoIt();

    // If we haven't already found the conversion position and there were
    // secondaries generated, then search for it
    if(!fEventAction->IsConvPosSet() && tN2ndariesTot > 0)
    {
      for(size_t lp1 = (*fSecondary).size() - tN2ndariesTot;
          lp1 < (*fSecondary).size(); ++lp1)
      {
        const G4VProcess* creator = (*fSecondary)[lp1]->GetCreatorProcess();
        if(creator)
        {
          G4String creatorName = creator->GetProcessName();
          if(creatorName == "phot" || creatorName == "compt" ||
             creatorName == "conv")
          {
            // since this is happening before the secondary is being tracked,
            // the vertex position has not been set yet (set in initial step)
            fEventAction->SetConvPos((*fSecondary)[lp1]->GetPosition());
          }
        }
      }
    }

    if(fOneStepPrimaries && thePrePV->GetName() == "scintillator")
      theTrack->SetTrackStatus(fStopAndKill);
  }

  if(!thePostPV)
  {  // out of world
    fExpectedNextStatus = Undefined;
    return;
  }

  if(theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  {
    // Optical photon only

    if(thePrePV->GetName() == "Slab")
      // force drawing of photons in WLS slab
      trackInformation->SetForceDrawTrajectory(true);
    else if(thePostPV->GetName() == "expHall")
      // Kill photons entering expHall from something other than Slab
      theTrack->SetTrackStatus(fStopAndKill);

    // Was the photon absorbed by the absorption process
    if(thePostPoint->GetProcessDefinedStep()->GetProcessName() ==
       "OpAbsorption")
    {
      fEventAction->IncAbsorption();
      trackInformation->AddTrackStatusFlag(absorbed);
    }

    boundaryStatus = boundary->GetStatus();

    if(thePostPoint->GetStepStatus() == fGeomBoundary)
    {
      // Check to see if the particle was actually at a boundary
      // Otherwise the boundary status may not be valid
      if(fExpectedNextStatus == StepTooSmall)
      {
        if(boundaryStatus != StepTooSmall)
        {
          G4ExceptionDescription ed;
          ed << "LXeSteppingAction::UserSteppingAction(): "
             << "No reallocation step after reflection!" << G4endl;
          G4Exception("LXeSteppingAction::UserSteppingAction()", "LXeExpl01",
                      FatalException, ed,
                      "Something is wrong with the surface normal or geometry");
        }
      }
      fExpectedNextStatus = Undefined;
      switch(boundaryStatus)
      {
        case Absorption:
          trackInformation->AddTrackStatusFlag(boundaryAbsorbed);
          fEventAction->IncBoundaryAbsorption();
          break;
        case Detection:  // Note, this assumes that the volume causing detection
                         // is the photocathode because it is the only one with
                         // non-zero efficiency
        {
          // Trigger sensitive detector manually since photon is
          // absorbed but status was Detection
          G4SDManager* SDman = G4SDManager::GetSDMpointer();
          G4String sdName    = "/LXeDet/pmtSD";
          LXePMTSD* pmtSD    = (LXePMTSD*) SDman->FindSensitiveDetector(sdName);
          if(pmtSD)
            pmtSD->ProcessHits_boundary(theStep, nullptr);
          trackInformation->AddTrackStatusFlag(hitPMT);
          break;
        }
        case FresnelReflection:
        case TotalInternalReflection:
        case LambertianReflection:
        case LobeReflection:
        case SpikeReflection:
        case BackScattering:
          trackInformation->IncReflections();
          fExpectedNextStatus = StepTooSmall;
          break;
        default:
          break;
      }
      if(thePostPV->GetName() == "sphere")
        trackInformation->AddTrackStatusFlag(hitSphere);
    }
  }
}
