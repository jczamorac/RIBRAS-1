#include "Reaction.hh"
#include "globals.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4GenericIon.hh"
#include <fstream>
#include <csignal>
#include "EventAction.hh"
#include "DetectorConstruction.hh"

using namespace CLHEP;
using namespace std;

Reaction::Reaction(const G4String& aName)
  : G4VProcess(aName){
  
  if (verboseLevel>1) {
    G4cout <<GetProcessName() << " is created "<< G4endl;
  }
  theProcessType=(G4ProcessType)6;       // Decay
  theProcessSubType=(G4ProcessType)231;  //DecayExt
}

Reaction::~Reaction(){
}

G4VParticleChange* Reaction::PostStepDoIt( const G4Track& aTrack,
					     const G4Step&){

  // Stop the current particle, if requested by G4UserLimits
  aParticleChange.Initialize(aTrack);
  
  if(reaction_here){// Here the physics of the Reaction is defined

    reaction_here=false;                                      // Set the flag so it does not do the reaction a second time
    gCESimulationManager->ThereWasAReaction();                // A reaction ocurred
    aParticleChange.ProposeTrackStatus(fStopAndKill);         // Kill the incident Particle

    G4double ThetaInCM =(5*CLHEP::RandFlat::shoot())*degree;  // 0 - 5 deg from z
    G4double randomPhiInCM = CLHEP::RandFlat::shoot()*2*pi;   // 0 - 2pi in transverse angle (azimuth)
    G4LorentzVector recoil4VecLab;
    G4LorentzVector ejectile4VecLab;

    // Calculate the 4Vector for Ejectile and Recoil particles
    gCESimulationManager->CalculateLab4Vectors(aTrack,ThetaInCM,randomPhiInCM,recoil4VecLab,ejectile4VecLab);     

    // Set the position of the interaction
    G4ThreeVector pos(aTrack.GetPosition().getX(),aTrack.GetPosition().getY(),aTrack.GetPosition().getZ());
    gCESimulationManager->SetInteractionPoint(pos);

    // Adding a the recoil particle as a secondary of this reaction
    aParticleChange.AddSecondary(gCESimulationManager->GetRecoilDynamicParticle(aTrack,recoil4VecLab),pos,true); 

    // Adding a the ejectile particle as a secondary of this reaction
    aParticleChange.AddSecondary(gCESimulationManager->GetEjectileDynamicParticle(aTrack,ejectile4VecLab),pos,true); 
  }
  return &aParticleChange;
}

G4double Reaction::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
							  G4double,
							  G4ForceCondition* condition){
  
  Inputs* Inputs = &Inputs::GetInputs();

  G4double eps= 0.1 *cm;//Small number saying that we have made it to the reaction point
                        //Ie when Zdiff < eps it is at the reaction point
  
  reaction_here=false;
  *condition=NotForced;

  G4String name = aTrack.GetVolume()->GetLogicalVolume()->GetName();

  if(name=="alvolog" && gCESimulationManager->GetWasThereACEReaction()==false){
    
    // Get z position of the target
    G4double ZCEReaction = /*gCESimulationManager->GetReactionZPoint()*/ Inputs->target_pos.getZ() + 301 *cm;

    // Get z position of the particle
    G4double ZCurrent = aTrack.GetPosition().getZ();

    G4double Zdiff = ZCEReaction-ZCurrent;

    if(Zdiff<0){
      //If the current Z position is past the reaction point
      //don't do the reaction.  Return DBL_MAX to ensure that this
      //processes isn't called
      return DBL_MAX;
    }
    if(Zdiff>eps){
      //if we are before the reaction point but not within eps of it
      //return the distance to the reaction point as the step length
      //NOTE:  unless the trajectory of the step is directly along the z-axis
      //it will have to make several jumps to get with eps of the reaction point.

      G4ThreeVector dir=aTrack.GetMomentumDirection();

      dir*=(ZCEReaction-ZCurrent);

      return dir.mag();
    }
    if(Zdiff<= eps && Zdiff>=0){
      //particle is now withing eps of reaction point.  Do reaction.
      //returns 0. to ensure that this process is called and defines the
      //length of the step.
      reaction_here=true;
      return 0.;
    }
  }
  else {
  //We are not within the target so return DBL_MAX to make sure
  //process isn't invoked.
  return DBL_MAX;
  }
}

