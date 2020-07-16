
// Local Headers
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "Reaction.hh"

// Default Headers
#include <fstream>
#include <csignal>

// Geant4 Headers
#include "globals.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4GenericIon.hh"

using namespace CLHEP;
using namespace std;

// --------------------------------------------------------------------------------------------------------------------------------- //
Reaction::Reaction(const G4String &aName)
    : G4VProcess(aName)
{
  theProcessType = (G4ProcessType)6;      // Decay
  theProcessSubType = (G4ProcessType)231; //DecayExt
}

// --------------------------------------------------------------------------------------------------------------------------------- //

Reaction::~Reaction()
{
}

// --------------------------------------------------------------------------------------------------------------------------------- //
// This method is called by Geant4 for every particle step, it tells the particle what to do next
// Inside this method contains the reaction process only called if flag "reaction_here" is true
// All used methods for calculate the 4Vectors can be found in EventAction.hh

G4VParticleChange *Reaction::PostStepDoIt(const G4Track &aTrack,
                                          const G4Step &)
{
  // Retrieving inputs
  Inputs *Inputs = &Inputs::GetInputs();

  // Stop the current particle, if requested by G4UserLimits
  aParticleChange.Initialize(aTrack);

  // Reaction process
  if (reaction_here)
  { // Here the physics of the Reaction is defined

    reaction_here = false;                            // Set the flag so it does not do the reaction a second time
    gCESimulationManager->ThereWasAReaction();        // A reaction ocurred

    aParticleChange.ProposeTrackStatus(fStopAndKill); // Kill the incident Particle

    G4double ThetaInCM = (40 * CLHEP::RandFlat::shoot() + 10) * degree; // 0 - 10 deg from z
    G4double randomPhiInCM = CLHEP::RandFlat::shoot() * 2 * pi;    // 0 - 2pi in transverse angle (azimuth)

    // Lorentz Vectors for each particle: ejectile, recoil, decay1 and decay2
    G4LorentzVector recoil4VecLab;
    G4LorentzVector ejectile4VecLab;
    bool newdecay = false;
    G4LorentzVector decayP1_4Vec;
    G4LorentzVector decayP2_4Vec;

    // Calculate the 4Vector for Ejectile and Recoil particles
    gCESimulationManager->CalculateLab4Vectors(aTrack, ThetaInCM, randomPhiInCM, recoil4VecLab, ejectile4VecLab);

    // Set the position of the interaction
    G4ThreeVector pos(aTrack.GetPosition().getX(), aTrack.GetPosition().getY(), aTrack.GetPosition().getZ());
    gCESimulationManager->SetInteractionPoint(pos);

    // Decay event
    gCESimulationManager->SetDecayThisEvent(newdecay);
    if (newdecay)
    {
      // This method calculates the 4vector for decay particles
      gCESimulationManager->DecayLab4Vectors(ejectile4VecLab, decayP1_4Vec, decayP2_4Vec);
      // Adding a the recoil particle as a secondary of this reaction
      aParticleChange.AddSecondary(gCESimulationManager->GetRecoilDynamicParticle(aTrack, recoil4VecLab), pos, true);
      // Adding a the decay product 1 as a secondary of this reaction
      aParticleChange.AddSecondary(gCESimulationManager->GetDecay1DynamicParticle(decayP1_4Vec), pos, true);
      // Adding a the decay product 2 as a secondary of this reaction
      aParticleChange.AddSecondary(gCESimulationManager->GetDecay2DynamicParticle(decayP2_4Vec), pos, true);
    }
    else
    {
      G4DynamicParticle *RecoilParticle = gCESimulationManager->GetRecoilDynamicParticle(aTrack, recoil4VecLab);
      G4DynamicParticle *EjectileParticle = gCESimulationManager->GetEjectileDynamicParticle(aTrack, ejectile4VecLab);

      // Adding a the recoil particle as a secondary of this reaction
      aParticleChange.AddSecondary(RecoilParticle , pos, true);

      // Adding a the ejectile particle as a secondary of this reaction
      aParticleChange.AddSecondary(EjectileParticle, pos, true);
    }
  }
  return &aParticleChange;
}
// --------------------------------------------------------------------------------------------------------------------------------- //
// This method calculates the particle position and tells Geant4 if any reaction should occur

G4double Reaction::PostStepGetPhysicalInteractionLength(const G4Track &aTrack,
                                                        G4double,
                                                        G4ForceCondition *condition)
{

  // Retrieving inputs
  Inputs *Inputs = &Inputs::GetInputs();

  //Small number saying that we have made it to the reaction point
  //Ie when Zdiff < eps it is at the reaction point
  G4double eps = 1e-5 * cm;

  reaction_here = false;
  *condition = NotForced;

  // Getting logical name of the volume which the particle is passing through
  G4String name = aTrack.GetVolume()->GetLogicalVolume()->GetName();

  // If the volume is "Log_Target" (Target logical name), it means the particle is inside the target and the reaction should occur
  if (name == "Log_Target" && gCESimulationManager->GetWasThereACEReaction() == false)
  {
    // Get z position of the target
    // 301 cm correspods to the center of the second solenoid
    G4double ZCEReaction = Inputs->target_pos.getZ() /* + 301 * cm */; 

    // Get z position of the particle
    G4double ZCurrent = aTrack.GetPosition().getZ();

    G4double Zdiff = ZCEReaction - ZCurrent;

    if (Zdiff < 0)
    {
      //If the current Z position is past the reaction point
      //don't do the reaction.  Return DBL_MAX to ensure that this
      //processes isn't called
      return DBL_MAX;
    }
    if (Zdiff > eps)
    {
      //if we are before the reaction point but not within eps of it
      //return the distance to the reaction point as the step length
      //NOTE:  unless the trajectory of the step is directly along the z-axis
      //it will have to make several jumps to get with eps of the reaction point.

      G4ThreeVector dir = aTrack.GetMomentumDirection();

      dir *= (ZCEReaction - ZCurrent);

      return dir.mag();
    }
    if (Zdiff <= eps && Zdiff >= 0)
    {
      //particle is now withing eps of reaction point.  Do reaction.
      //returns 0. to ensure that this process is called and defines the
      //length of the step.
      reaction_here = true;
      return 0.;
    }
  }
  else
  {
    //We are not within the target so return DBL_MAX to make sure
    //process isn't invoked.
    return DBL_MAX;
  }
}
// --------------------------------------------------------------------------------------------------------------------------------- //
