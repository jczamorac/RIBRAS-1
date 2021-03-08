#ifndef __REACTION_HH
#define __REACTION_HH

// Geant4 headers
#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

#include "TGraph.h"


class Reaction : public G4VProcess
{
public:
  // Flag for reaction
  G4bool reaction_here;

  //X section ionTable
  TGraph * xsectable;

  // Constructor
  Reaction(const G4String &processName = "Reaction");

  // Destructor
  virtual ~Reaction();

  // Methods used by Geant4 for every step
  virtual G4double PostStepGetPhysicalInteractionLength(
      const G4Track &track,
      G4double previousStepSize,
      G4ForceCondition *condition);

  virtual G4VParticleChange *PostStepDoIt(
      const G4Track &,
      const G4Step &);

  G4double GetTheta_Xsec();

  //  This methods are called by Geant4 when a particle has kinect enery 0
  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track &, G4ForceCondition *condition)
  {
    *condition = NotForced;

    return 9999999999.0;
  }

  virtual G4VParticleChange *AtRestDoIt(const G4Track &aTrack, const G4Step &)
  {
    aParticleChange.Initialize(aTrack);
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return &aParticleChange;
  };

  // No operation in AlongStepGPIL
  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track &, G4double, G4double, G4double &, G4GPILSelection *)
  {
    return -1.0;
  }

  // No operation in AlongStepDoIt
  virtual G4VParticleChange *AlongStepDoIt(const G4Track &, const G4Step &)
  {
    return NULL;
  }

private:
  // hide assignment operator as private
  Reaction &operator=(const Reaction &) { return *this; };
};

#endif
