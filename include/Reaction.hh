#ifndef __REACTION_HH
#define __REACTION_HH


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

class Reaction : public G4VProcess
{
public:
  G4bool reaction_here;

  Reaction(const G4String& processName ="Reaction" );

  virtual ~Reaction();

  virtual G4double PostStepGetPhysicalInteractionLength(
							const G4Track& track,
							G4double   previousStepSize,
							G4ForceCondition* condition
							);

  virtual G4VParticleChange* PostStepDoIt(
					  const G4Track& ,
					  const G4Step&
					  );

  //  no operation in  AtRestGPIL
  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track& ,G4ForceCondition* condition){
    //G4cout<<"IN AT REST GPIL"<<G4endl;
    *condition=NotForced;
    return 9999999999.0;
  }

  //  no operation in  AtRestDoIt
  virtual G4VParticleChange* AtRestDoIt(const G4Track& aTrack  ,const G4Step&){
    //G4cout<<"AT REST DO IT "<<G4endl;
    aParticleChange.Initialize(aTrack);


    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return &aParticleChange;
  };

  //  no operation in  AlongStepGPIL
  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track&,G4double,G4double,G4double&,G4GPILSelection*){
    G4cout<<"In AlongStep GPIL"<<G4endl;
    return -1.0;
  }

  //  no operation in  AlongStepDoIt
  virtual G4VParticleChange* AlongStepDoIt( const G4Track& , const G4Step&){
    return NULL;
  }

private:

  // hide assignment operator as private
  Reaction& operator=(const Reaction&){return *this;};

};

#endif /* __REACTION_HH */
