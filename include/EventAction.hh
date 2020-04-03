//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *2
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

 #ifndef EVENTACTION_HH_
 #define EVENTACTION_HH_
 
 
 #include "G4UserEventAction.hh"
 #include "G4String.hh"
 #include "G4Types.hh"
 #include "G4ClassificationOfNewTrack.hh"
 #include "G4TrackStatus.hh"
 #include "G4LorentzVector.hh"
 #include "globals.hh"
 #include <csignal>
 #include <fstream>
 #include <string>
 #include "G4ParticleDefinition.hh"
 #include "DetectorConstruction.hh"
 //#include "NoiseGenerator.hh"
 //#include "CrosstalkGenerator.hh"
// #include "MeV2ChargeConverter.hh"

using namespace CLHEP;
using namespace std;

class G4Event;
class RootSaver;
class G4VPhysicalVolume;
class G4Event;
class G4Run;
class G4Track;
class G4Step;
class G4ParticleGun;
class G4Ions;
class G4DynamicParticle;
class G4DecayProducts;
 
extern class EventAction *gCESimulationManager;
class EventAction : public G4UserEventAction
{
public:
//! Default constructor
EventAction();
//! Default destructor
virtual ~EventAction() {
    if (gCESimulationManager == this) {
    gCESimulationManager = (EventAction *)0;
}};
//! Beginning of event
void BeginOfEventAction(const G4Event* anEvent);
//! Digitize hits and store information
void EndOfEventAction(const G4Event* anEvent);
//! Set the RootSaver
inline void SetRootSaver( RootSaver* saver ) { rootSaver = saver; }

//Methods for Reaction
void CalculateLab4Vectors(const G4Track &BeamTrack,G4double CMScatteringAngle,G4double phi, G4LorentzVector & RecoilOut,G4LorentzVector & EjectileOut);
//!Calculate Recoil Particle
G4DynamicParticle * GetRecoilDynamicParticle(const G4Track & BeamTrack, const G4LorentzVector&);
//!Calculate Ejectile Particle
G4DynamicParticle * GetEjectileDynamicParticle(const G4Track & BeamTrack, const G4LorentzVector&);
//!Bool for Reaction
void ThereWasAReaction(){rThereWasACEReactionThisEvent=true;}
bool GetWasThereACEReaction(){return rThereWasACEReactionThisEvent;}
//!Get Z of the reaction point
inline G4double GetReactionZPoint(){return rZOfReactionInTarget;}
//!Set Interaction Point
inline void SetInteractionPoint(G4ThreeVector v){
    rInteractionPointX=v.x();
    rInteractionPointY=v.y();
    rInteractionPointZ=v.z();}

private:
         //! pointer to saver object
         RootSaver* rootSaver;
         //! hits collection name
         G4String hitsCollName;
         //! digits collection name
         //G4String digitsCollName;
         //! Hits collection ID
         G4int hitsCollID;

         G4Ions* RecoilParticle;
         G4Ions* EjectileParticle;

         G4bool rThereWasACEReactionThisEvent=false;

         double_t rZOfReactionInTarget;

         double_t rInteractionPointX;
         double_t rInteractionPointY;
         double_t rInteractionPointZ;

         double_t rInteractionPointBeamEnergy;

         double_t rRecoilKineticEnergy;

         double_t rRecoilTheta;
         double_t rEjectileTheta;

         void rReset();
 };
 
 #endif /* EVENTACTION_HH_ */
