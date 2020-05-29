 // Local headers
#include "RunAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

// Geant4 headers
#include "G4Run.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4Gamma.hh"
#include "G4DecayProducts.hh"
#include "G4GenericIon.hh"

using namespace std;

// -------------------------------------------------------------------------------//
RunAction::RunAction(EventAction *theEventAction) : eventAction(theEventAction)
{
        eventAction->SetRootSaver(&saver);
}

// -------------------------------------------------------------------------------//
void RunAction::BeginOfRunAction(const G4Run *aRun)
{
        MassMap = new TEnv("mass_table2.txt");
        eventAction->SetMassMap(MassMap);

        // For each run a new TTree is created, with default names
        saver.CreateTree();
}

// -------------------------------------------------------------------------------//
void RunAction::EndOfRunAction(const G4Run *aRun)
{
        saver.CloseTree();
        delete MassMap;
}
// -------------------------------------------------------------------------------//