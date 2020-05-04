// $Id: RunAction.cc 28 2010-01-12 10:24:06Z adotti $
/**
  * @file   RunAction.cc
  *
  * @date   17 Dec 2009
  * @author adotti
  *
  * @brief  Implements user class RunAction.
  */
#include "RunAction.hh"
#include "EventAction.hh"
#include "G4Run.hh"
#include "DetectorConstruction.hh"
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

RunAction::RunAction(EventAction *theEventAction) : eventAction(theEventAction)
{
        eventAction->SetRootSaver(&saver);
        MassMap = new TEnv("mass_table2.txt");
        eventAction->SetMassMap(MassMap);
}

void RunAction::BeginOfRunAction(const G4Run *aRun)
{
        G4cout << "Starting Run: " << aRun->GetRunID() << G4endl;
        // For each run a new TTree is created, with default names
        saver.CreateTree();

}

void RunAction::EndOfRunAction(const G4Run *aRun)
{
        saver.CloseTree();
        delete MassMap;
}
