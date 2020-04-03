 // $Id: RunAction.cc 28 2010-01-12 10:24:06Z adotti $
 /**
  * @file   RunAction.cc
  *
  * @date   17 Dec 2009
  * @author adotti
  *
  * @brief  Implements user class RunAction.
  */
 
#include <fstream>
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
 
ofstream sav("Inputs.txt");

RunAction::RunAction(EventAction* theEventAction ) :
        eventAction(theEventAction)
{
        eventAction->SetRootSaver( &saver );
}

void RunAction::BeginOfRunAction(const G4Run* aRun )
{
Inputs* Inputs = &Inputs::GetInputs();
sav << "Target Inputs: " << G4endl;
sav << "Z: " << Inputs->target_Z << G4endl;
sav << "A: " << Inputs->target_A << G4endl;
sav << "Mass: " << Inputs->target_mass /MeV << " MeV" << G4endl;
sav << " " << G4endl;

sav << "Recoil Inputs: " << G4endl;
sav << "Z: " << Inputs->recoil_Z << G4endl;
sav << "A: " << Inputs->recoil_A << G4endl;
sav << "Mass: " << Inputs->recoil_mass /MeV << " MeV" << G4endl;
sav << " " << G4endl;

sav << "Ejectile Inputs: " << G4endl;
sav << "Z: " << Inputs->ejectile_Z << G4endl;
sav << "A: " << Inputs->ejectile_A << G4endl;
sav << "Mass: " << Inputs->ejectile_mass /MeV << " MeV" << G4endl;
sav << " " << G4endl;
// G4cout<<"Starting Run: "<<aRun->GetRunID()<<G4endl;
//For each run a new TTree is created, with default names
saver.CreateTree();
}
 
 void RunAction::EndOfRunAction( const G4Run* aRun )
 {
        //  G4cout<<"Ending Run: "<<aRun->GetRunID()<<G4endl;
        //  G4cout<<"Number of events: "<<aRun->GetNumberOfEvent()<<G4endl;
         saver.CloseTree();
 }