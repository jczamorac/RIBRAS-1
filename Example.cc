// Default headers
#include <fstream>

// Geant4 headers
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4Version.hh"
#include "G4VisExecutive.hh"
#include "G4ScoringManager.hh"

#if  G4VERSION_NUMBER>=930
#include "G4UIExecutive.hh"
#else
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#endif

// Local headers
#include "DetectorConstruction.hh"
#include "BeamTestPrimaryGeneratorAction.hh"
#include "PhysicsList.hh"
#include "EventAction.hh"
#include "RunAction.hh"

// Physics list (Geant4 pre-made)
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"

using namespace std;

int main(int argc,char** argv)
{
  // Run manager
  G4RunManager * runManager = new G4RunManager();
  runManager->SetVerboseLevel(0);

  // mandatory Initialization classes
  G4VUserDetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

  // Local user Physics List
  G4VUserPhysicsList* physics = new PhysicsList();
  runManager->SetUserInitialization(physics);

  // mandatory User class
  G4VUserPrimaryGeneratorAction* gen_action = new BeamTestPrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);

  // Optional User Action classes
  EventAction* event_action = new EventAction();
  RunAction* run_action = new RunAction(event_action);
  runManager->SetUserAction( event_action );
  runManager->SetUserAction( run_action );

  // Initialize G4 kernel
  runManager->Initialize();

  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive();
  visManager->SetVerboseLevel(0);
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager * UImanager = G4UImanager::GetUIpointer();
  UImanager->SetVerboseLevel(0);

  if (argc!=1) {  // batch mode

      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
  }
  else {           // interactive mode : define UI session

#if  G4VERSION_NUMBER>=930
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute visQt.mac");
    else
      UImanager->ApplyCommand("/control/execute vis.mac");

#else
  #ifdef G4UI_USE_TCSH
    G4UIsession * ui = new G4UIterminal(new G4UItcsh);
  #else
    G4UIsession * ui = new G4UIterminal();
  #endif
    UImanager->ApplyCommand("/control/execute vis.mac");
  #endif
    ui->SessionStart();
    delete ui;

  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !


  delete runManager;

  return 0;
}
