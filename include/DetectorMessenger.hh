// $Id: DetectorMessenger.hh 33 2010-01-14 17:08:18Z adotti $

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

/**
  * @file
  * @brief defines class DetectorMessenger
  */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "TGraph.h"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

using namespace std;

class DetectorMessenger : public G4UImessenger
{
public:
  //! Constructor
  DetectorMessenger(DetectorConstruction *);
  //! Destructor
  ~DetectorMessenger();
  //! handle user commands
  void SetNewValue(G4UIcommand *command, G4String newValue);

private:
  DetectorConstruction *detector;

  // Creating directories
  G4UIdirectory *detDir;
  G4UIdirectory *secondSensorDir;
  G4UIdirectory *Dir;
  G4UIdirectory *targetDir;
  G4UIdirectory *recoilDir;
  G4UIdirectory *ejectileDir;
  G4UIdirectory *decayprod1Dir;
  G4UIdirectory *decayprod2Dir;
  G4UIdirectory *primaryDir;
  G4UIdirectory *detectorDir;
  G4UIdirectory *magneticDir;
  G4UIdirectory *xsectionDir;

  G4UIcmdWithoutParameter *updateCmd;

  G4UIcmdWithABool *setDUTsetupCmd;

  // Commands for target
  G4UIcmdWithAString *target_material;
  G4UIcmdWithADoubleAndUnit *target_width_cmd;
  G4UIcmdWithADoubleAndUnit *target_mass;
  G4UIcmdWithAnInteger *target_A;
  G4UIcmdWithAnInteger *target_Z;
  G4UIcmdWith3VectorAndUnit *target_pos;

  // Turn On/Off magnetic field
  G4UIcmdWithABool *magneticfieldon;

  // Commands for recoil particle
  G4UIcmdWithADoubleAndUnit *recoil_Ex;
  G4UIcmdWithADoubleAndUnit *recoil_mass;
  G4UIcmdWithAnInteger *recoil_A;
  G4UIcmdWithAnInteger *recoil_Z;

  // Commands for ejectile particle
  G4UIcmdWithADoubleAndUnit *ejectile_mass;
  G4UIcmdWithADoubleAndUnit *ejectile_Ex;
  G4UIcmdWithAnInteger *ejectile_A;
  G4UIcmdWithAnInteger *ejectile_Z;

  // Commands for decay particle 1
  G4UIcmdWithADoubleAndUnit *decayp1_mass;
  G4UIcmdWithADoubleAndUnit *decayp1_Ex;
  G4UIcmdWithAnInteger *decayp1_A;
  G4UIcmdWithAnInteger *decayp1_Z;

  // Commands for decay particle 2
  G4UIcmdWithADoubleAndUnit *decayp2_mass;
  G4UIcmdWithADoubleAndUnit *decayp2_Ex;
  G4UIcmdWithAnInteger *decayp2_A;
  G4UIcmdWithAnInteger *decayp2_Z;

  // Commands for primary beam
  G4UIcmdWithADoubleAndUnit *primary_energy;
  G4UIcmdWithAnInteger *primary_A;
  G4UIcmdWithAnInteger *primary_Z;
  G4UIcmdWith3VectorAndUnit *primary_pos;

  // Commands for detectors
  G4UIcmdWithABool *usingdetectors;

  // Commands for external Xsection
  G4UIcmdWithABool *usingXsec;
  G4UIcmdWithAString *xsecfile;

public:
  TGraph * graphtable;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
