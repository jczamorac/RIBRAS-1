// $Id: DetectorMessenger.cc 94 2010-01-26 13:18:30Z adotti $
/**
* @file
* @brief Implements class DetectorMessenger.
*/

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4ThreeVector.hh"
#include <sstream>
#include <stdlib.h>
#include <map>
#include <cstdio>
#include <iomanip>

using namespace CLHEP;
using namespace std;

DetectorMessenger::DetectorMessenger(DetectorConstruction * det)
:detector(det)
{
detDir = new G4UIdirectory("/det/");
detDir->SetGuidance("detector construction commands");

//Target Directory
targetDir = new G4UIdirectory("/det/target/");
targetDir->SetGuidance("target control");

//Recoil Directory
recoilDir = new G4UIdirectory("/det/recoil/");
recoilDir->SetGuidance("recoil control");

//Ejectile Directory
ejectileDir = new G4UIdirectory("/det/ejectile/");
ejectileDir->SetGuidance("ejectile control");

//Primary Directory
primaryDir = new G4UIdirectory("/det/primary/");
primaryDir->SetGuidance("Primary Beam control");

//Sensor Directory
secondSensorDir = new G4UIdirectory("/det/secondSensor/");
secondSensorDir->SetGuidance("comands related to the second sensor plane");

xShiftCmd = new G4UIcmdWithADoubleAndUnit("/det/secondSensor/xShift",this);
xShiftCmd->SetGuidance("Define x-shift of second sensor plane");
xShiftCmd->SetParameterName("xShift",true);
xShiftCmd->SetDefaultValue(0.);
xShiftCmd->SetUnitCategory("Length");
xShiftCmd->SetDefaultUnit("um");

yShiftCmd = new G4UIcmdWithADoubleAndUnit("/det/secondSensor/yShift",this);
yShiftCmd->SetGuidance("Define y-shift of second sensor plane");
yShiftCmd->SetParameterName("yShift",true);
yShiftCmd->SetDefaultValue(0.);
yShiftCmd->SetUnitCategory("Length");
yShiftCmd->SetDefaultUnit("um");

z_slitCmd = new G4UIcmdWithADoubleAndUnit("/det/slit/z",this);
z_slitCmd->SetGuidance("Define z pos of slit");
z_slitCmd->SetParameterName("z_slit",true);
z_slitCmd->SetDefaultValue(30.);
z_slitCmd->SetUnitCategory("Length");
z_slitCmd->SetDefaultUnit("mm");

x_slitCmd = new G4UIcmdWithADoubleAndUnit("/det/slit/x",this);
x_slitCmd->SetGuidance("Define x pos of slit");
x_slitCmd->SetParameterName("x_slit",true);
x_slitCmd->SetDefaultValue(0.);
x_slitCmd->SetUnitCategory("Length");
x_slitCmd->SetDefaultUnit("mm");

xShiftPocketCmd = new G4UIcmdWithADoubleAndUnit("/det/pocket/xShiftPocket",this);
xShiftPocketCmd->SetGuidance("Define x-shift of firtst pocket");
xShiftPocketCmd->SetParameterName("xShiftPocket",true);
xShiftPocketCmd->SetDefaultValue(0.);
xShiftPocketCmd->SetUnitCategory("Length");
xShiftPocketCmd->SetDefaultUnit("mm");

updateCmd = new G4UIcmdWithoutParameter("/det/update",this);
updateCmd->SetGuidance("force to recompute geometry.");
updateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
updateCmd->SetGuidance("if you changed geometrical value(s).");
updateCmd->AvailableForStates(G4State_Idle);

target_material = new G4UIcmdWithAString("/det/target/material", this);
target_material->SetGuidance("Specify the foil material name. This should be the G4_Material name (e.g. G4_lh2)");

target_A = new G4UIcmdWithAnInteger("/det/target/A", this);
target_A->SetGuidance("Specify a target mass number");

target_Z = new G4UIcmdWithAnInteger("/det/target/Z", this);
target_Z->SetGuidance("Specify a target atomic number");

target_mass = new G4UIcmdWithADoubleAndUnit("/det/target/mass", this);
target_mass->SetGuidance("Set target atomic mass (MeV)");
target_mass->SetDefaultUnit("MeV");

recoil_mass = new G4UIcmdWithADoubleAndUnit("/det/recoil/mass", this);
recoil_mass->SetGuidance("(Optional) Specify an atomic recoil mass in MeV");

recoil_Ex = new G4UIcmdWithADoubleAndUnit("/det/recoil/Ex", this);
recoil_Ex->SetGuidance("(Optional) Specify a recoil excitation energy in MeV");

recoil_A = new G4UIcmdWithAnInteger("/det/recoil/A", this);
recoil_A->SetGuidance("Specify a recoil mass number");

recoil_Z = new G4UIcmdWithAnInteger("/det/recoil/Z", this);
recoil_Z->SetGuidance("Specify a recoil atomic number");

ejectile_mass = new G4UIcmdWithADoubleAndUnit("/det/ejectile/mass", this);
ejectile_mass->SetGuidance("(Optional) Specify an atomic ejectile mass in MeV");

ejectile_Ex = new G4UIcmdWithADoubleAndUnit("/det/ejectile/Ex", this);
ejectile_Ex->SetGuidance("(Optional) Specify an ejectile excitation energy in MeV");

ejectile_A = new G4UIcmdWithAnInteger("/det/ejectile/A", this);
ejectile_A->SetGuidance("Specify a ejectile mass number");

ejectile_Z = new G4UIcmdWithAnInteger("/det/ejectile/Z", this);
ejectile_Z->SetGuidance("Specify a ejectile atomic number");

primary_energy = new G4UIcmdWithADoubleAndUnit("/det/primary/energy", this);
primary_energy->SetGuidance("Set primary beam energy in MeV");

primary_A = new G4UIcmdWithAnInteger("/det/primary/A", this);
primary_A->SetGuidance("Set mass of primary particle");

primary_Z = new G4UIcmdWithAnInteger("/det/primary/Z", this);
primary_Z->SetGuidance("Set atomic numer of primary particle");

primary_pos = new G4UIcmdWith3VectorAndUnit("/det/primary/pos", this);
primary_pos->SetGuidance("Set primary beam position");
primary_pos->SetParameterName("X","Y","Z",true,true);
primary_pos->SetDefaultUnit("cm");

magneticfieldon = new G4UIcmdWithABool("/det/field", this);
magneticfieldon->SetGuidance("Turn On (1), Turn Off (0)");

target_pos = new G4UIcmdWith3VectorAndUnit("/det/target/pos", this);
target_pos->SetGuidance("Set target position");
target_pos->SetParameterName("X","Y","Z",true,true);
target_pos->SetDefaultUnit("cm");

target_arial_density_cmd = new G4UIcmdWithADoubleAndUnit("/det/target/thickness", this);
target_arial_density_cmd->SetGuidance("Specify a target thickness (arial density) in units mg/cm2");

target_width_cmd = new G4UIcmdWithADoubleAndUnit("/det/target/width" ,this);
target_width_cmd->SetGuidance("Specify the target width in units mm");
target_width_cmd->SetDefaultUnit("mm");
}

DetectorMessenger::~DetectorMessenger()
{
  delete xShiftCmd;
  delete yShiftCmd;
  delete x_slitCmd;
  delete z_slitCmd;
  delete xShiftPocketCmd;
  delete secondSensorDir;
  delete target_material;
  delete target_Z;
  delete target_A;
  delete targetDir;
  delete recoilDir;
  delete recoil_mass;
  delete recoil_A;
  delete recoil_Z;
  delete recoil_Ex;
  delete ejectileDir;
  delete ejectile_Z;
  delete ejectile_mass;
  delete ejectile_Ex;
  delete ejectile_A;
  delete updateCmd;
  delete detDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  Inputs* Inputs = &Inputs::GetInputs();

  if ( command == xShiftCmd ) {
    G4ThreeVector pos= detector->SecondSensorPosition();
    pos.setX( xShiftCmd->GetNewDoubleValue(newValue) );
    detector->SetSecondSensorPosition(pos);
  }

  else if ( command == yShiftCmd ) {
    G4ThreeVector pos= detector->SecondSensorPosition();
    pos.setY( yShiftCmd->GetNewDoubleValue(newValue) );
    detector->SetSecondSensorPosition(pos);
  }

  else if ( command == z_slitCmd ){detector->SetSlit_z(z_slitCmd->GetNewDoubleValue(newValue));}

  else if ( command == x_slitCmd ){detector->SetSlit_x(x_slitCmd->GetNewDoubleValue(newValue));}

  else if ( command == xShiftPocketCmd ){detector->SetPocket1_x(xShiftPocketCmd->GetNewDoubleValue(newValue));}

  else if ( command == updateCmd ){detector->UpdateGeometry();}

  else if (command == target_A){Inputs->target_A = target_A->GetNewIntValue(newValue);}
  else if (command == target_Z){Inputs->target_Z = target_Z->GetNewIntValue(newValue);}
  else if (command == target_material){Inputs->g4_material_name = newValue;}
  else if (command == target_arial_density_cmd) {Inputs->arial_density = target_arial_density_cmd->GetNewDoubleValue(newValue);}
  else if (command == target_pos){Inputs->target_pos = target_pos->GetNew3VectorValue(newValue) - G4ThreeVector(0., 0., 301 *cm);}
  else if (command == target_width_cmd){Inputs->width = target_width_cmd->GetNewDoubleValue(newValue);}
  else if (command == target_mass){Inputs->target_mass = target_mass->GetNewDoubleValue(newValue);}

  else if (command == recoil_mass){Inputs->recoil_mass = recoil_mass->GetNewDoubleValue(newValue);}
  else if (command == recoil_Ex){Inputs->recoil_Ex = recoil_Ex->GetNewDoubleValue(newValue);}
  else if (command == recoil_A){Inputs->recoil_A = recoil_A->GetNewIntValue(newValue);}
  else if (command == recoil_Z){Inputs->recoil_Z = recoil_Z->GetNewIntValue(newValue);}

  else if (command == ejectile_mass){Inputs->ejectile_mass = ejectile_mass->GetNewDoubleValue(newValue);}
  else if (command == ejectile_Ex){Inputs->ejectile_Ex = ejectile_Ex->GetNewDoubleValue(newValue);}
  else if (command == ejectile_A){Inputs->ejectile_A = ejectile_A->GetNewIntValue(newValue);}
  else if (command == ejectile_Z){Inputs->ejectile_Z = ejectile_Z->GetNewIntValue(newValue);}

  else if (command == primary_energy){Inputs->primary_energy = primary_energy->GetNewDoubleValue(newValue);}
  else if (command == primary_A){Inputs->primary_A = primary_A->GetNewIntValue(newValue);}
  else if (command == primary_Z){Inputs->primary_Z = primary_Z->GetNewIntValue(newValue);}
  else if (command == primary_pos){Inputs->primary_pos = primary_pos->GetNew3VectorValue(newValue);}

  else if (command == magneticfieldon){Inputs->using_magneticfield = magneticfieldon->GetNewBoolValue(newValue);}
}
