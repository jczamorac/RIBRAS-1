// Geant4 headers
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4ThreeVector.hh"

// Local headers
#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

// Default headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <map>
#include <cstdio>
#include <iomanip>
#include <vector>



using namespace CLHEP;
using namespace std;

// ----------------------------------------------------------------------------- //

DetectorMessenger::DetectorMessenger(DetectorConstruction *det)
    : detector(det)
{

  // Main Directory
  detDir = new G4UIdirectory("/det/");
  detDir->SetGuidance("detector construction commands");

  // Target Directory
  targetDir = new G4UIdirectory("/det/target/");
  targetDir->SetGuidance("target control");

  // Recoil Directory
  recoilDir = new G4UIdirectory("/det/recoil/");
  recoilDir->SetGuidance("recoil control");

  // Ejectile Directory
  ejectileDir = new G4UIdirectory("/det/ejectile/");
  ejectileDir->SetGuidance("ejectile control");

  // Decay P1
  decayprod1Dir = new G4UIdirectory("/det/decayp1/");
  decayprod1Dir->SetGuidance("decayp1 control");

  // Decay P2
  decayprod2Dir = new G4UIdirectory("/det/decayp2/");
  decayprod2Dir->SetGuidance("decayp2 control");

  // Primary Directory
  primaryDir = new G4UIdirectory("/det/primary/");
  primaryDir->SetGuidance("Primary Beam control");

  // Sensor Directory
  secondSensorDir = new G4UIdirectory("/det/secondSensor/");
  secondSensorDir->SetGuidance("comands related to the second sensor plane");

  // Magnetic field Directory
  magneticDir = new G4UIdirectory("/det/magnetic/");
  magneticDir->SetGuidance("magnetic field control");

  // Update geometry command
  updateCmd = new G4UIcmdWithoutParameter("/det/update", this);
  updateCmd->SetGuidance("force to recompute geometry.");
  updateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  updateCmd->SetGuidance("if you changed geometrical value(s).");
  updateCmd->AvailableForStates(G4State_Idle);

  // Commands for target
  target_material = new G4UIcmdWithAString("/det/target/material", this);
  target_material->SetGuidance("Specify the foil material name. This should be the G4_Material name (e.g. G4_lh2)");

  target_A = new G4UIcmdWithAnInteger("/det/target/A", this);
  target_A->SetGuidance("Specify a target mass number");

  target_Z = new G4UIcmdWithAnInteger("/det/target/Z", this);
  target_Z->SetGuidance("Specify a target atomic number");

  target_mass = new G4UIcmdWithADoubleAndUnit("/det/target/mass", this);
  target_mass->SetGuidance("Set target atomic mass (MeV)");
  target_mass->SetDefaultUnit("MeV");

  target_pos = new G4UIcmdWith3VectorAndUnit("/det/target/pos", this);
  target_pos->SetGuidance("Set target position");
  target_pos->SetParameterName("X", "Y", "Z", true, true);
  target_pos->SetDefaultUnit("cm");

  target_width_cmd = new G4UIcmdWithADoubleAndUnit("/det/target/width", this);
  target_width_cmd->SetGuidance("Specify the target width in units mm");
  target_width_cmd->SetDefaultUnit("mm");

  // Commands for recoil
  recoil_mass = new G4UIcmdWithADoubleAndUnit("/det/recoil/mass", this);
  recoil_mass->SetGuidance("(Optional) Specify an atomic recoil mass in MeV");

  recoil_Ex = new G4UIcmdWithADoubleAndUnit("/det/recoil/Ex", this);
  recoil_Ex->SetGuidance("(Optional) Specify a recoil excitation energy in MeV");

  recoil_A = new G4UIcmdWithAnInteger("/det/recoil/A", this);
  recoil_A->SetGuidance("Specify a recoil mass number");

  recoil_Z = new G4UIcmdWithAnInteger("/det/recoil/Z", this);
  recoil_Z->SetGuidance("Specify a recoil atomic number");

  // Commands for ejectile
  ejectile_mass = new G4UIcmdWithADoubleAndUnit("/det/ejectile/mass", this);
  ejectile_mass->SetGuidance("(Optional) Specify an atomic ejectile mass in MeV");

  ejectile_Ex = new G4UIcmdWithADoubleAndUnit("/det/ejectile/Ex", this);
  ejectile_Ex->SetGuidance("(Optional) Specify an ejectile excitation energy in MeV");

  ejectile_A = new G4UIcmdWithAnInteger("/det/ejectile/A", this);
  ejectile_A->SetGuidance("Specify a ejectile mass number");

  ejectile_Z = new G4UIcmdWithAnInteger("/det/ejectile/Z", this);
  ejectile_Z->SetGuidance("Specify a ejectile atomic number");

  // Commands for decay particle 1
  decayp1_mass = new G4UIcmdWithADoubleAndUnit("/det/decayp1/mass", this);
  decayp1_mass->SetGuidance("(Optional) Specify an atomic  mass in MeV");

  decayp1_Ex = new G4UIcmdWithADoubleAndUnit("/det/decayp1/Ex", this);
  decayp1_Ex->SetGuidance("(Optional) Specify an excitation energy in MeV");

  decayp1_A = new G4UIcmdWithAnInteger("/det/decayp1/A", this);
  decayp1_A->SetGuidance("Specify a mass number");

  decayp1_Z = new G4UIcmdWithAnInteger("/det/decayp1/Z", this);
  decayp1_Z->SetGuidance("Specify an atomic number");

  // Commands for decay particle 2
  decayp2_mass = new G4UIcmdWithADoubleAndUnit("/det/decayp2/mass", this);
  decayp2_mass->SetGuidance("(Optional) Specify an atomic  mass in MeV");

  decayp2_Ex = new G4UIcmdWithADoubleAndUnit("/det/decayp2/Ex", this);
  decayp2_Ex->SetGuidance("(Optional) Specify an excitation energy in MeV");

  decayp2_A = new G4UIcmdWithAnInteger("/det/decayp2/A", this);
  decayp2_A->SetGuidance("Specify a mass number");

  decayp2_Z = new G4UIcmdWithAnInteger("/det/decayp2/Z", this);
  decayp2_Z->SetGuidance("Specify an atomic number");

  // Commands for primary beam
  primary_energy = new G4UIcmdWithADoubleAndUnit("/det/primary/energy", this);
  primary_energy->SetGuidance("Set primary beam energy in MeV");

  primary_A = new G4UIcmdWithAnInteger("/det/primary/A", this);
  primary_A->SetGuidance("Set mass of primary particle");

  primary_Z = new G4UIcmdWithAnInteger("/det/primary/Z", this);
  primary_Z->SetGuidance("Set atomic numer of primary particle");

  primary_pos = new G4UIcmdWith3VectorAndUnit("/det/primary/pos", this);
  primary_pos->SetGuidance("Set primary beam position");
  primary_pos->SetParameterName("X", "Y", "Z", true, true);
  primary_pos->SetDefaultUnit("cm");

  // Switch for magnetic field
  magneticfieldon = new G4UIcmdWithABool("/det/field", this);
  magneticfieldon->SetGuidance("Turn On (1), Turn Off (0)");

  // Define Xsection file
  usingXsec = new G4UIcmdWithABool("/det/xsec/ext_xsec", this);
  usingXsec->SetGuidance("Turn On (1), Turn Off (0)");
  xsecfile = new G4UIcmdWithAString("/det/xsec/filename", this);
  xsecfile->SetGuidance("Specify the filename of the Xsection located in the folder xsections");

}

// ----------------------------------------------------------------------------- //

DetectorMessenger::~DetectorMessenger()
{
  delete detDir;
  delete targetDir;
  delete recoilDir;
  delete ejectileDir;
  delete primaryDir;
  delete secondSensorDir;
  delete magneticDir;

  delete updateCmd;

  delete target_material;
  delete target_Z;
  delete target_A;
  delete target_mass;
  delete target_pos;
  delete target_width_cmd;
  delete recoil_mass;
  delete recoil_Ex;
  delete recoil_A;
  delete recoil_Z;
  delete ejectile_mass;
  delete ejectile_Ex;
  delete ejectile_A;
  delete ejectile_Z;
  delete decayp1_mass;
  delete decayp1_Ex;
  delete decayp1_A;
  delete decayp1_Z;
  delete decayp2_mass;
  delete decayp2_Ex;
  delete decayp2_A;
  delete decayp2_Z;
  delete primary_energy;
  delete primary_A;
  delete primary_Z;
  delete primary_pos;

  delete magneticfieldon;
  delete usingXsec;
  delete xsecfile;
}

// ----------------------------------------------------------------------------- //

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  // Pointer to inputs
  Inputs *Inputs = &Inputs::GetInputs();

  if (command == updateCmd)
  {
    detector->UpdateGeometry();
  }

  else if (command == target_A)
  {
    Inputs->target_A = target_A->GetNewIntValue(newValue);
  }
  else if (command == target_Z)
  {
    Inputs->target_Z = target_Z->GetNewIntValue(newValue);
  }
  else if (command == target_material)
  {
    Inputs->g4_material_name = newValue;
  }
  else if (command == target_pos)
  {
    Inputs->target_pos = target_pos->GetNew3VectorValue(newValue);
  }
  else if (command == target_width_cmd)
  {
    Inputs->width = target_width_cmd->GetNewDoubleValue(newValue);
  }
  else if (command == target_mass)
  {
    Inputs->target_mass = target_mass->GetNewDoubleValue(newValue);
  }

  else if (command == recoil_mass)
  {
    Inputs->recoil_mass = recoil_mass->GetNewDoubleValue(newValue);
  }
  else if (command == recoil_Ex)
  {
    Inputs->recoil_Ex = recoil_Ex->GetNewDoubleValue(newValue);
  }
  else if (command == recoil_A)
  {
    Inputs->recoil_A = recoil_A->GetNewIntValue(newValue);
  }
  else if (command == recoil_Z)
  {
    Inputs->recoil_Z = recoil_Z->GetNewIntValue(newValue);
  }

  else if (command == ejectile_mass)
  {
    Inputs->ejectile_mass = ejectile_mass->GetNewDoubleValue(newValue);
  }
  else if (command == ejectile_Ex)
  {
    Inputs->ejectile_Ex = ejectile_Ex->GetNewDoubleValue(newValue);
  }
  else if (command == ejectile_A)
  {
    Inputs->ejectile_A = ejectile_A->GetNewIntValue(newValue);
  }
  else if (command == ejectile_Z)
  {
    Inputs->ejectile_Z = ejectile_Z->GetNewIntValue(newValue);
  }

  else if (command == decayp1_mass)
  {
    Inputs->decayp1_mass = decayp1_mass->GetNewDoubleValue(newValue);
  }
  else if (command == decayp1_Ex)
  {
    Inputs->decayp1_Ex = decayp1_Ex->GetNewDoubleValue(newValue);
  }
  else if (command == decayp1_A)
  {
    Inputs->decayp1_A = decayp1_A->GetNewIntValue(newValue);
  }
  else if (command == decayp1_Z)
  {
    Inputs->decayp1_Z = decayp1_Z->GetNewIntValue(newValue);
  }

  else if (command == decayp2_mass)
  {
    Inputs->decayp2_mass = decayp2_mass->GetNewDoubleValue(newValue);
  }
  else if (command == decayp2_Ex)
  {
    Inputs->decayp2_Ex = decayp2_Ex->GetNewDoubleValue(newValue);
  }
  else if (command == decayp2_A)
  {
    Inputs->decayp2_A = decayp2_A->GetNewIntValue(newValue);
  }
  else if (command == decayp2_Z)
  {
    Inputs->decayp2_Z = decayp2_Z->GetNewIntValue(newValue);
  }

  else if (command == primary_energy)
  {
    Inputs->primary_energy = primary_energy->GetNewDoubleValue(newValue);
  }
  else if (command == primary_A)
  {
    Inputs->primary_A = primary_A->GetNewIntValue(newValue);
  }
  else if (command == primary_Z)
  {
    Inputs->primary_Z = primary_Z->GetNewIntValue(newValue);
  }
  else if (command == primary_pos)
  {
    Inputs->primary_pos = primary_pos->GetNew3VectorValue(newValue);
  }

  else if (command == magneticfieldon)
  {
    Inputs->using_magneticfield = magneticfieldon->GetNewBoolValue(newValue);
  }

  else if (command == usingXsec)
  {
    Inputs->using_xsection = usingXsec->GetNewBoolValue(newValue);
  }

  else if (command == xsecfile)
  {

    ifstream fXS(newValue);
		G4double theta=0, Xsecval=0;
		std::vector<double>  Theta;
    std::vector<double> Xsec;

		for (string line; getline(fXS, line);) {
			stringstream parse_xs(line);
			parse_xs >> theta >> Xsecval;
			Theta.push_back(theta);
			Xsec.push_back(Xsecval);
		}
		fXS.close();
		int v_size = Theta.size();
		float X[v_size];
		float Y[v_size];
    float Ymax = 0;
		for(int i=0; i<v_size; i++){
			X[i]=Theta.at(i); // in deg
			Y[i]=Xsec.at(i)*2.0*3.1415*sin(X[i]*3.1415/180.);
      if(Y[i] > Ymax) Ymax = Y[i];
			//cout<<X[i]<<" "<<Y[i]<<endl;
		}

		graphtable = new TGraph(v_size,X,Y);

    Inputs->xsection_graph = graphtable;
    Inputs->xmin = X[0];
    Inputs->xmax = X[v_size-1];
    Inputs->ymax = Ymax;

    //G4cout<<graphtable->GetN()<<"  "<<Ymax<<"  "<<Y[v_size-1]<<G4endl;
  }

}
