//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: BeamTestPrimaryGeneratorMessenger.cc,v 1.5 2006/06/29 16:33:59 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 

// ----------------------------------------------------------------------------- //

#include "BeamTestPrimaryGeneratorMessenger.hh"

#include "BeamTestPrimaryGeneratorAction.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

// ----------------------------------------------------------------------------- //

BeamTestPrimaryGeneratorMessenger::BeamTestPrimaryGeneratorMessenger(BeamTestPrimaryGeneratorAction* BeamTestGun)
:BeamTestAction(BeamTestGun)
{ 
   z_tarCmd = new G4UIcmdWithADoubleAndUnit("/gun/target/z",this);
   z_tarCmd->SetGuidance("pos z of target");
   z_tarCmd->SetParameterName("Z_target",true);
   z_tarCmd->SetDefaultValue(0.);
   z_tarCmd->SetUnitCategory("Length");
   z_tarCmd->SetDefaultUnit("mm");

}

// ----------------------------------------------------------------------------- //

BeamTestPrimaryGeneratorMessenger::~BeamTestPrimaryGeneratorMessenger()
{
  delete z_tarCmd;
}

// ----------------------------------------------------------------------------- //

void BeamTestPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == z_tarCmd )
   { 	BeamTestAction->SetZ_target(z_tarCmd->GetNewDoubleValue(newValue));}
}

// ----------------------------------------------------------------------------- //

