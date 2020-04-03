//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: MagneticFieldMessenger.cc,v 1.3 2002/12/13 11:34:34 gunter Exp $
// --------------------------------------------------------------
//
#include "MagneticFieldMessenger2.hh"
#include "MagneticField2.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

MagneticFieldMessenger2::MagneticFieldMessenger2(MagneticField2 * mpga)
:target2(mpga)
{
  fieldCmd2 = new G4UIcmdWithADoubleAndUnit("/mydet/currentValue2",this);
  fieldCmd2->SetGuidance("solenoid current");
  fieldCmd2->SetParameterName("current2",true);
  fieldCmd2->SetDefaultValue(22.);
  fieldCmd2->SetDefaultUnit("tesla");
}

MagneticFieldMessenger2::~MagneticFieldMessenger2()
{
  delete fieldCmd2;
}

void MagneticFieldMessenger2::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==fieldCmd2 )
  { target2->SetCurrent(fieldCmd2->GetNewDoubleValue(newValue)); }
}

G4String MagneticFieldMessenger2::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==fieldCmd2 )
  { cv = fieldCmd2->ConvertToString(target2->GetCurrent(),""); }

  return cv;
}

