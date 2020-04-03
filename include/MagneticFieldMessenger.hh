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
// $Id: MagneticFieldMessenger.hh,v 1.3 2002/12/13 11:34:29 gunter Exp $
// --------------------------------------------------------------
//
#ifndef MagneticFieldMessenger_h
#define MagneticFieldMessenger_h 1

class MagneticField;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class MagneticFieldMessenger: public G4UImessenger
{
  public:
    MagneticFieldMessenger(MagneticField* mpga);
    ~MagneticFieldMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    MagneticField * target;

    G4UIcmdWithADoubleAndUnit*  fieldCmd;
    G4UIcmdWithADoubleAndUnit*  fieldCmd2;

};

#endif


		


