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
// $Id: A01MagneticField.hh,v 1.3 2002/12/13 11:34:29 gunter Exp $
// --------------------------------------------------------------
//
// Original version: J. Zamora 20??
// Modified (October 2018) by: D. Flechas (dcflechasg@unal.edu.co)
//                             Andre (andreserra@ymail.com)
//              

#ifndef MagneticField2_H
#define MagneticField2_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
class MagneticFieldMessenger2;


class MagneticField2 : public G4MagneticField
{
public:
  MagneticField2();
  ~MagneticField2();

  virtual void GetFieldValue( const  double Point[3],
                               double *Bfield ) const;
  
private:
  MagneticFieldMessenger2* messenger;
  G4double Bz;
  G4double rmax;
  G4double l_med;
  G4double k;
  //G4double current;
  G4double current2;
  
public:
  inline void SetCurrent(G4double val) { current2 = val; }
  inline G4double GetCurrent() const { return current2; }
    
};

#endif

