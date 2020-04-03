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
// $Id: A01MagneticField.cc,v 1.4 2003/10/12 14:08:14 asaim Exp $
// --------------------------------------------------------------
//
// Original version: J. Zamora 20??
// Modified (October 2018) by: D. Flechas (dcflechasg@unal.edu.co)
//                             Andre (andreserra@ymail.com)
//

#include "MagneticField2.hh"
#include "MagneticFieldMessenger2.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

/*  units and constants */
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

MagneticField2::MagneticField2()
{
   messenger = new MagneticFieldMessenger2(this);
  k = 0.04078336599; //constant of RIBRAS solenoids 
  rmax = 19.2034*cm; // coil radiusBz
  l_med = 33.9979*cm; // half coil length
  current2 = 20.0*tesla; //electric current [Amp] (tesla units are included just for B)
  //current2 = 5.*tesla; //electric current [Amp] (tesla units are included just for B)

  
}

MagneticField2::~MagneticField2()
{
  delete messenger;
}

void MagneticField2::GetFieldValue(const double Point[3],double *Bfield) const
{
  

  if(Point[2]>191.*cm && Point[2]<411.*cm) 
  {
    
      G4double factor = 0.945;
      G4double r = sqrt(Point[0]*Point[0]+Point[1]*Point[1]);	
      G4double z_mas = (Point[2]-301*cm)+l_med;
      G4double z_menos = (Point[2]-301*cm)-l_med;
      G4double fr_mas = 1.0/(pow(z_mas*z_mas+rmax*rmax,1.5));
      G4double fr_menos = 1.0/(pow(z_menos*z_menos+rmax*rmax,1.5));
      G4double fz_mas = (z_mas)/(sqrt(z_mas*z_mas+rmax*rmax));
      G4double fz_menos = (z_menos)/(sqrt(z_menos*z_menos+rmax*rmax));
      G4double fr2_mas = (4*z_mas*z_mas-rmax*rmax)/(pow(z_mas*z_mas+rmax*rmax,3.5));
      G4double fr2_menos = (4*z_menos*z_menos-rmax*rmax)/(pow(z_menos*z_menos+rmax*rmax,3.5));
      G4double fz2_mas = (z_mas)/(pow(z_menos*z_menos+rmax*rmax,2.5));
      G4double fz2_menos = (z_menos)/(pow(z_mas*z_mas+rmax*rmax,2.5));
      G4double br = -0.5*r*k*current2*factor*rmax*rmax*(fr_menos-fr_mas)+0.375*r*r*r*rmax*rmax*k*current2*factor*(fr2_menos-fr2_mas);
  

      Bfield[2] = k*current2*factor*(fz_menos-fz_mas)-0.75*r*r*rmax*rmax*k*current2*factor*(fz2_mas-fz2_menos); 
      Bfield[0] = br*Point[0]/r;
      Bfield[1] = br*Point[1]/r;

      //G4cout <<Point[2]/cm<<"  "<< Bfield[2]<<G4endl;	
  }

  else
  { 
      Bfield[2] = 0.; 
      Bfield[0] = 0.;
      Bfield[1] = 0.;
      // G4cout << Bz<<G4endl;
  }

}

