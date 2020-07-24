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

#include "DetectorConstruction.hh"
#include "MagneticField.hh"
#include "MagneticFieldMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4Step.hh"

/*  units and constants */
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

using namespace CLHEP;
using namespace std;

// ----------------------------------------------------------------------------- //

MagneticField::MagneticField()
{
  messenger = new MagneticFieldMessenger(this);
  k = 0.04078336599;   //constant of RIBRAS solenoids
  rmax = 19.2034 * cm; // coil radius
  rmin = 0.1 * mm;
  l_med = 33.9979 * cm;   // half coil length
  current = 37. * tesla;  //electric current [Amp] (tesla units are included just for B)
  current2 = 22. * tesla; //electric current [Amp] (tesla units are included just for B)
}

// ----------------------------------------------------------------------------- //

MagneticField::~MagneticField()
{
  delete messenger;
}

// ----------------------------------------------------------------------------- //

void MagneticField::InicializaMag()
{
}

// ----------------------------------------------------------------------------- //

void MagneticField::GetFieldValue(const double Point[3], double *Bfield) const
{
  Inputs *Inputs = &Inputs::GetInputs();
  if (Point[2] > -110. * cm && Point[2] < 430. * cm && sqrt(Point[0] * Point[0] + Point[1] * Point[1]) < rmax && Inputs->using_magneticfield)
  {
      G4double z_mas = Point[2]+l_med;
      G4double z_menos = Point[2]-l_med;
      G4double r = sqrt(Point[0]*Point[0]+Point[1]*Point[1]);
      G4double fr_mas = 1.0/(pow(z_mas*z_mas+rmax*rmax,1.5));
      G4double fr_menos = 1.0/(pow(z_menos*z_menos+rmax*rmax,1.5));
      G4double fz_mas = (z_mas)/(sqrt(z_mas*z_mas+rmax*rmax));
      G4double fz_menos = (z_menos)/(sqrt(z_menos*z_menos+rmax*rmax));
      G4double fr2_mas = (4*z_mas*z_mas-rmax*rmax)/(pow(z_mas*z_mas+rmax*rmax,3.5));
      G4double fr2_menos = (4*z_menos*z_menos-rmax*rmax)/(pow(z_menos*z_menos+rmax*rmax,3.5));
      G4double fz2_mas = (z_mas)/(pow(z_menos*z_menos+rmax*rmax,2.5));
      G4double fz2_menos = (z_menos)/(pow(z_mas*z_mas+rmax*rmax,2.5));
      G4double factor = 0.945;
      G4double br = -0.5*r*k*current*factor*rmax*rmax*(fr_menos-fr_mas)+0.375*r*r*r*rmax*rmax*k*current*factor*(fr2_menos-fr2_mas);

      Bfield[2] = k*current*factor*(fz_menos-fz_mas)-0.75*r*r*rmax*rmax*k*current*factor*(fz2_mas-fz2_menos);
      Bfield[0] = br*Point[0]/r;
      Bfield[1] = br*Point[1]/r;

    if (Point[2] > 181. * cm)
    {
      r = sqrt(Point[0] * Point[0] + Point[1] * Point[1]);
      z_mas = (Point[2] - (301) * cm) + l_med;
      z_menos = (Point[2] - (301) * cm) - l_med;
      fr_mas = 1.0 / (pow(z_mas * z_mas + rmax * rmax, 1.5));
      fr_menos = 1.0 / (pow(z_menos * z_menos + rmax * rmax, 1.5));
      fz_mas = (z_mas) / (sqrt(z_mas * z_mas + rmax * rmax));
      fz_menos = (z_menos) / (sqrt(z_menos * z_menos + rmax * rmax));
      fr2_mas = (4 * z_mas * z_mas - rmax * rmax) / (pow(z_mas * z_mas + rmax * rmax, 3.5));
      fr2_menos = (4 * z_menos * z_menos - rmax * rmax) / (pow(z_menos * z_menos + rmax * rmax, 3.5));
      fz2_mas = (z_mas) / (pow(z_menos * z_menos + rmax * rmax, 2.5));
      fz2_menos = (z_menos) / (pow(z_mas * z_mas + rmax * rmax, 2.5));
      br = -0.5 * r * k * current2 * factor * rmax * rmax * (fr_menos - fr_mas) + 0.375 * r * r * r * rmax * rmax * k * current2 * factor * (fr2_menos - fr2_mas);

      Bfield[2] = Bfield[2] + k * current2 * factor * (fz_menos - fz_mas) - 0.75 * r * r * rmax * rmax * k * current2 * factor * (fz2_mas - fz2_menos);
      Bfield[0] = Bfield[0] + br * Point[0] / r;
      Bfield[1] = Bfield[1] + br * Point[1] / r;
    }
  }
  else
  {
    Bfield[2] = 0.;
    Bfield[0] = 0.;
    Bfield[1] = 0.;
  }
}
// ----------------------------------------------------------------------------- //