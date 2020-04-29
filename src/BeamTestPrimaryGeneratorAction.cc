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
// $Id:$
// GEANT4 tag $Name:$
//
// T. Aso        Original author
//
#include "BeamTestPrimaryGeneratorAction.hh"
#include "BeamTestPrimaryGeneratorMessenger.hh"

/*#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"*/

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "DetectorConstruction.hh"

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include <stdlib.h>
using namespace std;
#define pi 3.141592
#define N 465
#define N_inel 129
using namespace CLHEP;

//#include "kinematics.h"
//#include "kinematics2.h"

//#define N 258
//#define N_inel 396

G4ParticleGun* BeamTestPrimaryGeneratorAction::particleGun(0);

//SeccionEff* BeamTestPrimaryGeneratorAction::apuntador(0);

BeamTestPrimaryGeneratorAction::BeamTestPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  gunMessenger = new BeamTestPrimaryGeneratorMessenger(this);
}


BeamTestPrimaryGeneratorAction::~BeamTestPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}



void BeamTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  Inputs* Inputs = &Inputs::GetInputs();

  fposition = Inputs->primary_pos;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ionTable =  G4IonTable::GetIonTable();
  G4String particleName;

  particleGun->SetParticlePosition(fposition);

  G4ParticleDefinition* part = 0;
  G4double Energy = Inputs->primary_energy *MeV;

  part = ionTable->GetIon(Inputs->primary_Z,Inputs->primary_A,0.);
  particleGun->SetParticleDefinition(part);
  particleGun->SetParticleEnergy(Energy);

  G4double phi = 2*pi*G4UniformRand();
  G4double theta;
  do
  { theta = 
    // 0.066667*pi*(G4UniformRand()-0.5)
    // pi/2                                   //90°
    // pi/3                                   //60°
    // pi/6                                   //30°
    // (pi/6) + ( (pi/3) * G4UniformRand() )  //30° - 90°
    // (pi/6) + ( 0.261799 * G4UniformRand() )//30° - 45°
    // (pi/6) * G4UniformRand()               //0°  - 30°
    // (pi/48)*G4UniformRand()                //0°  - 15°
    // (pi/3) + ( (pi/6) * G4UniformRand() )  //60° - 90°
    // 0.0872665 * G4UniformRand()            //0°  - 5°
    (4 * G4UniformRand() + 2) * pi/180            //0°  - 2°

    ;}

    while(abs(theta) < 0.0349);

    // theta = 0.;
  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)));
  particleGun->GeneratePrimaryVertex(anEvent);


/*
//----------------------------------generador con densidad constante, cilindo vs gaussiana


 G4double x, y, z;
 G4double sigma_beam, radio_target;
 G4double x0, y0, z0;
 //G4double radio, sinomega;

  //seed();

	x=0.*CLHEP::mm;
	y=0.*CLHEP::mm;
	z=0.*CLHEP::mm;
	sigma_beam = 0.5*CLHEP::mm/2.35;  //alta intensidad 3mm-sigma
	radio_target = 1.0*CLHEP::mm;
	y0 = 0.*CLHEP::mm;
	z0 = 0*CLHEP::mm;
	x0 = 0.0*CLHEP::mm; //0.58
 G4double a = 0.5*CLHEP::mm, b=1.0*CLHEP::mm;


	for( ; ; ){
		//radio = radio_target*(2.*G4UniformRand() -1.);
		//sinomega = 2.*G4UniformRand() - 1.0;
		x = G4RandGauss::shoot(x0,0.002);
		y = G4RandGauss::shoot(y0,0.2);
		z = radio_target*(2.*G4UniformRand() -1.) + z0;

// Definindo a elipse usando que x²/a² + y²/b² = 1

	if( (((double)(x/a)*(x/a)) + (double)(y/b)*(y/b) <= 1) && abs(z) < (radio_target) ){break;}}

	if( abs(x) < 3*(abs(x0)+sigma_beam) && abs(y) < 3*(abs(y0)+sigma_beam) && abs(z) < (radio_target)){break;}
	}
	while(1==1);

	//x =0.;
	//y =0.;
	//z =0.;
	//particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
//--------------------------------------------------------------------------------------------
 //particleGun->GeneratePrimaryVertex(anEvent);*/
}
