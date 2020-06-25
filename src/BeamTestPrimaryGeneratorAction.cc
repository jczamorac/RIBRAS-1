// Local headers
#include "BeamTestPrimaryGeneratorAction.hh"
#include "BeamTestPrimaryGeneratorMessenger.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdlib.h>

// Geant4 headers
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "DetectorConstruction.hh"

#define pi 3.141592
#define N 465
#define N_inel 129
using namespace CLHEP;
using namespace std;

G4ParticleGun *BeamTestPrimaryGeneratorAction::particleGun(0);

// ----------------------------------------------------------------------------- //
BeamTestPrimaryGeneratorAction::BeamTestPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  gunMessenger = new BeamTestPrimaryGeneratorMessenger(this);
}

// ----------------------------------------------------------------------------- //

BeamTestPrimaryGeneratorAction::~BeamTestPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

// ----------------------------------------------------------------------------- //

void BeamTestPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  Inputs *Inputs = &Inputs::GetInputs();

  fposition = Inputs->primary_pos;
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4IonTable *ionTable = G4IonTable::GetIonTable();
  G4String particleName;

  G4ParticleDefinition *part = 0;
  G4double Energy = Inputs->primary_energy * MeV;

  part = ionTable->GetIon(Inputs->primary_Z, Inputs->primary_A, 0.);
  particleGun->SetParticleDefinition(part);
  particleGun->SetParticleEnergy(Energy);

  // ---------------------- Primary beam as a dot ------------------------- //

  /*   G4double phi = 2*pi*G4UniformRand();
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
    (1 * G4UniformRand() + 2) * pi/180            //2°  - 6°

    ;}

    while(abs(theta) < 0.0349 );

    theta = 0.00000001;
    
  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))); */

  // ----------------------------------------------------------------------- //

  //------------- Primary beam as a elipse -------------//

  G4double x, y, z;
  G4double sigma_beam, radio_target;
  G4double x0, y0, z0;

  x = 0. * mm;
  y = 0. * mm;
  z = 0. * mm;

  sigma_beam = 5 * mm / 2.35; // alta intensidad 3mm-sigma
  radio_target = 1.0 * mm;

  y0 = 0. * mm;
  x0 = 0. * mm;
  z0 = -115. *cm;

  G4double a = 0.5 * mm, b = 1.0 * mm;

  for (;;)
  {
    x = G4RandGauss::shoot(x0, sigma_beam);
    y = G4RandGauss::shoot(y0, sigma_beam);
    z = z0;

    if (abs(x) < 3 * (abs(x0) + sigma_beam) && abs(y) < 3 * (abs(y0) + sigma_beam))
    {
      break;
    }
  }
  //-----------------------------------------------------------------------------//

  G4ThreeVector sposition = G4ThreeVector(x ,y ,z);

  G4double theta = 0. , phi = 0.;

  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)));
  particleGun->SetParticlePosition(sposition);
  particleGun->GeneratePrimaryVertex(anEvent);
}
