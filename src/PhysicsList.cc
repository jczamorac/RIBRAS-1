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

//Local
#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "StepMax.hh"
#include "DetectorConstruction.hh"
#include <fstream>

//General
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4AutoDelete.hh"
#include "G4SystemOfUnits.hh"

//Hadron Physics
#include "G4HadronPhysicsINCLXX.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"

//Electromagnetic Physics
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmProcessOptions.hh"
#include "G4EmExtraPhysics.hh"

//Particles
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

using namespace std;

// ----------------------------------------------------------------------------- //

PhysicsList::PhysicsList() : G4VUserPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1. * mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForPositron = defaultCutValue;

  pMessenger = new PhysicsListMessenger(this);

  SetVerboseLevel(0);

  fEmPhysicsList = new G4EmStandardPhysics_option4();
}

// ----------------------------------------------------------------------------- //

PhysicsList::~PhysicsList()
{
}

// ----------------------------------------------------------------------------- //

void PhysicsList::ConstructParticle()
{
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor sLivedConstructor;
  sLivedConstructor.ConstructParticle();
}

// ----------------------------------------------------------------------------- //

// Geant4 classes need for adding physics processes to the simulation
#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "Reaction.hh"
#include "G4ionIonisation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ParticleTable.hh"

/////////Neutron Stuff////////////
// Neutron DATA
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4Neutron.hh"

// Processes for the neutrons
#include "G4HadronElasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"

// Models
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPInelastic.hh"

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  AddStepMax();

  // // For example Elastic scattering below 20 MeV
  // G4HadronElasticProcess* theNeutronElasticProcess = new G4HadronElasticProcess();

  // // Cross Section Data set
  // G4NeutronHPElasticData* theHPElasticData = new G4NeutronHPElasticData();
  // // Setting verbose 0
  // theHPElasticData->SetVerboseLevel(0);
  // /* theHPElasticData->DumpPhysicsTable(*G4Neutron::NeutronDefinition()); */
  // theNeutronElasticProcess->AddDataSet( theHPElasticData );

  // // Model
  // G4NeutronHPElastic* theNeutronElasticModel = new G4NeutronHPElastic();
  // theNeutronElasticModel->SetVerboseLevel(0);
  // theNeutronElasticProcess->RegisterMe(theNeutronElasticModel);

  // G4NeutronInelasticProcess * theNeutronInelasticProcess = new G4NeutronInelasticProcess();
  // G4NeutronHPInelasticData* theHPInelasticData = new G4NeutronHPInelasticData();
  // theHPInelasticData->SetVerboseLevel(0);
  // //  theHPInelasticData->DumpPhysicsTable(*G4Neutron::Neutron());
  // theNeutronInelasticProcess->AddDataSet( theHPInelasticData );

  // G4NeutronHPInelastic* theNeutronInelasticModel = new G4NeutronHPInelastic();
  // theNeutronInelasticModel->SetVerboseLevel(0);
  // theNeutronInelasticProcess->RegisterMe(theNeutronInelasticModel);

  // G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
  // pmanager->SetVerboseLevel(0);
  // pmanager->AddDiscreteProcess( theNeutronElasticProcess );
  // pmanager->AddDiscreteProcess( theNeutronInelasticProcess );

  return;
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::ConstructEM()
{
  fEmPhysicsList->ConstructProcess();

  auto theParticleIterator = GetParticleIterator();
  G4PhysicsListHelper *ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->CheckParticleList();
  ph->SetVerboseLevel(0);

  theParticleIterator->reset();
  while ((*theParticleIterator)())
  {

    G4ParticleDefinition *particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "GenericIon")
    { //Define physics for ions
      ph->RegisterProcess(new Reaction, particle);
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4ionIonisation, particle);
    }
    // else if (particleName == "gamma") {
    //   // gamma
    //   ph->RegisterProcess(new G4PhotoElectricEffect, particle);
    //   ph->RegisterProcess(new G4ComptonScattering, particle);
    //   ph->RegisterProcess(new G4GammaConversion, particle);

    // }
    // else if (particleName == "e-") {
    //   //electron
    //   ph->RegisterProcess(new G4eMultipleScattering, particle);
    //   ph->RegisterProcess(new G4eIonisation, particle);
    //   ph->RegisterProcess(new G4eBremsstrahlung, particle);

    // }
    // else if (particleName == "e+") {
    //   //positron
    //   ph->RegisterProcess(new G4eMultipleScattering, particle);
    //   ph->RegisterProcess(new G4eIonisation, particle);
    //   ph->RegisterProcess(new G4eBremsstrahlung, particle);
    //   ph->RegisterProcess(new G4eplusAnnihilation, particle);

    // }
    // else if (particleName == "proton" ||
    //          particleName == "pi-" ||
    //          particleName == "pi+")
    // {
    //   //proton
    //   ph->RegisterProcess(new Reaction, particle);
    //   ph->RegisterProcess(new G4hMultipleScattering, particle);
    //   ph->RegisterProcess(new G4hIonisation, particle);
    //   ph->RegisterProcess(new G4hBremsstrahlung, particle);
    //   ph->RegisterProcess(new G4hPairProduction, particle);
    // }
  }
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  // This process must exist in all threads.

  StepMax *stepMaxProcess = new StepMax();
  G4AutoDelete::Register(stepMaxProcess);

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
  {
    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();

    if (stepMaxProcess->IsApplicable(*particle) && pmanager)
    {
      pmanager->AddDiscreteProcess(stepMaxProcess);
    }
  }
}
