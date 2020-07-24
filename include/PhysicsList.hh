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
//  advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamples
//

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ParticleTypes.hh"

class G4VPhysicsConstructor;
class StepMax;
class PhysicsListMessenger;

class PhysicsList : public G4VUserPhysicsList
{
public:
    PhysicsList();
    virtual ~PhysicsList();

    void ConstructParticle();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);
    void SetDetectorCut(G4double cut);
    void AddPhysicsPhysics();
    void ConstructProcess();
    void AddStepMax();
    void AddPackage(const G4String &name);

private:
    void ConstructEM();

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;

    G4String emName;
    PhysicsListMessenger *pMessenger;

    G4VPhysicsConstructor *fEmPhysicsList;
};

#endif
