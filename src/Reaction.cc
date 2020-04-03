#include "Reaction.hh"

#include "globals.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4GenericIon.hh"
#include <fstream>
#include <csignal>
#include "EventAction.hh"
#include "DetectorConstruction.hh"
//#include <LorentzVector.h>

using namespace CLHEP;
using namespace std;

Reaction::Reaction(const G4String& aName)
  : G4VProcess(aName){
  
  if (verboseLevel>1) {
    G4cout <<GetProcessName() << " is created "<< G4endl;
  }
  theProcessType=(G4ProcessType)6;       // Decay
  theProcessSubType=(G4ProcessType)231;  //DecayExt
}

Reaction::~Reaction(){
}

G4VParticleChange* Reaction::PostStepDoIt( const G4Track& aTrack,
					     const G4Step&){
  // Stop the current particle, if requested by G4UserLimits

  //  G4cout<<"IN DO IT *********"<<G4endl;
  aParticleChange.Initialize(aTrack);
  
  if(reaction_here){
    // Here the physics of the Reaction is defined
    reaction_here=false;                               // Set the flag so it does not do the reaction a second time
    gCESimulationManager->ThereWasAReaction();
    aParticleChange.ProposeTrackStatus(fStopAndKill);         // Kill the incident Particle

    G4double ThetaInCM =(5*CLHEP::RandFlat::shoot())*degree;  // 0 - 5 deg from z
    G4double randomPhiInCM = CLHEP::RandFlat::shoot()*2*pi;   // 0 - 2pi in transverse angle (azimuth)
    G4LorentzVector recoil4VecLab;
    G4LorentzVector ejectile4VecLab;

    // G4cout << "Theta CM: " << ThetaInCM /degree << G4endl;

//	cin.get();

    gCESimulationManager->CalculateLab4Vectors(aTrack,ThetaInCM,randomPhiInCM,recoil4VecLab,ejectile4VecLab);     

    //THIS IS THE POSITION OF THE INTERACTION
    G4ThreeVector pos(aTrack.GetPosition().getX(),aTrack.GetPosition().getY(),aTrack.GetPosition().getZ());
    
    gCESimulationManager->SetInteractionPoint(pos);

    //Adding a the recoil particle as a secondary of this reaction
    aParticleChange.AddSecondary(gCESimulationManager->GetRecoilDynamicParticle(aTrack,recoil4VecLab),pos,true); 

    //Adding a the ejectile particle as a secondary of this reaction
    aParticleChange.AddSecondary(gCESimulationManager->GetEjectileDynamicParticle(aTrack,ejectile4VecLab),pos,true); 





    // aParticleChange.AddSecondary(gCESimulationManager->GetRecoilDynamicParticle(aTrack,recoil4Vec),pos,true); //Adding a the recoil particle as a secondary of this reaction

    // G4double LabNeutronT = gamma*NeutronEnergyCM + beta*beta*gamma*gamma*ProtonMass*cos(ThetaInCM)-NeutronMass;
    // G4cout<<"Lab Neturon T "<<LabNeutronT<<G4endl;

    //    aParticleChange.AddSecondary(par,pos,true);

    // G4cout<<"TEST*************"<<G4endl;

    
    // IncidentVectorInGlobalSystem.rotateZ(-phiForIncidentVector);
    // IncidentVectorInGlobalSystem.rotateY(-thetaForIncidentVector);

    // G4cout<<"New incident vec "<<IncidentVectorInGlobalSystem<<G4endl;

    // G4cout<<"TEST*************"<<G4endl;

    // G4cout<<"Mass From LabNeutronVector "<<sqrt( pow(LabNeutron.e(),2) - vec.mag2())<<G4endl;
    // G4cout<<"lab neutron vec "<<vec<<G4endl;
    // G4cout<<"lab neutron T "<<LabNeutron.e()-NeutronMass<<G4endl;

    //    G4cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<G4endl;
    
    // // G4cout<<"Recoil Lab Kinetic Energy "<<recoil.e() -NeutronMass<<G4endl;
    // // G4cout<<"Ejectile Lab Kinetic Energy "<<ejectile.e() -N16Mass<<G4endl;
    //  int t;
    //  G4cin>>t;
    // if(G4UniformRand()<BeamOut->getTFrac())
    // 	{
    // 	  aParticleChange.SetNumberOfSecondaries(2);
    // 	  aParticleChange.AddSecondary(BeamOut->ProjectileGS(),BeamOut->CEReactionPosition(),true);
    // 	  aParticleChange.AddSecondary(BeamOut->TargetExcitation(),BeamOut->CEReactionPosition(),true);

    // //  	  aParticleChange.DumpInfo();
    // // 	  getc(stdin);
    // //  	  aParticleChange.GetSecondary(0)->GetDynamicParticle()->DumpInfo();
    // //  	  aParticleChange.GetSecondary(1)->GetDynamicParticle()->DumpInfo();
    // //  	  getc(stdin);
    // 	}
    //       else
    // 	{
    // 	  aParticleChange.SetNumberOfSecondaries(1);
    // 	  aParticleChange.AddSecondary(BeamOut->CEReactionProduct(),BeamOut->CEReactionPosition(),true);
    // 	}


  }

  return &aParticleChange;
}

G4double Reaction::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
							  G4double,
							  G4ForceCondition* condition){
  
  Inputs* Inputs = &Inputs::GetInputs();

  G4double eps= 0.1 *cm;//Small number saying that we have made it to the reaction point
                       //Ie when Zdiff < eps it is at the reaction point
  
  reaction_here=false;
  *condition=NotForced;

  G4String name = aTrack.GetVolume()->GetLogicalVolume()->GetName();

  if(name=="alvolog" && gCESimulationManager->GetWasThereACEReaction()==false){
    
    G4double ZCEReaction = /*gCESimulationManager->GetReactionZPoint()*/ Inputs->target_pos.getZ() + 301 *cm;
    // G4cout << "Z Reaction " << ZCEReaction /cm << G4endl;
    G4double ZCurrent = aTrack.GetPosition().getZ();
    // G4int ID = aTrack.GetTrackID();
    // G4cout << "ID: " << ID << G4endl;

    // G4double K = aTrack.GetKineticEnergy();
    // G4cout << "Kinect Energy: " << K << G4endl;

    // G4cout << "Z Position " << ZCurrent/cm << G4endl;
    G4double Zdiff = ZCEReaction-ZCurrent;
    // G4cout << "ZDiff: " << Zdiff/cm << G4endl;

    if(Zdiff<0){
      //If the current Z position is past the reaction point
      //don't do the reaction.  Return DBL_MAX to ensure that this
      //processes isn't called
      return DBL_MAX;
    }
    if(Zdiff>eps){
      //if we are before the reaction point but not within eps of it
      //return the distance to the reaction point as the step length
      //NOTE:  unless the trajectory of the step is directly along the z-axis
      //it will have to make several jumps to get with eps of the reaction point.

      G4ThreeVector dir=aTrack.GetMomentumDirection();

      dir*=(ZCEReaction-ZCurrent);

      return dir.mag();
    }
    if(Zdiff<= eps && Zdiff>=0){
      //particle is now withing eps of reaction point.  Do reaction.
      //returns 0. to ensure that this process is called and defines the
      //length of the step.
      reaction_here=true;
      return 0.;
    }}
    else {
    //We are not within the target so return DBL_MAX to make sure
    //process isn't invoked.
    return DBL_MAX;
    }
}

