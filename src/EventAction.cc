// $Id: EventAction.cc 94 2010-01-26 13:18:30Z adotti $
/**
 * @file   EventAction.cc
 *
 * @date   17 Dec 2009
 * @author adotti
 * 
 * @brief  Implements user class EventAction.
 */

#include <fstream>
#include "EventAction.hh"
#include "RootSaver.hh"
//#include "SiDigi.hh"
#include "SiHit.hh"
//#include "SiDigitizer.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
// #include "G4DigiManager.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
//Include Geant4 Headers
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4Gamma.hh"
#include "G4DecayProducts.hh"
#include "G4GenericIon.hh"

//For debugging purposes
#include <csignal>
using namespace CLHEP;
using namespace std;

EventAction *gCESimulationManager = (EventAction *)0;

//////////////////////////////////////////////////////////////////////////////////////////
EventAction::EventAction() : rootSaver(0),
							 hitsCollName("SiHitCollection"),
							 //digitsCollName("SiDigitCollection"),
							 hitsCollID(-1)
{
	if (gCESimulationManager)
		delete gCESimulationManager;
	gCESimulationManager = this;
}
//////////////////////////////////////////////////////////////////////////////////////////
void EventAction::BeginOfEventAction(const G4Event *anEvent)
{
	Inputs *Inputs = &Inputs::GetInputs();

	if (anEvent->GetEventID() % 1000000 == 0)
	{
		G4cout << "Starting Event: " << anEvent->GetEventID() << G4endl;
	}
	//Retrieve the ID for the hit collection
	if (hitsCollID == -1)
	{
		G4SDManager *SDman = G4SDManager::GetSDMpointer();
		hitsCollID = SDman->GetCollectionID(hitsCollName);
	}

	//At the Beginning of each event set the flag to allow a Reaction to happen.
	//This ensures that the reaction happens exactly once in each event.
	//This is an issue because the Charge exchange process is added not to the beam
	//particle specifically but to an "G4GeneriIon".  In the case of 16C(p,n)16N
	//The 16N also has the CE process and will try to do this process.  Hence the
	//flag.
	rThereWasACEReactionThisEvent = false;

	//Setting particle definitions
	G4int ZOfEjectile = Inputs->ejectile_Z;
	G4int AOfEjectile = Inputs->ejectile_A;
	G4int ZOfRecoil = Inputs->recoil_Z;
	G4int AOfRecoil = Inputs->recoil_A;

	G4float ExOfEjectile = Inputs->ejectile_Ex;
	G4float ExOfRecoil = Inputs->recoil_Ex;

	//Check to see if the ejectile is a proton, neutron or ion
	//If the Ejectile is a proton
	if (ZOfEjectile == 1 && AOfEjectile == 1)
	{
		Inputs->EjectileParticle = G4Proton::Definition();
	}

	//If the Ejectile is a neutron
	else if (ZOfEjectile == 0 && AOfEjectile == 1)
	{
		Inputs->EjectileParticle = G4Neutron::Definition();
	}

	//If the Ejectile is a ion
	else
	{
		G4IonTable *ionTable = G4IonTable::GetIonTable();
		G4ParticleDefinition *part = 0;
		part = ionTable->GetIon(Inputs->ejectile_Z, Inputs->ejectile_A, Inputs->ejectile_Ex);
		Inputs->EjectileParticle = part;
	}

	//Check to see if the Recoil is a proton, neutron or ion
	//If the Recoil is a Proton
	if (ZOfRecoil == 1 && AOfRecoil == 1)
	{
		Inputs->RecoilParticle = G4Proton::Definition();
	}

	//If the Recoil is a neutron
	else if (ZOfRecoil == 0 && AOfRecoil == 1)
	{
		Inputs->RecoilParticle = G4Neutron::Definition();
	}

	//If the Recoil is a ion
	else
	{
		G4IonTable *ionTable = G4IonTable::GetIonTable();
		G4ParticleDefinition *part = 0;
		part = ionTable->GetIon(Inputs->recoil_Z, Inputs->recoil_A, Inputs->ejectile_Ex);
		Inputs->RecoilParticle = part;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
void EventAction::EndOfEventAction(const G4Event *anEvent)
{
	//Store information
	if (rootSaver)
	{
		//Retrieve digits collection
		//G4int digiCollID = digiManager->GetDigiCollectionID( digitsCollName );
		//const SiDigiCollection* digits = static_cast<const SiDigiCollection*>( digiManager->GetDigiCollection(digiCollID) );
		//Retrieve hits collections
		G4HCofThisEvent *hitsCollections = anEvent->GetHCofThisEvent();
		SiHitCollection *hits = 0;
		if (hitsCollections)
		{
			hits = static_cast<SiHitCollection *>(hitsCollections->GetHC(hitsCollID));
		}
		//Get Postion and Momentum of primary
		//This is needed to store in ntuple info @ z=0
		const G4ThreeVector &pos = anEvent->GetPrimaryVertex()->GetPosition();
		const G4ThreeVector &mom = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
		rootSaver->AddEvent(hits, pos, mom);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
void EventAction::CalculateLab4Vectors(const G4Track &BeamTrack,
									   G4double CMScatteringAngle,
									   G4double phi,
									   G4LorentzVector &RecoilOut,
									   G4LorentzVector &EjectileOut)
{

	// Retrieve inputs
	Inputs *Inputs = &Inputs::GetInputs();

	// Choose a random phi between 0 and 2 pi
	G4double randomPhiInCM = phi;
	G4double x = sin(CMScatteringAngle) * cos(randomPhiInCM);
	G4double y = sin(CMScatteringAngle) * sin(randomPhiInCM);
	G4double z = cos(CMScatteringAngle);
	G4ThreeVector RecoilDirectionInCM(-x, -y, -z); // Updated these so that the ejectile is moving forward
	G4ThreeVector EjectileDirectionInCM(x, y, z);  // Before the Ejectile was backward focused. I also
												   // had to adjust the boost, as it was boosting with -beta, and it should have been +beta

	G4ThreeVector BeamDirection = BeamTrack.GetMomentumDirection();
	G4ThreeVector BeamMomentum = BeamTrack.GetMomentum();

	G4double BeamKineticEnergy = BeamTrack.GetKineticEnergy();
	rInteractionPointBeamEnergy = BeamKineticEnergy;
	G4double BeamMass = BeamTrack.GetParticleDefinition()->GetPDGMass();

	G4NistManager *nist = G4NistManager::Instance();
	G4double TargetMass = (Inputs->target_mass == 0.0)
							  ? nist->GetAtomicMass(Inputs->target_Z, Inputs->target_A)
							  : Inputs->target_mass;

	G4double RecoilMass = (Inputs->recoil_mass == 0.0)
							  ? getMass(Inputs->recoil_Z, Inputs->recoil_A) 	  // Atomic mass (without electrons) 
							  : Inputs->recoil_mass;
	RecoilMass += Inputs->recoil_Ex;

	G4double EjectileMass = (Inputs->ejectile_mass == 0.0)
								? getMass(Inputs->ejectile_Z, Inputs->ejectile_A) // Isotopic mass (without electrons) */
								: Inputs->ejectile_mass;
	EjectileMass += Inputs->ejectile_Ex;

	if (!TargetMass || !EjectileMass || !RecoilMass)
	{
		G4cout << G4endl << "CE: Found zero mass in Target/Ejectile/Recoil" << G4endl << G4endl;
		G4Exception("CESimulationManager::CalculateLab4Vectors", 0, FatalErrorInArgument, "Error in masses");
	}

	// Kinematics in Center of Momentum frame
	G4double ELabBeam = BeamKineticEnergy + BeamMass;
	G4double ELabTarget = TargetMass;

	//  get COM parameters
	G4double beta_cm = BeamMomentum.z() / (ELabBeam + TargetMass);
	G4double gamma_cm = 1.0 / sqrt(1.0 - pow(beta_cm, 2.0));
	G4double S = 2. * ELabBeam * TargetMass + pow(BeamMass, 2.0) + pow(TargetMass, 2.0);
	G4double Pcm = 0.5 * sqrt(pow(S, 2.0) + pow(RecoilMass, 4.0) + pow(EjectileMass, 4.0) - 2 * S * pow(RecoilMass, 2.0) - 2 * pow(RecoilMass, 2.0) * pow(EjectileMass, 2.0) - 2 * S * pow(EjectileMass, 2.0)) / sqrt(S);

	//generate com  momenta (for now isotropic)
	RecoilDirectionInCM = RecoilDirectionInCM * Pcm;
	G4double RecoilCMEnergy = sqrt(pow(Pcm, 2.0) + pow(RecoilMass, 2.0));
	EjectileDirectionInCM = EjectileDirectionInCM * Pcm;
	G4double EjectileCMEnergy = sqrt(pow(Pcm, 2.0) + pow(EjectileMass, 2.0));
	G4LorentzVector CMRecoil4Vec(RecoilCMEnergy, RecoilDirectionInCM);
	G4LorentzVector CMEjectile4Vec(EjectileCMEnergy, EjectileDirectionInCM);

	// Boost to Lab frame
	RecoilOut = CMRecoil4Vec.boostZ(beta_cm);
	EjectileOut = CMEjectile4Vec.boostZ(beta_cm);

	// rotate back to the beam axis
	RecoilOut = RecoilOut.rotateUz(BeamDirection);
	EjectileOut = EjectileOut.rotateUz(BeamDirection);

	/*
	rRecoilKineticEnergy = RecoilOut.e() - RecoilOut.m();
	G4ThreeVector Vec = RecoilOut.getV();
	G4double recoangle = acos(Vec.z()/Vec.r());

	G4double rEjeKineticEnergy = EjectileOut.e() - EjectileOut.m();
	G4ThreeVector Vec2 = EjectileOut.getV();
	G4double ejeangle = acos(Vec2.z()/Vec2.r());
	G4double bangle = acos(BeamDirection.z()/BeamDirection.r());


	//G4cout<<recoangle/deg<<"  "<<ejeangle/deg<<"  "<<bangle/deg<<"  "<<CMScatteringAngle/deg<<G4endl;
	G4cout<<Vec2<<"  "<<EjectileOut<<G4endl;
	*/
}

//////////////////////////////////////////////////////////////////////////////////////////
G4DynamicParticle *EventAction::GetRecoilDynamicParticle(const G4Track &BeamTrack, const G4LorentzVector &Recoil4Vec)
{	// Creating Recoil Particle

	G4Event *anEvent;
	Inputs *Inputs = &Inputs::GetInputs();

	G4ThreeVector Vec = Recoil4Vec.getV().unit();

	rRecoilTheta = acos(Vec.z() / Vec.r());

	G4double Mass = Recoil4Vec.m();

	G4double LabKineticEnergy = Recoil4Vec.e() - Mass;

	G4DynamicParticle *par = new G4DynamicParticle(Inputs->RecoilParticle, Vec.unit(), LabKineticEnergy);

	return par;
}

//////////////////////////////////////////////////////////////////////////////////////////
G4DynamicParticle *EventAction::GetEjectileDynamicParticle(const G4Track &BeamTrack, const G4LorentzVector &Ejectile4Vec)
{	// Creating Ejectile Particle

	Inputs *Inputs = &Inputs::GetInputs();

	G4ThreeVector Vec = Ejectile4Vec.getV().unit();

	rEjectileTheta = acos(Vec.z() / Vec.r());

	G4double Mass = Ejectile4Vec.m();

	G4double LabKineticEnergy = Ejectile4Vec.e() - Mass;

	G4DynamicParticle *par = new G4DynamicParticle(Inputs->EjectileParticle, Vec.unit(), LabKineticEnergy);

	return par;
}

//////////////////////////////////////////////////////////////////////////////////////////
