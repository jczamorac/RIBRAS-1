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

	//G4cout<<ZOfEjectile<<" "<<AOfEjectile<<"  "<<runmassMap->GetValue(Form("%i_%i", ZOfEjectile, AOfEjectile),-1000.0)<<G4endl;

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

	//Decay particle 1
	if(Inputs->decayp1_A>1 && Inputs->decayp1_Z>1){
		G4IonTable *ionTable = G4IonTable::GetIonTable();
		G4ParticleDefinition *part = 0;
		part = ionTable->GetIon(Inputs->decayp1_Z, Inputs->decayp1_A, Inputs->decayp1_Ex);
		Inputs->DecayParticle1 = part;
	}
	else{
		if(Inputs->decayp1_A == 1 && Inputs->decayp1_Z == 1) Inputs->DecayParticle1 = G4Proton::Definition();
		else if(Inputs->decayp1_A == 1 && Inputs->decayp1_Z == 0) Inputs->DecayParticle1 = G4Neutron::Definition();
	}

	//Decay particle 2
	if(Inputs->decayp2_A>1 && Inputs->decayp2_Z>1){
		G4IonTable *ionTable = G4IonTable::GetIonTable();
		G4ParticleDefinition *part = 0;
		part = ionTable->GetIon(Inputs->decayp2_Z, Inputs->decayp2_A, Inputs->decayp2_Ex);
		Inputs->DecayParticle2 = part;
	}
	else{
		if(Inputs->decayp2_A == 1 && Inputs->decayp2_Z == 1) Inputs->DecayParticle2 = G4Proton::Definition();
		else if(Inputs->decayp2_A == 1 && Inputs->decayp2_Z == 0) Inputs->DecayParticle2 = G4Neutron::Definition();
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
												   // had to adjust the boost, as it was boosting with -beta, and it should have been +bet
	G4double ran_ex = 0;
	G4double rand_dis = 0;


	if(rDecayThisEvent){

		// MC to distribute the events with the cross section
		G4double Es = 0.136;
                do{
	                	rand_dis = 10.2*(CLHEP::RandFlat::shoot()) ;
                    ran_ex = 5*( CLHEP::RandFlat::shoot() ) + Es;

	        				}while( rand_dis > ExDistr(ran_ex, Es));
	}


	G4ThreeVector BeamDirection = BeamTrack.GetMomentumDirection();
	G4ThreeVector BeamMomentum = BeamTrack.GetMomentum();

	G4double BeamKineticEnergy = BeamTrack.GetKineticEnergy();
	rInteractionPointBeamEnergy = BeamKineticEnergy;
	G4double BeamMass = BeamTrack.GetParticleDefinition()->GetPDGMass();

	G4NistManager *nist = G4NistManager::Instance();
	G4double TargetMass = (Inputs->target_mass == 0.0)
							  //? nist->GetAtomicMass(Inputs->target_Z, Inputs->target_A)
								? 931.494*(runmassMap->GetValue(Form("%i_%i", Inputs->target_Z, Inputs->target_A),-1000.0)) - 0.511*Inputs->target_Z
							  : Inputs->target_mass;

	G4double RecoilMass = (Inputs->recoil_mass == 0.0)
								? 931.494*(runmassMap->GetValue(Form("%i_%i", Inputs->recoil_Z, Inputs->recoil_A),-1000.0)) - 0.511*Inputs->recoil_Z
							  : Inputs->recoil_mass;
	RecoilMass += Inputs->recoil_Ex;

	G4double EjectileMass = (Inputs->ejectile_mass == 0.0)
								? 931.494*(runmassMap->GetValue(Form("%i_%i", Inputs->ejectile_Z, Inputs->ejectile_A),-1000.0)) - 0.511*Inputs->ejectile_Z
								: Inputs->ejectile_mass;
	if(rDecayThisEvent) EjectileMass += ran_ex;
	else  EjectileMass += Inputs->ejectile_Ex;

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

	//G4Event *anEvent;
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

	//rEjectileTheta = acos(Vec.z() / Vec.r());

	G4double Mass = Ejectile4Vec.m();

	G4double LabKineticEnergy = Ejectile4Vec.e() - Mass;

	G4DynamicParticle *par = new G4DynamicParticle(Inputs->EjectileParticle, Vec.unit(), LabKineticEnergy);

	return par;
}

//////////////////////////////////////////////////////////////////////////////////////////


void EventAction::DecayLab4Vectors(const G4LorentzVector &ParentLV,  G4LorentzVector &DecayOut1, G4LorentzVector &DecayOut2 )
{
	Inputs *Inputs = &Inputs::GetInputs();
	G4double mass0 = ParentLV.m();
	G4double mass1 =  931.494*(runmassMap->GetValue(Form("%i_%i", Inputs->decayp1_Z, Inputs->decayp1_A),-1000.0)) - 0.511*Inputs->decayp1_Z;
	G4double mass2 =  931.494*(runmassMap->GetValue(Form("%i_%i", Inputs->decayp2_Z, Inputs->decayp2_A),-1000.0)) - 0.511*Inputs->decayp2_Z;
  G4double S = mass0*mass0;
	G4double Pcm = 0.5 * sqrt(pow(S, 2.0) + pow(mass1, 4.0) + pow(mass2, 4.0) - 2 * S * pow(mass1, 2.0) - 2 * pow(mass1, 2.0) * pow(mass2, 2.0) - 2 * S * pow(mass2, 2.0)) / sqrt(S);

	G4double ThetaCM = pi*(2.0*CLHEP::RandFlat::shoot() - 1.0);
	G4double PhiCM = CLHEP::RandFlat::shoot()*2*pi;
	G4double x = sin(ThetaCM) * cos(PhiCM);
	G4double y = sin(ThetaCM) * sin(PhiCM);
	G4double z = cos(ThetaCM);
	G4ThreeVector Decayp1DirectionInCM(-x, -y, -z);
	Decayp1DirectionInCM = Decayp1DirectionInCM * Pcm;
	G4ThreeVector Decayp2DirectionInCM(x, y, z);
	Decayp2DirectionInCM = Decayp2DirectionInCM * Pcm;
	G4double Decayp1CMEnergy = sqrt(pow(Pcm, 2.0) + pow(mass1, 2.0));
	G4double Decayp2CMEnergy = sqrt(pow(Pcm, 2.0) + pow(mass2, 2.0));

	G4LorentzVector CM_Decay1_4Vec(Decayp1CMEnergy, Decayp1DirectionInCM);
	G4LorentzVector CM_Decay2_4Vec(Decayp2CMEnergy, Decayp2DirectionInCM);
	G4double beta = ParentLV.beta();
	G4ThreeVector Vec = ParentLV.getV().unit();

	// Boost to Lab (parent direction)
	DecayOut1 = CM_Decay1_4Vec.boost(Vec,beta);
	DecayOut2 = CM_Decay2_4Vec.boost(Vec,beta);


}

G4DynamicParticle *EventAction::GetDecay1DynamicParticle(const G4LorentzVector &DecayP1_4Vec)
{	// Creating Ejectile Particle

	Inputs *Inputs = &Inputs::GetInputs();

	G4ThreeVector Vec = DecayP1_4Vec.getV().unit();

	//rEjectileTheta = acos(Vec.z() / Vec.r());

	G4double Mass = DecayP1_4Vec.m();

	G4double LabKineticEnergy = DecayP1_4Vec.e() - Mass;

	G4DynamicParticle *par = new G4DynamicParticle(Inputs->DecayParticle1, Vec.unit(), LabKineticEnergy);

	return par;
}

G4DynamicParticle *EventAction::GetDecay2DynamicParticle(const G4LorentzVector &DecayP2_4Vec)
{	// Creating Ejectile Particle

	Inputs *Inputs = &Inputs::GetInputs();

	G4ThreeVector Vec = DecayP2_4Vec.getV().unit();

	//rEjectileTheta = acos(Vec.z() / Vec.r());

	G4double Mass = DecayP2_4Vec.m();

	G4double LabKineticEnergy = DecayP2_4Vec.e() - Mass;

	G4DynamicParticle *par = new G4DynamicParticle(Inputs->DecayParticle2, Vec.unit(), LabKineticEnergy);

	return par;
}


G4double EventAction::ExDistr(G4double ex, G4double se){
	//Dissociation Energy Distribution based on  Nakamura, Physics Letters B 331 (1994) 296-301
	G4double f = pow(ex-se,1.5)/pow(ex,4.0);
	return f;
}
