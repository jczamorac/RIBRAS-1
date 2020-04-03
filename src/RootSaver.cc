 // $Id: RootSaver.cc 94 2010-01-26 13:18:30Z adotti $
 /**
  * @file   RootSaver.cc
  *
  * @date   17 Dec 2009
  * @author adotti
  * 
  * @brief  Implements class RootSaver.
  */
 
 #include "DetectorConstruction.hh"
 #include "RootSaver.hh"
 //#include "SiDigi.hh"
 #include "SiHit.hh"
 #include "TTree.h"
 #include "TFile.h"
 #include "TMath.h"
 #include "TRandom3.h"
 
 #include <fstream>
 #include <sstream>
 #include <iostream>
 #include <cassert>

 using namespace std;
 
 RootSaver::RootSaver() :
         rootTree(0),
         runCounter(0),
         nStrips(0),
         Signal0(0),
         Signal1(0),
         Signal2(0),	 
         Pos_x_det0(0),
         Pos_y_det0(0),
         Pos_z_det0(0),
         Pos_x_det1(0),
         Pos_y_det1(0),
         Pos_z_det1(0),
         Pos_x_det2(0),
	 Pos_y_det2(0),
	 Pos_z_det2(0),	 
         E_det0(0),
         E_det1(0),
         E_det2(0),	 
         TruthPosx(0),
	 TruthPosy(0),
	 TruthPosz(0),
         TruthAngle_theta(0),
	 TruthAngle_phi(0),
	 Px_dssd(0),
	 Py_dssd(0),
	 Pz_dssd(0),    
	 T_dssd(0),
	 T_sili1(0),
	 T_sili2(0),
	 T_dssd2(0),
	 Ekin_dssd2(0)
               
 {
 }
 
 RootSaver::~RootSaver()
 {
         //Close current file if needed
         if ( rootTree )
         {
                 CloseTree();
         }
 }
 
 void RootSaver::CreateTree( const std::string& fileName , const std::string& treeName )
 {
         if ( rootTree )
         {
                 std::cerr<<"TTree already created, first call CloseTree"<<std::endl;
                 return;
         }
         std::ostringstream fn;
         //fn << fileName << "_run" << runCounter++ << ".root";
	 fn << fileName << "_run_" << runCounter++ << ".root";
         //Create a new file and open it for writing, if the file already exists the file
         //is overwritten
         TFile* rootFile = TFile::Open( fn.str().data() , "recreate" );
         if ( rootFile == 0 || rootFile->IsZombie() )
         {
                 G4cerr<<"Error opening the file: "<<fn.str() <<" TTree will not be saved."<<G4endl;
                 return;
         }
         rootTree = new TTree( treeName.data() , treeName.data() );
         //TODO: Get detector strip numbers automatically
         nStrips = 128;
         Signal0 = new Float_t[nStrips];
         Signal1 = new Float_t[nStrips];
         Signal2 = new Float_t[nStrips];
	 
	 
         for ( Int_t strip = 0 ; strip < nStrips ; ++strip )
         {
                 Signal0[strip] = 0;
                 Signal1[strip] = 0;
                 Signal2[strip] = 0;
		 
         }
         //Digits variables
         rootTree->Branch( "E0", Signal0 ,   "E0[128]/F");
         rootTree->Branch( "E1", Signal1 , "E1[128]/F");
         rootTree->Branch( "E2", Signal2 ,   "E2[128]/F");
	 
         //Hits variables
         rootTree->Branch( "pos_x_det0" , &Pos_x_det0 );
	 rootTree->Branch( "pos_y_det0" , &Pos_y_det0 );
	 rootTree->Branch( "pos_z_det0" , &Pos_z_det0 );
         rootTree->Branch( "pos_x_det1" , &Pos_x_det1 );
	 rootTree->Branch( "pos_y_det1" , &Pos_y_det1 );
	 rootTree->Branch( "pos_z_det1" , &Pos_z_det1 );
         rootTree->Branch( "pos_x_det2" , &Pos_x_det2 );
	 rootTree->Branch( "pos_y_det2" , &Pos_y_det2 );
	 rootTree->Branch( "pos_z_det2" , &Pos_z_det2 );
	 
         rootTree->Branch( "ener_det0" , &E_det0 );
         rootTree->Branch( "ener_det1" , &E_det1 );
         rootTree->Branch( "ener_det2" , &E_det2 );
	 
         rootTree->Branch( "truthPosx" , &TruthPosx );
	 rootTree->Branch( "truthPosy" , &TruthPosy );
	 rootTree->Branch( "truthPosz" , &TruthPosz );
         rootTree->Branch( "truthAngle_theta" , &TruthAngle_theta );
	 rootTree->Branch( "truthAngle_phi" , &TruthAngle_phi );
	 rootTree->Branch( "px_dssd" , &Px_dssd );
	 rootTree->Branch( "py_dssd" , &Py_dssd );
	 rootTree->Branch( "pz_dssd" , &Pz_dssd );
	 rootTree->Branch( "t_dssd" , &T_dssd );
	 rootTree->Branch( "t_sili1" , &T_sili1 );
	 rootTree->Branch( "t_sili2" , &T_sili2 );
	 rootTree->Branch( "t_dssd2" , &T_dssd2 );
	 rootTree->Branch( "Ekin_dssd2" , &Ekin_dssd2 );
	 
 }	
 
 void RootSaver::CloseTree()
 {
         //Check if ROOT TTree exists,
         //in case get the associated file and close it.
         //Note that if a TFile goes above 2GB a new file
         //will be automatically opened. We have thus to get,
         //from the TTree the current opened file
         if ( rootTree )
         {
                 G4cout<<"Writing ROOT TTree: "<<rootTree->GetName()<<G4endl;
                 //rootTree->Print();
                 rootTree->Write();
                 TFile* currentFile = rootTree->GetCurrentFile();
                 if ( currentFile == 0 || currentFile->IsZombie() )
                 {
                         G4cerr<<"Error closing TFile "<<G4endl;
                         return;
                 }
                 currentFile->Close();
                 //The root is automatically deleted.
                 rootTree = 0;
                 delete[] Signal0;
                 delete[] Signal1;
                 delete[] Signal2;
		 
         }
 }
 
 void RootSaver::AddEvent( const SiHitCollection* const hits,
                                                   const G4ThreeVector& primPos, const G4ThreeVector& primMom )
 {
         //If root TTree is not created ends
         if ( rootTree == 0 )
         {
                 return;
         }
        
 
         //Store Hits infromation
         if ( hits )
         {
                 G4int nHits = hits->entries();
                 //Set defaults
                 E_det0 = 0;
                 E_det1 = 0;
                 E_det2 = 0;
		 
		 Ekin_dssd2 = 0;
                
		 Pos_x_det0 = -1000;
                 Pos_x_det1 = -1000;
                 Pos_x_det2 = -1000;
		 

		
	 	Pos_y_det0 = -1000;
	 	Pos_z_det0 = -1000;
         
	        Pos_y_det1 = -1000;
	 	Pos_z_det1 = -1000;

		
		Pos_y_det2 = -1000;
	 	Pos_z_det2 = -1000;
		
         
	 	Px_dssd = -1000;
	 	Py_dssd = -1000;
	 	Pz_dssd = -1000;    
	 	T_dssd = -1000;
	 	T_sili1 = -1000;
	 	T_sili2 = -1000;
		T_dssd2 = -1000;

	
                 //Loop on all hits, consider only the hits with isPrimary flag
                 //Position is weighted average of hit x()
                 for ( G4int h = 0 ; (h<nHits) ; ++h )
                 {
                         const SiHit* hit = static_cast<const SiHit*>( hits->GetHit( h ) );
                         //Uncomment this line if you want to record only
                         //primary energy depositions
                         //if ( hit->GetIsPrimary() == false ) continue;
			 G4int stripNum = hit->GetStripNumber();
                         G4int planeNum = hit->GetPlaneNumber();
                         G4ThreeVector pos = hit->GetPosition();
			 G4ThreeVector momentum = hit->GetIncidenceMomentumDirection();
                         G4double x = pos.x();
			 G4double y = pos.y();
			 G4double z = pos.z();
	
			 G4double momet_x = momentum.x();
			 G4double momet_y = momentum.y();
			 G4double momet_z = momentum.z();

			 G4double tiempo = hit->GetIncidenceTime();
			 G4double ekine = hit->GetIncidenceKineticEnergy();

                         //We save xyz in CLHEP::mm (detector coordinates)
                         x /= CLHEP::mm;
			 y /= CLHEP::mm;
			 z /= CLHEP::mm;
			 tiempo /= CLHEP::ns;
			 ekine /= CLHEP::MeV;
                         //We save energy in CLHEP::MeV
                         Float_t edep = static_cast<Float_t>(hit->GetEdep());
                         edep /= CLHEP::MeV;

			 G4cout<<"Strip "<<stripNum<<" detector "<<planeNum<<" x "<<x<<"  y  "<<y<<"  hits  "<<nHits<<G4endl;

                         if ( planeNum == 0 )
                         {
                                 if ( hit->GetIsPrimary() == true ){ 
							Pos_x_det0 = x;
							Pos_y_det0 = y;
							Pos_z_det0 = z;
							T_sili1 = tiempo;
							}
				 double Econv = Digital(edep);
                                 E_det0 += Econv;
				 Signal0[stripNum] += Econv;
                         }
                         else if ( planeNum == 1)
                         {
                                 if ( hit->GetIsPrimary() == true ){ 
							Pos_x_det1 = x;
							Pos_y_det1 = y;
							Pos_z_det1 = z;
							Px_dssd = momet_x;
							Py_dssd = momet_y;
							Pz_dssd = momet_z;
							T_dssd = tiempo;
							//G4cout<<"energia:  "<<E_dssd<<G4endl;
							
							}
                                 double Econv = Digital(edep);
                                 E_det1 += Econv;
				 Signal1[stripNum] += Econv;
                         }
	                         else if ( planeNum == 2)
                         {
                                 if ( hit->GetIsPrimary() == true ){
					 		Pos_x_det2 = x;
							Pos_y_det2 = y;
							Pos_z_det2 = z;
							T_sili2 = tiempo;
							}
                                 double Econv = Digital(edep);
                                 E_det2 += Econv;
				 Signal2[stripNum] += Econv;
                         }

			   
                         else
                         {
                                 G4cerr<<"Hit Error: Plane number "<<planeNum<<" expected max value: 9"<<G4endl;
                                 continue;
                         }
 
                 }
         }
         else
         {
                 G4cerr<<"Error: No hits collection passed to RootSaver"<<G4endl;
         }
         TruthPosx = static_cast<Float_t>( primPos.x() );
	 TruthPosy = static_cast<Float_t>( primPos.y() );
	 TruthPosz = static_cast<Float_t>( primPos.z() );
         //Measure angle of the beam in xz plane measured from z+ direction
         // -pi<Angle<=pi (positive when close to x positiove direction)
         Float_t sign_z = ( primMom.z()>= 0 ) ? +1 : -1;
         Float_t sign_x = ( primMom.x()>= 0 ) ? +1 : -1;
         TruthAngle_theta = ( primMom.z() != 0 ) ?
                         TMath::PiOver2()*sign_x*(1-sign_z)+std::atan( primMom.x()/primMom.z() )
                         : sign_x*TMath::PiOver2(); //beam perpendicular to z
         TruthAngle_theta /= CLHEP::deg;

	  //Measure angle of the beam in xy. Phi
         // -pi/2<Phi<=pi/2 
	
	 Float_t sign_y = ( primMom.y()>= 0 ) ? +1 : -1;
	 TruthAngle_phi = ( primMom.x() != 0 ) ?
                         sign_y*std::abs(atan( primMom.y()/primMom.x()) )
                         : sign_y*TMath::PiOver2(); //beam perpendicular to x
         TruthAngle_phi /= CLHEP::deg;
	//G4cout<<"x mom "<<primMom.x()<<" y mom "<<primMom.y()<<G4endl;
	//if(TruthAngle_phi>7 && E_dssd2 > 0)G4cout<<"atan "<<std::abs(atan( primMom.y()/primMom.x() ))<<" en degree "<<TruthAngle_phi<<G4endl;
	//if(TruthAngle_phi>9 && E_dssd2 > 2)G4cout<<"x mom "<<primMom.x()<<" y mom "<<primMom.y()<<" z mom "<< primMom.z()<<" pos x  "<<TruthPosx<<" pos y  "<<TruthPosy<<" pos z "<<TruthPosz<<" posy "<< Pos_y_dssd2<<G4endl;

         rootTree->Fill();
 }


double RootSaver::Digital(double Eraw){

	//TRandom3* myRandom = new TRandom3();
	double ion_pot = 3.6e-6; //eV
	double fanofactor = 0.1;
	double sigmaElnoise = 11.75e-3; //keV
	int Nelectrons = floor(Eraw/ion_pot);
	Nelectrons = gRandom -> Gaus(Nelectrons, sqrt(Nelectrons*fanofactor));
	double energypoint = Nelectrons*ion_pot;
	energypoint = gRandom -> Gaus(energypoint, sigmaElnoise); //electronic noise

	//delete myRandom;

	return energypoint;

}


