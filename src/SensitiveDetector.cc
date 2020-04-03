 // $Id: SensitiveDetector.cc 100 2010-01-26 16:12:59Z adotti $
 /**
  * @file   SensitiveDetector.cc
  * @brief  Implements sensitive part of simulation.
  *
  * @date   10 Dec 2009
  * @author Andreas Schaelicke
  */
 
 #include "SensitiveDetector.hh"
 
 #include "G4Step.hh"
 #include "Randomize.hh"
 
 #include "G4HCofThisEvent.hh"
 
 #include "G4HCtable.hh"
 #include "G4SDManager.hh"
 
 
 SensitiveDetector::SensitiveDetector(G4String SDname)
   : G4VSensitiveDetector(SDname)
 {
   // 'collectionName' is a protected data member of base class G4VSensitiveDetector.
   // Here we declare the name of the collection we will be using.
   collectionName.insert("SiHitCollection");
  
   // Note that we may add as many collection names we would wish: ie
   // a sensitive detector can have many collections.
 }
 
 SensitiveDetector::~SensitiveDetector()
 {}
 
 G4bool SensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *)
 {
   // step is guaranteed to be in Strip volume : no need to check for volume
   
   G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
   // energy deposit in this step 
   G4double edep = step->GetTotalEnergyDeposit();
   //check if step is due to primary particle: it has track ID 1 and parent 0
   // The primary is the track with ID 1 and with no parent
   G4double tiempo = step->GetPreStepPoint()->GetGlobalTime();

   G4bool isPrimary = ( step->GetTrack()->GetTrackID() == 1 && step->GetTrack()->GetParentID() == 0 ) ? true : false;
 
   if (edep <= 0.) return false;
 
   // get step points in world coordinate system
   G4ThreeVector point1 = step->GetPreStepPoint()->GetPosition();
   G4ThreeVector point2 = step->GetPostStepPoint()->GetPosition();

 


   // randomize point of energy deposition
   G4ThreeVector pointE = point1 + G4UniformRand()*(point2 - point1);      
 
   G4int stripCopyNo = touchable->GetReplicaNumber();
   G4int planeCopyNo = touchable->GetReplicaNumber(1);
 
   SiHit* hit = new SiHit(stripCopyNo,planeCopyNo,isPrimary);
   hitCollection->insert(hit);

  // ojo, a ver si trasforma las coordenadas del detetor
//--------------------------------------------------------------------------------------------

G4AffineTransform const& toLocal = step->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform();
// G4AffineTransform        toWorld = toLocal.Inverse();
G4ThreeVector pointELocal  = toLocal.TransformPoint(pointE);
 G4ThreeVector momentumD = toLocal.NetRotation() * step->GetPreStepPoint()->GetMomentumDirection();
hit->SetIncidenceMomentumDirection(momentumD);
   hit->SetPosition(pointELocal);

        G4double ke = step->GetPreStepPoint()->GetKineticEnergy();
	hit->SetIncidenceKineticEnergy(ke);

/*
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
   hit->SetLogicalVolume(thePhysical->GetLogicalVolume());
   G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
   aTrans.Invert();
   hit->SetRotation(aTrans.NetRotation());
   hit->SetPosition(aTrans.NetTranslation());
// Get Transportaion Matrix
 G4ThreeVector position = aTrans.NetRotation() * ( aTrans.NetTranslation() + pointE ); 
 G4ThreeVector momentumD = aTrans.NetRotation() * step->GetPreStepPoint()->GetMomentumDirection();
*/
 //-------------------------------------------------------------------------------------------


   // set energy deposition
   hit->AddEdep(edep);
   // store position of energy deposition
  
   hit->SetIncidenceTime(tiempo);
  
   return true;
 }
 
 void SensitiveDetector::Initialize(G4HCofThisEvent* HCE)
 {
   // ------------------------------
   // -- Creation of the collection
   // ------------------------------
   // -- collectionName[0] is "SiHitCollection", as declared in constructor
   hitCollection = new SiHitCollection(GetName(), collectionName[0]);
 
   // ----------------------------------------------------------------------------
   // -- and attachment of this collection to the "Hits Collection of this Event":
   // ----------------------------------------------------------------------------
   // -- To insert the collection, we need to get an index for it. This index
   // -- is unique to the collection. It is provided by the GetCollectionID(...)
   // -- method (which calls what is needed in the kernel to get this index).
   static G4int HCID = -1;
   if (HCID<0) HCID = GetCollectionID(0); // <<-- this is to get an ID for collectionName[0]
   HCE->AddHitsCollection(HCID, hitCollection);
 }
 
 void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
 {
 /*      
   // test output of hits
   G4cout << "EndOfEvent method of SD `" << GetName() << "' called." << G4endl;
   for (size_t i = 0; i < hitCollection->GetSize(); ++i ) {
     G4cout<<i<<G4endl;
         ( *hitCollection)[i]->Print();      
   }
 
 */
   // -- we could have attached the collection to the G4HCofThisEvent in this
   // -- method as well (instead of Initialize).
 }
