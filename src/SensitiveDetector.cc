// Local headers
#include "SensitiveDetector.hh"

// Geant4 headers
#include "G4Step.hh"
#include "Randomize.hh"
#include "G4HCofThisEvent.hh"
#include "G4HCtable.hh"
#include "G4SDManager.hh"

// ------------------------------------------------------------------------------------- //

SensitiveDetector::SensitiveDetector(G4String SDname) : G4VSensitiveDetector(SDname)
{
  // 'collectionName' is a protected data member of base class G4VSensitiveDetector.
  // Here we declare the name of the collection we will be using.
  collectionName.insert("SiHitCollection");

  // Note that we may add as many collection names we would wish: ie
  // a sensitive detector can have many collections.
}

// ------------------------------------------------------------------------------------- //

SensitiveDetector::~SensitiveDetector()
{
}

// ------------------------------------------------------------------------------------- //

G4bool SensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *)
{
  // Step is guaranteed to be in Strip volume : no need to check for volume
  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();

  // Energy deposit in this step
  G4double edep = step->GetTotalEnergyDeposit();

  // If edep <= 0, return false, it means there wasn't a hit
  if (edep <= 0.)
    return false;

  // Check if step is due to primary particle: it has track ID 1 and parent 0
  // The primary is the track with ID 1 and with no parent
  G4bool isPrimary = ((step->GetTrack()->GetParentID()) == 1) ? false : true;

  // Get hit time
  G4double tiempo = step->GetPreStepPoint()->GetGlobalTime();

  // Getting logical volume where occured the hit
  G4LogicalVolume *LogicalVolume = step->GetTrack()->GetVolume()->GetLogicalVolume();
  G4String LogicalName = LogicalVolume->GetName();

  // Get step points in world coordinate system
  G4ThreeVector point1 = step->GetPreStepPoint()->GetPosition();
  G4ThreeVector point2 = step->GetPostStepPoint()->GetPosition();

  // randomize point of energy deposition
  G4ThreeVector pointE = point1 + G4UniformRand() * (point2 - point1);

  G4int stripCopyNo = touchable->GetReplicaNumber();
  G4int planeCopyNo = touchable->GetReplicaNumber(1);

  // Saving all data retrieved from step
  SiHit *hit = new SiHit(stripCopyNo, planeCopyNo, isPrimary);
  hitCollection->insert(hit);

  // Get momentum, kinect energy, hit coordinates
  G4AffineTransform const &toLocal = step->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform();
  G4ThreeVector pointELocal = toLocal.TransformPoint(pointE);
  G4ThreeVector momentumD = toLocal.NetRotation() * step->GetPreStepPoint()->GetMomentumDirection();
  G4double ke = step->GetPreStepPoint()->GetKineticEnergy();

  // Saving everything
  hit->SetIncidenceMomentumDirection(momentumD);
  hit->SetPosition(pointELocal);
  hit->SetIncidenceKineticEnergy(ke);
  hit->SetLogicalVolume(LogicalVolume);
  hit->AddEdep(edep);
  hit->SetIncidenceTime(tiempo);

  // Setting hit bools
  if (LogicalName == "Log_Target")            // If this happen, it means that the hit was on target
  {
    hit->SetHitOnTarget();
  }
  else if (LogicalName == "SensorStripD00" ||
           LogicalName == "SensorStripD01" ||
           LogicalName == "SensorStripD02" ||
           LogicalName == "SensorStripD03" ||
           LogicalName == "SensorStripD04" ||
           LogicalName == "SensorStripD05" ||
           LogicalName == "SensorStripD06" ||
           LogicalName == "SensorStripD07")   // Else, hit on detector
  {
    hit->SetHitOnDetector();
  }

  // Return true means that a hit ocurred
  return true;
}

// ------------------------------------------------------------------------------------- //

void SensitiveDetector::Initialize(G4HCofThisEvent *HCE)
{
  // -------------------------------- //
  // -- Creation of the collection -- //
  // -------------------------------- //

  // -- collectionName[0] is "SiHitCollection", as declared in constructor
  hitCollection = new SiHitCollection(GetName(), collectionName[0]);

  // ----------------------------------------------------------------------------
  // -- and attachment of this collection to the "Hits Collection of this Event":
  // ----------------------------------------------------------------------------
  // -- To insert the collection, we need to get an index for it. This index
  // -- is unique to the collection. It is provided by the GetCollectionID(...)
  // -- method (which calls what is needed in the kernel to get this index).
  static G4int HCID = -1;
  if (HCID < 0)
    HCID = GetCollectionID(0); // <<-- this is to get an ID for collectionName[0]
  HCE->AddHitsCollection(HCID, hitCollection);
}

// ------------------------------------------------------------------------------------- //

void SensitiveDetector::EndOfEvent(G4HCofThisEvent *)
{
  /*   // test output of hits
  G4cout << "EndOfEvent method of SD `" << GetName() << "' called." << G4endl;
  for (size_t i = 0; i < hitCollection->GetSize(); ++i)
  {
    (*hitCollection)[i]->Print();
    G4cout << " " << G4endl;
  } */
}

// ------------------------------------------------------------------------------------- //