#ifndef SiHit_h
#define SiHit_h 1

// Geant4 headers
#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4THitsCollection.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4ParticleDefinition.hh"

class SiHit : public G4VHit
{
public:
  // Constructor
  SiHit(const G4int strip, const G4int plane, const G4bool isPrimary);

  // Destructor
  ~SiHit();

  // Print on screen a Hit
  void Print();

public:
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

public:
  // Deposited energy
  void AddEdep(const double e) { eDep += e; }
  G4double GetEdep() const { return eDep; }

  // Hit Coordinates
  void SetPosition(const G4ThreeVector &pos) { position = pos; }
  G4ThreeVector GetPosition() const { return position; }

  // Rotation matrix
  inline void SetRotation(G4RotationMatrix rotation) { fRotation = rotation; }
  inline G4RotationMatrix GetRotation() const { return fRotation; }

  // Logical volume
  inline void SetLogicalVolume(G4LogicalVolume *volume) { pLogicalVolume = volume; }
  inline const G4LogicalVolume *GetLogicalVolume() const { return pLogicalVolume; }

  // Get time
  inline void SetIncidenceTime(G4double time) { fITime = time; }
  inline G4double GetIncidenceTime() const { return fITime; }

  // Momentum
  inline void SetIncidenceMomentumDirection(G4ThreeVector momentum) { fIMomentumD = momentum; }
  inline G4ThreeVector GetIncidenceMomentumDirection() const { return fIMomentumD; }

  // Kinect Energy
  inline void SetIncidenceKineticEnergy(G4double ekin) { fIKEnergy = ekin; }
  inline G4double GetIncidenceKineticEnergy() const { return fIKEnergy; }

  // Particle ID
  inline void SetParticleID(int i) { ParticleID = i; }
  int GetParticleID() const { return ParticleID; }

  // Hit on Detector
  inline void SetHitOnDetector() { HitOnDetector = true; }
  G4bool GetHitOnDetector() const { return HitOnDetector; }

  // Hit on Target
  inline void SetHitOnTarget() { HitOnTarget = true; }
  G4bool GetHitOnTarget() const { return HitOnTarget; }

  // Recoil Theta CM
  inline void SetRecoilThetaCM( G4double T ){ RecoilThetaCM = T; }
  G4double GetRecoilThetaCM() const { return RecoilThetaCM; }
  
  // This variables are initialized and saved in the constructor
  G4int GetStripNumber() const { return stripNumber; }
  G4int GetPlaneNumber() const { return planeNumber; }
  G4bool GetIsPrimary() const { return isPrimary; }

private:
  const G4int stripNumber, planeNumber;
  G4double eDep;
  G4ThreeVector position;
  const G4bool isPrimary;
  G4bool HitOnDetector = false;
  G4bool HitOnTarget = false;
  G4double fITime;
  G4ThreeVector fIMomentumD;
  G4double fIKEnergy;
  G4int ParticleID;
  G4double RecoilThetaCM;

  G4RotationMatrix fRotation;
  const G4LogicalVolume *pLogicalVolume;
};

// Define the "hit collection" using the template class G4THitsCollection:
typedef G4THitsCollection<SiHit> SiHitCollection;

// New and delete overloaded operators:
extern G4Allocator<SiHit> SiHitAllocator;

inline void *SiHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *)SiHitAllocator.MallocSingle();
  return aHit;
}
inline void SiHit::operator delete(void *aHit)
{
  SiHitAllocator.FreeSingle((SiHit *)aHit);
}

#endif
