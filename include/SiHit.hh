 // $Id: SiHit.hh 100 2010-01-26 16:12:59Z adotti $
 #ifndef SiHit_h
 #define SiHit_h 1
 
 /**
  * @file
  * @brief Define user class SiHit.
  *
  * We define "our" hit format : it is caracterized by its plane and strip
  * numbers, and an energy value, the accumulated energy in this strip
  */
 
 #include "G4VHit.hh"
 #include "G4Allocator.hh"
 #include "G4ThreeVector.hh"
 #include "G4THitsCollection.hh"

#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4ParticleDefinition.hh" 


 /*!
  * \brief This class stores information of a hit.
  *
  * It contains
  *  - strip and plane number
  *  - deposited energy
  *  - position information
  */
 class SiHit : public G4VHit {
 public:
   /// Constructor
   SiHit(const G4int strip, const G4int plane, const G4bool isPrimary);
   /// Destructor
   ~SiHit();
   //! Print on screen a Hit
   void Print();
   
 public:
   //! \name The new and delete operators are overloaded for performances reasons:
   /*! -- Tricky business here... :-(, but provided for you below */
   //@{
   inline void *operator    new(size_t);
   inline void  operator delete(void *aHit);
   //@}
 
 public:
   //! \name  simple set and get methods
   //@{
   void          AddEdep(const double e)                { eDep += e; }
   void          SetPosition(const G4ThreeVector & pos) { position = pos; }

  // Rotation matrix
  inline void SetRotation(G4RotationMatrix rotation) {fRotation = rotation;}
  inline G4RotationMatrix GetRotation() const {return fRotation;}

  // Logical volume
  inline void SetLogicalVolume(G4LogicalVolume* volume) {pLogicalVolume = volume;}
  inline const G4LogicalVolume* GetLogicalVolume() const {return pLogicalVolume;}
// get time
  inline void SetIncidenceTime(G4double time) {fITime = time;}
  inline G4double GetIncidenceTime() const {return fITime;}
//momentum
inline void SetIncidenceMomentumDirection(G4ThreeVector momentum) {fIMomentumD = momentum;}
  inline G4ThreeVector GetIncidenceMomentumDirection() const {return fIMomentumD;}
 

  inline void SetIncidenceKineticEnergy(G4double ekin) {fIKEnergy = ekin;}
  inline G4double GetIncidenceKineticEnergy() const {return fIKEnergy;}
 
   G4double      GetEdep()        const { return eDep;}
   G4ThreeVector GetPosition()    const { return position; }
   G4int         GetStripNumber() const { return stripNumber; }
   G4int         GetPlaneNumber() const { return planeNumber; }
   G4bool            GetIsPrimary()   const { return isPrimary; }
   //@}
 
 private:
   const G4int   stripNumber, planeNumber;
   G4double      eDep;
   G4ThreeVector position;
   const G4bool  isPrimary;
   G4double fITime;
   G4ThreeVector fIMomentumD;
   G4double fIKEnergy;

  G4RotationMatrix fRotation;
  const G4LogicalVolume* pLogicalVolume;
 };
 
 // Define the "hit collection" using the template class G4THitsCollection:
 typedef G4THitsCollection<SiHit> SiHitCollection;
 
 
 // -- new and delete overloaded operators:
 extern G4Allocator<SiHit> SiHitAllocator;
 
 inline void* SiHit::operator new(size_t)
 {
   void *aHit;
   aHit = (void *) SiHitAllocator.MallocSingle();
   return aHit;
 }
 inline void SiHit::operator delete(void *aHit)
 {
   SiHitAllocator.FreeSingle((SiHit*) aHit);
 }
 
 #endif
