 // $Id: SiHit.cc 100 2010-01-26 16:12:59Z adotti $
 /**
  * @file
  * @brief Implementation of user class SiHit.
 */
 
 #include "SiHit.hh"
 
 // -- one more nasty trick for new and delete operator overloading:
 G4Allocator<SiHit> SiHitAllocator;
 
 SiHit::SiHit(const G4int strip, const G4int plane, const G4bool primary )
   : stripNumber(strip), planeNumber(plane) , isPrimary(primary)// <<-- note BTW this is the only way to initialize a "const" member
 {
   eDep     = 0.0;
 }
 
 SiHit::~SiHit()
 {
 }
 
 void SiHit::Print()
 {
         G4cout<<"Hit: Plane= "<<planeNumber<<" Strip= "<<stripNumber<<" E= "<<eDep/CLHEP::MeV<<" CLHEP::MeV isPrimary="<<(isPrimary?"true":"false")<<G4endl;
 }
