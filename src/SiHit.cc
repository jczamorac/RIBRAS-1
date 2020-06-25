// Local headers
#include "SiHit.hh"

// One more nasty trick for new and delete operator overloading:
G4Allocator<SiHit> SiHitAllocator;

// ----------------------------------------------------------------------------------------------------------------------------------//
SiHit::SiHit(const G4int strip, const G4int plane, const G4bool primary)
    : stripNumber(strip), planeNumber(plane), isPrimary(primary) // <<-- note BTW this is the only way to initialize a "const" member
{
  eDep = 0.0;
}
// ----------------------------------------------------------------------------------------------------------------------------------//

SiHit::~SiHit()
{
}

// ----------------------------------------------------------------------------------------------------------------------------------//

void SiHit::Print()
{
}

// ----------------------------------------------------------------------------------------------------------------------------------//
