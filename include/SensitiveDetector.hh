 // $Id: SensitiveDetector.hh 22 2009-12-22 12:36:46Z schaelic $
 
 #ifndef SensitiveDetector_h
 #define SensitiveDetector_h 1
 
 /**
  * @file
  * @brief Defines SensitiveDetector class.
  */
 
 #include "G4VSensitiveDetector.hh"
 class DetectorConstruction;
 class RunAction;
 
 #include "SiHit.hh"              // <<- the hit "format" we define
 class G4HCofThisEvent;           // <<- means "H(it) C(ollections) of This Event"
 
 
 /*!
  * \brief Defines sensitve part of detector geometry.
  *
  * Stores Hits with 
  *  * deposited energy
  *  * position
  * in <i>Hit Collections of This Event</i>
  *
  * /sa ProcessHits()
  */
 class SensitiveDetector : public G4VSensitiveDetector
 {
 public:
   /// Constructor
   SensitiveDetector(G4String SDname);
   /// Destructor
   ~SensitiveDetector();
 
 public:
   /// @name methods from base class G4VSensitiveDetector
   //@{
   /// Mandatory base class method : it must to be overloaded:
    G4bool ProcessHits(G4Step *step, G4TouchableHistory *ROhist);
 
   /// (optional) method of base class G4VSensitiveDetector
   void Initialize(G4HCofThisEvent* HCE);
   /// (optional) method of base class G4VSensitiveDetector
   void EndOfEvent(G4HCofThisEvent* HCE);
   //@}
 
 
 private:
   SiHitCollection*      hitCollection;
 };
 
 #endif