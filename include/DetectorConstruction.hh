#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

// Local headers
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "MagneticField.hh"
#include "SensitiveDetector.hh"
#include "globals.hh"

// Geant4 headers
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

class MagneticField;
class MagneticField2;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
class G4Ions;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  // Constructor
  DetectorConstruction();
  // Destructor
  ~DetectorConstruction();
  // Construct geometry of the setup
  G4VPhysicalVolume *Construct();
  // Update geometry
  void UpdateGeometry();

private:
  // Define needed materials
  void DefineMaterials();
  // Initialize geometry parameters
  void ComputeParameters();
  // Construct geometry of the Beam Telescope
  void ConstructDetectors();
  // Construct geometry of the Device-under-test
  void ConstructSetup(void);

  // Function to construct detectors
  G4VPhysicalVolume *Construct_D_0_0();
  G4VPhysicalVolume *Construct_D_0_1();
  G4VPhysicalVolume *Construct_D_0_2();
  G4VPhysicalVolume *Construct_D_0_3();
  G4VPhysicalVolume *Construct_D_0_4();
  G4VPhysicalVolume *Construct_D_0_5();
  G4VPhysicalVolume *Construct_D_0_6();
  G4VPhysicalVolume *Construct_D_0_7();

private:
  // -------------------------------------------------------------------------  //

  // Messenger
  DetectorMessenger *messenger;

  // -------------------------------------------------------------------------  //

  // Materials
  G4Material *air;
  G4Material *silicon;
  G4Material *vacuum;
  G4Material *H2;
  G4Material *AlN;
  G4Material *steel;
  G4Material *tungsten;
  G4Material *lead;
  G4Material *tantalum;
  G4Material *sio2;
  G4Material *CH4;

  // -------------------------------------------------------------------------  //

  // Logical Volumes
  G4LogicalVolume *logicWorld;          // World

  G4LogicalVolume *Log_Solenoid1;       // Solenoid 1
  G4LogicalVolume *Log_Solenoid2;       // Solenoid 2

  G4LogicalVolume *Log_Magnet1;         // Magnetic Field 1
  G4LogicalVolume *Log_Magnet2;         // Magnetic Field 2

  G4LogicalVolume *Log_Target;          // Target

  G4LogicalVolume *logicSensorStripD00; // Strips
  G4LogicalVolume *logicSensorStripD01;
  G4LogicalVolume *logicSensorStripD02;
  G4LogicalVolume *logicSensorStripD03;
  G4LogicalVolume *logicSensorStripD04;
  G4LogicalVolume *logicSensorStripD05;
  G4LogicalVolume *logicSensorStripD06;
  G4LogicalVolume *logicSensorStripD07;

  // -------------------------------------------------------------------------  //

  // Physical Volumes
  G4VPhysicalVolume *Phys_Solenoid1;    // Solenoid 1
  G4VPhysicalVolume *Phys_Solenoid2;    // Solenoid 2
  G4VPhysicalVolume *Phys_Magnet1;      // Magnetic Field 1
  G4VPhysicalVolume *Phys_Magnet2;      // Magnetic Field 2
  G4VPhysicalVolume *Phys_Target;       // Target

  // -------------------------------------------------------------------------  //

  // Solid Volumes
  G4Tubs *Sol_Solenoid1;                // Solenoid 1
  G4Tubs *Sol_Solenoid2;                // Solenoid 2
  G4Tubs *Sol_Magnet1;                  // Magnetic Field 1
  G4Tubs *Sol_Magnet2;                  // Magnetic Field 2
  G4Box  *Sol_Target;                   // Target

  // -------------------------------------------------------------------------  //

  // Magnetic Field
  MagneticField *magneticField;

  // -------------------------------------------------------------------------  //

  // Parameters
  G4double halfWorldLength;             // World half length
  G4int noOfSensorStrips;               // Number of Strips in each detector

  // -------------------------------------------------------------------------  //

  // Detectors parameters
  G4double Lengthy_dssd_t1;             // Detectors height
  G4double Thickness_dssd_t1;           // Detectors thickness
  G4double Lengthx_dssd_t1;             // Detectors length

  // -------------------------------------------------------------------------  //
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// This class contain all the inputs needed by the simulation
class Inputs
{
public:
  static Inputs &GetInputs()
  {
    static Inputs instance;
    return instance;
  }

public: // Flags:
  G4bool initialized;
  G4bool using_magneticfield;

public: // Target
  G4Material *TargetMaterial;                 // Target material
  G4ThreeVector target_pos;                   // Target Position 
  G4double width;                             // Target thickness
  G4String g4_material_name;      
  G4double target_mass;                       // Target Mass
  G4int target_A, target_Z;                   // Target A , Z

public: // Recoil / Ejectile
  G4double recoil_mass, recoil_Ex;            // Recoil particle mass, excitation energy
  G4int recoil_A, recoil_Z;                   // Recoil particle A , Z

  G4double ejectile_mass, ejectile_Ex;        // Ejectile particle mass, excitation energy
  G4int ejectile_A, ejectile_Z;               // Ejectile particle A , Z

  G4ParticleDefinition *RecoilParticle;
  G4ParticleDefinition *EjectileParticle;

public: // Primary beam
  G4double primary_energy;                    // Primary beam energy
  G4int primary_Z, primary_A;                 // Primary beam Z , A
  G4ThreeVector primary_pos;                  // Primary beam vertex position

public: // Detectors
  G4ThreeVector detector1_pos;                
  G4ThreeVector detector2_pos;

private:
  Inputs() : // Initializing parameters
             initialized(false),
             using_magneticfield(false), TargetMaterial(nullptr),
             width(1), g4_material_name(""),
             recoil_mass(0.0), recoil_Ex(0.0),
             recoil_A(0), recoil_Z(0), target_mass(0),
             target_A(0), target_Z(0), target_pos(0),
             ejectile_mass(0.0), ejectile_Ex(0.0),
             ejectile_A(0), ejectile_Z(0),
             primary_energy(0), primary_Z(0), primary_A(0),
             primary_pos(0), RecoilParticle(0), EjectileParticle(0),
             detector1_pos(0), detector2_pos(0){};
  Inputs(Inputs const &) = delete;
  void operator=(Inputs const &) = delete;
};
#endif
