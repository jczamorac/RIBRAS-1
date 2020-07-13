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
#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"

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
  G4Material *CD2;
  G4Material *CH2;

  // -------------------------------------------------------------------------  //

  // Logical Volumes
  G4LogicalVolume *logicWorld; // World

  G4LogicalVolume *Log_Solenoid1; // Solenoid 1
  G4LogicalVolume *Log_Solenoid2; // Solenoid 2

  G4LogicalVolume *Log_Magnet1; // Magnetic Field 1
  G4LogicalVolume *Log_Magnet2; // Magnetic Field 2

  G4LogicalVolume *Log_Target; // Target

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
  G4VPhysicalVolume *Phys_Solenoid1; // Solenoid 1
  G4VPhysicalVolume *Phys_Solenoid2; // Solenoid 2
  G4VPhysicalVolume *Phys_Magnet1;   // Magnetic Field 1
  G4VPhysicalVolume *Phys_Magnet2;   // Magnetic Field 2
  G4VPhysicalVolume *Phys_Target;    // Target

  // -------------------------------------------------------------------------  //

  // Solid Volumes
  G4Tubs *Sol_Solenoid1; // Solenoid 1
  G4Tubs *Sol_Solenoid2; // Solenoid 2
  G4Tubs *Sol_Magnet1;   // Magnetic Field 1
  G4Tubs *Sol_Magnet2;   // Magnetic Field 2
  G4Box *Sol_Target;     // Target

  // -------------------------------------------------------------------------  //

  // Magnetic Field
  MagneticField *magneticField;

  // -------------------------------------------------------------------------  //

  // Parameters
  G4double halfWorldLength; // World half length
  G4int noOfSensorStrips;   // Number of Strips in each detector

  // -------------------------------------------------------------------------  //

  // Detectors parameters
  G4double Lengthy_dssd_t1;   // Detectors height
  G4double Thickness_dssd_t1; // Detectors thickness
  G4double Lengthx_dssd_t1;   // Detectors length

  // Detector position
  G4ThreeVector DetectorPosition[8];
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

public:                       // Target
  G4Material *TargetMaterial; // Target material
  G4ThreeVector target_pos;   // Target Position
  G4double width;             // Target thickness
  G4String g4_material_name;
  G4double target_mass;     // Target Mass
  G4int target_A, target_Z; // Target A , Z

public:                            // Recoil / Ejectile
  G4double recoil_mass, recoil_Ex; // Recoil particle mass, excitation energy
  G4int recoil_A, recoil_Z;        // Recoil particle A , Z

  G4double ejectile_mass, ejectile_Ex; // Ejectile particle mass, excitation energy
  G4int ejectile_A, ejectile_Z;        // Ejectile particle A , Z

  G4ParticleDefinition *RecoilParticle;
  G4double rTheta;
  G4double rKinectEnergy;
  G4ParticleDefinition *EjectileParticle;

public: // Decay products 1,2
  G4double decayp1_mass, decayp1_Ex;
  G4int decayp1_A, decayp1_Z;

  G4double decayp2_mass, decayp2_Ex;
  G4int decayp2_A, decayp2_Z;

  G4ParticleDefinition *DecayParticle1;
  G4ParticleDefinition *DecayParticle2;

public:                       // Primary beam
  G4double primary_energy;    // Primary beam energy
  G4int primary_Z, primary_A; // Primary beam Z , A
  G4ThreeVector primary_pos;  // Primary beam vertex position

public: // Detectors
  G4bool using_detectors;

public: // Magnetic Field
  G4double Current1, Current2;

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
             decayp1_mass(0), decayp1_Ex(0), decayp1_A(0), decayp1_Z(0),
             decayp2_mass(0), decayp2_Ex(0), decayp2_A(0), decayp2_Z(0),
             DecayParticle1(0), DecayParticle2(0),
             primary_energy(0), primary_Z(0), primary_A(0),
             primary_pos(0), RecoilParticle(0), EjectileParticle(0),
             rTheta(0), rKinectEnergy(0),
             Current1(0), Current2(0){};
  Inputs(Inputs const &) = delete;
  void operator=(Inputs const &) = delete;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// This class builds all detectors, you must create an object, define its parameters, rotate or not, set color and position, and them construct it.
class Detector
{
public:
  // Constructor
  Detector(const G4int DetectorNumber, const G4double Height, const G4double Lenght, const G4double Width, G4int nStrips);

  // Destructor
  ~Detector();

public:
  // Construct Detector
  G4VPhysicalVolume *Construct(G4LogicalVolume *LogicalMother);

  void SetPosition(G4ThreeVector Pos) { Position = Pos; }

  void SetColor(G4Color Colour) { Color = Colour; }

  void Rotate(G4double RotX, G4double RotY, G4double RotZ);

private:
  const G4int rDetector;
  const G4double rHeight, rLenght, rWidth;
  G4int rStrips;
  G4Color Color;
  G4ThreeVector Position;
  G4String DetectorName, StripsName;
  G4RotationMatrix *RotationMatrix = new G4RotationMatrix;
};

#endif
