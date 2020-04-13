// $Id: DetectorConstruction.hh 33 2010-01-14 17:08:18Z adotti $
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

/**
* @file
* @brief Defines mandatory user class DetectorConstruction.
*/

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"

class MagneticField;
class MagneticField2;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
class G4Ions;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*!
\brief This mandatory user class defines the geometry.

It is responsible for
- Definition of material, and
- Construction of geometry

\sa Construct()
*/
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  //! Constructor
  DetectorConstruction();
  //! Destructor
  ~DetectorConstruction();
public:
  //! Construct geometry of the setup
  G4VPhysicalVolume* Construct();

  //! Update geometry
  void UpdateGeometry();

  //! \name some simple set & get functions
  //@{
  G4ThreeVector FirstSensorPosition() const  { return vector_D_0_0; }
  G4ThreeVector SecondSensorPosition() const { return vector_D_0_1; }
  //G4ThreeVector ThirdSensorPosition() const { return vector_D_0_2; }
  //G4ThreeVector CuartoSensorPosition() const { return posCuartoSensor; }

  G4ThreeVector SetFirstSensorPosition(const G4ThreeVector & pos) { return vector_D_0_0=pos; }
  G4ThreeVector SetSecondSensorPosition(const G4ThreeVector & pos) { return vector_D_0_1=pos; }
  // G4ThreeVector SetThirdSensorPosition(const G4ThreeVector & pos) { return vector_D_0_2=pos; }
  //G4ThreeVector SetCuartoSensorPosition(const G4ThreeVector & pos) { return posCuartoSensor=pos; }



//-------------------------algunos mensajes------------------------------
  G4double GetSlit_z() const { return slit_z; }
  G4double SetSlit_z(const G4double el_z)  { return slit_z=el_z; }

  G4double GetSlit_x() const { return slit_x; }
  G4double SetSlit_x(const G4double el_x)  { return slit_x=el_x; }

  G4double GetPocket1_x() const { return pocket_x; }
  G4double SetPocket1_x(const G4double el_x_pocket)  { return pocket_x=el_x_pocket; }
//-----------------------------------------------------------------------


  //@}
private:
  //! define needed materials
  void DefineMaterials();
  //! initialize geometry parameters
  void ComputeParameters();
  //! Construct geometry of the Beam Telescope
  G4VPhysicalVolume* ConstructDetectors();
  //! Construct geometry of the Device-under-test
  void ConstructSetup(void);


  //G4VPhysicalVolume* ConstructRings();
  //! Construct geometry of the Device-under-test
  G4VPhysicalVolume* ConstructRing_D_0_0();
  //! Construct geometry of the 2Pocket
  G4VPhysicalVolume* ConstructRing_D_0_1();

/* G4VPhysicalVolume* ConstructRing_D_0_2(); */

  G4double* rotacion(G4double dx , G4double dy, G4double dz, G4double angulo);

private:

  //! \name Materials
  //@{
  G4Material* air;
  G4Material* silicon;
  G4Material* vacuum;
  G4Material* H2;
  G4Material* AlN;
  G4Material* steel;
  G4Material* tungsten;
  G4Material* lead;
  G4Material* tantalum;
  G4Material* sio2;
  G4Material* CH4;

  //@}

  //! \name Geometry
  //@{

  //! global mother volume
  G4LogicalVolume * logicWorld;


  MagneticField* magneticField;
  MagneticField2* magneticField2;


  //G4VPhysicalVolume * physiSensorslit;

  G4VPhysicalVolume * physiTarget;
  G4VPhysicalVolume * physiCilindro;
  G4VPhysicalVolume* solenoidef;
  G4VPhysicalVolume* Phy_Solenoid;

  //  Região do campo magnético


  //G4VPhysicalVolume* Phy_Magnet;
  //G4VPhysicalVolume* magneticof;




//@}

  //! \name Parameters
  //@{
  G4double halfWorldLength;

  G4int noOfSensorStrips;
  G4double Lengthy_sili;
  G4double Thickness_sili;

  G4double Lengthx_sili;

  G4double Lengthy_dssd;
  G4double Thickness_dssd;
  G4double Thickness_dssd_gr;
  G4double Lengthx_dssd;

  G4ThreeVector posFirstSensor;
  G4ThreeVector posSecondSensor;
  G4ThreeVector posSecondSensor1a;
  G4ThreeVector posThirdSensor;

  G4double slit_x;
  G4double slit_z;
  G4double pocket_x;

  G4double fwhm_beam;
  G4double fwhm_target;
  //@}

  //------------ring D_0_0
  G4VPhysicalVolume * physiSensorStripD_0_0;
  G4VPhysicalVolume* detector1phys;
  G4ThreeVector vector_D_0_0;
  G4double phi_D_0_0;
  G4double coorx_D_0_0;
  G4double coory_D_0_0;
  G4double coorz_D_0_0;
  G4double coorx_D_0_0_e;
  G4double coory_D_0_0_e;
  G4double coorz_D_0_0_e;
  G4double theta_D_0_0;


  //------------ring D_0_1
  G4VPhysicalVolume* physiSensorStripD_0_1;
  G4VPhysicalVolume* physiSensorRing_D_0_1;
  G4ThreeVector vector_D_0_1;
  G4double phi_D_0_1;
  G4double coorx_D_0_1;
  G4double coory_D_0_1;
  G4double coorz_D_0_1;
  G4double coorx_D_0_1_e;
  G4double coory_D_0_1_e;
  G4double coorz_D_0_1_e;
  G4double theta_D_0_1;

  G4double Lengthy_dssd_t1;
  G4double Thickness_dssd_t1;
  G4double Lengthx_dssd_t1;

  //! \name UI Messenger
  //@{
  DetectorMessenger * messenger;
  //@}
  //
 };

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Static singleton to handle inputs for CESimulations across translation units (*.hh/*.cc file pairs)
 class Inputs {
 public:
 	static Inputs& GetInputs() {
 		static Inputs    instance;
 		return instance;
 	}
 public: // Flags:
 	G4bool initialized;
 	G4bool using_lh2;
 	G4bool using_recfoil;
 	G4bool using_cylfoil;
 	G4bool using_lenda;
 	G4bool using_cagra;
 	G4bool using_source;
  G4bool using_magneticfield;
 	G4double source_energy;
 public: // Target
 	// Uppercase -> set internally (may derive from input file down the line)
 	G4Material* TargetMaterial;
  G4ThreeVector target_pos;
 	// Lowercase -> set externally (mac file)
 	G4double radius;
 	G4double height, width;
 	G4double arial_density;
 	G4String g4_material_name;
 public: // Recoil / Ejectile
 	G4double recoil_mass, recoil_Ex;
 	G4int recoil_A, recoil_Z;
 	G4double target_mass;
 	G4int target_A, target_Z;
 	G4double ejectile_mass, ejectile_Ex;
 	G4int ejectile_A, ejectile_Z;
  G4ParticleDefinition* RecoilParticle;
  G4ParticleDefinition* EjectileParticle;  
  public: //Primary beam
  G4double primary_energy;
  G4int primary_Z, primary_A;
  G4ThreeVector primary_pos;
  
 private:
 	Inputs() :
 		initialized(false),
 		using_lh2(false), using_recfoil(false),
 		using_cylfoil(false), using_lenda(false),
 		using_cagra(false), using_source(false), using_magneticfield(false),
 		source_energy(0.0), TargetMaterial(nullptr),
 		radius(0.0), height(0.0), width(1),
 		arial_density(0.0),	g4_material_name(""),
 		recoil_mass(0.0), recoil_Ex(0.0),
 		recoil_A(0), recoil_Z(0), target_mass(0),
 		target_A(0), target_Z(0), target_pos(0),
 		ejectile_mass(0.0), ejectile_Ex(0.0),
 		ejectile_A(0), ejectile_Z(0),
    primary_energy(0), primary_Z(0), primary_A(0),
    primary_pos(0), RecoilParticle(0), EjectileParticle(0) {};
 	Inputs(Inputs const&) = delete;
 	void operator=(Inputs const&) = delete;
 };
 #endif
