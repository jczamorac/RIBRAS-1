// $Id: DetectorConstruction.cc 94 2010-01-26 13:18:30Z adotti $
/**
* @file
* @brief Implements mandatory user class DetectorConstruction.
*/
// Local headers
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "MagneticField.hh"
#include "SensitiveDetector.hh"

// Default headers
#include <sstream>

// Geant4 headers
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Transform3D.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4UserLimits.hh"
#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

using namespace std;
using namespace CLHEP;

SensitiveDetector *Sensitive = (SensitiveDetector *)0;

//--------------------------------------------------------------------------------------------------------------//

DetectorConstruction::DetectorConstruction()
{
  //Create a messanger (defines custom UI commands)
  messenger = new DetectorMessenger(this);

  //--------- Material definition ---------//
  DefineMaterials();

  //--------- Sizes of the principal geometrical components (solids)  ---------//
  ComputeParameters();
  ConstructSetup();
  magneticField = new MagneticField();
}

//--------------------------------------------------------------------------------------------------------------//

DetectorConstruction::~DetectorConstruction()
{
  delete messenger;
}

//--------------------------------------------------------------------------------------------------------------//

Detector::Detector(const G4int DetectorNumber, const G4double Height, const G4double Lenght, const G4double Width, G4int nStrips)
    : rDetector(DetectorNumber), rHeight(Height), rLenght(Lenght), rWidth(Width), rStrips(nStrips)
{
  // Getting Detector Name
  // e.g
  // Detector number is X so its name is going to be Detector_X and its strips SensorStripD0X
  ostringstream DetectorN, StripsN;

  DetectorN << "Detector " << rDetector;
  StripsN << "SensorStripD0" << rDetector;

  DetectorName = DetectorN.str().data();
  StripsName = StripsN.str().data();
}

//--------------------------------------------------------------------------------------------------------------//

Detector::~Detector()
{
}

//--------------------------------------------------------------------------------------------------------------//

void Detector::Rotate(G4double RotX, G4double RotY, G4double RotZ)
{
  // Angles in degree
  RotationMatrix->rotateX(RotX * deg);
  RotationMatrix->rotateY(RotY * deg);
  RotationMatrix->rotateZ(RotZ * deg);
}

//--------------------------------------------------------------------------------------------------------------//

G4VPhysicalVolume *Detector::Construct(G4LogicalVolume *LogicalMother)
{
  // Detector Material
  G4Material *Material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  // Creating a Physical Volume
  G4Box *DetectorBox = new G4Box(DetectorName, rLenght / 2, rHeight / 2, rWidth / 2);

  // Creating a Logical Volume
  G4LogicalVolume *LogVolume = new G4LogicalVolume(DetectorBox, Material, DetectorName);

  // Placing Detector
  G4VPhysicalVolume *PhysicalVolume = new G4PVPlacement(RotationMatrix, // Rotation Matrix
                                                        Position,       // Position
                                                        LogVolume,      // Respective Logical Volume
                                                        DetectorName,   // Name
                                                        LogicalMother,  // Logical Mother Volume
                                                        false,          // No Boolean Operation
                                                        rDetector,      // Copy Number
                                                        true);          // IDK

  // Strips parameters
  G4double StripHeight = rHeight / 2;
  G4double StripLenght = rLenght / (2 * rStrips);
  G4double StripWidth = rWidth / 2;

  // Creating Physical Volume
  G4Box *StripsBox = new G4Box(StripsName, StripLenght, StripHeight, StripWidth);

  // Creating LogicalVolume
  G4LogicalVolume *LogStripVolume = new G4LogicalVolume(StripsBox, Material, StripsName);

  //Placing Strips
  G4VPhysicalVolume *PhysicalStripsVolume = new G4PVReplica(StripsName,
                                                            LogStripVolume,
                                                            LogVolume,
                                                            kXAxis,
                                                            rStrips,
                                                            2.0 * StripLenght);

  // Set "User Limits"
  G4UserLimits *userLimits = new G4UserLimits(0.1 * mm);
  LogStripVolume->SetUserLimits(userLimits);

  // Setting Detector Sensitive
  if (!Sensitive)
  {
    Sensitive = new SensitiveDetector("/myDet/SiStripSD");
    //We register now the SD with the manager
    G4SDManager::GetSDMpointer()->AddNewDetector(Sensitive);
  }
  LogStripVolume->SetSensitiveDetector(Sensitive);

  // Color
  LogStripVolume->SetVisAttributes(new G4VisAttributes(Color));

  // Return Your Beaultiful Detector
  return PhysicalVolume;
}

//--------------------------------------------------------------------------------------------------------------//

void DetectorConstruction::DefineMaterials()
{

  G4double a, z, density; // z=mean number of protons;
  G4double temperature, pressure;
  G4int ncomponents, natoms;
  G4String name, symbol; // a=mass of a CLHEP::mole;

  // Define Elements
  a = 1.01 * CLHEP::g / CLHEP::mole;
  G4Element *elH = new G4Element(name = "Hydrogen", symbol = "H", z = 1., a);
  density = 1.33e-11 * CLHEP::g / CLHEP::cm3;
  pressure = 1.0913e-10 * CLHEP::atmosphere;
  temperature = 200. * CLHEP::kelvin;
  H2 = new G4Material(name = "Hydrogen gas", density, ncomponents = 1,
                      kStateGas, temperature, pressure);
  H2->AddElement(elH, natoms = 2);

  // Define the AlN
  G4Element *elC = new G4Element(name = "Carbono", symbol = "C", z = 6., a = 12.01 * CLHEP::g / CLHEP::mole);
  G4Element *elHi = new G4Element(name = "Hidrogênio", symbol = "H", z = 1., a = 1.01 * CLHEP::g / CLHEP::mole);
  density = 0.98 * CLHEP::kg / CLHEP::m3;
  temperature = 200. * CLHEP::kelvin;
  pressure = 1. * CLHEP::atmosphere;
  CH4 = new G4Material(name = "Metano", density, ncomponents = 2, kStateGas, temperature, pressure);
  CH4->AddElement(elC, natoms = 1);
  CH4->AddElement(elHi, natoms = 4);

  // Define POLYETHYLENE deuterado
  G4Element *deuteron = new G4Element(name = "Deuteron", symbol = "D", z = 1., a = 2. * CLHEP::g / CLHEP::mole);
  CD2 = new G4Material(name = "CD2", 0.94 * g / cm3, ncomponents = 2, kStateSolid); /* , temperature, pressure); */
  CD2->AddElement(elC, natoms = 1);
  CD2->AddElement(elHi, natoms = 2);

  // Get Materials from NIST database
  G4NistManager *man = G4NistManager::Instance();
  man->SetVerbose(0);

  // Define NIST materials
  air = man->FindOrBuildMaterial("G4_AIR");
  silicon = man->FindOrBuildMaterial("G4_Si");
  vacuum = man->FindOrBuildMaterial("G4_Galactic");
  steel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  tungsten = man->FindOrBuildMaterial("G4_W");
  tantalum = man->FindOrBuildMaterial("G4_Ta");
  sio2 = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  CH2 = man->FindOrBuildMaterial("G4_POLYETHYLENE");
}

//--------------------------------------------------------------------------------------------------------------//

void DetectorConstruction::ComputeParameters()
{
  // This function defines the defaults of the geometry construction
  // World size
  halfWorldLength = 20.0 * CLHEP::m;

  // Number of strips in the detector
  noOfSensorStrips = 60;

  // Detector Parameters
  Lengthy_dssd_t1 = 14 * CLHEP::cm;
  Lengthx_dssd_t1 = 25 * CLHEP::cm;
  Thickness_dssd_t1 = 300. * CLHEP::um;

  // Detector position

  G4double rPosition_z = -12.5;

  // Rear detectors
  DetectorPosition[0] = G4ThreeVector(-8., 0., rPosition_z) * cm;
  DetectorPosition[1] = G4ThreeVector(8., 0., rPosition_z) * cm;
  DetectorPosition[2] = G4ThreeVector(0., -8., rPosition_z) * cm;
  DetectorPosition[3] = G4ThreeVector(0., 8., rPosition_z) * cm;

  // Front detectors
  DetectorPosition[4] = G4ThreeVector(0., 10, 17.) * cm;
  DetectorPosition[5] = G4ThreeVector(0., -10, 17.) * cm;
  DetectorPosition[6] = G4ThreeVector(10, 0., 17.) * cm;
  DetectorPosition[7] = G4ThreeVector(-10, 0., 17.) * cm;
}

//--------------------------------------------------------------------------------------------------------------//

void DetectorConstruction::ConstructSetup(void)
{
}

//--------------------------------------------------------------------------------------------------------------//

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  // Retrieving inputs
  Inputs *Inputs = &Inputs::GetInputs();

  G4NistManager *man = G4NistManager::Instance();
  man->SetVerbose(0);

  // Available colors
  G4Color green(0.0, 1.0, 0.0),
      red(1.0, 0.0, 0.0),
      yellow(1.0, 1.0, 0.0),
      cyan(0.2, 0.5, 0.5),
      blue(0.0, 0.0, 1.0);

  //---- This function is called by G4 when the detector has to be created -----//
  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------//

  //----------------------------------------------------//
  // World - Aqui o "mundo" é criado. O mundo é Lógico. //
  //----------------------------------------------------//

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.0 * halfWorldLength);

  G4Box *solidWorld = new G4Box("world", halfWorldLength, halfWorldLength, halfWorldLength);
  logicWorld = new G4LogicalVolume(solidWorld, vacuum, "World", 0, 0, 0);

  //  Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume *physiWorld = new G4PVPlacement(0,               // no rotation
                                                    G4ThreeVector(), // at (0,0,0)
                                                    logicWorld,      // its logical volume
                                                    "World",         // its name
                                                    0,               // its mother  volume
                                                    false,           // no boolean operations
                                                    0);              // copy number

  logicWorld->SetVisAttributes(G4VisAttributes::Invisible); // Turn the world invisible

  // ---------------------------------------------------------------------------------------------------- //

  //-----Here, I'm building the 2 solenoids and their magnetic field------//

  //////////////////
  //  Solenoid 1  //
  //////////////////

  // Solenoid 1 parameters
  G4double Solenoid_length = 100.0 * CLHEP::cm;
  G4double Solenoid_diameter_inner = 30.05 * CLHEP::cm;
  G4double Solenoid_diameter_outer = 100.0 * CLHEP::cm;

  // Position
  G4ThreeVector Pos_Solenoid = G4ThreeVector();

  // Creating Solenoid 1
  Sol_Solenoid1 = new G4Tubs("Sol_Solenoid1", Solenoid_diameter_inner / 2.0, Solenoid_diameter_outer / 2.0, Solenoid_length / 2.0, 0., 360. * CLHEP::deg);
  Log_Solenoid1 = new G4LogicalVolume(Sol_Solenoid1, G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe"), "Log_Solenoid1");

  // Placing Solenoid 1
  Phys_Solenoid1 = new G4PVPlacement(0, Pos_Solenoid, Log_Solenoid1, "Solenoid 1", logicWorld, false, 0, true);

  // Setting Solenoid 1 color
  Log_Solenoid1->SetVisAttributes(new G4VisAttributes(green));

  ///////////////////////////////
  /// Magnetic field region 1 ///
  ///////////////////////////////

  // Creating Solenoid 1 magnetic field
  G4double Mag_diameter = 30.0 * cm;
  G4double Mag_length = 68.0 * cm; //coil length

  Sol_Magnet1 = new G4Tubs("Sol_Magnet1", 0., Mag_diameter / 2.0, Mag_length / 2.0, 0., 360. * deg);

  Log_Magnet1 = new G4LogicalVolume(Sol_Magnet1, G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"), "Log_Magnet1", 0, 0, 0);

  Phys_Magnet1 = new G4PVPlacement(0, Pos_Solenoid, Log_Magnet1, "Magnetic Field 1", logicWorld, false, 0, true);

  // Set "user limits" for drawing smooth curve
  G4UserLimits *userLimits = new G4UserLimits(1.e-5 * mm);
  Log_Magnet1->SetUserLimits(userLimits);

  //------- Building Target -------//

  //////////////
  //  Target  //
  //////////////

  // Target position
  G4ThreeVector Target_pos = Inputs->target_pos;

  // Target Material
  G4Material *TargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Inputs->g4_material_name);
  Inputs->TargetMaterial = TargetMaterial;

  // Rotating target
  G4RotationMatrix *rotation = new G4RotationMatrix;
  rotation->rotateZ(180. * CLHEP::deg);

  // Creating Target
  Sol_Target = new G4Box("Sol_Target", 1.5 * CLHEP::cm, 1.5 * CLHEP::cm, Inputs->width / 2 * mm);
  Log_Target = new G4LogicalVolume(Sol_Target, man->FindOrBuildMaterial("G4_LITHIUM_FLUORIDE"), "Log_Target");

  // Placing Target
  Phys_Target = new G4PVPlacement(rotation, Target_pos, Log_Target, "Target", Log_Magnet1, false, 8, true);

  // Color
  Log_Target->SetVisAttributes(new G4VisAttributes(red));

  // Setting target sensitive
  Log_Target->SetSensitiveDetector(Sensitive);

  //--------- Setting magnetic field ---------//

  ////////////////////
  // Magnetic field //
  ////////////////////

  G4bool fieldIsInitialized = false;
  if (!fieldIsInitialized)
  {
    G4MagIntegratorStepper *fStepper;

    G4FieldManager *fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magneticField);
    fieldMgr->CreateChordFinder(magneticField);

    G4double minEps = 1.0e-10 * cm; //   Minimum & value for smallest steps
    G4double maxEps = 1.0e-8 * cm;  //   Maximum & value for largest steps

    fieldMgr->SetMinimumEpsilonStep(minEps);
    fieldMgr->SetMaximumEpsilonStep(maxEps);

    fieldMgr->GetChordFinder()->SetDeltaChord(0.001 * mm);

    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(magneticField);

    // Note that for magnetic field that do not vary with time,
    fStepper = new G4HelixExplicitEuler(fEquation);

    fieldMgr->SetDeltaIntersection(0.1 * mm);
    fieldMgr->SetAccuraciesWithDeltaOneStep(0.01 * mm);
    fieldMgr->SetDeltaOneStep(0.01 * mm); // 0.5 micrometer

    fieldIsInitialized = true;
  }

  //-------------------------------------- Building All Detectors ----------------------------------------//

  ///////////////
  // Detectors //
  ///////////////

  Detector Detector_0(0, Lengthy_dssd_t1, Lengthx_dssd_t1, Thickness_dssd_t1, noOfSensorStrips);
  Detector_0.Rotate(0., 90., 0.);
  Detector_0.SetPosition(DetectorPosition[0]);
  Detector_0.SetColor(cyan);
  Detector_0.Construct(Log_Magnet1); // This method requires a Logical Mother Volume

  Detector Detector_1(1, Lengthy_dssd_t1, Lengthx_dssd_t1, Thickness_dssd_t1, noOfSensorStrips);
  Detector_1.Rotate(0., 90., 0.);
  Detector_1.SetPosition(DetectorPosition[1]);
  Detector_1.SetColor(cyan);
  Detector_1.Construct(Log_Magnet1); // This method requires a Logical Mother Volume

  Detector Detector_2(2, Lengthy_dssd_t1, Lengthx_dssd_t1, Thickness_dssd_t1, noOfSensorStrips);
  Detector_2.Rotate(90., 0., 90.);
  Detector_2.SetPosition(DetectorPosition[2]);
  Detector_2.SetColor(cyan);
  Detector_2.Construct(Log_Magnet1); // This method requires a Logical Mother Volume

  Detector Detector_3(3, Lengthy_dssd_t1, Lengthx_dssd_t1, Thickness_dssd_t1, noOfSensorStrips);
  Detector_3.Rotate(90., 0., 90.);
  Detector_3.SetPosition(DetectorPosition[3]);
  Detector_3.SetColor(cyan);
  Detector_3.Construct(Log_Magnet1); // This method requires a Logical Mother Volume

  Detector Detector_4(4, Lengthy_dssd_t1, Lengthx_dssd_t1, Thickness_dssd_t1, noOfSensorStrips);
  Detector_4.Rotate(90., 0., 90.);
  Detector_4.SetPosition(DetectorPosition[4]);
  Detector_4.SetColor(cyan);
  Detector_4.Construct(Log_Magnet1); // This method requires a Logical Mother Volume

  Detector Detector_5(5, Lengthy_dssd_t1, Lengthx_dssd_t1, Thickness_dssd_t1, noOfSensorStrips);
  Detector_5.Rotate(90., 0., 90.);
  Detector_5.SetPosition(DetectorPosition[5]);
  Detector_5.SetColor(cyan);
  Detector_5.Construct(Log_Magnet1); // This method requires a Logical Mother Volume

  Detector Detector_6(6, Lengthy_dssd_t1, Lengthx_dssd_t1, Thickness_dssd_t1, noOfSensorStrips);
  Detector_6.Rotate(0., 90., 0.);
  Detector_6.SetPosition(DetectorPosition[6]);
  Detector_6.SetColor(cyan);
  Detector_6.Construct(Log_Magnet1); // This method requires a Logical Mother Volume

  Detector Detector_7(7, Lengthy_dssd_t1, Lengthx_dssd_t1, Thickness_dssd_t1, noOfSensorStrips);
  Detector_7.Rotate(0., 90., 0.);
  Detector_7.SetPosition(DetectorPosition[7]);
  Detector_7.SetColor(cyan);
  Detector_7.Construct(Log_Magnet1); // This method requires a Logical Mother Volume

  //--------------------------------------------------------------------------------------------------------------//

  // Return World!
  return physiWorld;
}

//--------------------------------------------------------------------------------------------------------------//

#include "G4RunManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

void DetectorConstruction::UpdateGeometry()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Reconstruct world with new parameters
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
//--------------------------------------------------------------------------------------------------------------//