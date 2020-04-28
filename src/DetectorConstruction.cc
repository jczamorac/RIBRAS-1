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
  CD2 = new G4Material(name = "CD2", 0.94 *g/cm3, ncomponents = 2, kStateSolid);/* , temperature, pressure); */
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
  //lead  = man->FindOrBuildMaterial("G4_Pb");
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
  Lengthy_dssd_t1 = 10.0 * CLHEP::cm;
  Thickness_dssd_t1 = 300. * CLHEP::um;
  Lengthx_dssd_t1 = 30.0 * CLHEP::cm;
}

//--------------------------------------------------------------------------------------------------------------//

void DetectorConstruction::ConstructSetup(void)
{
}

//--------------------------------------------------------------------------------------------------------------//

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  Inputs *Inputs = &Inputs::GetInputs();

  // Available colors
  G4Color green(0.0, 1.0, 0.0),
      red(1.0, 0.0, 0.0),
      yellow(1.0, 1.0, 0.0),
      orange(1.0, 1.0, 0.5),
      blue(0.0, 0.0, 1.0);

  //---- This function is called by G4 when the detector has to be created -----//
  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------//

  //----------------------------------------------------//
  // World - Aqui o "mundo" é criado. O mundo é Lógico. //
  //----------------------------------------------------//

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2. * halfWorldLength);

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

  //------ Building Second Solenoid ------//

  ///////////////////
  //  Solenoid 2  //
  ///////////////////

  // Solenoid 2 parameters
  G4double raiointerno = 30.05 * CLHEP::cm;
  G4double raioexterno = 100.0 * CLHEP::cm;
  G4double comprimento = 100.0 * CLHEP::cm;

  // Solenoid 2 Position
  G4ThreeVector position = G4ThreeVector(0., 0., 301. * CLHEP::cm);

  // Creating Solenoid 2
  Sol_Solenoid2 = new G4Tubs("Sol_Solenoid2", raiointerno / 2.0, raioexterno / 2.0, comprimento / 2.0, 0., 360. * CLHEP::deg);
  Log_Solenoid2 = new G4LogicalVolume(Sol_Solenoid2, G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe"), "Log_Solenoid2");

  // Placing Splenoid 2
  Phys_Solenoid2 = new G4PVPlacement(0, position, Log_Solenoid2, "Solenoid 2", logicWorld, false, 0, true);

  // Color
  Log_Solenoid2->SetVisAttributes(new G4VisAttributes(green));

  ///////////////////////////////
  /// Magnetic field region 2 ///
  ///////////////////////////////

  // Creating Solenoid 2 magnetic field
  G4double diametromag = 30.0 * cm;
  G4double comprimentomag = 68.0 * cm; //coil length

  Sol_Magnet2 = new G4Tubs("Sol_Magnet2", 0., diametromag / 2.0, comprimentomag / 2.0, 0., 360. * deg);

  Log_Magnet2 = new G4LogicalVolume(Sol_Magnet2, G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"), "Log_Magnet2");

  Phys_Magnet2 = new G4PVPlacement(0,
                                   position,
                                   Log_Magnet2,
                                   "Magnetic Field 2",
                                   logicWorld,
                                   false,
                                   0,
                                   true);

  // Set "user limits" for drawing smooth curve
  G4UserLimits *limiteuser = new G4UserLimits(1.e-5 * mm);
  Log_Magnet2->SetUserLimits(limiteuser);

  //------- Building Target -------//

  //////////////
  //  Target  //
  //////////////

  // Target position
  G4ThreeVector Target_pos = Inputs->target_pos;

  // Target Material
  G4Material *TargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Inputs->g4_material_name);
  Inputs->TargetMaterial = TargetMaterial;

  // Creating Target
  Sol_Target = new G4Box("Sol_Target", 7.5 * CLHEP::cm, 7.5 * CLHEP::cm, Inputs->width * CLHEP::mm);
  Log_Target = new G4LogicalVolume(Sol_Target, CD2, "Log_Target");

  // Placing Target
  Phys_Target = new G4PVPlacement(0, Target_pos, Log_Target, "Target", Log_Magnet2, false, 0, true);

  // Color
  Log_Target->SetVisAttributes(new G4VisAttributes(red));

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
    fStepper = new G4HelixExplicitEuler(fEquation); // mais indicado para campos intensos e não suaves
    //fStepper = new G4ClassicalRK4( fEquation );       // integrador mais usado

    fieldMgr->SetDeltaIntersection(0.1 * mm);
    fieldMgr->SetAccuraciesWithDeltaOneStep(0.01 * mm);
    fieldMgr->SetDeltaOneStep(0.01 * mm); // 0.5 micrometer

    fieldIsInitialized = true;
  }

  //--------- Building All Detectors -----------//

  ///////////////
  // Detectors //
  ///////////////

  // This method create the detectors
  ConstructDetectors();

  return physiWorld;
}

//--------------------------------------------------------------------------------------------------------------//

///////////////
// Detectors //
///////////////

/* G4VPhysicalVolume * */ void DetectorConstruction::ConstructDetectors()
{
  // Constructing detectors
  Construct_D_0_0();
  Construct_D_0_1();
  Construct_D_0_2();
  Construct_D_0_3();
  Construct_D_0_4();
  Construct_D_0_5();
  Construct_D_0_6();
  Construct_D_0_7();

  static SensitiveDetector *sensitive = 0;
  if (!sensitive)
  {
    sensitive = new SensitiveDetector("/myDet/SiStripSD");
    //We register now the SD with the manager
    G4SDManager::GetSDMpointer()->AddNewDetector(sensitive);
  }

  // These commands sets the detectors sensitive
  logicSensorStripD00->SetSensitiveDetector(sensitive);
  logicSensorStripD01->SetSensitiveDetector(sensitive);
  logicSensorStripD02->SetSensitiveDetector(sensitive);
  logicSensorStripD03->SetSensitiveDetector(sensitive);
  logicSensorStripD04->SetSensitiveDetector(sensitive);
  logicSensorStripD05->SetSensitiveDetector(sensitive);
  logicSensorStripD06->SetSensitiveDetector(sensitive);
  logicSensorStripD07->SetSensitiveDetector(sensitive);
  Log_Target->SetSensitiveDetector(sensitive);
}

//--------------------------------------------------------------------------------------------------------------//

////////////////
/////////////////////////////// Detector 0 ///////////////////////////////////
////////////////

G4VPhysicalVolume *DetectorConstruction::Construct_D_0_0()
{

  // Retrieving Inputs
  Inputs *Inputs = &Inputs::GetInputs();

  // Detector 1 position
  G4ThreeVector posicao_detector1 = Inputs->detector1_pos;

  // Creating the detector
  G4Box *Detector_0 = new G4Box("Detector0", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *DetectorLog_0 = new G4LogicalVolume(Detector_0, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector1");

  // Rotating detector
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateY(90. * CLHEP::deg);

  // Placing detector 1
  G4VPhysicalVolume *DetectorPhys_0 = new G4PVPlacement(rotacaox,
                                                        posicao_detector1,
                                                        DetectorLog_0,
                                                        "Detector_0",
                                                        Log_Magnet2,
                                                        false,
                                                        0,
                                                        true);
  //////////////
  //  Strips  //
  //////////////

  // Creating Strips
  // Strips parameters
  G4double halfSensorStripSizeX = Lengthx_dssd_t1 / (2.0 * noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1 / 2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1 / 2.;

  // Creating physical strips
  G4Box *solidSensorStripD00 = new G4Box("SensorStripD00", halfSensorStripSizeX, halfSensorStripSizeY, halfSensorStripSizeZ);
  logicSensorStripD00 = new G4LogicalVolume(solidSensorStripD00, silicon, "SensorStripD00");

  // Placing Strips
  G4VPhysicalVolume *physiSensorStripD_0_0 = new G4PVReplica("SensorStripD00",    //its name
                                                             logicSensorStripD00, //its logical volume
                                                             DetectorLog_0,       //its mother
                                                             kXAxis,              //axis of replication
                                                             noOfSensorStrips,    //number of replica
                                                             2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);          //witdth of replica

  // Set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD00->SetUserLimits(userLimits);

  // Colors
  G4Color yellow(1.0, 1.0, 0.0), red(1.0, 0.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  DetectorLog_0->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD00->SetVisAttributes(new G4VisAttributes(green));

  return DetectorPhys_0;
}

//--------------------------------------------------------------------------------------------------------------//

////////////////
/////////////////////////////// Detector 1 ///////////////////////////////////
////////////////

G4VPhysicalVolume *DetectorConstruction::Construct_D_0_1()
{

  // Retrieving Inputs
  Inputs *Inputs = &Inputs::GetInputs();

  // Detector position
  G4ThreeVector posicao_detector1 = Inputs->detector2_pos;

  // Creating Detector
  G4Box *solidSensor_D_0_1 = new G4Box("SensorD_0_1", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *logicSensorPlane = new G4LogicalVolume(solidSensor_D_0_1, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "SensorL_D_0_1");

  // Rotating detectors
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateY(90.0 * CLHEP::deg);

  // Placing Detector 2
  G4VPhysicalVolume *DetectorPhys_1 = new G4PVPlacement(rotacaox,
                                                        posicao_detector1,
                                                        logicSensorPlane,
                                                        "Detector_1",
                                                        Log_Magnet2,
                                                        false,
                                                        1,
                                                        true);

  //////////////
  //  Strips  //
  //////////////

  // Creating Strips for the second detector
  G4double halfSensorStripSizeX = Lengthx_dssd_t1 / (2.0 * noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1 / 2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1 / 2.;

  G4Box *solidSensorStripD01 = new G4Box("SensorStripD01", halfSensorStripSizeX, halfSensorStripSizeY, halfSensorStripSizeZ);

  logicSensorStripD01 = new G4LogicalVolume(solidSensorStripD01, silicon, "SensorStripD01");

  G4VPhysicalVolume *PhysiSensorStripD_0_1 = new G4PVReplica("SensorStripD01",    //its name
                                                             logicSensorStripD01, //its logical volume
                                                             logicSensorPlane,    //its mother
                                                             kXAxis,              //axis of replication
                                                             noOfSensorStrips,    //number of replica
                                                             2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);              //witdth of replica

  // set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD01->SetUserLimits(userLimits);

  G4Color red(1.0, 0.0, 0.0), yellow(1.0, 1.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  logicSensorPlane->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD01->SetVisAttributes(new G4VisAttributes(green));

  return DetectorPhys_1;
}

////////////////
/////////////////////////////// Detector 2 ///////////////////////////////////
////////////////

G4VPhysicalVolume *DetectorConstruction::Construct_D_0_2()
{

  // Detector Position
  G4ThreeVector posicao_detector2 = G4ThreeVector(7.5 * CLHEP::cm, 0., 18. * CLHEP::cm);

  // Building detectors
  G4Box *detector2 = new G4Box("detector3", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *detector2log = new G4LogicalVolume(detector2, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector2");

  // Rotating
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateY(90. * CLHEP::deg);

  // Placing detector
  G4VPhysicalVolume *DetectorPhys_2 = new G4PVPlacement(rotacaox,
                                                        posicao_detector2,
                                                        detector2log,
                                                        "Detector_2",
                                                        Log_Magnet2,
                                                        false,
                                                        2,
                                                        true);
  //////////////
  //  Strips  //
  //////////////

  G4double halfSensorStripSizeX = Lengthx_dssd_t1 / (2.0 * noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1 / 2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1 / 2.;

  G4Box *solidSensorStripD02 = new G4Box("SensorStripD02", halfSensorStripSizeX, halfSensorStripSizeY, halfSensorStripSizeZ);

  logicSensorStripD02 = new G4LogicalVolume(solidSensorStripD02, silicon, "SensorStripD02");

  G4VPhysicalVolume *physiSensorStripD_0_2 = new G4PVReplica("SensorStripD02",    //its name
                                                             logicSensorStripD02, //its logical volume
                                                             detector2log,        //its mother
                                                             kXAxis,              //axis of replication
                                                             noOfSensorStrips,    //number of replica
                                                             2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);          //witdth of replica

  // set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD02->SetUserLimits(userLimits);

  G4Color yellow(1.0, 1.0, 0.0), red(1.0, 0.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  detector2log->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD02->SetVisAttributes(new G4VisAttributes(green));

  return DetectorPhys_2;
}

////////////////
/////////////////////////////// Detector 3 ///////////////////////////////////
////////////////

G4VPhysicalVolume *DetectorConstruction::Construct_D_0_3()
{
  // Detector 1 position
  G4ThreeVector posicao_detector1 = G4ThreeVector(7.5 * cm, 0., -18. * cm);

  // Creating the detector
  G4Box *Detector_3 = new G4Box("Detector3", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *DetectorLog_3 = new G4LogicalVolume(Detector_3, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector3");

  // Rotating detector
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateY(90. * CLHEP::deg);

  // Placing detector 1
  G4VPhysicalVolume *DetectorPhys_3 = new G4PVPlacement(rotacaox,
                                                        posicao_detector1,
                                                        DetectorLog_3,
                                                        "Detector_3",
                                                        Log_Magnet2,
                                                        false,
                                                        3,
                                                        true);
  //////////////
  //  Strips  //
  //////////////

  // Creating Strips
  // Strips parameters
  G4double halfSensorStripSizeX = Lengthx_dssd_t1 / (2.0 * noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1 / 2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1 / 2.;

  // Creating physical strips
  G4Box *solidSensorStripD03 = new G4Box("SensorStripD03", halfSensorStripSizeX, halfSensorStripSizeY, halfSensorStripSizeZ);
  logicSensorStripD03 = new G4LogicalVolume(solidSensorStripD03, silicon, "SensorStripD03");

  // Placing Strips
  G4VPhysicalVolume *physiSensorStripD_0_3 = new G4PVReplica("SensorStripD03",    //its name
                                                             logicSensorStripD03, //its logical volume
                                                             DetectorLog_3,       //its mother
                                                             kXAxis,              //axis of replication
                                                             noOfSensorStrips,    //number of replica
                                                             2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);          //witdth of replica

  // Set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD03->SetUserLimits(userLimits);

  // Colors
  G4Color yellow(1.0, 1.0, 0.0), red(1.0, 0.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  DetectorLog_3->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD03->SetVisAttributes(new G4VisAttributes(green));

  return DetectorPhys_3;
}

////////////////
/////////////////////////////// Detector 4 ///////////////////////////////////
////////////////

G4VPhysicalVolume *DetectorConstruction::Construct_D_0_4()
{
  // Detector 1 position
  G4ThreeVector posicao_detector1 = G4ThreeVector(0., 7.5 * cm, -18. * cm);

  // Creating the detector
  G4Box *Detector_4 = new G4Box("Detector4", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *DetectorLog_4 = new G4LogicalVolume(Detector_4, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector4");

  // Rotating detector
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateX(90. * deg);
  rotacaox->rotateZ(90. * deg);

  // Placing detector 1
  G4VPhysicalVolume *DetectorPhys_4 = new G4PVPlacement(rotacaox,
                                                        posicao_detector1,
                                                        DetectorLog_4,
                                                        "Detector_4",
                                                        Log_Magnet2,
                                                        false,
                                                        4,
                                                        true);
  //////////////
  //  Strips  //
  //////////////

  // Creating Strips
  // Strips parameters
  G4double halfSensorStripSizeX = Lengthx_dssd_t1 / (2.0 * noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1 / 2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1 / 2.;

  // Creating physical strips
  G4Box *solidSensorStripD04 = new G4Box("SensorStripD04", halfSensorStripSizeX, halfSensorStripSizeY, halfSensorStripSizeZ);
  logicSensorStripD04 = new G4LogicalVolume(solidSensorStripD04, silicon, "SensorStripD04");

  // Placing Strips
  G4VPhysicalVolume *physiSensorStripD_0_4 = new G4PVReplica("SensorStripD04",    //its name
                                                             logicSensorStripD04, //its logical volume
                                                             DetectorLog_4,       //its mother
                                                             kXAxis,              //axis of replication
                                                             noOfSensorStrips,    //number of replica
                                                             2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);          //witdth of replica

  // Set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD04->SetUserLimits(userLimits);

  // Colors
  G4Color yellow(1.0, 1.0, 0.0), red(1.0, 0.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  DetectorLog_4->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD04->SetVisAttributes(new G4VisAttributes(green));

  return DetectorPhys_4;
}

////////////////
/////////////////////////////// Detector 5 ///////////////////////////////////
////////////////

G4VPhysicalVolume *DetectorConstruction::Construct_D_0_5()
{
  // Detector 1 position
  G4ThreeVector posicao_detector1 = G4ThreeVector(0., -7.5 * cm, -18. * cm);

  // Creating the detector
  G4Box *Detector_5 = new G4Box("Detector5", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *DetectorLog_5 = new G4LogicalVolume(Detector_5, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector5");

  // Rotating detector
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateX(90. * deg);
  rotacaox->rotateZ(90. * deg);

  // Placing detector 1
  G4VPhysicalVolume *DetectorPhys_5 = new G4PVPlacement(rotacaox,
                                                        posicao_detector1,
                                                        DetectorLog_5,
                                                        "Detector_5",
                                                        Log_Magnet2,
                                                        false,
                                                        5,
                                                        true);
  //////////////
  //  Strips  //
  //////////////

  // Creating Strips
  // Strips parameters
  G4double halfSensorStripSizeX = Lengthx_dssd_t1 / (2.0 * noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1 / 2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1 / 2.;

  // Creating physical strips
  G4Box *solidSensorStripD05 = new G4Box("SensorStripD05", halfSensorStripSizeX, halfSensorStripSizeY, halfSensorStripSizeZ);
  logicSensorStripD05 = new G4LogicalVolume(solidSensorStripD05, silicon, "SensorStripD05");

  // Placing Strips
  G4VPhysicalVolume *physiSensorStripD_0_5 = new G4PVReplica("SensorStripD05",    //its name
                                                             logicSensorStripD05, //its logical volume
                                                             DetectorLog_5,       //its mother
                                                             kXAxis,              //axis of replication
                                                             noOfSensorStrips,    //number of replica
                                                             2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);          //witdth of replica

  // Set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD05->SetUserLimits(userLimits);

  // Colors
  G4Color yellow(1.0, 1.0, 0.0), red(1.0, 0.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  DetectorLog_5->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD05->SetVisAttributes(new G4VisAttributes(green));

  return DetectorPhys_5;
}

////////////////
/////////////////////////////// Detector 6 ///////////////////////////////////
////////////////

G4VPhysicalVolume *DetectorConstruction::Construct_D_0_6()
{
  // Detector 1 position
  G4ThreeVector posicao_detector1 = G4ThreeVector(0., -7.5 * cm, 18. * cm);

  // Creating the detector
  G4Box *Detector_6 = new G4Box("Detector6", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *DetectorLog_6 = new G4LogicalVolume(Detector_6, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector6");

  // Rotating detector
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateX(90. * deg);
  rotacaox->rotateZ(90. * deg);

  // Placing detector 1
  G4VPhysicalVolume *DetectorPhys_6 = new G4PVPlacement(rotacaox,
                                                        posicao_detector1,
                                                        DetectorLog_6,
                                                        "Detector_6",
                                                        Log_Magnet2,
                                                        false,
                                                        6,
                                                        true);
  //////////////
  //  Strips  //
  //////////////

  // Creating Strips
  // Strips parameters
  G4double halfSensorStripSizeX = Lengthx_dssd_t1 / (2.0 * noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1 / 2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1 / 2.;

  // Creating physical strips
  G4Box *solidSensorStripD06 = new G4Box("SensorStripD06", halfSensorStripSizeX, halfSensorStripSizeY, halfSensorStripSizeZ);
  logicSensorStripD06 = new G4LogicalVolume(solidSensorStripD06, silicon, "SensorStripD06");

  // Placing Strips
  G4VPhysicalVolume *physiSensorStripD_0_6 = new G4PVReplica("SensorStripD06",    //its name
                                                             logicSensorStripD06, //its logical volume
                                                             DetectorLog_6,       //its mother
                                                             kXAxis,              //axis of replication
                                                             noOfSensorStrips,    //number of replica
                                                             2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);          //witdth of replica

  // Set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD06->SetUserLimits(userLimits);

  // Colors
  G4Color yellow(1.0, 1.0, 0.0), red(1.0, 0.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  DetectorLog_6->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD06->SetVisAttributes(new G4VisAttributes(green));

  return DetectorPhys_6;
}

////////////////
/////////////////////////////// Detector 7 ///////////////////////////////////
////////////////

G4VPhysicalVolume *DetectorConstruction::Construct_D_0_7()
{
  // Detector 1 position
  G4ThreeVector posicao_detector1 = G4ThreeVector(0., 7.5 * cm, 18. * cm);

  // Creating the detector
  G4Box *Detector_7 = new G4Box("Detector7", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *DetectorLog_7 = new G4LogicalVolume(Detector_7, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector7");

  // Rotating detector
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateX(90. * deg);
  rotacaox->rotateZ(90. * deg);

  // Placing detector 1
  G4VPhysicalVolume *DetectorPhys_7 = new G4PVPlacement(rotacaox,
                                                        posicao_detector1,
                                                        DetectorLog_7,
                                                        "Detector_7",
                                                        Log_Magnet2,
                                                        false,
                                                        7,
                                                        true);
  //////////////
  //  Strips  //
  //////////////

  // Creating Strips
  // Strips parameters
  G4double halfSensorStripSizeX = Lengthx_dssd_t1 / (2.0 * noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1 / 2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1 / 2.;

  // Creating physical strips
  G4Box *solidSensorStripD07 = new G4Box("SensorStripD07", halfSensorStripSizeX, halfSensorStripSizeY, halfSensorStripSizeZ);
  logicSensorStripD07 = new G4LogicalVolume(solidSensorStripD07, silicon, "SensorStripD07");

  // Placing Strips
  G4VPhysicalVolume *physiSensorStripD_0_7 = new G4PVReplica("SensorStripD07",    //its name
                                                             logicSensorStripD07, //its logical volume
                                                             DetectorLog_7,       //its mother
                                                             kXAxis,              //axis of replication
                                                             noOfSensorStrips,    //number of replica
                                                             2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);          //witdth of replica

  // Set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD07->SetUserLimits(userLimits);

  // Colors
  G4Color yellow(1.0, 1.0, 0.0), red(1.0, 0.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  DetectorLog_7->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD07->SetVisAttributes(new G4VisAttributes(green));

  return DetectorPhys_7;
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
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
