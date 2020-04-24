// $Id: DetectorConstruction.cc 94 2010-01-26 13:18:30Z adotti $
/**
* @file
* @brief Implements mandatory user class DetectorConstruction.
*/

// Default libraries
#include <fstream>

// Local headers
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "MagneticField.hh"
#include "SensitiveDetector.hh"

// Geant4 headers
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

DetectorConstruction::~DetectorConstruction()
{
  delete messenger;
}

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
}

void DetectorConstruction::ComputeParameters()
{
  // This function defines the defaults of the geometry construction
  // World size
  halfWorldLength = 20.0 * CLHEP::m;

  // Number of strips in the detector
  noOfSensorStrips = 60;

  // Detector Parameters
  Lengthy_dssd_t1 = 5.0 * CLHEP::cm;
  Thickness_dssd_t1 = 300. * CLHEP::um;
  Lengthx_dssd_t1 = 30.0 * CLHEP::cm;

}

void DetectorConstruction::ConstructSetup(void)
{
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  Inputs *Inputs = &Inputs::GetInputs();

  //This function is called by G4 when the detector has to be created
  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

  //------------------------------
  // World - Aqui o "mundo" é criado. O mundo é Lógico.
  //------------------------------

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
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);          //Turn the world invisible

  ///////////////////
  //  Solenoide 1  //
  ///////////////////

  // Solenoid 1 parameters
  G4double Solenoid_length = 100.0 * CLHEP::cm;
  G4double Solenoid_diameter_inner = 30.05 * CLHEP::cm;
  G4double Solenoid_diameter_outer = 100.0 * CLHEP::cm;

  // Position
  G4ThreeVector Pos_Solenoid = G4ThreeVector();

  // Creating Solenoid 1
  G4VSolid *Sol_Solenoid = new G4Tubs("magneticexterior", Solenoid_diameter_inner / 2.0, Solenoid_diameter_outer / 2.0, Solenoid_length / 2.0, 0., 360. * CLHEP::deg);

  G4LogicalVolume *Log_Solenoid = new G4LogicalVolume(Sol_Solenoid, G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe"), "Log_Solenoid");

  // Placing Solenoid 1
  Phy_Solenoid = new G4PVPlacement(0, Pos_Solenoid, Log_Solenoid, "P_Solenoid", logicWorld, false, 0, true);

  // Available colors
  G4Color green(0.0, 1.0, 0.0),
          red(1.0, 0.0, 0.0),
          yellow(1.0, 1.0, 0.0),
          orange(1.0, 1.0, 0.5),
          blue(0.0, 0.0, 1.0);

  // Setting Solenoid 1 color
  Log_Solenoid->SetVisAttributes(new G4VisAttributes(green));

  ///////////////////////////////
  /// Magnetic field region 1 ///
  ///////////////////////////////

  // Creating Solenoid 1 magnetic field
  G4double Mag_diameter = 30.0 * cm;
  G4double Mag_length = 68.0 * cm; //coil length

  G4VSolid *Sol_Magnet = new G4Tubs("S_Magnet", 0., Mag_diameter / 2.0, Mag_length / 2.0, 0., 360. * deg);

  G4LogicalVolume *Log_Magnet = new G4LogicalVolume(Sol_Magnet, G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"), "Log_Magnet", 0, 0, 0);

  G4VPhysicalVolume *Phy_Magnet = new G4PVPlacement(0, Pos_Solenoid, Log_Magnet, "P_Magnet", logicWorld, false, 0, true);

  // Set "user limits" for drawing smooth curve
  G4UserLimits *userLimits = new G4UserLimits(1.e-5 * mm);
  Log_Magnet->SetUserLimits(userLimits);

  ///////////////////
  //  Solenoide 2  //
  ///////////////////

  // Creating Solenoid 2
  G4ThreeVector position = G4ThreeVector(0., 0., 301. * CLHEP::cm);
  G4double raiointerno = 30.05 * CLHEP::cm;
  G4double raioexterno = 100.0 * CLHEP::cm;
  G4double comprimento = 100.0 * CLHEP::cm;

  G4VSolid *solenoide = new G4Tubs("solenoide", raiointerno / 2.0, raioexterno / 2.0, comprimento / 2.0, 0., 360. * CLHEP::deg);

  G4LogicalVolume *logsolenoide = new G4LogicalVolume(solenoide, G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe"), "logsolenoide");

  solenoidef = new G4PVPlacement(0, position, logsolenoide, "solenoidef", logicWorld, false, 0, true);

  logsolenoide->SetVisAttributes(new G4VisAttributes(green));

  ///////////////////////////////
  /// Magnetic field region 2 ///
  ///////////////////////////////

  //Creating Solenoid 2 magnetic field
  G4double diametromag = 30.0 * cm;
  G4double comprimentomag = 68.0 * cm; //coil length

  G4VSolid *magsolenoide = new G4Tubs("magsolenoide", 0., diametromag / 2.0, comprimentomag / 2.0, 0., 360. * deg);

  logmagnetico = new G4LogicalVolume(magsolenoide, G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"), "logmagnetico");

  G4VPhysicalVolume *magneticof = new G4PVPlacement(0,
                                                    position,
                                                    logmagnetico,
                                                    "magneticof",
                                                    logicWorld,
                                                    false,
                                                    0,
                                                    true);

  // Set "user limits" for drawing smooth curve
  G4UserLimits *limiteuser = new G4UserLimits(1.e-5 * mm);
  logmagnetico->SetUserLimits(limiteuser);

  //////////////
  //  Target  //
  //////////////

  // Target position
  G4ThreeVector Target_pos = Inputs->target_pos;

  // Target Material
  G4Material *TargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Inputs->g4_material_name);
  Inputs->TargetMaterial = TargetMaterial;

  // Creating Target
  G4Box *alvo = new G4Box("alvo", 7.5 * CLHEP::cm, 7.5 * CLHEP::cm, Inputs->width * CLHEP::mm);
  G4LogicalVolume *alvolog = new G4LogicalVolume(alvo, TargetMaterial, "alvolog");

  // Placing Target
  G4VPhysicalVolume *alvophys = new G4PVPlacement(0, Target_pos, alvolog, "alvof", logmagnetico, false, 0, true);

  // Color
  alvolog->SetVisAttributes(new G4VisAttributes(red));

  // This method create the detectors
  ConstructDetectors();

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


  return physiWorld;
}

////////////////
// Detectores //
////////////////

G4VPhysicalVolume *DetectorConstruction::ConstructDetectors()
{

  // Constructing detectors
  ConstructRing_D_0_0();
  ConstructRing_D_0_1();
  ConstructRing_D_0_2();

  static SensitiveDetector *sensitive = 0;
  if (!sensitive)
  {

    sensitive = new SensitiveDetector("/myDet/SiStripSD");
    //We register now the SD with the manager
    G4SDManager::GetSDMpointer()->AddNewDetector(sensitive);
  }

  // dos lineas para anexar un detector diferente en una replica de detectores

  const G4LogicalVolume *log = detector1phys->GetLogicalVolume();
  log->GetDaughter(0)->GetLogicalVolume()->SetSensitiveDetector(sensitive);

  const G4LogicalVolume *ring_sen_D01 = physiSensorRing_D_0_1->GetLogicalVolume();
  ring_sen_D01->GetDaughter(0)->GetLogicalVolume()->SetSensitiveDetector(sensitive);

  const G4LogicalVolume *ring_sen_D02 = detector2phys->GetLogicalVolume();
  ring_sen_D02->GetDaughter(0)->GetLogicalVolume()->SetSensitiveDetector(sensitive);

  return detector1phys;
}

                             ////////////////
/////////////////////////////// Detector 1 ///////////////////////////////////
                             ////////////////

G4VPhysicalVolume *DetectorConstruction::ConstructRing_D_0_0()
{

  // Retrieving Inputs
  Inputs *Inputs = &Inputs::GetInputs();

  // Detector 1 position
  G4ThreeVector posicao_detector1 = Inputs->detector1_pos;

  // Creating the detector
  G4Box *detector1 = new G4Box("detector1", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *detector1log = new G4LogicalVolume(detector1, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector1");

  // Rotating detector
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateY(90. * CLHEP::deg);

  // Placing detector 1
  detector1phys = new G4PVPlacement(rotacaox,
                                    posicao_detector1,
                                    detector1log,
                                    "Detector_1",
                                    logmagnetico,
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
  G4LogicalVolume *logicSensorStripD00 = new G4LogicalVolume(solidSensorStripD00, silicon, "SensorStripD00");

  // Placing Strips
  physiSensorStripD_0_0 = new G4PVReplica("SensorStripD00",    //its name
                                          logicSensorStripD00, //its logical volume
                                          detector1log,        //its mother
                                          kXAxis,              //axis of replication
                                          noOfSensorStrips,    //number of replica
                                          2.0 * halfSensorStripSizeX);
  //                                      Lengthx_sili);          //witdth of replica

  // Set "user limits" StepSize
  G4UserLimits *userLimits = new G4UserLimits(0.1 * CLHEP::mm);
  logicSensorStripD00->SetUserLimits(userLimits);

  // Colors
  G4Color yellow(1.0, 1.0, 0.0), red(1.0, 0.0, 0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  detector1log->SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD00->SetVisAttributes(new G4VisAttributes(green));

  return detector1phys;
}

                             ////////////////
/////////////////////////////// Detector 2 ///////////////////////////////////
                             ////////////////

G4VPhysicalVolume *DetectorConstruction::ConstructRing_D_0_1()
{

  // Retrieving Inputs
  Inputs *Inputs = &Inputs::GetInputs();

  // Detector position
  G4ThreeVector posicao_detector2 = Inputs->detector2_pos;

  // Creating Detector
  G4Box *solidSensor_D_0_1 = new G4Box("SensorD_0_1", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *logicSensorPlane = new G4LogicalVolume(solidSensor_D_0_1, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "SensorL_D_0_1");

  // Rotating detectors
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateY(90.0 * CLHEP::deg);

  // Placing Detector 2
  physiSensorRing_D_0_1 = new G4PVPlacement(rotacaox,
                                            posicao_detector2,
                                            logicSensorPlane,
                                            "Detector_2",
                                            logmagnetico,
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

  G4LogicalVolume *logicSensorStripD01 = new G4LogicalVolume(solidSensorStripD01, silicon, "SensorStripD01");

  physiSensorStripD_0_1 = new G4PVReplica("SensorStripD01",    //its name
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
  logicSensorStripD01->SetVisAttributes(new G4VisAttributes(red));

  return physiSensorRing_D_0_1;
}

<<<<<<< HEAD
=======
G4VPhysicalVolume *DetectorConstruction::ConstructRing_D_0_2()
{

  // Retrieving Inputs
  Inputs *Inputs = &Inputs::GetInputs();

  // Coordenadas do detector 2
  G4ThreeVector posicao_detector2 = G4ThreeVector(7.5, 0., -18);

  // Criando o detector
  G4Box *detector2 = new G4Box("detector3", Lengthx_dssd_t1 / 2., Lengthy_dssd_t1 / 2., Thickness_dssd_t1 / 2.);
  G4LogicalVolume *detector2log = new G4LogicalVolume(detector2, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector2");

  // Rodando os detectores em 90 graus
  G4RotationMatrix *rotacaox = new G4RotationMatrix;
  rotacaox->rotateY(0. * CLHEP::deg);

  // Alocando o detector 1
  detector2phys = new G4PVPlacement(rotacaox,
                                    posicao_detector2,
                                    detector2log,
                                    "Detector_3",
                                    logmagnetico,
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

  G4LogicalVolume *logicSensorStripD02 = new G4LogicalVolume(solidSensorStripD02, silicon, "SensorStripD02");

  physiSensorStripD_0_0 = new G4PVReplica("SensorStripD02",    //its name
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

  return detector2phys;
}


G4double *DetectorConstruction::rotacion(G4double dx, G4double dy, G4double dz, G4double angulo)
{
  static G4double valores[3];
  valores[0] = dx * cos(angulo * 3.141562 / 180.) - dy * sin(angulo * 3.141562 / 180.);
  valores[1] = dx * sin(angulo * 3.141562 / 180.) + dy * cos(angulo * 3.141562 / 180.);
  valores[2] = dz;

  return valores;
}

>>>>>>> cb07d2e02789bdab2dbb7f4d731366b95fb38191
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
