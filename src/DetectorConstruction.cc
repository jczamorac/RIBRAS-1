// $Id: DetectorConstruction.cc 94 2010-01-26 13:18:30Z adotti $
/**
* @file
* @brief Implements mandatory user class DetectorConstruction.
*/

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

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

#include "SensitiveDetector.hh"
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
#include <fstream>

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

#include "MagneticField.hh"
#include "MagneticField2.hh"

using namespace std;
using namespace CLHEP;

//  Inputs* Inputs = &Inputs::GetInputs();

DetectorConstruction::DetectorConstruction()
{
  //Create a messanger (defines custom UI commands)
  messenger = new DetectorMessenger(this);

  //--------- Material definition ---------
  DefineMaterials();

  //--------- Sizes of the principal geometrical components (solids)  ---------
  ComputeParameters();
  ConstructSetup();
  magneticField = new MagneticField();
  magneticField2 = new MagneticField2();
}

DetectorConstruction::~DetectorConstruction()
{
  delete messenger;
}

void DetectorConstruction::DefineMaterials()
{

  G4double a, z, density;            // z=mean number of protons;
  G4double temperature, pressure;
  G4int ncomponents, natoms;
  G4String name, symbol;             // a=mass of a CLHEP::mole;

  // define Elements
  a = 1.01*CLHEP::g/CLHEP::mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  density     = 1.33e-11*CLHEP::g/CLHEP::cm3;
  pressure    = 1.0913e-10*CLHEP::atmosphere;
  temperature = 200.*CLHEP::kelvin;
  H2 = new G4Material(name="Hydrogen gas", density, ncomponents=1,
                                    kStateGas,temperature,pressure);
  H2->AddElement(elH, natoms=2);


  //define the AlN
  G4Element* elC  = new G4Element(name="Carbono",symbol="C" , z= 6., a=12.01*CLHEP::g/CLHEP::mole);
  G4Element* elHi  = new G4Element(name="Hidrogênio",symbol="H" , z= 1., a=1.01*CLHEP::g/CLHEP::mole);
  density = 0.98*CLHEP::kg/CLHEP::m3;
  temperature = 200.*CLHEP::kelvin;
  pressure =  1.*CLHEP::atmosphere;
  CH4 = new G4Material(name="Metano", density, ncomponents=2, kStateGas, temperature, pressure);
  CH4->AddElement(elC, natoms=1);
  CH4->AddElement(elHi, natoms=4);

  //Get Materials from NIST database
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(0);

  //Define NIST materials
  air     = man->FindOrBuildMaterial("G4_AIR");
  silicon = man->FindOrBuildMaterial("G4_Si");
  vacuum  = man->FindOrBuildMaterial("G4_Galactic");
  steel  = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  tungsten  = man->FindOrBuildMaterial("G4_W");
  //lead  = man->FindOrBuildMaterial("G4_Pb");
  tantalum  = man->FindOrBuildMaterial("G4_Ta");
  sio2  = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
}

void DetectorConstruction::ComputeParameters()
{
  //This function defines the defaults
  //of the geometry construction

  // ** World **
  halfWorldLength = 20.0* CLHEP::m;

  // ** Esta simulacion ya no tiene strips (sub-det dentro de los volumenes) **
  // por eso el tamanho de strip es el mismo del detector
  noOfSensorStrips = 60;
  Lengthy_sili = 100.*CLHEP::mm;
  Thickness_sili = 6.5*CLHEP::mm;
  Lengthx_sili  = 56.*CLHEP::mm;

  Lengthy_dssd = 64.*CLHEP::mm;
  Thickness_dssd = 285.0*CLHEP::um;
  Thickness_dssd_gr = 1.0*CLHEP::mm;
  //Thickness_dssd = 0.005*CLHEP::um;
  Lengthx_dssd  = 64. * CLHEP::mm;

  fwhm_beam = 3.1*CLHEP::mm;
  fwhm_target = 5.5*CLHEP::mm;

  //--------------------ponemos algunos mensajes----------------------------

  //slit_z = GetSlit_z();
  //slit_x = GetSlit_x();

  //pocket_x = GetPocket1_x();

//-------------------------------------------------------------------------

  Lengthy_dssd_t1 = 5.0*CLHEP::cm;
  Thickness_dssd_t1 = 300.*CLHEP::um;
  Lengthx_dssd_t1  = 30.0*CLHEP::cm;

  theta_D_0_0 = /*15.71059*/90.*CLHEP::deg;
  theta_D_0_1 = -15.71059*CLHEP::deg;
  //theta_D_0_2 = -45.71059*CLHEP::deg;

  // ** RingD00
  coorx_D_0_0 = -(450.*CLHEP::mm)*cos(-(90.- theta_D_0_0/CLHEP::deg)*3.141562/180.);
  coory_D_0_0 = 0.*CLHEP::mm;
  coorz_D_0_0 = -(450.*CLHEP::mm)*sin(-(90.- theta_D_0_0/CLHEP::deg)*3.141562/180.);

  // ** RingD01
  coorx_D_0_1 = -(550.*CLHEP::mm)*cos(-(90.- theta_D_0_1/CLHEP::deg)*3.141562/180.);
  coory_D_0_1 = 0.*CLHEP::mm;
  coorz_D_0_1 = -(550.*CLHEP::mm)*sin(-(90.- theta_D_0_1/CLHEP::deg)*3.141562/180.);

  // ** RingD02
/* coorx_D_0_2 = -(560.*CLHEP::mm)*cos(-(90.- theta_D_0_2/CLHEP::deg)*3.141562/180.);
  coory_D_0_2 = 0.*CLHEP::mm;
  coorz_D_0_2 = -(560.*CLHEP::mm)*sin(-(90.- theta_D_0_2/CLHEP::deg)*3.141562/180.);*/

}



void DetectorConstruction::ConstructSetup(void)
{
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  Inputs* Inputs = &Inputs::GetInputs();

  //This function is called by G4 when the detector has to be created
  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------


  //------------------------------
  // World - Aqui o "mundo" é criado. O mundo é Lógico.
  //------------------------------

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*halfWorldLength);
  G4cout << "Computed tolerance = "
        << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/CLHEP::mm
        << " CLHEP::mm" << G4endl;

  G4Box * solidWorld= new G4Box("world",halfWorldLength,halfWorldLength,halfWorldLength);
  logicWorld= new G4LogicalVolume( solidWorld, vacuum, "World", 0, 0, 0);

  //  Must place the World Physical volume unrotated at (0,0,0).
  //
  G4VPhysicalVolume * physiWorld = new G4PVPlacement(0,               // no rotation
                                                    G4ThreeVector(), // at (0,0,0)
                                                    logicWorld,      // its logical volume
                                                    "World",         // its name
                                                    0,               // its mother  volume
                                                    false,           // no boolean operations
                                                    0);              // copy number
  logicWorld -> SetVisAttributes(G4VisAttributes::Invisible);         //Turn the world invisible


  //
  //  Solenoide 1
  //


G4double Solenoid_length = 100.0*CLHEP::cm;
G4double Solenoid_diameter_inner = 30.05*CLHEP::cm;
G4double Solenoid_diameter_outer = 100.0*CLHEP::cm;

G4ThreeVector Pos_Solenoid = G4ThreeVector();

G4VSolid* Sol_Solenoid = new G4Tubs("magneticexterior", Solenoid_diameter_inner/2.0, Solenoid_diameter_outer/2.0, Solenoid_length/2.0,0.,360.*CLHEP::deg);

G4LogicalVolume* Log_Solenoid = new G4LogicalVolume(Sol_Solenoid, G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe"), "Log_Solenoid");

Phy_Solenoid = new G4PVPlacement(0, Pos_Solenoid, Log_Solenoid, "P_Solenoid", logicWorld, false, 0, true);

G4Color green(0.0, 1.0, 0.0),
        red(1.0,0.0,0.0),
        yellow(1.0,1.0,0.0),
        orange(1.0, 1.0, 0.5),
        blue(0.0, 0.0, 1.0);

Log_Solenoid -> SetVisAttributes(new G4VisAttributes(green));



/////////////////////////////
/// Magnetic field region ///
/////////////////////////////

G4double Mag_diameter = 30.0*cm;
G4double Mag_length = 68.0*cm; //coil length


G4VSolid* Sol_Magnet = new G4Tubs("S_Magnet",0.,Mag_diameter/2.0,Mag_length/2.0,0.,360.*deg);

G4LogicalVolume* Log_Magnet  = new G4LogicalVolume(Sol_Magnet,G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),"Log_Magnet",0,0,0);

G4VPhysicalVolume* Phy_Magnet = new G4PVPlacement(0,Pos_Solenoid, Log_Magnet, "P_Magnet", logicWorld, false, 0, true);

// set "user limits" for drawing smooth curve
G4UserLimits* userLimits = new G4UserLimits(1.e-5*mm);
Log_Magnet->SetUserLimits(userLimits);







    //
    // Solenoide 2
    //

G4ThreeVector position = G4ThreeVector(0., 0., 301.*CLHEP::cm);
G4double raiointerno = 30.05*CLHEP::cm;
G4double raioexterno = 100.0*CLHEP::cm;
G4double comprimento = 100.0*CLHEP::cm;

G4VSolid* solenoide = new G4Tubs("solenoide", raiointerno/2.0, raioexterno/2.0, comprimento/2.0, 0. ,360.*CLHEP::deg);

G4LogicalVolume* logsolenoide = new G4LogicalVolume(solenoide, G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe"), "logsolenoide");

solenoidef = new G4PVPlacement(0, position, logsolenoide, "solenoidef", logicWorld, false, 0, true);

logsolenoide -> SetVisAttributes(new G4VisAttributes(green));

    //
    // Região do campo magnético do solenóide 2
    //

G4double diametromag = 30.0*cm;
G4double comprimentomag = 68.0*cm; //coil length
//G4int zero = 0;

G4VSolid* magsolenoide = new G4Tubs("magsolenoide",0.,diametromag/2.0,comprimentomag/2.0,0.,360.*deg);

G4LogicalVolume* logmagnetico  = new G4LogicalVolume(magsolenoide, G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"), "logmagnetico");

G4VPhysicalVolume* magneticof = new G4PVPlacement(0,
              position,
              logmagnetico,
              "magneticof",
              logicWorld,
              false,
              0,
              true);

// set "user limits" for drawing smooth curve
G4UserLimits* limiteuser = new G4UserLimits(1.e-5*mm);
logmagnetico->SetUserLimits(limiteuser);
//logicWorld->SetUserLimits(limiteuser);

//
// Alvo
//

G4ThreeVector Target_pos = Inputs->target_pos;
G4Material* TargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Inputs->g4_material_name);
Inputs->TargetMaterial = TargetMaterial;

G4Box* alvo = new G4Box("alvo",7.5*CLHEP::cm,7.5*CLHEP::cm,Inputs->width *CLHEP::mm);
//G4VSolid* alvo = new G4Tubs("alvo", 0., raiointerno/2.0, 0.5*CLHEP::mm, 0., 360.*deg);

G4LogicalVolume* alvolog = new G4LogicalVolume(alvo, TargetMaterial, "alvolog");

G4VPhysicalVolume* alvophys = new G4PVPlacement(0, Target_pos, alvolog, "alvof", logmagnetico, false, 0, true);

alvolog -> SetVisAttributes(new G4VisAttributes(red));
  
////////////////////
// Magnetic field //
////////////////////

G4bool fieldIsInitialized = false;
if(!fieldIsInitialized)
{

  //magneticField->SetCurrent(21.336); //electric current of solenoid 1 (Amp)

  // G4cout<<"Creating Magnetic Field 1: I = "<<magneticField->GetCurrent()<<" Amp"<<G4endl;
  // G4cout<<"Creating Magnetic Field 2: I = "<<magneticField2->GetCurrent()<<" Amp"<<G4endl;


  G4MagIntegratorStepper* fStepper;


  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetDetectorField(magneticField);
  fieldMgr->CreateChordFinder(magneticField);


  G4double minEps= 1.0e-10*cm;  //   Minimum & value for smallest steps
  G4double maxEps= 1.0e-8*cm;  //   Maximum & value for largest steps

  fieldMgr->SetMinimumEpsilonStep( minEps );
  fieldMgr->SetMaximumEpsilonStep( maxEps );


  fieldMgr->GetChordFinder()->SetDeltaChord( 0.001*mm);

  G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(magneticField);
  //G4EqMagElectricField* fEquation = new G4EqMagElectricField(magneticField);

  // Note that for magnetic field that do not vary with time,
  //  the default number of variables suffices.
  // or ..
  fStepper = new G4HelixExplicitEuler( fEquation );   // mais indicado para campos intensos e não suaves
  //fStepper = new G4HelixSimpleRunge( fEquation );
  //fStepper = new G4HelixImplicitEuler( fEquation );
  //fStepper = new G4CashKarpRKF45( fEquation );
  // fStepper = new G4ClassicalRK4( fEquation );   // integrador mais usado

  fieldMgr->SetDeltaIntersection(0.1*mm);
  fieldMgr->SetAccuraciesWithDeltaOneStep(0.01*mm );
  fieldMgr->SetDeltaOneStep( 0.01*mm );  // 0.5 micrometer

  /* G4ChordFinder( magneticField,
                1.0e-2 * mm,
              0 );
*/
  fieldIsInitialized = true;
}

  ConstructDetectors();

  return physiWorld;
}

G4VPhysicalVolume* DetectorConstruction::ConstructDetectors()
{

  ConstructRing_D_0_0();
  ConstructRing_D_0_1();
  //ConstructRing_D_0_2();




static SensitiveDetector* sensitive = 0;
if (!sensitive)
{

  sensitive = new SensitiveDetector("/myDet/SiStripSD");
  //We register now the SD with the manager
  G4SDManager::GetSDMpointer()->AddNewDetector(sensitive);

}


//******* dos lineas para anexar un detector diferente en una replica de detectores

    const G4LogicalVolume* log = detector1phys->GetLogicalVolume();
    log->GetDaughter(0)->GetLogicalVolume()->SetSensitiveDetector(sensitive);

    const G4LogicalVolume* ring_sen_D01 = physiSensorRing_D_0_1->GetLogicalVolume();
    ring_sen_D01->GetDaughter(0)->GetLogicalVolume()->SetSensitiveDetector(sensitive);

  /*const G4LogicalVolume* ring_sen_D02 = physiSensorRing_D_0_2->GetLogicalVolume();
  ring_sen_D02->GetDaughter(0)->GetLogicalVolume()->SetSensitiveDetector(sensitive);*/


  return detector1phys;
}


G4VPhysicalVolume* DetectorConstruction::ConstructRing_D_0_0()
{

//
// Detectores
//

// Coordenadas detectores

G4ThreeVector posicao_detector1 = G4ThreeVector(-6.5*CLHEP::cm , 0., 276.*CLHEP::cm);
G4ThreeVector posicao_detector1teste = G4ThreeVector( -65.*CLHEP::cm, 0., 276.*CLHEP::cm);

// Criando os detectores físicos

G4Box* detector1 = new G4Box("detector1",/* 300.*CLHEP::um, 2.5*CLHEP::cm, 15.*CLHEP::cm*/Lengthx_dssd_t1/2.,Lengthy_dssd_t1/2.,Thickness_dssd_t1/2.);
G4LogicalVolume* detector1log = new G4LogicalVolume(detector1, G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"), "logdetector1");
//G4VPhysicalVolume* detector1phys = new G4PVPlacement(0, posicao_detector1teste, detector1log, "detectorf1", logicWorld, false, 0, true);

/*  G4Box * solidSensor_D_0_0 = new G4Box("SensorD_0_0",
                                  Lengthx_dssd_t1/2.,Lengthy_dssd_t1/2.,Thickness_dssd_t1/2.);

  G4LogicalVolume * logicSensorPlane = new G4LogicalVolume(solidSensor_D_0_0, // its solid
                                                          silicon,     //its material
                                                          "SensorL_D_0_0");   //its name*/


  //for(int jd00=0;jd00<1;jd00++){
int jd00 = 0;
G4RotationMatrix * rot_D_0_0 = new G4RotationMatrix;
//phi_D_0_0 = (-(360./12.)*jd00);
//rot_D_0_0->rotateZ(phi_D_0_0*CLHEP::deg);
    rot_D_0_0->rotateY(theta_D_0_0);
G4RotationMatrix * rotacaox = new G4RotationMatrix;
rotacaox -> rotateY(90.0*CLHEP::deg);


//if(jd00 % 2 == 0){ ring_D_0_0 = rotacion(coorx_D_0_0, coory_D_0_0, coorz_D_0_0, -phi_D_0_0);}
//else{ring_D_0_0 = rotacion(coorx_D_0_0_e, coory_D_0_0_e, coorz_D_0_0_e, -phi_D_0_0);}

//vector_D_0_0 = G4ThreeVector(*(ring_D_0_0), *(ring_D_0_0+1), *(ring_D_0_0+2));
//vector_D_0_0 = G4ThreeVector(coorx_D_0_0, coory_D_0_0, coorz_D_0_0);

std::ostringstream nombre_d_0_0;
  nombre_d_0_0 << "D_0_0_" << jd00 ;

//G4cout<<phi_D_0_0*deg/deg<<"  "<<vector_D_0_0<<G4endl;
//G4cout<<jd00<<"  "<<theta_D_0_0/CLHEP::deg<<"  "<<-phi_D_0_0*CLHEP::deg/CLHEP::deg<<"  "<<sqrt(coorx_D_0_0*coorx_D_0_0 + coorz_D_0_0*coorz_D_0_0)<<G4endl;
  detector1phys =
      new G4PVPlacement(rotacaox,
                        posicao_detector1,
                        detector1log,
                        nombre_d_0_0.str().data(),
                        logicWorld,
                        false,
                        jd00);
//}

  //
  // Strips
  //

  G4double halfSensorStripSizeX = Lengthx_dssd_t1/(2.0*noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1/2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1/2.;

  G4Box * solidSensorStripD00 =
    new G4Box("SensorStripD00",
              halfSensorStripSizeX,halfSensorStripSizeY,halfSensorStripSizeZ);

  G4LogicalVolume * logicSensorStripD00 =
    new G4LogicalVolume(solidSensorStripD00,silicon,"SensorStripD00");

  physiSensorStripD_0_0 =
    new G4PVReplica("SensorStripD00",       //its name
                    logicSensorStripD00,    //its logical volume
                    detector1log,           //its mother
                    kXAxis,                 //axis of replication
                    noOfSensorStrips,       //number of replica
                    2.0*halfSensorStripSizeX);
  //                Lengthx_sili);            //witdth of replica



  // set "user limits" StepSize
G4UserLimits* userLimits = new G4UserLimits(0.1*CLHEP::mm);
logicSensorStripD00->SetUserLimits(userLimits);


G4Color yellow(1.0,1.0,0.0),red(1.0,0.0,0.0);/*, orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);*/
  detector1log -> SetVisAttributes(new G4VisAttributes(red));
  logicSensorStripD00 -> SetVisAttributes(new G4VisAttributes(red));

  return detector1phys;

}

G4VPhysicalVolume* DetectorConstruction::ConstructRing_D_0_1()
{


  G4ThreeVector posicao_detector2 = G4ThreeVector(-6.5*CLHEP::cm,0.*CLHEP::cm , 326.*CLHEP::cm);


  G4Box * solidSensor_D_0_1 = new G4Box("SensorD_0_1",
                                  Lengthx_dssd_t1/2.,Lengthy_dssd_t1/2.,Thickness_dssd_t1/2.);

  G4LogicalVolume * logicSensorPlane = new G4LogicalVolume(solidSensor_D_0_1, // its solid
                                                          G4NistManager::Instance()->FindOrBuildMaterial("G4_Si"),     //its material
                                                          "SensorL_D_0_1");   //its name


G4RotationMatrix * rotacaox = new G4RotationMatrix;
rotacaox -> rotateY(90.0*CLHEP::deg);



  //for(int jd00=0;jd00<1;jd00++){
int jd01 = 1;
G4RotationMatrix * rot_D_0_1 = new G4RotationMatrix;
//phi_D_0_0 = (-(360./12.)*jd00);
//rot_D_0_0->rotateZ(phi_D_0_0*CLHEP::deg);
    rot_D_0_1->rotateY(theta_D_0_1);


//if(jd00 % 2 == 0){ ring_D_0_0 = rotacion(coorx_D_0_0, coory_D_0_0, coorz_D_0_0, -phi_D_0_0);}
//else{ring_D_0_0 = rotacion(coorx_D_0_0_e, coory_D_0_0_e, coorz_D_0_0_e, -phi_D_0_0);}

//vector_D_0_0 = G4ThreeVector(*(ring_D_0_0), *(ring_D_0_0+1), *(ring_D_0_0+2));
vector_D_0_1 = G4ThreeVector(coorx_D_0_1, coory_D_0_1, coorz_D_0_1);


std::ostringstream nombre_d_0_1;
  nombre_d_0_1 << "D_0_1_" << jd01 ;

//G4cout<<phi_D_0_0*deg/deg<<"  "<<vector_D_0_0<<G4endl;
//G4cout<<jd00<<"  "<<theta_D_0_0/CLHEP::deg<<"  "<<-phi_D_0_0*CLHEP::deg/CLHEP::deg<<"  "<<sqrt(coorx_D_0_0*coorx_D_0_0 + coorz_D_0_0*coorz_D_0_0)<<G4endl;
    physiSensorRing_D_0_1 =
      new G4PVPlacement(rotacaox,
                        posicao_detector2,
                        logicSensorPlane,
                        nombre_d_0_1.str().data(),
                        logicWorld,
                        false,
                        jd01);
//}

//
// Strips
//

  G4double halfSensorStripSizeX = Lengthx_dssd_t1/(2.0*noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1/2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1/2.;


  G4Box * solidSensorStripD01 =
    new G4Box("SensorStripD01",
              halfSensorStripSizeX,halfSensorStripSizeY,halfSensorStripSizeZ);

  G4LogicalVolume * logicSensorStripD01 =
    new G4LogicalVolume(solidSensorStripD01,silicon,"SensorStripD01");

  physiSensorStripD_0_1 =
    new G4PVReplica("SensorStripD01",           //its name
                    logicSensorStripD01,        //its logical volume
                    logicSensorPlane,           //its mother
                    kXAxis,                     //axis of replication
                    noOfSensorStrips,           //number of replica
                    2.0*halfSensorStripSizeX);
  //                Lengthx_sili);            //witdth of replica

  // set "user limits" StepSize
G4UserLimits* userLimits = new G4UserLimits(0.1*CLHEP::mm);
logicSensorStripD01->SetUserLimits(userLimits);


  G4Color red(1.0,0.0,0.0),yellow(1.0,1.0,0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  logicSensorPlane -> SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD01 -> SetVisAttributes(new G4VisAttributes(red));

  return physiSensorRing_D_0_1;
}


/* G4VPhysicalVolume* DetectorConstruction::ConstructRing_D_0_2()
{


  G4Box * solidSensor_D_0_2 = new G4Box("SensorD_0_2",
                                  Lengthx_dssd_t1/2.,Lengthy_dssd_t1/2.,Thickness_dssd_t1/2.);

  G4LogicalVolume * logicSensorPlane = new G4LogicalVolume(solidSensor_D_0_2, // its solid
                                                          silicon,     //its material
                                                          "SensorL_D_0_2");   //its name



  //for(int jd00=0;jd00<1;jd00++){
int jd02 = 2;
G4RotationMatrix * rot_D_0_2 = new G4RotationMatrix;
//phi_D_0_0 = (-(360./12.)*jd00);
//rot_D_0_0->rotateZ(phi_D_0_0*CLHEP::deg);
    rot_D_0_2->rotateY(theta_D_0_2);


//if(jd00 % 2 == 0){ ring_D_0_0 = rotacion(coorx_D_0_0, coory_D_0_0, coorz_D_0_0, -phi_D_0_0);}
//else{ring_D_0_0 = rotacion(coorx_D_0_0_e, coory_D_0_0_e, coorz_D_0_0_e, -phi_D_0_0);}

//vector_D_0_0 = G4ThreeVector(*(ring_D_0_0), *(ring_D_0_0+1), *(ring_D_0_0+2));
vector_D_0_2 = G4ThreeVector(coorx_D_0_2, coory_D_0_2, coorz_D_0_2);


std::ostringstream nombre_d_0_2;
  nombre_d_0_2 << "D_0_2_" << jd02 ;

//G4cout<<phi_D_0_0*deg/deg<<"  "<<vector_D_0_0<<G4endl;
//G4cout<<jd00<<"  "<<theta_D_0_0/CLHEP::deg<<"  "<<-phi_D_0_0*CLHEP::deg/CLHEP::deg<<"  "<<sqrt(coorx_D_0_0*coorx_D_0_0 + coorz_D_0_0*coorz_D_0_0)<<G4endl;
    physiSensorRing_D_0_2 =
      new G4PVPlacement(rot_D_0_2,
                        vector_D_0_2,
                        logicSensorPlane,
                        nombre_d_0_2.str().data(),
                        logicWorld,
                        false,
                        jd02);
//}

//
  // Strips

  G4double halfSensorStripSizeX = Lengthx_dssd_t1/(2.0*noOfSensorStrips);
  G4double halfSensorStripSizeY = Lengthy_dssd_t1/2.;
  G4double halfSensorStripSizeZ = Thickness_dssd_t1/2.;


  G4Box * solidSensorStripD02 =
    new G4Box("SensorStripD02",
              halfSensorStripSizeX,halfSensorStripSizeY,halfSensorStripSizeZ);

  G4LogicalVolume * logicSensorStripD02 =
    new G4LogicalVolume(solidSensorStripD02,silicon,"SensorStripD02");

  physiSensorStripD_0_2 =
    new G4PVReplica("SensorStripD02",           //its name
                    logicSensorStripD02,           //its logical volume
                    logicSensorPlane,           //its mother
                    kXAxis,                     //axis of replication
                    noOfSensorStrips,           //number of replica
                    2.0*halfSensorStripSizeX);
  //                Lengthx_sili);            //witdth of replica

  // set "user limits" StepSize
G4UserLimits* userLimits = new G4UserLimits(0.1*CLHEP::mm);
logicSensorStripD02->SetUserLimits(userLimits);


  G4Color red(1.0,0.0,0.0),yellow(1.0,1.0,0.0), orange(1.0, 1.0, 0.5), blue(0.0, 0.0, 1.0), green(0.0, 1.0, 0.0);
  logicSensorPlane -> SetVisAttributes(new G4VisAttributes(yellow));
  logicSensorStripD02 -> SetVisAttributes(new G4VisAttributes(red));

  return physiSensorRing_D_0_20;
} */

G4double* DetectorConstruction::rotacion(G4double dx , G4double dy, G4double dz, G4double angulo)
{
static G4double valores[3];
    valores[0] = dx*cos(angulo*3.141562/180.) - dy*sin(angulo*3.141562/180.);
    valores[1] = dx*sin(angulo*3.141562/180.) + dy*cos(angulo*3.141562/180.);
    valores[2] = dz;

return valores;
}



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
