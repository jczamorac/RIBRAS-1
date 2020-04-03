 // $Id: DetectorMessenger.hh 33 2010-01-14 17:08:18Z adotti $

 #ifndef DetectorMessenger_h
 #define DetectorMessenger_h 1

 /**
  * @file
  * @brief defines class DetectorMessenger
  */

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 #include "globals.hh"
 #include "G4UImessenger.hh"
 #include "G4UIcmdWithAnInteger.hh"


 class DetectorConstruction;
 class G4UIdirectory;
 class G4UIcmdWithAString;
 class G4UIcmdWithADoubleAndUnit;
 class G4UIcmdWith3VectorAndUnit;
 class G4UIcmdWithoutParameter;
 class G4UIcmdWithABool;

 using namespace std;

 class DetectorMessenger: public G4UImessenger
 {
 public:
   //! Constructor
   DetectorMessenger(DetectorConstruction* );
   //! Destructor
   ~DetectorMessenger();

   //! handle user commands
   void SetNewValue(G4UIcommand* command, G4String newValue);

 private:

   DetectorConstruction*      detector;

   G4UIdirectory*             detDir;
   G4UIdirectory*             secondSensorDir;
   G4UIdirectory*             Dir;
   G4UIdirectory*             targetDir;
   G4UIdirectory*             recoilDir;
   G4UIdirectory*             ejectileDir;
   G4UIdirectory*             primaryDir;

   G4UIcmdWithADoubleAndUnit* xShiftCmd;
   G4UIcmdWithADoubleAndUnit* yShiftCmd;
   G4UIcmdWithADoubleAndUnit* thetaCmd;

   G4UIcmdWithADoubleAndUnit* z_slitCmd;
   G4UIcmdWithADoubleAndUnit* x_slitCmd;
   G4UIcmdWithADoubleAndUnit* xShiftPocketCmd;

   G4UIcmdWithoutParameter*   updateCmd;

   G4UIcmdWithABool*                      setDUTsetupCmd;

   G4UIcmdWithAString*        target_material;
 	 G4UIcmdWithADoubleAndUnit* target_width_cmd;
 	 G4UIcmdWithADoubleAndUnit* target_height_cmd;
 	 G4UIcmdWithADoubleAndUnit* target_radius_cmd;
 	 G4UIcmdWithADoubleAndUnit* target_arial_density_cmd;
 	 G4UIcmdWithABool*          target_recfoil;
 	 G4UIcmdWithABool*          target_cylfoil;
 	 G4UIcmdWithABool*          target_lh2;
   G4UIcmdWithABool*          magneticfieldon;
 	 G4UIcmdWithADoubleAndUnit* target_mass;
 	 G4UIcmdWithADoubleAndUnit* recoil_mass;
 	 G4UIcmdWithADoubleAndUnit* ejectile_mass;
 	 G4UIcmdWithADoubleAndUnit* recoil_Ex;
 	 G4UIcmdWithADoubleAndUnit* ejectile_Ex;
 	 G4UIcmdWithAnInteger*      target_A;
 	 G4UIcmdWithAnInteger*      recoil_A;
 	 G4UIcmdWithAnInteger*      ejectile_A;
 	 G4UIcmdWithAnInteger*      target_Z;
 	 G4UIcmdWithAnInteger*      recoil_Z;
 	 G4UIcmdWithAnInteger*      ejectile_Z;
   G4UIcmdWithADoubleAndUnit* primary_energy;
   G4UIcmdWithAnInteger*      primary_A;
 	 G4UIcmdWithAnInteger*      primary_Z;
   G4UIcmdWith3VectorAndUnit* primary_pos;
   G4UIcmdWith3VectorAndUnit* target_pos;
 };

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 #endif
