// $Id: RootSaver.hh 47 2010-01-19 16:20:46Z adotti $
/**
  * @file   RootSaver.hh
  *
  * @date   17 Dec 2009
  * @author adotti
  *
  * @brief  Save hits and digits in a ROOT TTree
  */

#ifndef ROOTSAVER_HH_
#define ROOTSAVER_HH_

// Local headers
#include "SiHit.hh"

// Default headers
#include <string>
#include <fstream>
#include <iostream>

// ROOT headers
#include "TH2F.h"
#include "TVectorD.h"
#include <TTree.h>

using namespace std;

class TFile;
class RootSaver
{
public:
  //! Default constructor,
  RootSaver();

  //! Default destructor.
  virtual ~RootSaver();

  virtual void CreateTree(const std::string &fileName = "tree", const std::string &treeName = "SiTelescope");

  virtual void CloseTree();

  virtual void AddEvent(const SiHitCollection *const hits,
                        const G4ThreeVector &primaryPos, const G4ThreeVector &primaryMom);

  virtual double Digital(double Eraw);

private:
  // Pointer to the ROOT TTree
  TTree *rootTree;

  // 2D histogram
  TH2F *hist;

  unsigned int runCounter; // Run counter to uniquely identify ROOT file

  Double_t HitsOnDetector = 0;
  Double_t HitsOnTarget = 0;

  // Number of strips of each module
  Int_t nStrips;

  // Number of hits in each detector
  G4int RecoilHit[8] = {0};
  G4int EjectileHit[8] = {0};

  // Signal (energy) in each strip for first module
  Float_t *Signal0;
  Float_t *Signal1;
  Float_t *Signal2;
  Float_t *Signal3;
  Float_t *Signal4;
  Float_t *Signal5;
  Float_t *Signal6;
  Float_t *Signal7;

  // Kinect Energy per strip
  Float_t *Ekin0;
  Float_t *Ekin1;
  Float_t *Ekin2;
  Float_t *Ekin3;
  Float_t *Ekin4;
  Float_t *Ekin5;
  Float_t *Ekin6;
  Float_t *Ekin7;

  // Recoil Theta CM per strip
  Float_t *rThetaCM0;
  Float_t *rThetaCM1;
  Float_t *rThetaCM2;
  Float_t *rThetaCM3;
  Float_t *rThetaCM4;
  Float_t *rThetaCM5;
  Float_t *rThetaCM6;
  Float_t *rThetaCM7;

  // "Truth" position of module
  Float_t Pos_x_det[8] = {0};
  Float_t Pos_y_det[8] = {0};
  Float_t Pos_z_det[8] = {0};

  Float_t Pos_x_target = 0;
  Float_t Pos_y_target = 0;
  Float_t Pos_z_target = 0;

  // Sum of Hits Edep in module
  Float_t E_det[8] = {0};

  // X of the primary at origin
  Float_t TruthPosx;
  Float_t TruthPosy;
  Float_t TruthPosz;

  // Angle in the xz plane (measured from z-axis) of primary at origin
  Float_t TruthAngle_theta;
  Float_t TruthAngle_phi;

  // Momentum
  Float_t Px_dssd;
  Float_t Py_dssd;
  Float_t Pz_dssd;

  //tof
  Float_t T_dssd;
  Float_t T_sili[8] = {0};

  // Etot
  Float_t Etot;

  Int_t StripNumber;
};

#endif
