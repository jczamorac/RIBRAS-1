//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 13 15:34:09 2020 by ROOT version 6.21/01
// from TChain SiTelescope/
//////////////////////////////////////////////////////////

#ifndef Analysis_h
#define Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include "TF1.h"
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TH2Poly.h>
#include "TEnv.h"

#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <Math/Vector3D.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <numeric> // std::iota

// Header file for the classes stored in the TTree if any.

class Analysis
{
public:
   TTree *fChain;  //!pointer to the analyzed TTree or TChain
   Int_t fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t E0[60];
   Float_t E1[60];
   Float_t E2[60];
   Float_t E3[60];
   Float_t E4[60];
   Float_t E5[60];
   Float_t E6[60];
   Float_t E7[60];
   Float_t EKin0[60];
   Float_t EKin1[60];
   Float_t EKin2[60];
   Float_t EKin3[60];
   Float_t EKin4[60];
   Float_t EKin5[60];
   Float_t EKin6[60];
   Float_t EKin7[60];
   Float_t rThetaCM0[60];
   Float_t rThetaCM1[60];
   Float_t rThetaCM2[60];
   Float_t rThetaCM3[60];
   Float_t rThetaCM4[60];
   Float_t rThetaCM5[60];
   Float_t rThetaCM6[60];
   Float_t rThetaCM7[60];
   Float_t pos_x_det0;
   Float_t pos_y_det0;
   Float_t pos_z_det0;
   Float_t pos_x_det1;
   Float_t pos_y_det1;
   Float_t pos_z_det1;
   Float_t pos_x_det2;
   Float_t pos_y_det2;
   Float_t pos_z_det2;
   Float_t pos_x_det3;
   Float_t pos_y_det3;
   Float_t pos_z_det3;
   Float_t pos_x_det4;
   Float_t pos_y_det4;
   Float_t pos_z_det4;
   Float_t pos_x_det5;
   Float_t pos_y_det5;
   Float_t pos_z_det5;
   Float_t pos_x_det6;
   Float_t pos_y_det6;
   Float_t pos_z_det6;
   Float_t pos_x_det7;
   Float_t pos_y_det7;
   Float_t pos_z_det7;
   Float_t pos_x_target;
   Float_t pos_y_target;
   Float_t pos_z_target;
   Float_t ener_det0;
   Float_t ener_det1;
   Float_t ener_det2;
   Float_t ener_det3;
   Float_t ener_det4;
   Float_t ener_det5;
   Float_t ener_det6;
   Float_t ener_det7;
   Float_t truthPosx;
   Float_t truthPosy;
   Float_t truthPosz;
   Float_t truthAngle_theta;
   Float_t truthAngle_phi;
   Float_t px_dssd;
   Float_t py_dssd;
   Float_t pz_dssd;
   Float_t t_dssd;
   Float_t t_sili0;
   Float_t t_sili1;
   Float_t t_sili2;
   Float_t t_sili3;
   Float_t t_sili4;
   Float_t t_sili5;
   Float_t t_sili6;
   Float_t t_sili7;
   Int_t Strip_Number0;
   Int_t Strip_Number1;
   Int_t Strip_Number2;
   Int_t Strip_Number3;
   Int_t Strip_Number4;
   Int_t Strip_Number5;
   Int_t Strip_Number6;
   Int_t Strip_Number7;

   Float_t t_dssd2;
   Float_t ETot;

   // List of branches
   TBranch *b_E0;               //!
   TBranch *b_E1;               //!
   TBranch *b_E2;               //!
   TBranch *b_E3;               //!
   TBranch *b_E4;               //!
   TBranch *b_E5;               //!
   TBranch *b_E6;               //!
   TBranch *b_E7;               //!
   TBranch *b_EKin0;            //!
   TBranch *b_EKin1;            //!
   TBranch *b_EKin2;            //!
   TBranch *b_EKin3;            //!
   TBranch *b_EKin4;            //!
   TBranch *b_EKin5;            //!
   TBranch *b_EKin6;            //!
   TBranch *b_EKin7;            //!
   TBranch *b_rThetaCM0;        //!
   TBranch *b_rThetaCM1;        //!
   TBranch *b_rThetaCM2;        //!
   TBranch *b_rThetaCM3;        //!
   TBranch *b_rThetaCM4;        //!
   TBranch *b_rThetaCM5;        //!
   TBranch *b_rThetaCM6;        //!
   TBranch *b_rThetaCM7;        //!
   TBranch *b_pos_x_det0;       //!
   TBranch *b_pos_y_det0;       //!
   TBranch *b_pos_z_det0;       //!
   TBranch *b_pos_x_det1;       //!
   TBranch *b_pos_y_det1;       //!
   TBranch *b_pos_z_det1;       //!
   TBranch *b_pos_x_det2;       //!
   TBranch *b_pos_y_det2;       //!
   TBranch *b_pos_z_det2;       //!
   TBranch *b_pos_x_det3;       //!
   TBranch *b_pos_y_det3;       //!
   TBranch *b_pos_z_det3;       //!
   TBranch *b_pos_x_det4;       //!
   TBranch *b_pos_y_det4;       //!
   TBranch *b_pos_z_det4;       //!
   TBranch *b_pos_x_det5;       //!
   TBranch *b_pos_y_det5;       //!
   TBranch *b_pos_z_det5;       //!
   TBranch *b_pos_x_det6;       //!
   TBranch *b_pos_y_det6;       //!
   TBranch *b_pos_z_det6;       //!
   TBranch *b_pos_x_det7;       //!
   TBranch *b_pos_y_det7;       //!
   TBranch *b_pos_z_det7;       //!
   TBranch *b_pos_x_target;     //!
   TBranch *b_pos_y_target;     //!
   TBranch *b_pos_z_target;     //!
   TBranch *b_ener_det0;        //!
   TBranch *b_ener_det1;        //!
   TBranch *b_ener_det2;        //!
   TBranch *b_ener_det3;        //!
   TBranch *b_ener_det4;        //!
   TBranch *b_ener_det5;        //!
   TBranch *b_ener_det6;        //!
   TBranch *b_ener_det7;        //!
   TBranch *b_truthPosx;        //!
   TBranch *b_truthPosy;        //!
   TBranch *b_truthPosz;        //!
   TBranch *b_truthAngle_theta; //!
   TBranch *b_truthAngle_phi;   //!
   TBranch *b_px_dssd;          //!
   TBranch *b_py_dssd;          //!
   TBranch *b_pz_dssd;          //!
   TBranch *b_t_dssd;           //!
   TBranch *b_t_sili0;          //!
   TBranch *b_t_sili1;          //!
   TBranch *b_t_sili2;          //!
   TBranch *b_t_sili3;          //!
   TBranch *b_t_sili4;          //!
   TBranch *b_t_sili5;          //!
   TBranch *b_t_sili6;          //!
   TBranch *b_t_sili7;          //!
   TBranch *b_t_dssd2;          //!
   TBranch *b_ETot;             //!
   TBranch *b_Strip_Number0;    //!
   TBranch *b_Strip_Number1;    //!
   TBranch *b_Strip_Number2;    //!
   TBranch *b_Strip_Number3;    //!
   TBranch *b_Strip_Number4;    //!
   TBranch *b_Strip_Number5;    //!
   TBranch *b_Strip_Number6;    //!
   TBranch *b_Strip_Number7;    //!

   // 2D Histograms
   TH2F *EvS0 = new TH2F("Energy vs Strip 0", "Energy vs Strip", 60, 0, 60, 1000, 0, 20);
   TH2F *EvS1 = new TH2F("Energy vs Strip 1", "Energy vs Strip", 60, 0, 60, 1000, 0, 20);
   TH2F *EvS2 = new TH2F("Energy vs Strip 2", "Energy vs Strip", 60, 0, 60, 1000, 0, 20);
   TH2F *EvS3 = new TH2F("Energy vs Strip 3", "Energy vs Strip", 60, 0, 60, 1000, 0, 20);
   TH2F *EvS4 = new TH2F("Energy vs Strip 4", "Energy vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EvS5 = new TH2F("Energy vs Strip 5", "Energy vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EvS6 = new TH2F("Energy vs Strip 6", "Energy vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EvS7 = new TH2F("Energy vs Strip 7", "Energy vs Strip", 60, 0, 60, 1000, 0, 45);

   TH2F *EvThetaCM0 = new TH2F("Energy vs ThetaCM0", "Energy vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EvThetaCM1 = new TH2F("Energy vs ThetaCM1", "Energy vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EvThetaCM2 = new TH2F("Energy vs ThetaCM2", "Energy vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EvThetaCM3 = new TH2F("Energy vs ThetaCM3", "Energy vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EvThetaCM4 = new TH2F("Energy vs ThetaCM4", "Energy vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EvThetaCM5 = new TH2F("Energy vs ThetaCM5", "Energy vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EvThetaCM6 = new TH2F("Energy vs ThetaCM6", "Energy vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EvThetaCM7 = new TH2F("Energy vs ThetaCM7", "Energy vs ThetaCM", 180, 0, 180, 1000, 0, 40);

   TH2F *EkinvTheta0 = new TH2F("Ekin vs ThetaCM 0", "Ekin vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EkinvTheta1 = new TH2F("Ekin vs ThetaCM 1", "Ekin vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EkinvTheta2 = new TH2F("Ekin vs ThetaCM 2", "Ekin vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EkinvTheta3 = new TH2F("Ekin vs ThetaCM 3", "Ekin vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EkinvTheta4 = new TH2F("Ekin vs ThetaCM 4", "Ekin vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EkinvTheta5 = new TH2F("Ekin vs ThetaCM 5", "Ekin vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EkinvTheta6 = new TH2F("Ekin vs ThetaCM 6", "Ekin vs ThetaCM", 180, 0, 180, 1000, 0, 40);
   TH2F *EkinvTheta7 = new TH2F("Ekin vs ThetaCM 7", "Ekin vs ThetaCM", 180, 0, 180, 1000, 0, 40);

   TH2F *EKinvS0 = new TH2F("Ekin vs Strip 0", "Ekin vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EKinvS1 = new TH2F("Ekin vs Strip 1", "Ekin vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EKinvS2 = new TH2F("Ekin vs Strip 2", "Ekin vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EKinvS3 = new TH2F("Ekin vs Strip 3", "Ekin vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EKinvS4 = new TH2F("Ekin vs Strip 4", "Ekin vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EKinvS5 = new TH2F("Ekin vs Strip 5", "Ekin vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EKinvS6 = new TH2F("Ekin vs Strip 6", "Ekin vs Strip", 60, 0, 60, 1000, 0, 45);
   TH2F *EKinvS7 = new TH2F("Ekin vs Strip 7", "Ekin vs Strip", 60, 0, 60, 1000, 0, 45);

   TH2F *TimevS0 = new TH2F("Tempo de voo vs Strip 0", "Tempo de voo vs Strip", 60, 0, 60, 1000, 0, 400);
   TH2F *TimevS1 = new TH2F("Tempo de voo vs Strip 1", "Tempo de voo vs Strip", 60, 0, 60, 1000, 0, 400);
   TH2F *TimevS2 = new TH2F("Tempo de voo vs Strip 2", "Tempo de voo vs Strip", 60, 0, 60, 1000, 0, 400);
   TH2F *TimevS3 = new TH2F("Tempo de voo vs Strip 3", "Tempo de voo vs Strip", 60, 0, 60, 1000, 0, 400);
   TH2F *TimevS4 = new TH2F("Tempo de voo vs Strip 4", "Tempo de voo vs Strip", 60, 0, 60, 1000, 0, 400);
   TH2F *TimevS5 = new TH2F("Tempo de voo vs Strip 5", "Tempo de voo vs Strip", 60, 0, 60, 1000, 0, 400);
   TH2F *TimevS6 = new TH2F("Tempo de voo vs Strip 6", "Tempo de voo vs Strip", 60, 0, 60, 1000, 0, 400);
   TH2F *TimevS7 = new TH2F("Tempo de voo vs Strip 7", "Tempo de voo vs Strip", 60, 0, 60, 1000, 0, 400);

   TH2F *ThetavS0 = new TH2F("Strip vs Theta 0", "Strip vs Theta", 60, 0, 60, 180, 0, 180);
   TH2F *ThetavS1 = new TH2F("Strip vs Theta 1", "Strip vs Theta", 60, 0, 60, 180, 0, 180);
   TH2F *ThetavS2 = new TH2F("Strip vs Theta 2", "Strip vs Theta", 60, 0, 60, 180, 0, 180);
   TH2F *ThetavS3 = new TH2F("Strip vs Theta 3", "Strip vs Theta", 60, 0, 60, 180, 0, 180);
   TH2F *ThetavS4 = new TH2F("Strip vs Theta 4", "Strip vs Theta", 60, 0, 60, 180, 0, 180);
   TH2F *ThetavS5 = new TH2F("Strip vs Theta 5", "Strip vs Theta", 60, 0, 60, 180, 0, 180);
   TH2F *ThetavS6 = new TH2F("Strip vs Theta 6", "Strip vs Theta", 60, 0, 60, 180, 0, 180);
   TH2F *ThetavS7 = new TH2F("Strip vs Theta 7", "Strip vs Theta", 60, 0, 60, 180, 0, 180);

   TH2F *TvEkin0 = new TH2F("Tempo de voo vs Energia cinetica 0", "Tempo de voo vs Energia cinetica", 60, 0, 45, 1000, 0, 400);
   TH2F *TvEkin1 = new TH2F("Tempo de voo vs Energia cinetica 1", "Tempo de voo vs Energia cinetica", 60, 0, 45, 1000, 0, 400);
   TH2F *TvEkin2 = new TH2F("Tempo de voo vs Energia cinetica 2", "Tempo de voo vs Energia cinetica", 60, 0, 45, 1000, 0, 400);
   TH2F *TvEkin3 = new TH2F("Tempo de voo vs Energia cinetica 3", "Tempo de voo vs Energia cinetica", 60, 0, 45, 1000, 0, 400);
   TH2F *TvEkin4 = new TH2F("Tempo de voo vs Energia cinetica 4", "Tempo de voo vs Energia cinetica", 60, 0, 45, 1000, 0, 400);
   TH2F *TvEkin5 = new TH2F("Tempo de voo vs Energia cinetica 5", "Tempo de voo vs Energia cinetica", 60, 0, 45, 1000, 0, 400);
   TH2F *TvEkin6 = new TH2F("Tempo de voo vs Energia cinetica 6", "Tempo de voo vs Energia cinetica", 60, 0, 45, 1000, 0, 400);
   TH2F *TvEkin7 = new TH2F("Tempo de voo vs Energia cinetica 7", "Tempo de voo vs Energia cinetica", 60, 0, 45, 1000, 0, 400);

   TH2F *TvE0 = new TH2F("Tempo de voo vs Energia depositada 0", "Tempo de voo vs Energia Depositada", 60, 0, 45, 1000, 0, 400);
   TH2F *TvE1 = new TH2F("Tempo de voo vs Energia depositada 1", "Tempo de voo vs Energia Depositada", 60, 0, 45, 1000, 0, 400);
   TH2F *TvE2 = new TH2F("Tempo de voo vs Energia depositada 2", "Tempo de voo vs Energia Depositada", 60, 0, 45, 1000, 0, 400);
   TH2F *TvE3 = new TH2F("Tempo de voo vs Energia depositada 3", "Tempo de voo vs Energia Depositada", 60, 0, 45, 1000, 0, 400);
   TH2F *TvE4 = new TH2F("Tempo de voo vs Energia depositada 4", "Tempo de voo vs Energia Depositada", 60, 0, 45, 1000, 0, 400);
   TH2F *TvE5 = new TH2F("Tempo de voo vs Energia depositada 5", "Tempo de voo vs Energia Depositada", 60, 0, 45, 1000, 0, 400);
   TH2F *TvE6 = new TH2F("Tempo de voo vs Energia depositada 6", "Tempo de voo vs Energia Depositada", 60, 0, 45, 1000, 0, 400);
   TH2F *TvE7 = new TH2F("Tempo de voo vs Energia depositada 7", "Tempo de voo vs Energia Depositada", 60, 0, 45, 1000, 0, 400);

   Analysis(TTree *tree = 0);
   virtual ~Analysis();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop();
   virtual Bool_t Notify();
   virtual void Show(Long64_t entry = -1);
   void Soma();
   void MergeRootfile(TDirectory *target, TList *sourcelist);
   double Current1 = 16.80;
   double Current2 = 60.40;
};

#endif

#ifdef Analysis_cxx
Analysis::Analysis(TTree *tree) : fChain(0)
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   if (tree == 0)
   {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen())
      {
         f = new TFile("Memory Directory");
      }
      f->GetObject("SiTelescope", tree);

#else  // SINGLE_TREE
      stringstream filename;

      // Go to where are the root files
      string root = "/home/leo/Desktop/RIBRAS/ROOT/tree_";

      // File
      filename.precision(2);
      filename.setf(std::ios::fixed, std::ios::floatfield);
      filename << root << Analysis::Current1 << "_" << Analysis::Current2 << ".root";
      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain *chain = new TChain("SiTelescope", "");
      chain->Add(filename.str().data());
      tree = chain;
#endif // SINGLE_TREE
   }
   Init(tree);
}

Analysis::~Analysis()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis::LoadTree(Long64_t entry)
{
   // Set the environment to read one entry
   if (!fChain)
      return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0)
      return centry;
   if (fChain->GetTreeNumber() != fCurrent)
   {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree)
      return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("E0", E0, &b_E0);
   fChain->SetBranchAddress("E1", E1, &b_E1);
   fChain->SetBranchAddress("E2", E2, &b_E2);
   fChain->SetBranchAddress("E3", E3, &b_E3);
   fChain->SetBranchAddress("E4", E4, &b_E4);
   fChain->SetBranchAddress("E5", E5, &b_E5);
   fChain->SetBranchAddress("E6", E6, &b_E6);
   fChain->SetBranchAddress("E7", E7, &b_E7);
   fChain->SetBranchAddress("EKin0", EKin0, &b_EKin0);
   fChain->SetBranchAddress("EKin1", EKin1, &b_EKin1);
   fChain->SetBranchAddress("EKin2", EKin2, &b_EKin2);
   fChain->SetBranchAddress("EKin3", EKin3, &b_EKin3);
   fChain->SetBranchAddress("EKin4", EKin4, &b_EKin4);
   fChain->SetBranchAddress("EKin5", EKin5, &b_EKin5);
   fChain->SetBranchAddress("EKin6", EKin6, &b_EKin5);
   fChain->SetBranchAddress("EKin7", EKin7, &b_EKin7);
   fChain->SetBranchAddress("rThetaCM0", &rThetaCM0, &b_rThetaCM0);
   fChain->SetBranchAddress("rThetaCM1", &rThetaCM1, &b_rThetaCM1);
   fChain->SetBranchAddress("rThetaCM2", &rThetaCM2, &b_rThetaCM2);
   fChain->SetBranchAddress("rThetaCM3", &rThetaCM3, &b_rThetaCM3);
   fChain->SetBranchAddress("rThetaCM4", &rThetaCM4, &b_rThetaCM4);
   fChain->SetBranchAddress("rThetaCM5", &rThetaCM5, &b_rThetaCM5);
   fChain->SetBranchAddress("rThetaCM6", &rThetaCM6, &b_rThetaCM6);
   fChain->SetBranchAddress("rThetaCM7", &rThetaCM7, &b_rThetaCM7);
   fChain->SetBranchAddress("pos_x_det0", &pos_x_det0, &b_pos_x_det0);
   fChain->SetBranchAddress("pos_y_det0", &pos_y_det0, &b_pos_y_det0);
   fChain->SetBranchAddress("pos_z_det0", &pos_z_det0, &b_pos_z_det0);
   fChain->SetBranchAddress("pos_x_det1", &pos_x_det1, &b_pos_x_det1);
   fChain->SetBranchAddress("pos_y_det1", &pos_y_det1, &b_pos_y_det1);
   fChain->SetBranchAddress("pos_z_det1", &pos_z_det1, &b_pos_z_det1);
   fChain->SetBranchAddress("pos_x_det2", &pos_x_det2, &b_pos_x_det2);
   fChain->SetBranchAddress("pos_y_det2", &pos_y_det2, &b_pos_y_det2);
   fChain->SetBranchAddress("pos_z_det2", &pos_z_det2, &b_pos_z_det2);
   fChain->SetBranchAddress("pos_x_det3", &pos_x_det3, &b_pos_x_det3);
   fChain->SetBranchAddress("pos_y_det3", &pos_y_det3, &b_pos_y_det3);
   fChain->SetBranchAddress("pos_z_det3", &pos_z_det3, &b_pos_z_det3);
   fChain->SetBranchAddress("pos_x_det4", &pos_x_det4, &b_pos_x_det4);
   fChain->SetBranchAddress("pos_y_det4", &pos_y_det4, &b_pos_y_det4);
   fChain->SetBranchAddress("pos_z_det4", &pos_z_det4, &b_pos_z_det4);
   fChain->SetBranchAddress("pos_x_det5", &pos_x_det5, &b_pos_x_det5);
   fChain->SetBranchAddress("pos_y_det5", &pos_y_det5, &b_pos_y_det5);
   fChain->SetBranchAddress("pos_z_det5", &pos_z_det5, &b_pos_z_det5);
   fChain->SetBranchAddress("pos_x_det6", &pos_x_det6, &b_pos_x_det6);
   fChain->SetBranchAddress("pos_y_det6", &pos_y_det6, &b_pos_y_det6);
   fChain->SetBranchAddress("pos_z_det6", &pos_z_det6, &b_pos_z_det6);
   fChain->SetBranchAddress("pos_x_det7", &pos_x_det7, &b_pos_x_det7);
   fChain->SetBranchAddress("pos_y_det7", &pos_y_det7, &b_pos_y_det7);
   fChain->SetBranchAddress("pos_z_det7", &pos_z_det7, &b_pos_z_det7);
   fChain->SetBranchAddress("pos_x_target", &pos_x_target, &b_pos_x_target);
   fChain->SetBranchAddress("pos_y_target", &pos_y_target, &b_pos_y_target);
   fChain->SetBranchAddress("pos_z_target", &pos_z_target, &b_pos_z_target);
   fChain->SetBranchAddress("ener_det0", &ener_det0, &b_ener_det0);
   fChain->SetBranchAddress("ener_det1", &ener_det1, &b_ener_det1);
   fChain->SetBranchAddress("ener_det2", &ener_det2, &b_ener_det2);
   fChain->SetBranchAddress("ener_det3", &ener_det3, &b_ener_det3);
   fChain->SetBranchAddress("ener_det4", &ener_det4, &b_ener_det4);
   fChain->SetBranchAddress("ener_det5", &ener_det5, &b_ener_det5);
   fChain->SetBranchAddress("ener_det6", &ener_det6, &b_ener_det6);
   fChain->SetBranchAddress("ener_det7", &ener_det7, &b_ener_det7);
   fChain->SetBranchAddress("truthPosx", &truthPosx, &b_truthPosx);
   fChain->SetBranchAddress("truthPosy", &truthPosy, &b_truthPosy);
   fChain->SetBranchAddress("truthPosz", &truthPosz, &b_truthPosz);
   fChain->SetBranchAddress("truthAngle_theta", &truthAngle_theta, &b_truthAngle_theta);
   fChain->SetBranchAddress("truthAngle_phi", &truthAngle_phi, &b_truthAngle_phi);
   fChain->SetBranchAddress("px_dssd", &px_dssd, &b_px_dssd);
   fChain->SetBranchAddress("py_dssd", &py_dssd, &b_py_dssd);
   fChain->SetBranchAddress("pz_dssd", &pz_dssd, &b_pz_dssd);
   fChain->SetBranchAddress("t_dssd", &t_dssd, &b_t_dssd);
   fChain->SetBranchAddress("t_sili0", &t_sili0, &b_t_sili0);
   fChain->SetBranchAddress("t_sili1", &t_sili1, &b_t_sili1);
   fChain->SetBranchAddress("t_sili2", &t_sili2, &b_t_sili2);
   fChain->SetBranchAddress("t_sili3", &t_sili3, &b_t_sili3);
   fChain->SetBranchAddress("t_sili4", &t_sili4, &b_t_sili4);
   fChain->SetBranchAddress("t_sili5", &t_sili5, &b_t_sili5);
   fChain->SetBranchAddress("t_sili6", &t_sili6, &b_t_sili6);
   fChain->SetBranchAddress("t_sili7", &t_sili7, &b_t_sili7);
   fChain->SetBranchAddress("t_dssd2", &t_dssd2, &b_t_dssd2);
   fChain->SetBranchAddress("Etot", &ETot, &b_ETot);
   fChain->SetBranchAddress("Strip_Number0", &Strip_Number0, &b_Strip_Number0);
   fChain->SetBranchAddress("Strip_Number1", &Strip_Number1, &b_Strip_Number1);
   fChain->SetBranchAddress("Strip_Number2", &Strip_Number2, &b_Strip_Number2);
   fChain->SetBranchAddress("Strip_Number3", &Strip_Number3, &b_Strip_Number3);
   fChain->SetBranchAddress("Strip_Number4", &Strip_Number4, &b_Strip_Number4);
   fChain->SetBranchAddress("Strip_Number5", &Strip_Number5, &b_Strip_Number5);
   fChain->SetBranchAddress("Strip_Number6", &Strip_Number6, &b_Strip_Number6);
   fChain->SetBranchAddress("Strip_Number7", &Strip_Number7, &b_Strip_Number7);
   Notify();
}

Bool_t Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t Analysis::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_cxx
