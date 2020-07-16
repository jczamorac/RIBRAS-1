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
#include <numeric>   // std::iota
#include <algorithm> // std::sort

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
   Float_t t_sili1;
   Float_t t_sili2;
   Float_t t_dssd2;
   Float_t Ekin_dssd2;
   Int_t Strip_Number;

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
   TBranch *b_t_sili1;          //!
   TBranch *b_t_sili2;          //!
   TBranch *b_t_dssd2;          //!
   TBranch *b_Ekin_dssd2;       //!
   TBranch *b_Strip_Number;     //!

   // 2D Histograms
   TH2F *EvS0 = new TH2F("Energy vs Strip 0", // Histogram's name
                         "Energy vs Strip",   // Histogram's title
                         60,                  // Bins on X
                         0,                   // Initial X
                         60,                  // Final X
                         60,                  // Bins on Y
                         0,                   // Initial Y
                         0);                  // Final Y
   TH2F *EvS1 = new TH2F("Energy vs Strip 1", "Energy vs Strip", 60, 0, 60, 60, 0, 0);
   TH2F *EvS2 = new TH2F("Energy vs Strip 2", "Energy vs Strip", 60, 0, 60, 60, 0, 0);
   TH2F *EvS3 = new TH2F("Energy vs Strip 3", "Energy vs Strip", 60, 0, 60, 60, 0, 0);
   TH2F *Ekin55 = new TH2F("Energy vs Ekin", "Energy vs Ekin", 60, 0, 60, 60, 0, 0);
   TH3F *Teste = new TH3F("Energy vs Ekin", "Energy vs Ekin", 60, 0, 60, 60, 0, 0, 60, 0, 0);

   Analysis(TTree *tree = 0);
   virtual ~Analysis();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop();
   virtual Bool_t Notify();
   virtual void Show(Long64_t entry = -1);
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

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain *chain = new TChain("SiTelescope", "");
      chain->Add("tree_40.00.root/SiTelescope");
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
   fChain->SetBranchAddress("t_sili1", &t_sili1, &b_t_sili1);
   fChain->SetBranchAddress("t_sili2", &t_sili2, &b_t_sili2);
   fChain->SetBranchAddress("t_dssd2", &t_dssd2, &b_t_dssd2);
   fChain->SetBranchAddress("Ekin_dssd2", &Ekin_dssd2, &b_Ekin_dssd2);
   fChain->SetBranchAddress("Strip_Number", &Strip_Number, &b_Strip_Number);
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
