#define Analysis_cxx
#include "Analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

///=============Mass definition============
#define aum 931.494043         // in MeV
#define e_mass 548.57990945e-6 // in aum

#define H1 1.00782503223 // in aum
#define H1_ch 1
#define H1_mass (H1 - H1_ch * e_mass) * aum //in MeV

#define B8 8.024607326 // in aum
#define B8_ch 5
#define B8_mass (B8 - B8_ch * e_mass) * aum //in MeV

#define C12 12.0000000000 // in aum
#define C12_ch 6
#define C12_mass (C12 - C12_ch * e_mass) * aum //in MeV

#define Ar40 39.96238312370 // in aum
#define Ar40_ch 18
#define Ar40_mass (Ar40 - Ar40_ch * e_mass) * aum //in MeV

///=============Two Body Kinematics===========
//Kinematics based on the Mandelstam invariant variables
double omega(double x, double y, double z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

double *kine_2b(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject)
{

   //in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u; //---Mandelstam variables
   double p1, p3;
   double J_LtoCM; //jacobian Lab to CM

   s = pow(m1, 2) + pow(m2, 2) + 2 * m2 * Et1;
   u = pow(m2, 2) + pow(m3, 2) - 2 * m2 * Et3;

   m4_ex = sqrt((cos(thetalab) * omega(s, pow(m1, 2), pow(m2, 2)) * omega(u, pow(m2, 2), pow(m3, 2)) - (s - pow(m1, 2) - pow(m2, 2)) * (pow(m2, 2) + pow(m3, 2) - u)) / (2 * pow(m2, 2)) + s + u - pow(m2, 2));
   Ex = m4_ex - m4;

   t = pow(m2, 2) + pow(m4_ex, 2) - 2 * m2 * Et4;

   //for normal kinematics
   theta_cm = acos((pow(s, 2) + s * (2 * t - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) + (pow(m1, 2) - pow(m2, 2)) * (pow(m3, 2) - pow(m4_ex, 2))) / (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m3, 2), pow(m4_ex, 2))));

   //for inverse kinematics Note: this angle corresponds to the recoil
   //theta_cm = TMath::Pi() - acos( ( pow(s,2) +s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) + (pow(m1,2) - pow(m2,2))*(pow(m3,2) - pow(m4_ex,2)) )/( omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4_ex,2))) ) ;

   p1 = sqrt(pow(Et1, 2) - pow(m1, 2));
   p3 = sqrt(pow(Et3, 2) - pow(m3, 2));

   J_LtoCM = abs(((omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m3, 2), pow(m4, 2))) / (4 * s * p1 * p3)) * (1. + Et1 / m2 - cos(thetalab) * (Et3 * p1) / (m2 * p3)));

   static double output[3];
   output[0] = theta_cm;
   output[1] = Ex;
   output[2] = J_LtoCM;
   return output;
}

void Analysis::Loop()
{
   // Number of entries
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //   In a ROOT session, you can do:
   //      root> .L Analysis.C
   //      root> Analysis t
   //      root> t.GetEntry(12); // Fill t data members with entry number 12
   //      root> t.Show();       // Show values of entry 12
   //      root> t.Show(16);     // Read and show values of entry 16
   //      root> t.Loop();       // Loop on all entries
   //

   //     This is the loop skeleton where:
   //    jentry is the global entry number in the chain
   //    ientry is the entry number in the current Tree
   //  Note that the argument to GetEntry must be:
   //    jentry for TChain::GetEntry
   //    ientry for TTree::GetEntry and TBranch::GetEntry
   //
   //       To read only selected branches, Insert statements like:
   // METHOD1:
   //    fChain->SetBranchStatus("*",0);  // disable all branches
   //    fChain->SetBranchStatus("branchname",1);  // activate branchname
   // METHOD2: replace line
   //    fChain->GetEntry(jentry);       //read all branches
   //by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0)
      return;
   // Creating TFile
   TFile *file = new TFile("histograms2B.root", "recreate");

   // Maker Style
   EvS0->SetMarkerStyle(8);
   EvS1->SetMarkerStyle(8);
   EvS2->SetMarkerStyle(8);
   EvS3->SetMarkerStyle(8);
   Ekin55->SetMarkerStyle(8);
   Teste->SetZTitle("Energia Armazenada");
   Teste->SetXTitle("Energia cinetica");
   Teste->SetYTitle("Strip");

   // Total Energy per strip
   Double_t E0Total[100] = {0};
   Double_t E1Total[100] = {0};
   Double_t E2Total[100] = {0};
   Double_t E3Total[100] = {0};

   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;

      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      for (int i = 0; i < 60; i++)
      {
         Teste->Fill(EKin0[i], i , E0[i]);
      }
      Ekin55->Fill(EKin0[55], E0[55] );

      for (int i = 0; i < 60; i++)
      {
         E0Total[i] += E0[i];
         E1Total[i] += E1[i];
         E2Total[i] += E2[i];
         E3Total[i] += E3[i];
      }
   }

   for (int n = 0; n < 60; n++)
   {
      EvS0->Fill(n, E0Total[n]);
      EvS1->Fill(n, E1Total[n]);
      EvS2->Fill(n, E2Total[n]);
      EvS3->Fill(n, E3Total[n]);
   }
   EvS0->Write();
   EvS1->Write();
   EvS2->Write();
   EvS3->Write();
   Ekin55->Write();
   Teste->Write();

   file->Close();
}
