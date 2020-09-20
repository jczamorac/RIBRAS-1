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

void Analysis::Soma()
{
   TFile *file = new TFile("0.root", "open");

   TH2F *EvStrip0 = (TH2F *)file->Get("Energy vs Strip 0");
   TH2F *EvStrip1 = (TH2F *)file->Get("Energy vs Strip 1");
   TH2F *EvStrip2 = (TH2F *)file->Get("Energy vs Strip 2");
   TH2F *EvStrip3 = (TH2F *)file->Get("Energy vs Strip 3");
   TH2F *EvStrip4 = (TH2F *)file->Get("Energy vs Strip 4");
   TH2F *EvStrip5 = (TH2F *)file->Get("Energy vs Strip 5");
   TH2F *EvStrip6 = (TH2F *)file->Get("Energy vs Strip 6");
   TH2F *EvStrip7 = (TH2F *)file->Get("Energy vs Strip 7");

   TH2F *EvThetaCM0 = (TH2F *)file->Get("Energy vs ThetaCM0");
   TH2F *EvThetaCM1 = (TH2F *)file->Get("Energy vs ThetaCM1");
   TH2F *EvThetaCM2 = (TH2F *)file->Get("Energy vs ThetaCM2");
   TH2F *EvThetaCM3 = (TH2F *)file->Get("Energy vs ThetaCM3");
   TH2F *EvThetaCM4 = (TH2F *)file->Get("Energy vs ThetaCM4");
   TH2F *EvThetaCM5 = (TH2F *)file->Get("Energy vs ThetaCM5");
   TH2F *EvThetaCM6 = (TH2F *)file->Get("Energy vs ThetaCM6");
   TH2F *EvThetaCM7 = (TH2F *)file->Get("Energy vs ThetaCM7");

   TH2F *EkinvThetaCM0 = (TH2F *)file->Get("Ekin vs ThetaCM 0");
   TH2F *EkinvThetaCM1 = (TH2F *)file->Get("Ekin vs ThetaCM 1");
   TH2F *EkinvThetaCM2 = (TH2F *)file->Get("Ekin vs ThetaCM 2");
   TH2F *EkinvThetaCM3 = (TH2F *)file->Get("Ekin vs ThetaCM 3");
   TH2F *EkinvThetaCM4 = (TH2F *)file->Get("Ekin vs ThetaCM 4");
   TH2F *EkinvThetaCM5 = (TH2F *)file->Get("Ekin vs ThetaCM 5");
   TH2F *EkinvThetaCM6 = (TH2F *)file->Get("Ekin vs ThetaCM 6");
   TH2F *EkinvThetaCM7 = (TH2F *)file->Get("Ekin vs ThetaCM 7");

   TH2F *EkinvStrip0 = (TH2F *)file->Get("Ekin vs Strip 0");
   TH2F *EkinvStrip1 = (TH2F *)file->Get("Ekin vs Strip 1");
   TH2F *EkinvStrip2 = (TH2F *)file->Get("Ekin vs Strip 2");
   TH2F *EkinvStrip3 = (TH2F *)file->Get("Ekin vs Strip 3");
   TH2F *EkinvStrip4 = (TH2F *)file->Get("Ekin vs Strip 4");
   TH2F *EkinvStrip5 = (TH2F *)file->Get("Ekin vs Strip 5");
   TH2F *EkinvStrip6 = (TH2F *)file->Get("Ekin vs Strip 6");
   TH2F *EkinvStrip7 = (TH2F *)file->Get("Ekin vs Strip 7");

   TH2F *TvStrip0 = (TH2F *)file->Get("Tempo de voo vs Strip 0");
   TH2F *TvStrip1 = (TH2F *)file->Get("Tempo de voo vs Strip 1");
   TH2F *TvStrip2 = (TH2F *)file->Get("Tempo de voo vs Strip 2");
   TH2F *TvStrip3 = (TH2F *)file->Get("Tempo de voo vs Strip 3");
   TH2F *TvStrip4 = (TH2F *)file->Get("Tempo de voo vs Strip 4");
   TH2F *TvStrip5 = (TH2F *)file->Get("Tempo de voo vs Strip 5");
   TH2F *TvStrip6 = (TH2F *)file->Get("Tempo de voo vs Strip 6");
   TH2F *TvStrip7 = (TH2F *)file->Get("Tempo de voo vs Strip 7");

   TH2F *ThetaCMvStrip0 = (TH2F *)file->Get("Strip vs Theta 0");
   TH2F *ThetaCMvStrip1 = (TH2F *)file->Get("Strip vs Theta 1");
   TH2F *ThetaCMvStrip2 = (TH2F *)file->Get("Strip vs Theta 2");
   TH2F *ThetaCMvStrip3 = (TH2F *)file->Get("Strip vs Theta 3");
   TH2F *ThetaCMvStrip4 = (TH2F *)file->Get("Strip vs Theta 4");
   TH2F *ThetaCMvStrip5 = (TH2F *)file->Get("Strip vs Theta 5");
   TH2F *ThetaCMvStrip6 = (TH2F *)file->Get("Strip vs Theta 6");
   TH2F *ThetaCMvStrip7 = (TH2F *)file->Get("Strip vs Theta 7");

   // TH2F *TvEdep0 = (TH2F *)file->Get("Tempo de voo vs Energia depositada 0");
   // TH2F *TvEdep1 = (TH2F *)file->Get("Tempo de voo vs Energia depositada 1");
   // TH2F *TvEdep2 = (TH2F *)file->Get("Tempo de voo vs Energia depositada 2");
   // TH2F *TvEdep3 = (TH2F *)file->Get("Tempo de voo vs Energia depositada 3");
   // TH2F *TvEdep4 = (TH2F *)file->Get("Tempo de voo vs Energia depositada 4");
   // TH2F *TvEdep5 = (TH2F *)file->Get("Tempo de voo vs Energia depositada 5");
   // TH2F *TvEdep6 = (TH2F *)file->Get("Tempo de voo vs Energia depositada 6");
   // TH2F *TvEdep7 = (TH2F *)file->Get("Tempo de voo vs Energia depositada 7");

   // TH2F *TvEKin0 = (TH2F *)file->Get("Tempo de voo vs Energia cinetica 0");
   // TH2F *TvEKin1 = (TH2F *)file->Get("Tempo de voo vs Energia cinetica 1");
   // TH2F *TvEKin2 = (TH2F *)file->Get("Tempo de voo vs Energia cinetica 2");
   // TH2F *TvEKin3 = (TH2F *)file->Get("Tempo de voo vs Energia cinetica 3");
   // TH2F *TvEKin4 = (TH2F *)file->Get("Tempo de voo vs Energia cinetica 4");
   // TH2F *TvEKin5 = (TH2F *)file->Get("Tempo de voo vs Energia cinetica 5");
   // TH2F *TvEKin6 = (TH2F *)file->Get("Tempo de voo vs Energia cinetica 6");
   // TH2F *TvEKin7 = (TH2F *)file->Get("Tempo de voo vs Energia cinetica 7");

   TFile *file1 = new TFile("Histogramas.root", "recreate");

   TH2F *EvStripTraseiro = (TH2F *)EvStrip0->Clone("EvStrip Traseiros");
   EvStripTraseiro->Add(EvStrip1);
   EvStripTraseiro->Add(EvStrip2);
   EvStripTraseiro->Add(EvStrip3);
   EvStripTraseiro->Write();

   TH2F *EvStripFrontais = (TH2F *)EvStrip4->Clone("EvStrip Frontais");
   EvStripFrontais->Add(EvStrip5);
   EvStripFrontais->Add(EvStrip6);
   EvStripFrontais->Add(EvStrip7);
   EvStripFrontais->Write();

   TH2F *EvThetaCMTraseiro = (TH2F *)EvThetaCM0->Clone("EvThetaCM Traseiros");
   EvThetaCMTraseiro->Add(EvThetaCM1);
   EvThetaCMTraseiro->Add(EvThetaCM2);
   EvThetaCMTraseiro->Add(EvThetaCM3);
   EvThetaCMTraseiro->Write();

   TH2F *EvThetaCMFrontais = (TH2F *)EvThetaCM4->Clone("EvThetaCM Frontais");
   EvThetaCMFrontais->Add(EvThetaCM5);
   EvThetaCMFrontais->Add(EvThetaCM6);
   EvThetaCMFrontais->Add(EvThetaCM7);
   EvThetaCMFrontais->Write();

   TH2F *EkinvThetaCMTraseiro = (TH2F *)EkinvThetaCM0->Clone("EkinvThetaCM Traseiros");
   EkinvThetaCMTraseiro->Add(EkinvThetaCM1);
   EkinvThetaCMTraseiro->Add(EkinvThetaCM2);
   EkinvThetaCMTraseiro->Add(EkinvThetaCM3);
   EkinvThetaCMTraseiro->Write();

   TH2F *EkinvThetaCMFrontais = (TH2F *)EkinvThetaCM4->Clone("EkinvThetaCM Frontais");
   EkinvThetaCMFrontais->Add(EkinvThetaCM5);
   EkinvThetaCMFrontais->Add(EkinvThetaCM6);
   EkinvThetaCMFrontais->Add(EkinvThetaCM7);
   EkinvThetaCMFrontais->Write();

   TH2F *EkinvStripTraseiro = (TH2F *)EkinvStrip0->Clone("EkinvStrip Traseiros");
   EkinvStripTraseiro->Add(EkinvStrip1);
   EkinvStripTraseiro->Add(EkinvStrip2);
   EkinvStripTraseiro->Add(EkinvStrip3);
   EkinvStripTraseiro->Write();

   TH2F *EkinvStripFrontais = (TH2F *)EkinvStrip4->Clone("EkinvStrip Frontais");
   EkinvStripFrontais->Add(EkinvStrip5);
   EkinvStripFrontais->Add(EkinvStrip6);
   EkinvStripFrontais->Add(EkinvStrip7);
   EkinvStripFrontais->Write();

   TH2F *TvStripTraseiro = (TH2F *)TvStrip0->Clone("TvStrip Traseiros");
   TvStripTraseiro->Add(TvStrip1);
   TvStripTraseiro->Add(TvStrip2);
   TvStripTraseiro->Add(TvStrip3);
   TvStripTraseiro->Write();

   TH2F *TvStripFrontais = (TH2F *)TvStrip4->Clone("TvStrip Frontais");
   TvStripFrontais->Add(TvStrip5);
   TvStripFrontais->Add(TvStrip6);
   TvStripFrontais->Add(TvStrip7);
   TvStripFrontais->Write();

   TH2F *ThetaCMvStripTraseiro = (TH2F *)ThetaCMvStrip0->Clone("ThetaCMvStrip Traseiros");
   ThetaCMvStripTraseiro->Add(ThetaCMvStrip1);
   ThetaCMvStripTraseiro->Add(ThetaCMvStrip2);
   ThetaCMvStripTraseiro->Add(ThetaCMvStrip3);
   ThetaCMvStripTraseiro->Write();

   TH2F *ThetaCMvStripFrontais = (TH2F *)ThetaCMvStrip4->Clone("ThetaCMvStrip Frontais");
   ThetaCMvStripFrontais->Add(ThetaCMvStrip5);
   ThetaCMvStripFrontais->Add(ThetaCMvStrip6);
   ThetaCMvStripFrontais->Add(ThetaCMvStrip7);
   ThetaCMvStripFrontais->Write();

   TH2F *TvEdepTraseiros = (TH2F *)TvEdep0->Clone("TvEdep Traseiros");
   TvEdepTraseiros->Add(TvEdep1);
   TvEdepTraseiros->Add(TvEdep2);
   TvEdepTraseiros->Add(TvEdep3);
   TvEdepTraseiros->Write();

   TH2F *TvEdepFrontais = (TH2F *)TvEdep4->Clone("TvEdep Frontais");
   TvEdepFrontais->Add(TvEdep5);
   TvEdepFrontais->Add(TvEdep6);
   TvEdepFrontais->Add(TvEdep7);
   TvEdepFrontais->Write();

   // TH2F *TvEkinTraseiros = (TH2F *)TvEKin0->Clone("TvEkin Traseiros");
   // TvEkinTraseiros->Add(TvEKin1);
   // TvEkinTraseiros->Add(TvEKin2);
   // TvEkinTraseiros->Add(TvEKin3);
   // TvEkinTraseiros->Write();

   // TH2F *TvEkinFrontais = (TH2F *)TvEKin4->Clone("TvEdep Frontais");
   // TvEkinFrontais->Add(TvEKin5);
   // TvEkinFrontais->Add(TvEKin6);
   // TvEkinFrontais->Add(TvEKin7);
   // TvEkinFrontais->Write();
}

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
   // File name
   stringstream f;
   f << "Histograms2D - " << Analysis::Current1 << "_" << Analysis::Current2 << "A.root";

   // Number of entries
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   if (fChain == 0)
      return;
   // Creating TFile
   TFile *file = new TFile(f.str().data(), "recreate");

   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      // Strip vs Edep
      // Strip vs Recoil Theta CM
      // Strip vs Time
      if (E0[Strip_Number0] > 0)
      {
         EvS0->Fill(Strip_Number0, E0[Strip_Number0]);
         EvThetaCM0->Fill(rThetaCM0[Strip_Number0], E0[Strip_Number0]);
         ThetavS0->Fill(Strip_Number0, rThetaCM0[Strip_Number0]);
         TimevS0->Fill(Strip_Number0, t_sili0);
         TvE0->Fill(E0[Strip_Number0], t_sili0);
      }
      if (E1[Strip_Number1] > 0)
      {
         EvS1->Fill(Strip_Number1, E1[Strip_Number1]);
         EvThetaCM1->Fill(rThetaCM1[Strip_Number1], E1[Strip_Number1]);
         ThetavS1->Fill(Strip_Number1, rThetaCM1[Strip_Number1]);
         TimevS1->Fill(Strip_Number1, t_sili1);
         TvE1->Fill(E1[Strip_Number1], t_sili1);
      }
      if (E2[Strip_Number2] > 0)
      {
         EvS2->Fill(Strip_Number2, E2[Strip_Number2]);
         EvThetaCM2->Fill(rThetaCM2[Strip_Number2], E2[Strip_Number2]);
         ThetavS2->Fill(Strip_Number2, rThetaCM2[Strip_Number2]);
         TimevS2->Fill(Strip_Number2, t_sili2);
         TvE2->Fill(E2[Strip_Number2], t_sili2);
      }
      if (E3[Strip_Number3] > 0)
      {
         EvS3->Fill(Strip_Number3, E3[Strip_Number3]);
         EvThetaCM3->Fill(rThetaCM3[Strip_Number3], E3[Strip_Number3]);
         ThetavS3->Fill(Strip_Number3, rThetaCM3[Strip_Number3]);
         TimevS3->Fill(Strip_Number3, t_sili3);
         TvE3->Fill(E3[Strip_Number3], t_sili3);
      }
      if (E4[Strip_Number4] > 0)
      {
         EvS4->Fill(Strip_Number4, E4[Strip_Number4]);
         EvThetaCM4->Fill(rThetaCM4[Strip_Number4], E4[Strip_Number4]);
         ThetavS4->Fill(Strip_Number4, rThetaCM4[Strip_Number4]);
         TimevS4->Fill(Strip_Number4, t_sili4);
         TvE4->Fill(E4[Strip_Number4], t_sili4);
      }
      if (E5[Strip_Number5] > 0)
      {
         EvS5->Fill(Strip_Number5, E5[Strip_Number5]);
         EvThetaCM5->Fill(rThetaCM5[Strip_Number5], E5[Strip_Number5]);
         ThetavS5->Fill(Strip_Number5, rThetaCM5[Strip_Number5]);
         TimevS5->Fill(Strip_Number5, t_sili5);
         TvE5->Fill(E5[Strip_Number5], t_sili5);
      }
      if (E6[Strip_Number6] > 0)
      {
         EvS6->Fill(Strip_Number6, E6[Strip_Number6]);
         EvThetaCM6->Fill(rThetaCM6[Strip_Number6], E6[Strip_Number6]);
         ThetavS6->Fill(Strip_Number6, rThetaCM6[Strip_Number6]);
         TimevS6->Fill(Strip_Number6, t_sili6);
         TvE6->Fill(E6[Strip_Number6], t_sili6);
      }
      if (E7[Strip_Number7] > 0)
      {
         EvS7->Fill(Strip_Number7, E7[Strip_Number7]);
         EvThetaCM7->Fill(rThetaCM7[Strip_Number7], E7[Strip_Number7]);
         ThetavS7->Fill(Strip_Number7, rThetaCM7[Strip_Number7]);
         TimevS7->Fill(Strip_Number7, t_sili7);
         TvE7->Fill(E7[Strip_Number7], t_sili7);
      }

      // Strip vs Kinect Energy
      // Recoil Theta CM vs Kinect Energy
      if (EKin0[Strip_Number0] > 0)
      {
         EKinvS0->Fill(Strip_Number0, EKin0[Strip_Number0]);
         EkinvTheta0->Fill(rThetaCM0[Strip_Number0], EKin0[Strip_Number0]);
         TvEkin0->Fill(EKin0[Strip_Number0], t_sili0);
      }

      if (EKin1[Strip_Number1] > 0)
      {
         EKinvS1->Fill(Strip_Number1, EKin1[Strip_Number1]);
         EkinvTheta1->Fill(rThetaCM1[Strip_Number1], EKin1[Strip_Number1]);
         TvEkin1->Fill(EKin1[Strip_Number1], t_sili1);
      }

      if (EKin2[Strip_Number2] > 0)
      {
         EKinvS2->Fill(Strip_Number2, EKin2[Strip_Number2]);
         EkinvTheta2->Fill(rThetaCM2[Strip_Number2], EKin2[Strip_Number2]);
         TvEkin2->Fill(EKin2[Strip_Number2], t_sili2);
      }

      if (EKin3[Strip_Number3] > 0)
      {
         EKinvS3->Fill(Strip_Number3, EKin3[Strip_Number3]);
         EkinvTheta3->Fill(rThetaCM3[Strip_Number3], EKin3[Strip_Number3]);
         TvEkin3->Fill(EKin3[Strip_Number3], t_sili3);
      }

      if (EKin4[Strip_Number4] > 0)
      {
         EKinvS4->Fill(Strip_Number4, EKin4[Strip_Number4]);
         EkinvTheta4->Fill(rThetaCM4[Strip_Number4], EKin4[Strip_Number4]);
         TvEkin4->Fill(EKin4[Strip_Number4], t_sili4);
      }

      if (EKin5[Strip_Number5] > 0)
      {
         EKinvS5->Fill(Strip_Number5, EKin5[Strip_Number5]);
         EkinvTheta5->Fill(rThetaCM5[Strip_Number5], EKin5[Strip_Number5]);
         TvEkin5->Fill(EKin5[Strip_Number5], t_sili5);
      }

      if (EKin6[Strip_Number6] > 0)
      {
         EKinvS6->Fill(Strip_Number6, EKin6[Strip_Number6]);
         EkinvTheta6->Fill(rThetaCM6[Strip_Number6], EKin6[Strip_Number6]);
         TvEkin6->Fill(EKin6[Strip_Number6], t_sili6);
      }

      if (EKin7[Strip_Number7] > 0)
      {
         EKinvS7->Fill(Strip_Number7, EKin7[Strip_Number7]);
         EkinvTheta7->Fill(rThetaCM7[Strip_Number7], EKin7[Strip_Number7]);
         TvEkin7->Fill(EKin7[Strip_Number7], t_sili7);
      }
   }

   EvS0->Write();
   EvS1->Write();
   EvS2->Write();
   EvS3->Write();
   EvS4->Write();
   EvS5->Write();
   EvS6->Write();
   EvS7->Write();
   EvThetaCM0->Write();
   EvThetaCM1->Write();
   EvThetaCM2->Write();
   EvThetaCM3->Write();
   EvThetaCM4->Write();
   EvThetaCM5->Write();
   EvThetaCM6->Write();
   EvThetaCM7->Write();
   EkinvTheta0->Write();
   EkinvTheta1->Write();
   EkinvTheta2->Write();
   EkinvTheta3->Write();
   EkinvTheta4->Write();
   EkinvTheta5->Write();
   EkinvTheta6->Write();
   EkinvTheta7->Write();
   EKinvS0->Write();
   EKinvS1->Write();
   EKinvS2->Write();
   EKinvS3->Write();
   EKinvS4->Write();
   EKinvS5->Write();
   EKinvS6->Write();
   EKinvS7->Write();
   TimevS0->Write();
   TimevS1->Write();
   TimevS2->Write();
   TimevS3->Write();
   TimevS4->Write();
   TimevS5->Write();
   TimevS6->Write();
   TimevS7->Write();
   ThetavS0->Write();
   ThetavS1->Write();
   ThetavS2->Write();
   ThetavS3->Write();
   ThetavS4->Write();
   ThetavS5->Write();
   ThetavS6->Write();
   ThetavS7->Write();
   TvEkin0->Write();
   TvEkin1->Write();
   TvEkin2->Write();
   TvEkin3->Write();
   TvEkin4->Write();
   TvEkin5->Write();
   TvEkin6->Write();
   TvEkin7->Write();
   TvE0->Write();
   TvE1->Write();
   TvE2->Write();
   TvE3->Write();
   TvE4->Write();
   TvE5->Write();
   TvE6->Write();
   TvE7->Write();

   file->Close();
}
