// $Id: RootSaver.cc 94 2010-01-26 13:18:30Z adotti $
/**
  * @file   RootSaver.cc
  *
  * @date   17 Dec 2009
  * @author adotti
  *
  * @brief  Implements class RootSaver.
  */

// Local Headers
#include "DetectorConstruction.hh"
#include "RootSaver.hh"
//#include "SiDigi.hh"
#include "SiHit.hh"
#include "G4DigiManager.hh"

// Root headers
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TGraph.h"

// Default headers
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>

using namespace std;
using namespace CLHEP;

//-------------------------------------------------------------------------//

RootSaver::RootSaver() : //Initializing parameters
                         rootTree(0), runCounter(0), nStrips(0),
                         Signal0(0), Signal1(0), Signal2(0),
                         Signal3(0), Signal5(0), Signal6(0),
                         Signal7(0),
                         TruthPosx(0), TruthPosy(0), TruthPosz(0),
                         TruthAngle_theta(0), TruthAngle_phi(0),
                         Px_dssd(0), Py_dssd(0), Pz_dssd(0),
                         T_dssd(0),
                         Ekin_dssd2(0), StripNumber(0)
{
}

//-------------------------------------------------------------------------//

RootSaver::~RootSaver()
{
        //Close current file if needed
        if (rootTree)
        {
                CloseTree();
        }
}

//-------------------------------------------------------------------------//

void RootSaver::CreateTree(const std::string &fileName, const std::string &treeName)
{
        if (rootTree)
        {
                std::cerr << "TTree already created, first call CloseTree" << std::endl;
                return;
        }

        // Path to where ROOT should save the files
        //G4String Path = "/home/leo/Desktop/RIBRAS/ROOT/";
        G4String Path = "";

        // Creating ROOT file
        std::ostringstream fn;
        fn << Path << fileName << "_run_" << runCounter++ << ".root";

        // Create a new file and open it for writing, if the file already exists the file is overwritten
        TFile *rootFile = TFile::Open(fn.str().data(), "recreate");

        if (rootFile == 0 || rootFile->IsZombie())
        {
                G4cerr << "Error opening the file: " << fn.str() << " TTree will not be saved." << G4endl;
                return;
        }

        rootTree = new TTree(treeName.data(), treeName.data());

        nStrips = 60;
        Signal0 = new Float_t[nStrips];
        Signal1 = new Float_t[nStrips];
        Signal2 = new Float_t[nStrips];
        Signal3 = new Float_t[nStrips];
        Signal4 = new Float_t[nStrips];
        Signal5 = new Float_t[nStrips];
        Signal6 = new Float_t[nStrips];
        Signal7 = new Float_t[nStrips];

        for (Int_t strip = 0; strip < nStrips; ++strip)
        {
                Signal0[strip] = 0;
                Signal1[strip] = 0;
                Signal2[strip] = 0;
                Signal3[strip] = 0;
                Signal4[strip] = 0;
                Signal5[strip] = 0;
                Signal6[strip] = 0;
                Signal7[strip] = 0;
        }

        // Digits variables
        //-------- Total energy per strip ----------//
        rootTree->Branch("E0", Signal0, "E0[60]/F");
        rootTree->Branch("E1", Signal1, "E1[60]/F");
        rootTree->Branch("E2", Signal2, "E2[60]/F");
        rootTree->Branch("E3", Signal3, "E3[60]/F");
        rootTree->Branch("E4", Signal4, "E4[60]/F");
        rootTree->Branch("E5", Signal5, "E5[60]/F");
        rootTree->Branch("E6", Signal6, "E6[60]/F");
        rootTree->Branch("E7", Signal7, "E7[60]/F");
        //------------------------------------------//

        //--- Hit Position ( Detector Reference ) ---//
        rootTree->Branch("pos_x_det0", &Pos_x_det[0]);
        rootTree->Branch("pos_y_det0", &Pos_y_det[0]);
        rootTree->Branch("pos_z_det0", &Pos_z_det[0]);

        rootTree->Branch("pos_x_det1", &Pos_x_det[1]);
        rootTree->Branch("pos_y_det1", &Pos_y_det[1]);
        rootTree->Branch("pos_z_det1", &Pos_z_det[1]);

        rootTree->Branch("pos_x_det2", &Pos_x_det[2]);
        rootTree->Branch("pos_y_det2", &Pos_y_det[2]);
        rootTree->Branch("pos_z_det2", &Pos_z_det[2]);

        rootTree->Branch("pos_x_det3", &Pos_x_det[3]);
        rootTree->Branch("pos_y_det3", &Pos_y_det[3]);
        rootTree->Branch("pos_z_det3", &Pos_z_det[3]);

        rootTree->Branch("pos_x_det4", &Pos_x_det[4]);
        rootTree->Branch("pos_y_det4", &Pos_y_det[4]);
        rootTree->Branch("pos_z_det4", &Pos_z_det[4]);

        rootTree->Branch("pos_x_det5", &Pos_x_det[5]);
        rootTree->Branch("pos_y_det5", &Pos_y_det[5]);
        rootTree->Branch("pos_z_det5", &Pos_z_det[5]);

        rootTree->Branch("pos_x_det6", &Pos_x_det[6]);
        rootTree->Branch("pos_y_det6", &Pos_y_det[6]);
        rootTree->Branch("pos_z_det6", &Pos_z_det[6]);

        rootTree->Branch("pos_x_det7", &Pos_x_det[7]);
        rootTree->Branch("pos_y_det7", &Pos_y_det[7]);
        rootTree->Branch("pos_z_det7", &Pos_z_det[7]);
        //-------------------------------------------//

        //------- Total energy in the detector-------//
        rootTree->Branch("ener_det0", &E_det[1]);
        rootTree->Branch("ener_det1", &E_det[2]);
        rootTree->Branch("ener_det2", &E_det[3]);
        rootTree->Branch("ener_det3", &E_det[4]);
        rootTree->Branch("ener_det4", &E_det[5]);
        rootTree->Branch("ener_det5", &E_det[6]);
        rootTree->Branch("ener_det6", &E_det[7]);
        rootTree->Branch("ener_det7", &E_det[8]);
        //-------------------------------------------//

        rootTree->Branch("truthPosx", &TruthPosx);
        rootTree->Branch("truthPosy", &TruthPosy);
        rootTree->Branch("truthPosz", &TruthPosz);

        rootTree->Branch("truthAngle_theta", &TruthAngle_theta);
        rootTree->Branch("truthAngle_phi", &TruthAngle_phi);

        rootTree->Branch("px_dssd", &Px_dssd);
        rootTree->Branch("py_dssd", &Py_dssd);
        rootTree->Branch("pz_dssd", &Pz_dssd);

        rootTree->Branch("t_dssd", &T_dssd);
        rootTree->Branch("t_sili1", &T_sili[1]);
        rootTree->Branch("t_sili2", &T_sili[2]);
        rootTree->Branch("t_dssd2", &T_dssd);
        rootTree->Branch("Ekin_dssd2", &Ekin_dssd2);

        rootTree->Branch("Strip_Number", &StripNumber);
}

//-------------------------------------------------------------------------//

void RootSaver::CloseTree()
{
        // Check if ROOT TTree exists,
        // in case get the associated file and close it.
        // Note that if a TFile goes above 2GB a new file
        // will be automatically opened. We have thus to get,
        // from the TTree the current opened file
        if (rootTree)
        {

                G4cout << "Writing ROOT TTree: " << rootTree->GetName() << G4endl;
                //rootTree->Print();
                rootTree->Write();
                TFile *currentFile = rootTree->GetCurrentFile();
                if (currentFile == 0 || currentFile->IsZombie())
                {
                        G4cerr << "Error closing TFile " << G4endl;
                        return;
                }
                currentFile->Close();
                //The root is automatically deleted.

                // n = 0;
                // auto canvas = new TCanvas();
                // G4String root = "/home/leo/Desktop/RIBRAS/ROOT/tree_run_";
                // std::ostringstream f;
                // f << root << runCounter - 1 << ".root";
                // TFile *File = TFile::Open(f.str().data());
                // TTreeReader myReader("SiTelescope", File);

                // TTreeReaderValue<Float_t> posx(myReader, "pos_x_det0");
                // TTreeReaderValue<Float_t> posy(myReader, "pos_y_det0");

                // G4double x[100], y[100];

                // while (myReader.Next())
                // {
                //         x[n] = *posx;
                //         y[n] = *posy;
                //         n++;
                // }

                // TGraph graph(n, x, y);

                // graph.SetTitle("Position");
                // graph.SetMarkerColor(kRed);
                // graph.SetMarkerStyle(21);
                // graph.DrawClone("APE");

                // std::ostringstream fn;
                // G4String Path = "/home/leo/Desktop/RIBRAS/Graficos_Posicao/";
                // fn << Path << "posicao"
                //    << "_run_" << runCounter << ".png";
                // canvas->Print(fn.str().data());

                rootTree = 0;
                delete[] Signal0;
                delete[] Signal1;
                delete[] Signal2;
                delete[] Signal3;
                delete[] Signal4;
                delete[] Signal5;
                delete[] Signal6;
                delete[] Signal7;
        }

        G4cout << "Total Hits: " << TotalHits << G4endl;
}

//-------------------------------------------------------------------------//

void RootSaver::AddEvent(const SiHitCollection *const hits, const G4ThreeVector &primPos, const G4ThreeVector &primMom)
{
        //If root TTree is not created ends
        if (rootTree == 0)
        {
                return;
        }

        // Store Hits information
        if (hits->entries())
        {
                hits;
                G4int nHits = hits->entries();
                // Set defaults values
                Ekin_dssd2 = 0;

                Px_dssd = -1000;
                Py_dssd = -1000;
                Pz_dssd = -1000;

                T_dssd = -1000;

                StripNumber = -1000;

                TotalHits++;

                // Loop on all hits, consider only the hits with isPrimary flag
                // Position is weighted average of hit x()

                for (G4int h = 0; (h < nHits); ++h)
                {
                        const SiHit *hit = static_cast<const SiHit *>(hits->GetHit(h));

                        // Getting what strip ocurred a hit
                        G4int stripNum = hit->GetStripNumber();
                        StripNumber = stripNum;

                        // Same for detector
                        G4int planeNum = hit->GetPlaneNumber();

                        // Getting position of hit (detector reference)
                        G4ThreeVector pos = hit->GetPosition();
                        G4double x = pos.x();
                        G4double y = pos.y();
                        G4double z = pos.z();

                        // Getting momentum
                        G4ThreeVector momentum = hit->GetIncidenceMomentumDirection();
                        G4double momet_x = momentum.x();
                        G4double momet_y = momentum.y();
                        G4double momet_z = momentum.z();

                        // Hit time
                        G4double tiempo = hit->GetIncidenceTime();

                        // We save xyz in mm (detector coordinates)
                        x /= CLHEP::mm;
                        y /= CLHEP::mm;
                        z /= CLHEP::mm;

                        // Time in nanoseconds
                        tiempo /= CLHEP::ns;

                        // We save energy in MeV
                        Float_t edep = static_cast<Float_t>(hit->GetEdep());
                        edep /= CLHEP::MeV;

                        // Saving information for each detector
                        if (planeNum == 0)
                        {
                                Pos_x_det[planeNum] = x;
                                Pos_y_det[planeNum] = y;
                                Pos_z_det[planeNum] = z;
                                T_sili[planeNum] = tiempo;
                                double Econv = Digital(edep);
                                E_det[planeNum] += Econv;
                                Signal0[stripNum] += Econv;
                        }
                        else if (planeNum == 1)
                        {
                                Pos_x_det[planeNum] = x;
                                Pos_y_det[planeNum] = y;
                                Pos_z_det[planeNum] = z;
                                T_sili[planeNum] = tiempo;
                                double Econv = Digital(edep);
                                E_det[planeNum] += Econv;
                                Signal1[stripNum] += Econv;
                        }
                        else if (planeNum == 2)
                        {
                                Pos_x_det[planeNum] = x;
                                Pos_y_det[planeNum] = y;
                                Pos_z_det[planeNum] = z;
                                T_sili[planeNum] = tiempo;
                                double Econv = Digital(edep);
                                E_det[planeNum] += Econv;
                                Signal2[stripNum] += Econv;
                        }
                        else if (planeNum == 3)
                        {
                                Pos_x_det[planeNum] = x;
                                Pos_y_det[planeNum] = y;
                                Pos_z_det[planeNum] = z;
                                T_sili[planeNum] = tiempo;
                                double Econv = Digital(edep);
                                E_det[planeNum] += Econv;
                                Signal3[stripNum] += Econv;
                        }
                        else if (planeNum == 4)
                        {
                                Pos_x_det[planeNum] = x;
                                Pos_y_det[planeNum] = y;
                                Pos_z_det[planeNum] = z;
                                T_sili[planeNum] = tiempo;
                                double Econv = Digital(edep);
                                E_det[planeNum] += Econv;
                                Signal4[stripNum] += Econv;
                        }
                        else if (planeNum == 5)
                        {
                                Pos_x_det[planeNum] = x;
                                Pos_y_det[planeNum] = y;
                                Pos_z_det[planeNum] = z;
                                T_sili[planeNum] = tiempo;
                                double Econv = Digital(edep);
                                E_det[planeNum] += Econv;
                                Signal5[stripNum] += Econv;
                        }
                        else if (planeNum == 6)
                        {
                                Pos_x_det[planeNum] = x;
                                Pos_y_det[planeNum] = y;
                                Pos_z_det[planeNum] = z;
                                T_sili[planeNum] = tiempo;
                                double Econv = Digital(edep);
                                E_det[planeNum] += Econv;
                                Signal6[stripNum] += Econv;
                        }
                        else if (planeNum == 7)
                        {
                                Pos_x_det[planeNum] = x;
                                Pos_y_det[planeNum] = y;
                                Pos_z_det[planeNum] = z;
                                T_sili[planeNum] = tiempo;
                                double Econv = Digital(edep);
                                E_det[planeNum] += Econv;
                                Signal7[stripNum] += Econv;
                        }
                        else
                        {
                                G4cerr << "Hit Error: Plane number " << planeNum << " expected max value: 8" << G4endl;
                                continue;
                        }
                        // G4cout.precision(5);
                        // G4cout << "Strip: " << stripNum << " ||"
                        //        << " Detector " << planeNum << " ||"
                        //        << " x: " << x << " ||"
                        //        << " y: " << y << " ||"
                        //        << " Hits: " << nHits << G4endl;
                }
                rootTree->Fill();
        }
        else
        {
                // G4cerr << "Error: No hits collection passed to RootSaver" << G4endl;
        }
        TruthPosx = static_cast<Float_t>(primPos.x());
        TruthPosy = static_cast<Float_t>(primPos.y());
        TruthPosz = static_cast<Float_t>(primPos.z());
        //Measure angle of the beam in xz plane measured from z+ direction
        // -pi<Angle<=pi (positive when close to x positiove direction)
        Float_t sign_z = (primMom.z() >= 0) ? +1 : -1;
        Float_t sign_x = (primMom.x() >= 0) ? +1 : -1;
        TruthAngle_theta = (primMom.z() != 0) ? TMath::PiOver2() * sign_x * (1 - sign_z) + std::atan(primMom.x() / primMom.z())
                                              : sign_x * TMath::PiOver2(); //beam perpendicular to z
        TruthAngle_theta /= CLHEP::deg;

        //Measure angle of the beam in xy. Phi
        // -pi/2<Phi<=pi/2

        Float_t sign_y = (primMom.y() >= 0) ? +1 : -1;
        TruthAngle_phi = (primMom.x() != 0) ? sign_y * std::abs(atan(primMom.y() / primMom.x()))
                                            : sign_y * TMath::PiOver2(); //beam perpendicular to x
        TruthAngle_phi /= CLHEP::deg;
        //G4cout<<"x mom "<<primMom.x()<<" y mom "<<primMom.y()<<G4endl;
        //if(TruthAngle_phi>7 && E_dssd2 > 0)G4cout<<"atan "<<std::abs(atan( primMom.y()/primMom.x() ))<<" en degree "<<TruthAngle_phi<<G4endl;
        //if(TruthAngle_phi>9 && E_dssd2 > 2)G4cout<<"x mom "<<primMom.x()<<" y mom "<<primMom.y()<<" z mom "<< primMom.z()<<" pos x  "<<TruthPosx<<" pos y  "<<TruthPosy<<" pos z "<<TruthPosz<<" posy "<< Pos_y_dssd2<<G4endl;
}

//-------------------------------------------------------------------------//

double RootSaver::Digital(double Eraw)
{

        //TRandom3* myRandom = new TRandom3();
        double ion_pot = 3.6e-6; //eV
        double fanofactor = 0.1;
        double sigmaElnoise = 11.75e-3; //keV
        int Nelectrons = floor(Eraw / ion_pot);
        Nelectrons = gRandom->Gaus(Nelectrons, sqrt(Nelectrons * fanofactor));
        double energypoint = Nelectrons * ion_pot;
        energypoint = gRandom->Gaus(energypoint, sigmaElnoise); //electronic noise

        //delete myRandom;

        return energypoint;
}
