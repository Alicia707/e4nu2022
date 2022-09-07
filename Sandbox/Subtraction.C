#ifndef SUBTRACTION_CXX
#define SUBTRACTION_CXX

#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TGraph.h>
#include <vector>
#include "Subtraction.h"

void Subtraction::prot2_pi1_rot_func(TLorentzVector V4_2prot_uncorr[2], TLorentzVector V4_el, TLorentzVector V4_1pi, TVector3 V3q, int q_pi, double Ecal[2], double p_miss_perp[2], double P_tot[2])
{
  //8.30.22 Removing arguments that don't even get used (V3_2prot_corr, )
  /*
    8.30.22 Changed corrected to uncorrected 
     - Removed V3 dependence -> We only want 4-vectors 
     - Implemented good programming practices aka consistent naming schema
  */
 //Check charge of pion matches what we want to fill :)
 const int num_prot=2;
 TVector3 V3_2p_rotated[2],V3_1pi_rotated;
 bool status_pi=true;
 bool status_prot[2] = {true};
 double N_all=0,N_1p_1pi[2]={0};
 double P_2pto1p[2]={0},N_2p_det=0;
 Float_t pimi_phimin, pimi_phimax, cphil, cphir;
 double rot_angle = 0;

 for(int g = 0; g < N_tot; g++)
 {
     rot_angle = gRandom->Uniform(0,2*TMath::Pi());

     for(int i = 0; i < num_prot; i++)
     {//get rid of corrections (FOR NOW)
       V3_2p_rotated[i] = V4_2prot_uncorr[i].Vect();
       V3_2p_rotated[i].Rotate(rot_angle,V3q);
       status_prot[i] = PFiducialCut(fbeam_en, V3_2p_rotated[i]);
     }

     V3_1pi_rotated = V4_1pi.Vect();
     //We rotate the particle around the momentum transfer vector (q)
     V3_1pi_rotated.Rotate(rot_angle,V3q);
     status_pi = Pi_phot_fid_united(fbeam_en, V3_1pi_rotated, q_pi);

     if( status_prot[0] && !status_prot[1] && status_pi) N_1p_1pi[0]++;
     if(!status_prot[0] && status_prot[1] && status_pi) N_1p_1pi[1]++;
     if( status_prot[0] && status_prot[1] && status_pi) N_all++;
  }

  for(int h = 0; h < 2; h++)
  {
     prot1_pi1_en_calc(V4_2prot_uncorr[h], V4_1pi, q_pi, V4_el, &Ecal[h], &p_miss_perp[h]);
     if(N_all > 0)
        P_tot[h] = -N_1p_1pi[h]/N_all;
      else
        P_tot[h] = 0;
  }
}

void Subtraction::prot1_pi1_en_calc(TLorentzVector V4prot, TLorentzVector V4pi, int q_pi, TLorentzVector V4_el, double *Ecal, double *p_miss_perp)
{
    double m_prot=0.9382720813;
    TVector3 V3_total = V4prot.Vect() + V4pi.Vect() + V4_el.Vect();
    *Ecal = V4_el.E() + V4prot.E() - m_prot + V4pi.E();
    *p_miss_perp = TMath::Sqrt(V3_total.Px()*V3_total.Px()+V3_total.Py()*V3_total.Py());
}
#endif
