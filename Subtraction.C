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
#include "Subtraction.h"

void Subtraction::prot1_pi2_rot_func(TVector3  V3prot, TVector3 V3pi[2], TLorentzVector V4prot, TLorentzVector V4pi[2], int q_pi[2], TLorentzVector V4_el, double Ecal[2], double p_miss_perp[2], double P_1p1pi[]){

    const int N_pi=2;
    double rotation_ang;
    TVector3 V3_rot_pi[2], V3_p_rot;
    bool status_pi[2]={true};

    double N_all = 0;
    double N_1p1pi[2]={0};

       for(int g=0; g<N_tot; g++){

         rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
         V3_p_rot= V3prot;

         V3_p_rot.Rotate(rotation_ang,V3q);

         for(int i=0;i<N_pi;i++){

  	        V3_rot_pi[i]=V3pi[i];
  	        V3_rot_pi[i].Rotate(rotation_ang,V3q);
  	        status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);

         }

         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] ) N_1p1pi[0]++;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] ) N_1p1pi[1]++;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] ) N_all++;

       }
         //----------------------1p2pi->1p1pi
         for(int h=0;h<N_pi;h++)
         {
           prot1_pi1_en_calc(V4prot, V4pi[h], q_pi[h], V4_el, &Ecal[h], &p_miss_perp[h]);
           if(N_all!=0)
            P_1p1pi[h] = -(N_1p1pi[h]/N_all);
           else
            P_1p1pi[h]=0; //N_all!=0 statement
         }
       }

void Subtraction::prot1_pi3_rot_func(TVector3  V3prot, TVector3 V3pi[3], TLorentzVector V4prot, TLorentzVector V4pi[3], int q_pi[3], TLorentzVector V4_el, double Ecal[3], double p_miss_perp[3], double P_tot[3]){
    const int N_pi=3;
    double P_1p3pito1p1pi[3] = {0};
    double rotation_ang;
    TVector3 V3_rot_pi[N_pi], V3_p_rot;
    bool status_pi[N_pi]={true};
    double N_all = 0;
    double N_1p1pi[3]={0},N_1p2pi[3]={0};

    for(int g=0; g<N_tot; g++){

         rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
         V3_p_rot= V3prot;

         V3_p_rot.Rotate(rotation_ang,V3q);

         for(int i=0;i<N_pi;i++){

  	        V3_rot_pi[i]=V3pi[i];
  	        V3_rot_pi[i].Rotate(rotation_ang,V3q);
  	        status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
         }

         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] && !status_pi[2]) N_1p1pi[0]=N_1p1pi[0]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] && !status_pi[2]) N_1p1pi[1]=N_1p1pi[1]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && !status_pi[1] && status_pi[2]) N_1p1pi[2]=N_1p1pi[2]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] && !status_pi[2]) N_1p2pi[0]=N_1p2pi[0]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] && status_pi[2]) N_1p2pi[1]=N_1p2pi[1]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] && status_pi[2]) N_1p2pi[2]=N_1p2pi[2]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] && status_pi[2])  N_all=N_all+1;
    }
    for(int i=0;i<3;i++)
    {
      if(N_all!=0)
      {
        P_1p3pito1p1pi[i] = -(N_1p1pi[i]/N_all);
      }
      else
      {
        P_1p3pito1p1pi[i] = 0;
      }
    }
    //--------------1p2pi->1p1pi------
    double P_1p2pito1p1pi[6] = {0};
    int count = 0;
    TVector3 V3pi2[2];
    TLorentzVector V4pi2[2];
    int q_pi2[2] = {0};
    double P_1p1pi[2] = {0};
    double Ecal2[2] = {0};
    double p_miss_perp2[2] = {0};
    if(N_all!=0)
    {
      for(int i=0; i<3; i++)
        {
          for(int j=0; j<3; j++)
          {
            if(i<j)
            {
              P_1p1pi[0] = P_1p1pi[1] = 0;
              V3pi2[0] = V3pi[i];
              V3pi2[1] = V3pi[j];
              V4pi2[0] = V4pi[i];
              V4pi2[1] = V4pi[j];
              q_pi2[0] = q_pi[i];
              q_pi2[1] = q_pi[j];
              prot1_pi2_rot_func(V3prot, V3pi2, V4prot, V4pi2, q_pi2, V4_el, Ecal2, p_miss_perp2, P_1p1pi);
              Ecal[i] = Ecal2[0];
              Ecal[j] = Ecal2[1];
              p_miss_perp[i] = p_miss_perp2[0];
              p_miss_perp[j] = p_miss_perp2[1];
              P_1p3pito1p1pi[i] += -P_1p1pi[0]*N_1p2pi[count]/N_all;
              P_1p3pito1p1pi[j] += -P_1p1pi[1]*N_1p2pi[count]/N_all;
              count = count+1;
            }
          }
        }
        //Ask about this
        P_tot[0] = P_1p3pito1p1pi[0] + P_1p3pito1p1pi[0] + P_1p3pito1p1pi[1];
        P_tot[1] = P_1p3pito1p1pi[1] + P_1p3pito1p1pi[0] + P_1p3pito1p1pi[2];
        P_tot[2] = P_1p3pito1p1pi[2] + P_1p3pito1p1pi[1] + P_1p3pito1p1pi[2];
      }
}


void Subtraction::prot2_pi1_rot_func(TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_1pi, TLorentzVector V4_2prot_corr[2], TLorentzVector V4_1pi, int q_pi, TLorentzVector V4_el, double Ecal[2], double p_miss_perp[2], double P_tot[2]){

    const int N_2prot=2;
    TVector3 V3_2p_rotated[2],V3_1pirot;
    bool pi1_stat=true;
    double N_all=0,N_1p_1pi[2]={0};
    double P_2pto1p[2]={0},N_2p_det=0;
    Float_t pimi_phimin, pimi_phimax, cphil, cphir;
    double rot_angle = 0;

    prot1_pi1_en_calc(V4_2prot_corr[0], V4_1pi, q_pi, V4_el, &Ecal[0], &p_miss_perp[0]);
    prot1_pi1_en_calc(V4_2prot_corr[1], V4_1pi, q_pi, V4_el, &Ecal[1], &p_miss_perp[1]);

       for(int g=0; g<N_tot; g++)
       {
           rot_angle=gRandom->Uniform(0,2*TMath::Pi());


           V3_2p_rotated[0]=V3_2prot_uncorr[0];
           V3_2p_rotated[1]=V3_2prot_uncorr[1];
           V3_2p_rotated[0].Rotate(rot_angle,V3q);
           V3_2p_rotated[1].Rotate(rot_angle,V3q);

           V3_1pirot=V3_1pi;
           V3_1pirot.Rotate(rot_angle,V3q);
           pi1_stat=Pi_phot_fid_united(fbeam_en, V3_1pirot, q_pi);


           if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_1p_1pi[0]++;
           if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_1p_1pi[1]++;
           if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_all++;
        }

  for(int h = 0; h < 2; h++)
     {
       if(N_all != 0)
          P_tot[h] = -N_1p_1pi[h]/N_all;
        else
          P_tot[h] = 0;
     }
  }


void Subtraction::prot2_pi2_rot_func(TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_2pi[2], TLorentzVector V4_2prot_corr[2], TLorentzVector V4_2pi[2], int q_pi[2], TLorentzVector V4_el, double Ecal[2][2], double p_miss_perp[2][2], double P_tot_2p[2][2])
{

    const int N_2prot=2,N_2pi=2;
    TVector3 V3_2p_rotated[N_2prot],V3_2pirot[N_2pi];
    bool pi2_stat[N_2pi]={true};
    double rot_angle;
    double N_2p_1pi[N_2pi]={0},N_1p_2pi[N_2prot]={0},N_all=0,N_1p_1pi[N_2prot][N_2pi]={0};
    double N_pidet=0,N_piundet=0;
    double P_2pto1p[N_2prot]={0},N_2p_det=0;
    double P_1p1pi[N_2pi]={0};
    double P_2p1pito1p1pi[2]={0},Ptot=0;
    double P_2p2pito1p1pi[N_2prot]={0},P_2p2pito1p2pi[N_2prot]={0},P_2p2pito2p1pi[N_2prot]={0};
    P_tot_2p[0][0] = P_tot_2p[0][1] = P_tot_2p[1][0] = P_tot_2p[1][1] = 0;

    for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());

       for(int k=0; k<N_2pi; k++){

         V3_2p_rotated[k]=V3_2prot_uncorr[k];
         V3_2p_rotated[k].Rotate(rot_angle,V3q);


  	     V3_2pirot[k]=V3_2pi[k];
  	     V3_2pirot[k].Rotate(rot_angle,V3q);
  	     pi2_stat[k]=Pi_phot_fid_united(fbeam_en, V3_2pirot[k], q_pi[k]);
       }

       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_2p_1pi[0]++;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_2p_1pi[1]++;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[0]++;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[1]++;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[0][0]++;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[0][1]++;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[1][0]++;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[1][1]++;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_all++;
     }
     if(N_all!=0)
     {
       double prob2p2pito1p1pi[2][2] = {0};
     //---------------------------------------------------2p2pi->1p2pi-------------------------------------------------------
        double P_tot[2] = {0};
        double Ecal2[2] = {0};
        double p_miss_perp2[2] = {0};
      for(int i=0;i<2;i++)
      {
        prot1_pi2_rot_func(V3_2prot_uncorr[i], V3_2pi, V4_2prot_corr[i], V4_2pi, q_pi, V4_el, Ecal2, p_miss_perp2, P_tot);
        Ecal[i][0] = Ecal2[0];
        Ecal[i][1] = Ecal2[1];
        p_miss_perp[i][0] = p_miss_perp2[0];
        p_miss_perp[i][1] = p_miss_perp2[1];
        //Changed below to reflect the rotation function being performed. From N_2p_1pi -> N_1p_2pi
        prob2p2pito1p1pi[i][0]= -(N_1p_2pi[i]/N_all)*P_tot[0];
        prob2p2pito1p1pi[i][1]= -(N_1p_2pi[i]/N_all)*P_tot[1];
      }
      //---------------------------------------------------2p2pi->2p1pi->1p1pi-------------------------------------------------------
      P_tot[0] = P_tot[1] = 0;
    for(int i=0;i<2;i++)
    {
      prot2_pi1_rot_func(V3_2prot_corr, V3_2prot_uncorr, V3_2pi[i], V4_2prot_corr, V4_2pi[i], q_pi[i], V4_el, Ecal2, p_miss_perp2, P_tot);
      //Below is redundant
      //Ecal[0][i] = Ecal2[0];
      //Ecal[1][i] = Ecal2[1];
      p_miss_perp[0][i] = p_miss_perp2[0];
      p_miss_perp[1][i] = p_miss_perp2[1];
      prob2p2pito1p1pi[0][i] = -(N_2p_1pi[i]/N_all)*P_tot[0];
      prob2p2pito1p1pi[1][i] = -(N_2p_1pi[i]/N_all)*P_tot[1];
      //Positive bc it is doubly subtracted ^
    }
    //----------------------------------------------------2p2pi->1p1pi---------------------------------------------------------
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
          prob2p2pito1p1pi[i][j] += N_1p_1pi[i][j]/N_all; //probability = 2stepWeight - 1stepWeight (prob = prob - N1...)
          P_tot_2p[i][j] = prob2p2pito1p1pi[i][j];
        }
    }
   // P_tot_2p = prob2p2pito1p1pi;
  }
  else
  {
    //else (N_all == 0)
    P_tot_2p[0][0] = P_tot_2p[0][1] = P_tot_2p[1][0] = P_tot_2p[1][1] = 0;
  }

}

void Subtraction::prot3_pi1_rot_func(TVector3 V3_3prot_corr[3],TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, TLorentzVector V4_3prot_corr[3], TLorentzVector V4_pi, int q_pi, TLorentzVector V4_el, double Ecal[3], double p_miss_perp[3], double P_tot_3p[3])
{
    const int N_3prot=3;
    TVector3 V3_3p_rotated[N_3prot],V3_pirot;
    bool pi_stat=true;
    double rot_angle;

    double N_all=0,N_1p1pi[N_3prot]={0},N_2p1pi[N_3prot]={0};
    double  P_3p1pito1p1pi[N_3prot]={0};
    double   N_pidet=0,N_piundet=0;
    TVector3 V3_2p_corr[N_3prot],V3_2p_uncorr[N_3prot];
    double P_2pto1p[2]={0},N_2p_det=0;
    int count=0;
    double N_p1[N_3prot]={0},N_p_three=0;
    double P_3pto2p[N_3prot][2]={0};
    double P_3p1pito2p1pi[N_3prot]={0};
    double P_2p1pito1p1pi[2]={0},Ptot=0;

       for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());
       for(int k=0; k<N_3prot; k++){

         V3_3p_rotated[k]=V3_3prot_uncorr[k];
         V3_3p_rotated[k].Rotate(rot_angle,V3q);
       }

       V3_pirot=V3_pi;
       V3_pirot.Rotate(rot_angle,V3q);
       pi_stat=Pi_phot_fid_united(fbeam_en, V3_pirot, q_pi);

       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[0]=N_1p1pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[1]=N_1p1pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[2]=N_1p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[0]=N_2p1pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[1]=N_2p1pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[2]=N_2p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_all=N_all+1;
     }
    if(N_all!=0)
    {
      for(int z=0;z<N_3prot;z++){

   //---------------------------------- 3p 1pi ->1p 1pi   ----------------------------------------------
          P_3p1pito1p1pi[z] = -(N_1p1pi[z]/N_all);

   //---------------------------------- 3p 1pi ->2p 1pi   ----------------------------------------------
   TVector3 V3_prot_corr[2];
   TVector3 V3_prot_uncorr[2];
   TLorentzVector V4_prot_corr[2];
   double Ecal2[2] = {0};
   double p_miss_perp2[2] = {0};
   for(int i=0;i<N_3prot;i++){       //looping through 2p combinations  out of 3p
     if(z!=i && z<i)
     {               // 3 pairs of 2proton combinations with z, i indexes(z<i)
          P_2p1pito1p1pi[0]=P_2p1pito1p1pi[1]=0;
          Ptot=0;
          V3_prot_corr[0] = V3_3prot_corr[z];
          V3_prot_corr[1] = V3_3prot_corr[i];
          V3_prot_uncorr[0] = V3_3prot_uncorr[z];
          V3_prot_uncorr[1] = V3_3prot_uncorr[i];
          V4_prot_corr[0] = V4_3prot_corr[z];
          V4_prot_corr[1] = V4_3prot_corr[z];
          prot2_pi1_rot_func(V3_prot_corr,V3_prot_uncorr,V3_pi,V4_prot_corr, V4_pi, q_pi, V4_el, Ecal2, p_miss_perp2, P_2p1pito1p1pi);
          Ecal[z] = Ecal2[0];
          Ecal[i] = Ecal2[1];
          p_miss_perp[z] = p_miss_perp2[0];
          p_miss_perp[i] = p_miss_perp2[1];
          P_3p1pito2p1pi[z] += -(N_2p1pi[count]/N_all)*(P_2p1pito1p1pi[0]); //Probability will automatically be + in case that it is 2p2pi->1p1pi
          P_3p1pito2p1pi[i] += -(N_2p1pi[count]/N_all)*(P_2p1pito1p1pi[1]);

  	      count=count+1;
  		  	}
  		  }
      }//looping through 3p
      P_tot_3p[0]=P_3p1pito2p1pi[0]+P_3p1pito1p1pi[0];
      P_tot_3p[1]=P_3p1pito2p1pi[1]+P_3p1pito1p1pi[1];
      P_tot_3p[2]=P_3p1pito2p1pi[2]+P_3p1pito1p1pi[2];

      }
      else
      {
        P_tot_3p[0]= P_tot_3p[1]=P_tot_3p[2]=0;
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
