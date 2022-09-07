#ifndef SUBTRACTION_H
#define SUBTRACTION_H

#include <TVectorT.h>
#include <TVector3.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include "Fiducial.h"


struct Subtraction {
  //initialize Beam Energy String, Target Name String, Map of binding energies and Fiducial cuts
  //N_tot of rotations, q vector is set to (0,0,0)
  void InitSubtraction(std::string in_beam_en, std::string in_target_name, std::map<std::string,double> in_bind_en, int  in_nrot, Fiducial *in_fiducial) {
      fbeam_en = in_beam_en;
      target_name = in_target_name;
      bind_en = in_bind_en;
      fiducialcut = in_fiducial;
      N_tot = in_nrot;
      V3q.SetX(0);
      V3q.SetY(0);
      V3q.SetZ(0);
  }


  void prot2_pi1_rot_func(TLorentzVector V4_2prot_uncorr[2], TLorentzVector V4_el, TLorentzVector V4_1pi, TVector3 V3q, int q_pi, double Ecal[2], double p_miss_perp[2], double P_tot[2]);
  void prot1_pi1_en_calc(TLorentzVector V4prot, TLorentzVector V4pi, int q_pi, TLorentzVector V4_el, double *Ecal, double *p_miss_perp);

  void  SetQVector(TVector3 qin) {
    V3q.SetX(qin.X());
    V3q.SetY(qin.Y());
    V3q.SetZ(qin.Z());
  }

  void  ResetQVector() {
    V3q.SetX(0);
    V3q.SetY(0);
    V3q.SetZ(0);
  }

  void  PrintQVector() {
    std::cout << "Subtraction Class stored q vector: ( " << V3q.X() << " , " << V3q.Y() << " , " << V3q.Z() << " ) " << std::endl;
  }

  Bool_t EFiducialCut(std::string beam_en, TVector3 momentum) {
    return fiducialcut->EFiducialCut(beam_en, momentum);
  }
  Bool_t PFiducialCut(std::string beam_en, TVector3 momentum) {
    return fiducialcut->PFiducialCut(beam_en, momentum);
  }
  bool Pi_phot_fid_united(std::string beam_en, TVector3 V3_pi_phot, int q_pi_phot) {
    return fiducialcut->Pi_phot_fid_united(beam_en,V3_pi_phot, q_pi_phot);
  }


  Fiducial *fiducialcut;
  TVector3 V3q;
  std::string fbeam_en;
  std::string target_name;
  std::map<std::string,double> bind_en;
  int N_tot;



};

#endif
