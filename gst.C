#define gst_cxx
#include "gst.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TVector3.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGaxis.h>

#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "Constants.h"

using namespace std;

void gst::Loop() {

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;

	// --------------------------------------------------------------

	TString FileName = "/pnfs/genie/persistent/users/apapadop/e4v_FilteredSamples/FilteredGenie_56Fe_4461GeV_1p0pi_1p1pi.root";
	TFile* FilteredFile = new TFile(FileName,"recreate");
	TTree *maketreenew = fChain->GetTree()->CloneTree(0);

	// --------------------------------------------------------------

	int SelectedEvents = 0;

	// --------------------------------------------------------------

	TF1 *myElectronFit = new TF1("myElectronFit","[0]+[1]/x",0.,5.);
	myElectronFit->SetParameters(13.5,15);

	TF1* myPiMinusFit = new TF1("myPiMinusFit","(x<0.35)*(25.+7./TMath::Power(x,1.)) + (x>0.35)*(16.+10/TMath::Power(x,1.))",0,5.);
	TH1F *h1_E_cal_1p1pi_pimi_TRUE = new TH1F("h1_E_cal_1p1pi_pimi_TRUE", "", 0., 0.01);
	// --------------------------------------------------------------

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
//	for (Long64_t jentry=0; jentry<10000000;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		// --------------------------------------------------------------

		if (jentry%1000 == 0) {std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/fChain->GetEntries()*100. << " %"<< std::endl;}

		// --------------------------------------------------------------

		// Electron threshods for 4.4 GeV

		if (El < 1.1) { continue; }
		if (Q2 < 0.8) { continue; }
		if (W > 2) { continue; }

		TVector3 ElectronV3(pxl,pyl,pzl);
		double el_momentum = ElectronV3.Mag();
		double el_theta = ElectronV3.Theta();

		double theta_min = myElectronFit->Eval(el_momentum);
		if (el_theta*180./TMath::Pi() < theta_min) { continue; }

		// --------------------------------------------------------------

		int ProtonCounter = 0;
		int ChargedPionCounter = 0;
		int PiPlusCounter = 0;
		int PiMinusCounter = 0;
		int GammaCounter = 0;

		//Loop for Hadrons

		for (int i = 0; i < nf; i++) {

			double theta = TMath::ACos(cthf[i]) * 180./TMath::Pi();

			if (pdgf[i] == 2212  && pf[i] > 0.3 && theta > MinThetaProton) { ProtonCounter++; }
			if (pdgf[i] == 211  && pf[i] > 0.15 && theta > MinThetaPiPlus) { PiPlusCounter++; ChargedPionCounter++; }
			if (pdgf[i] == -211  && pf[i] > 0.15) {

				TVector3 PiMinusV3(pxf[i],pyf[i],pzf[i]);
				double piminus_momentum = PiMinusV3.Mag();
				double piminus_theta = PiMinusV3.Theta();

				double piminus_theta_min = myPiMinusFit->Eval(piminus_momentum);
				if (piminus_theta*180./TMath::Pi() > piminus_theta_min) {

					PiMinusCounter++;
					ChargedPionCounter++;

				}

			}

			if (pdgf[i] == 22  && pf[i] > 0.3 && theta > MinThetaGamma) { GammaCounter++; }


		}

		// --------------------------------------------------------------

                // Chosen topology

                if (
                    !((ProtonCounter == 1 && ChargedPionCounter == 1))

                ) { continue; }



		// --------------------------------------------------------------

		// 1p0pi + 2p0pi
//		if ( !( (ProtonCounter == 1 && ChargedPionCounter == 0) || (ProtonCounter == 2 && ChargedPionCounter == 0) ) ) { continue; }
//		if (GammaCounter != 0) { continue; }

		// --------------------------------------------------------------

		SelectedEvents++;
		maketreenew->Fill();

		// --------------------------------------------------------------


   } // End of the loop over the events

   std::cout << "Selected events = " << SelectedEvents << std::endl;
   std::cout << "Created file " << FileName << std::endl;

   maketreenew->Write();
   FilteredFile->Close();

} // End of the Loop()
