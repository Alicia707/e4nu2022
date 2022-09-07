#define GENIE_ANALYSIS_C

#include "genie_analysis.h"
#include "Constants.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TH1D.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TMath.h>
#include <exception>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH3.h>
#include <TGraph.h>

#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

vector<double> CalculateCalKineVars(double ECal,TLorentzVector FSElectron) 
{

	vector<double> CalKineVars; CalKineVars.clear();

	TLorentzVector V4_beam_Cal(0,0,ECal,ECal);
	double nu_Cal = -(FSElectron-V4_beam_Cal).E();
	double Q2_Cal = -(FSElectron-V4_beam_Cal).Mag2();
	double x_bjk_Cal = Q2_Cal/(2*m_prot*nu_Cal);
	TVector3 V3_q_Cal = (V4_beam_Cal-FSElectron).Vect();
	double W_var_Cal = TMath::Sqrt((m_prot+nu_Cal)*(m_prot+nu_Cal)-V3_q_Cal*V3_q_Cal);


	CalKineVars.push_back(nu_Cal); // 0-th element: energy transfer using Ecal
	CalKineVars.push_back(Q2_Cal); // 1st element: Q2 using Ecal
	CalKineVars.push_back(x_bjk_Cal); // 2nd element: xB using Ecal
	CalKineVars.push_back(W_var_Cal); // 3rd element: invariant mass using Ecal

	return CalKineVars;

}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Loading all the constants from Constant.h (e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565,
// H3_bind_en, He4_bind_en, C12_bind_en, B_bind_en, He3_bind_en, D2_bind_en, Fe_bind_en, Mn_bind_en

void genie_analysis::Loop(Int_t choice) 
{
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	double histoweight = 1.;
	//Choice = 0 is for analysis of CLAS data while choice = 1 is for the analysis of GENIE Simulation
	if (choice != 1 && choice != 0)
	{
		std::cout << "This parameter value is not implemented in genie_analysis::Loop(). It should be either 0 or 1. The given value is " << choice << std::endl;
		std::exit(0);
	}

	std::map<std::string,double>bind_en;
	std::map<std::string,double>target_mass;
	std::map<std::string,double>residual_target_mass;
	std::map<std::string, double> Ecal_offset; //that might not be necessary for simulation data

	target_name = ftarget; //std string for target name
	en_beam["1161"]=1.161;
	en_beam["2261"]=2.261;
	en_beam["4461"]=4.461;

	en_beam_Ecal["1161"]=1.161;
	en_beam_Ecal["2261"]=2.261;
	en_beam_Ecal["4461"]=4.461;

	en_beam_Eqe["1161"]=1.161;
	en_beam_Eqe["2261"]=2.261;
	en_beam_Eqe["4461"]=4.461;

	if (fChain == 0) return;
	// 8.31.22 -> GetEntriesFast is returning garbage, unclear where defined, GetEntries appears to work
	//Long64_t nentries = fChain->GetEntriesFast();
	Long64_t const nentries = fChain->GetEntries();
	//nentries =8000000;

	//Resolutions for Smearing for GENIE simulation data
	double reso_p = 0.01; // smearing for the proton
	double reso_e = 0.005; // smearing for the electrons
	double reso_pi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)

	// Resolution defined above seems to be insufficient at 1.1 GeV -> tripled it for all particles
	if(fbeam_en == "1161") { reso_p = 3*reso_p; reso_e = 3*reso_e; reso_pi = 3*reso_pi; }

	double Wcut = 2; //cut for all beam energies < 2
	double Q2cut = 0; // cut for 1.1 GeV > 0.1, for 2.2 GeV > 0.4 and 4.4 GeV > 0.8

	const int n_slice=3; // Stick to the 3 slices
	const double pperp_min[n_slice]={0.,0.2,0.4};
	const double pperp_max[n_slice]={0.2,0.4,10.};

	TVector3 V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;

	TString E_acc_file;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) //1.1 GeV Configuration parameters and cuts
	{
		E_acc_file="1_161";
		Q2cut = 0.1;
	}


	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) //2.2 GeV Configuration parameters and cuts
	{
		E_acc_file="2_261";
		Q2cut = 0.4;
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) //4.4 GeV Configuration parameters and cuts
	{
		E_acc_file="4_461";
		Q2cut = 0.8;
	}

	//Further constants for binding energies and target masses
	Ecal_offset["3He"]=0.004;
	Ecal_offset["4He"]=0.005;
	Ecal_offset["C12"]=0.005;
	Ecal_offset["56Fe"]=0.011;

	bind_en["3He"] = He3_bind_en-D2_bind_en + Ecal_offset["3He"]; //the offset is used to shift the peak to be at 0
	bind_en["4He"] = He4_bind_en-H3_bind_en + Ecal_offset["4He"];
	bind_en["C12"] = C12_bind_en-B_bind_en	+ Ecal_offset["C12"];
	bind_en["56Fe"]= Fe_bind_en-Mn_bind_en	+ Ecal_offset["56Fe"];
	bind_en["CH2"] = C12_bind_en-B_bind_en;

	target_mass["3He"] = 2*m_prot+m_neut-He3_bind_en;
	target_mass["4He"] = 2*m_prot+2*m_neut-He4_bind_en;
	target_mass["C12"] = 6*m_prot+6*m_neut-C12_bind_en;
	target_mass["56Fe"]= 26*m_prot+30*m_neut-Fe_bind_en;
	target_mass["CH2"] = 6*m_prot+6*m_neut-C12_bind_en;

	residual_target_mass["3He"] = m_prot+m_neut-D2_bind_en;
	residual_target_mass["4He"] = m_prot+2*m_neut-H3_bind_en;
	residual_target_mass["C12"] = 5*m_prot+6*m_neut-B_bind_en;
	residual_target_mass["56Fe"]= 25*m_prot+30*m_neut-Mn_bind_en;
	residual_target_mass["CH2"] = 25*m_prot+30*m_neut-Mn_bind_en;

	gRandom = new TRandom3();
	gRandom->SetSeed(10);

	TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);
	TLorentzVector V4_target(0,0,0,target_mass[ftarget]);

	//Acceptance Maps

	TString WhichMap = "e2a_maps";
	TFile* file_acceptance;
	TFile* file_acceptance_p;
	TFile* file_acceptance_pip;

	TString Target = "12C";
	if (ftarget.c_str() == "3He") { Target = "3He"; }
	if (ftarget.c_str() == "4He") { Target = "4He"; }

	if (choice == 1)
	{ 
		//Only need acceptance maps for GENIE simulation data
		file_acceptance 	= TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+".root");
		file_acceptance_p 	= TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_p.root");
		file_acceptance_pip = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_pip.root");
	}

	double XSecScale = 1.;
	TFile* XSecFile = TFile::Open("/uboone/app/users/apapadop/R-3_0_6/mySplines/xsec_gxspl-FNALbig.root");

	TGraph* gr = NULL;

	if (XSecFile)
	{
		TDirectory* dir = (TDirectory*)(XSecFile->Get("nu_mu_C12"));
		gr = (TGraph*)dir->Get("tot_cc");
	}

	// ---------------------------------------------------------------------------------------------------------------
	//Output file definition

	TFile *file_out;
	if(choice == 1) { file_out = new TFile(Form("genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");}

	// ---------------------------------------------------------------------------------------------------------------

	//initialize Fiducial functions for EC limits
	fiducialcut->InitEClimits();
	std::cout << " Test InitEClimits Loop " << fiducialcut->up_lim1_ec->Eval(60) << std::endl;

	// -------------------------------------------------------------------------------------------------------

	//Binning for energy reconstruction histograms
	int n_bins;
	double *x_values;
	double *x_qe;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.)
	{
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=17;i++) { x_values[i]=0.4+i*0.04; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=20;i++) { x_values[i+18]=1.08+(i+1)*0.02; x_qe[i+18] = (x_values[i+18] - en_beam[fbeam_en]) / en_beam[fbeam_en]; }
	}

	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.)
	{
		n_bins=54;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=23;i++) { x_values[i]=i*0.09; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=30;i++) { x_values[i+24]=2.07+(i+1)*0.03; x_qe[i+24] = (x_values[i+24] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.)
	{
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=21;i++)	{ x_values[i]=i*0.2; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=16;i++)	{ x_values[i+22]=4.2+(i+1)*0.05; x_qe[i+22] = (x_values[i+22] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}

	// -------------------------------------------------------------------------------------------------------

	//Definitions of further Histograms

	int NBinsNu = 300, NBinsQ2 = 300;
	double MinNu = 0., MaxNu = 4.; double MinQ2 = 0., MaxQ2 = 6.;

	//Definitions of further Histograms
	TH1F *h1_E_rec_pipl = new TH1F("h1_E_rec_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_pimi = new TH1F("h1_E_rec_pimi","",n_bins,x_values);
	TH1F *h1_E_tot_pipl = new TH1F("h1_E_tot_pipl","",n_bins,x_values); //Total events
	TH1F *h1_E_tot_pimi = new TH1F("h1_E_tot_pimi","",n_bins,x_values);

	TH1F *h1_E_tot_2p1pi_1p1pi_pipl	= new TH1F("h1_E_tot_2p1pi_1p1pi_pipl","",n_bins,x_values);
	TH1F *h1_E_tot_2p1pi_1p1pi_pimi	= new TH1F("h1_E_tot_2p1pi_1p1pi_pimi","",n_bins,x_values);
	TH1F *h1_E_rec_2p1pi_1p1pi_pipl	= new TH1F("h1_E_rec_2p1pi_1p1pi_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_2p1pi_1p1pi_pimi= new TH1F("h1_E_rec_2p1pi_1p1pi_pimi","",n_bins,x_values);


	//Plots for E_Cal
	//11.3.21 EDIT: All events with at least 1p1pi
	TH1F *h1_E_cal_1p1pi_pimi_tot 	= new TH1F("h1_E_cal_1p1pi_pimi_tot","",n_bins,x_values);
	TH1F *h1_E_cal_1p1pi_pipl_tot 	= new TH1F("h1_E_cal_1p1pi_pipl_tot","",n_bins,x_values);
	TH1F *h1_E_cal_1p1pi_pimi_ONLY 	= new TH1F("h1_E_cal_1p1pi_pimi_ONLY","",n_bins,x_values);
	TH1F *h1_E_cal_1p1pi_pipl_ONLY 	= new TH1F("h1_E_cal_1p1pi_pipl_ONLY","",n_bins,x_values);

	//genie_truth analysis histograms and true undetected histograms
	TH1F* h1_E_cal_1p1pi_pimi_TRUE = new TH1F("h1_E_cal_1p1pi_pimi_TRUE","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_TRUE = new TH1F("h1_E_cal_1p1pi_pipl_TRUE","",n_bins, x_values);
	//true undetected histograms
	TH1F* h1_E_cal_1p1pi_pimi_2pNotDetected 		= new TH1F("h1_E_cal_1p1pi_pimi_2pNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pimi_2piNotDetected 		= new TH1F("h1_E_cal_1p1pi_pimi_2piNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pimi_2p2piNotDetected 		= new TH1F("h1_E_cal_1p1pi_pimi_2p2piNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_2pNotDetected 		= new TH1F("h1_E_cal_1p1pi_pipl_2pNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_2piNotDetected 		= new TH1F("h1_E_cal_1p1pi_pipl_2piNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_2p2piNotDetected 		= new TH1F("h1_E_cal_1p1pi_pipl_2p2piNotDetected", "", n_bins, x_values);

	// Vector containing kinematic variables using Ecal
	vector<double> CalKineVars{};
	// Weight to fill the plots mentioned above
	double LocalWeight;

	// Signal Event Counter -> 1e1p0pi events (everything lese is bkg)
	int SignalEvents 		= 0;
	int QESignalEvents 		= 0;
	int MECSignalEvents 	= 0;
	int RESSignalEvents 	= 0;
	int DISSignalEvents 	= 0;
	int OtherSignalEvents 	= 0;

	//----- Counts for events to match up -----//
	//----- 1p1pi -----//
	int C1p1piALL 	= 0;
	int C1p1piPIMI 	= 0;
	int C1p1piPIPL 	= 0;
	//----- 2p1pi -----//
	int C2p1piALL 	= 0;
	int C2p1piPIMI 	= 0;
	int C2p1piPIPL 	= 0;
	//----- 3p1pi -----//
	int C3p1piALL 	= 0;
	int C3p1piPIMI 	= 0;
	int C3p1piPIPL 	= 0;
	//----- 1p2pi -----//
	int C1p2piALL 	= 0;
	int C1p2piPIMI 	= 0;
	int C1p2piPIPL 	= 0;
	//----- 2p2pi -----//
	int C2p2piALL 	= 0;
	int C2p2piPIMI 	= 0;
	int C2p2piPIPL 	= 0;
	//----- 1p3pi -----//
	int C1p3piALL 	= 0;
	int C1p3piPIMI 	= 0;
	int C1p3piPIPL 	= 0;
	//----- Counter for Cut Check ----///EDIT
	int CpCut = 0; //prot cut
	int CeCut = 0; //e cut 
	int CpiCut = 0; //pi cut
	int CpFidCut = 0; 
	int CeFidCut = 0; 
	int CpiFidCut = 0; 
	int CpMomCut = 0; 
	int CeMomCut = 0; 
	int CpiMomCut = 0; 

	int EQESignalEventsWithin5Perc = 0, EQESignalEventsWithin5Perc_FirstSlice = 0, EQESignalEventsWithin5Perc_SecondSlice = 0, EQESignalEventsWithin5Perc_ThirdSlice = 0;
	int ECalSignalEventsWithin5Perc = 0, ECalSignalEventsWithin5Perc_FirstSlice = 0, ECalSignalEventsWithin5Perc_SecondSlice = 0, ECalSignalEventsWithin5Perc_ThirdSlice = 0;
	int PMiss_FirstSlice = 0, PMiss_SecondSlice = 0, PMiss_ThirdSlice = 0;

	/** Beginning of Event Loop **/
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
//	for(Long64_t jentry = 0; jentry < 1000000; jentry++)
	{
//	for (Long64_t jentry=0; jentry<Nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		//Read Entry
		int nb = GetEntry(jentry);
		if (nb == 0) { std::cout <<"Event loop: 0 byte read for entry " << jentry << ". Indicate failure in reading the file" <<	std::endl;}

		if (jentry%10000 == 0) {std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/fChain->GetEntries()*100. << " %"<< std::endl;}

		if( jentry%200000 == 0 )
		{
			gDirectory->Write("hist_Files", TObject::kOverwrite);
			//cout<<jentry<<endl;
		}

		if(jentry == 0)
		{ //first entry to initialize TorusCurrent, Fiducials and Subtraction classes

			//The TorusField has to be set before the Fiducialcut parameters are initialized
			if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. ) //1.1 GeV, we are not using the 1.1 GeV data with 1500 current field
			{
				 fTorusCurrent = 750;
			}
			else if( (en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) || (en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) ) //2.2 GeV	or 4.4 GeV
			{
				 fTorusCurrent = 2250;
			}
			else { std::cout << "genie_analysis::Loop(): fTorusCurrent could not be assigned" << std::endl;}

			fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
			fiducialcut->SetFiducialCutParameters(fbeam_en);
			std::cout << " EventLoop: Finished setting up fiducial cut class " << std::endl;
			rotation->InitSubtraction(fbeam_en, target_name, bind_en, N_tot, fiducialcut);
			std::cout << " EventLoop: Finished setting up rotation initialize " << std::endl;
		}

		//Resets q vector to (0,0,0)
		rotation->ResetQVector();

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		double SmearedPe;
		double SmearedEe;
		double e_acc_ratio = 1.;	//will be 1 for CLAS data

		// Outgoing e',	Uncorr and corrected are the same read from root file.
		//V4_el and V3_el will be changed by smearing for GENIE simulation data
		TLorentzVector 	V4_el(pxl,pyl,pzl,El);
		TLorentzVector 	V4_el_uncorr(pxl,pyl,pzl,El);
		TVector3 		V3_el(pxl,pyl,pzl);

		double el_momentum 	= V3_el.Mag();
		double el_theta 	= V3_el.Theta();
		double el_theta_deg;

		if (choice == 1)
		{ //smearing, fiducials and acceptance ratio for GENIE simulation data
			//Smearing of Electron Vector from Simulation
			SmearedPe = gRandom->Gaus(pl,reso_e*pl);
			SmearedEe = sqrt( SmearedPe*SmearedPe + e_mass * e_mass );
			//Ali Look here -> Maybe just set V4 and then make V3 from there?
			V3_el.SetXYZ(SmearedPe/pl * pxl,SmearedPe/pl * pyl,SmearedPe/pl * pzl);
			double phi_ElectronOut = V3_el.Phi(); //in Radians

			V3_el.SetPhi(phi_ElectronOut + TMath::Pi() ); // Vec.Phi() is between (-180,180), GENIE coordinate system flipped with respect to CLAS
			V4_el.SetPxPyPzE(V3_el.X(),V3_el.Y(),V3_el.Z(), SmearedEe);
			V4_el.SetPhi(V3_el.Phi());
			//Fiducial Cuts with the smeared values
			if ( !EFiducialCut(fbeam_en,V3_el) )
			{
				CeFidCut++;
				continue; // Electron theta & phi fiducial cuts
			} 

			phi_ElectronOut += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
			el_momentum = V3_el.Mag(); //Momentum after smearing
			el_theta = V3_el.Theta(); //Angle after smearing

			//acceptance_c takes phi in radians and here unmodified by 30 degree.
			//11 - might be particle type (electron)
			/*e_acc_ratio = acceptance_c(el_momentum, cos(el_theta), phi_ElectronOut, 11,file_acceptance);
			if ( fabs(e_acc_ratio) != e_acc_ratio ) 
			{
				CeCut++;
				continue; 
			}*/
		}

		// Explicit cuts on electron momentum
		if (fbeam_en=="1161" && el_momentum < 0.4) 
		{
			CeMomCut++;
			continue; 
		}
		if (fbeam_en=="2261" && el_momentum < 0.55) 
		{
			CeMomCut++;
			continue; 
		}
		if (fbeam_en=="4461" && el_momentum < 1.1) 
		{
			CeMomCut++;
			continue; 
		}

		//Definition as for data. It is also correct for GENIE simulation data since V3_el is rotated above by 180 degree in phi
		double el_phi_mod = V3_el.Phi()*TMath::RadToDeg()  + 30; //Add 30 degree for plotting and photon phi cut
		if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree


		//Calculated Mott Cross Section and Weights for Inclusive Histograms
		//Wght and e_acc_ratio is 1 for CLAS data
		//double Mott_cross_sec = ( pow(fine_struc_const,2.)*(cos(el_theta)+1))/(2*pow(El,2.)*pow((1-cos(el_theta)),2.));

		double reco_Q2 = -(V4_el-V4_beam).Mag2();
		double Q4 = reco_Q2 * reco_Q2;
		double Mott_cross_sec = (1./Q4) * XSecScale;  //not really the Mott_cross_sec

		// ---------------------------------------------------------------------------------------------------------------------

		// For neutrino scattering
		// switch to true for nu scattering to account for the difference in the propagator

		bool neutrino = false;

		if (neutrino && XSecFile)
		{
			XSecScale = gr->Eval(Ev);
			Mott_cross_sec = XSecScale;
		}

		// ---------------------------------------------------------------------------------------------------------------------

		double WeightIncl = wght*e_acc_ratio / Mott_cross_sec;//wght -> genie weight

		// Securing ourselves against infinities
		if ( fabs(WeightIncl) != WeightIncl ) { continue; }

		//Calculation of Reconstructed Energy from ELectron only
		//using the same value of single nucleon separation E Ecal and Eqe
		double m_delta = 1.232;
		double E_rec = (m_delta*m_delta-(m_prot-bind_en[ftarget])*(m_prot-bind_en[ftarget])+2*(m_prot-bind_en[ftarget])*V4_el.E())/(2*(m_prot-bind_en[ftarget]-V4_el.E()+V4_el.Rho()*cos(el_theta)));

		//Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
		double nu = -(V4_el-V4_beam).E();
		double x_bjk = reco_Q2/(2*m_prot*nu);

		// QE selection
		//if ( fabs(x_bjk - 1.) > 0.2) { continue; }

		// ---------------------------------------------------------------------------------------------------------------------

		TVector3 V3_q = (V4_beam-V4_el).Vect();
		double V3_q_theta_deg = V3_q.Theta()*TMath::RadToDeg();
		double V3_q_phi_deg = V3_q.Phi()*TMath::RadToDeg();//Changed both to same method 
		//Ensure that our degrees lie within the unit circle 
		if (V3_q_phi_deg > 360) 
		{ //adjusted to ensure that we actually get put in our selection 
			while(V3_q_phi_deg > 360)
			{
				V3_q_phi_deg = V3_q_phi_deg - 360.; 
			}
		}
		if (V3_q_phi_deg < 0) 
		{
			while(V3_q_phi_deg < 0)
			{
				V3_q_phi_deg = V3_q_phi_deg + 360.; 
			}
		}
		double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

		//converting theta to degrees
		el_theta_deg = el_theta*TMath::RadToDeg();

		//Cuts on Q2 and W, only keep events with Q2 > Q2cut and W < Wcut
		if ( reco_Q2 < Q2cut || W_var > Wcut) continue;

		//Set q vector for the following rotations for the subtraction procedure
		rotation->SetQVector(V3_q);
//		rotation->PrintQVector();

		//Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

		//Index variables for hadrons (p and pions)
		int index_p[20]; //index for each proton
		int index_pi[20]; //index for each pion
		int ind_pi_phot[20]; //index for pions and photons
		int index_pipl[20]; //index for each pi plus
		int index_pimi[20]; //index for each pi minus

		int charge_pi[20]; //Charge for the pions and photons
		//Smeared Momentum and Energy values for GENIE (simulation) data
		double Smeared_Pp[20]; //smeared momentum values for protons
		double Smeared_Ep[20]; //smeared energy values for protons
		double Smeared_Ppi[20]; //smeared momentum values for pions
		double Smeared_Epi[20]; //smeared energy values for pions

		//Number of hadrons
		//Changed notation from num_p, num_pi, etc. to num_p_det -> Look here: Carry thru with all other programs
		int num_p_det = 0;
		int num_pi_det = 0;
		int num_pi_phot = 0; //couting all pions and photons
		int num_pimi_det = 0;
		int num_pipl_det = 0;
		int num_pi_phot_nonrad = 0; //counting all pions and non-radiation photons
		int num_phot_rad = 0; //counting radiation photons
		//Index and number variables for neutral particles
		int ec_index_n[20];
		int ec_num_n = 0;
		bool ec_radstat_n[20];

		//Array initialize to -1 or false
		for (int i = 0; i < 20; i++)
		{
			index_p[i] = -1;   index_pi[i] = -1;   index_pipl[i] = -1;   index_pimi[i] = -1;   ind_pi_phot[i] = -1;
			ec_index_n[i] = -1;   ec_radstat_n[i] = false;
			charge_pi[i] = -2; //default number should be not a possible real charge
			Smeared_Pp[i]  = 0; Smeared_Ep[i]  = 0;  //default 0 momentum and energy after smearing
			Smeared_Ppi[i] = 0; Smeared_Epi[i] = 0;  //default 0 momentum and energy after smearing
		}

		const double phot_rad_cut = 40;
		const double phot_e_phidiffcut = 30; //electron - photon phi difference cut

		// Creating vectors to store id of particles in the array
		vector <int> ProtonID; vector <int> PiPlusID; vector <int> PiMinusID; vector <int> PhotonID;
		ProtonID.clear(); PiPlusID.clear(); PiMinusID.clear();  PhotonID.clear();

		//counters for genie_truth analysis
		//- - - - - Creation of Vectors - - - - -//
		vector<int> true_PionIndexCounter;
		vector<int> true_ProtonIndexCounter;
		vector<int> true_piplIndexCounter;
		vector<int> true_pimiIndexCounter;
		vector<int> true_ChargedPionIndexCounter;
		int true_ProtonCounter = 0;
		int true_ChargedPionCounter = 0;
		int true_PiPlusCounter = 0;
		int true_PiMinusCounter = 0;
		int true_GammaCounter = 0;

		//Loop for Hadrons
		for (int i = 0; i < nf; i++)
		{
			//Start of proton selection
			if (pdgf[i] == 2212  && pf[i] > 0.3)
			{
				if ( choice == 1)
				{ //GENIE data

					//genie_truth analysis
					true_ProtonCounter++;
					true_ProtonIndexCounter.push_back(i);

					//Smearing of proton
					double temp_smear_P = gRandom->Gaus(pf[i],reso_p*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_prot * m_prot );

					TVector3 V3_prot_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_prot = V3_prot_corr.Phi();
					V3_prot_corr.SetPhi(phi_prot + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS
					if (!PFiducialCut(fbeam_en, V3_prot_corr) ) 
					{
						CpFidCut++;
						continue; 
					} // Proton theta & phi fiducial cuts

					num_p_det = num_p_det + 1;
					index_p[num_p_det - 1] = i;
					ProtonID.push_back(i);
					Smeared_Pp[num_p_det - 1] = temp_smear_P;
					Smeared_Ep[num_p_det - 1] = temp_smear_E;
				}
				else
				{ //CLAS data does not need Fiducial Cut again
						num_p_det = num_p_det + 1;
						index_p[num_p_det - 1] = i;
						ProtonID.push_back(i);
				}
			}

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			if (pdgf[i] == -211  && pf[i] > 0.15)
			{ //PI minus
				if ( choice == 1)
				{ //GENIE data
					//genie_truth analysis
					true_PiMinusCounter++;
					true_ChargedPionCounter++;
					true_pimiIndexCounter.push_back(i);
					true_PionIndexCounter.push_back(i);
					true_ChargedPionIndexCounter.push_back(i);

					//Smearing of pi minus
					double temp_smear_P = gRandom->Gaus(pf[i],reso_pi*pf[i]);
					//8.29.22 Changed sqrt to Tmath::Sqrt EDIT
					double temp_smear_E = TMath::Sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

					TVector3 V3_pi_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_pion = V3_pi_corr.Phi();
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					// Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
					if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, -1) )
					{
						CpiFidCut++;
						continue;
					}
					num_pimi_det = num_pimi_det + 1;
					num_pi_det = num_pi_det + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pimi[num_pi_det - 1] = i;
					index_pi[num_pi_det - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiMinusID.push_back(i);
					charge_pi[num_pi_det - 1] = -1;
					Smeared_Ppi[num_pi_det - 1] = temp_smear_P;
					Smeared_Epi[num_pi_det - 1] = temp_smear_E;
				}
			}//end of pimi

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			if ( pdgf[i] == 211  && pf[i] > 0.15)
			{
				if (choice == 1)
				{ //GENIE data
					//genie_truth analysis
					true_PiPlusCounter++;
					true_ChargedPionCounter++;
					true_piplIndexCounter.push_back(i);
					true_PionIndexCounter.push_back(i);
					true_ChargedPionIndexCounter.push_back(i);

					//Smearing of pi plus
					double temp_smear_P = gRandom->Gaus(pf[i],reso_pi*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

					TVector3 V3_pi_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_pion = V3_pi_corr.Phi();
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					// Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
					if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, 1) )     
					{
						CpiFidCut++;
						continue; 
					}
					//Implemented better programming practices EDIT
					num_pipl_det 		+= 1;
					num_pi_det  		+= 1;
					num_pi_phot 		+= 1;
					num_pi_phot_nonrad 	+= 1;
					index_pipl[num_pi_det - 1] = i;
					index_pi[num_pi_det - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;	

					charge_pi[num_pi_det - 1] = 1;
					Smeared_Ppi[num_pi_det - 1] = temp_smear_P;
					Smeared_Epi[num_pi_det - 1] = temp_smear_E;
					PiPlusID.push_back(i);
				}
			} //end of pipl
			if (pdgf[i] == 22  && pf[i] > 0.3)
			{
				continue; //6.29.22 No photons!
			}//end of gamma
		} //end of hadron loop

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------

		//Skip event if there is at least one radiation photon
		if (num_phot_rad > 0) 
		{
		  continue;
		}
		// -------------------------------------------------------------------------------------------------------------------------------------------------------------
		// For GENIE samples, identify the interaction type

		int Interaction = -1;
		if (choice == 1) 
		{
			if (qel) { Interaction = 1; }
			if (mec) { Interaction = 2; }
			if (res) { Interaction = 3; }
			if (dis) { Interaction = 4; }

		}

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------

		//genie_truth analysis:
		if(num_p_det == 1 && num_pi_det == 1) //detected 1p and detected 1pi //do sep for det pipl and det pimi (LOOK HERE) -> globally
		{ //define E_cal up here
			double E_cal = Ef[index_p[0]] + El + Ef[index_pi[0]] - m_prot;
			if(true_ProtonCounter == 2 && true_ChargedPionCounter == 1) //2nd proton undetected
			{
				if(true_PiPlusCounter == 1 && true_PiMinusCounter == 0)
				{
					h1_E_cal_1p1pi_pipl_2pNotDetected->Fill(E_cal, histoweight);
				}
				else if(true_PiPlusCounter == 0 && true_PiMinusCounter == 1)
				{
					h1_E_cal_1p1pi_pimi_2pNotDetected->Fill(E_cal, histoweight);
				}
				else
					cout<<"This should not happen! 2p1pi undetect with extra pion! ln 1252" << endl;
			}
			else if(true_ProtonCounter == 1 && true_ChargedPionCounter == 2) //2nd pi undetected
			{
				//check if the first pion is positively charged
				if(charge_pi[0] > 0) //detected pion is pipl
				{
					if(true_PiPlusCounter == 2 || true_PiMinusCounter == 1) //extra charged pion is pipl or pimi
					{
						h1_E_cal_1p1pi_pipl_2piNotDetected->Fill(E_cal, histoweight);
					}
					else
					{
						cout << "This should not happen. Line 1270" << endl;
					}
				}
				else if(charge_pi[0] < 0) //detected pion is pimi
				{
					if(true_PiPlusCounter == 1 || true_PiMinusCounter == 2)//extra charged pion is pipl or PiMinusID
					{
						h1_E_cal_1p1pi_pimi_2piNotDetected->Fill(E_cal, histoweight);
					}
					else
					{
						cout << "This should not happen. Line 1282" << endl;
					}
				}
			}
			else if(true_ProtonCounter == 2 && true_ChargedPionCounter == 2)
			{
				if(charge_pi[0] > 0 ) //detected pion is pipl
				{
					if(true_PiPlusCounter == 2 || true_PiMinusCounter == 1) //extra pi is either pipl or pimi
					{
						h1_E_cal_1p1pi_pipl_2p2piNotDetected->Fill(E_cal, histoweight);
					}
					else
					{
						cout << "This should not happen. Ln 1302. " << endl;
					}
				}
				else if(charge_pi[0] < 0) //detected pi is pimi
				{
					if(true_PiPlusCounter == 1 || true_PiMinusCounter == 2)
					{
						h1_E_cal_1p1pi_pimi_2p2piNotDetected->Fill(E_cal, histoweight);
					}
					else
					{
						cout << "This should not happen. Ln 1313." << endl;
					}
				}
				else
				{
					cout << "Somehow the charge of the first pion is 0, this is not possible! Line 1318" << endl;
					cout << "The charge of the first pion is: " << charge_pi[0] << endl;
				}
			}
		}

		if(choice == 1 && true_PiPlusCounter == 1 && true_ProtonCounter == 1)
		{
			double E_cal = Ef[true_ProtonIndexCounter[0]] + El + Ef[true_piplIndexCounter[0]]-m_prot;
			if(true_ProtonCounter == 1 && true_ChargedPionCounter == 1 && true_PiPlusCounter == 1)
			{
				h1_E_cal_1p1pi_pipl_TRUE->Fill(E_cal, histoweight);
			}
		}

		//genie_truth analysis
		if(choice == 1 && true_ProtonCounter == 1 && true_PiMinusCounter == 1)
		{
			double E_cal = Ef[true_ProtonIndexCounter[0]] + El + Ef[true_pimiIndexCounter[0]]-m_prot;
			if(true_ProtonCounter == 1 && true_ChargedPionCounter == 1 && true_PiMinusCounter == 1)
			{
				h1_E_cal_1p1pi_pimi_TRUE->Fill(E_cal, histoweight);
			}
		}

		//Events with exactly 2 protons
		if(num_p_det == 2)
		{
			//8.29.22 Implemented good programming practices (arrays) EDIT
			//START EDIT
			TLorentzVector V4_2prot_uncorr[2]; 
			TLorentzVector V4_2prot_corr[2];
			double p_acc_ratio[2] = {1.}; 

			for(int i = 0; i < 2; i++)
			{
				V4_2prot_uncorr[i].SetPxPyPzE(pxf[index_p[i]], pyf[index_p[i]], pzf[index_p[i]], TMath::Sqrt(m_prot*m_prot+pf[index_p[i]]*pf[index_p[i]]));
			}
			//Calculate corrected Lorentz Vectors for each proton
			if (choice == 1)
			{ //GENIE data, fiducials are done in hadron loop
				double prot_phi[2]; 
				double prot_theta[2]; 
				double prot_mom_corr[2];
				bool prot_continue_flag[2] = {false};
				//Loop over two protons in hadron loop 
				for(int i = 0; i < 2; i++)
				{
					V4_2prot_corr[i].SetPxPyPzE(Smeared_Pp[i]/pf[index_p[i]] * pxf[index_p[i]],Smeared_Pp[i]/pf[index_p[i]] * pyf[index_p[i]],Smeared_Pp[i]/pf[index_p[i]] * pzf[index_p[i]], Smeared_Ep[i]);
					prot_phi[i] = (V4_2prot_corr[i].Vect()).Phi(); 
					prot_phi[i] += TMath::Pi(); 
					V4_2prot_corr[i].SetPhi(prot_phi[i]);

					prot_theta[i] = (V4_2prot_corr[i].Vect()).Theta(); 
					prot_mom_corr[i] = (V4_2prot_corr[i].Vect()).Mag(); //3Vector magnitude is diff from LorentzVector mag 

					/*p_acc_ratio[i] = acceptance_c(prot_mom_corr[i], cos(prot_theta[i]), prot_phi[i], 2212, file_acceptance_p);
					//If the acceptance ratio isn't positive, continue to next event
					if(fabs(p_acc_ratio[i]) != p_acc_ratio[i]) {prot_continue_flag[i] = true;}*/
				}
				if(prot_continue_flag[0] == true || prot_continue_flag[1] == true)
				{
					CpCut++;
					continue;
				}
			}

			//Total proton weight EDIT
			double weight_protons = p_acc_ratio[0] * p_acc_ratio[1];

			TVector3 V3_2prot_uncorr[2];
			TVector3 V3_2prot_corr[2];
			for(int i = 0; i < 2; i++)
			{
				V3_2prot_uncorr[i] = V4_2prot_uncorr[i].Vect();
				V3_2prot_corr[i] = V4_2prot_corr[i].Vect();
			}
			//8.29.22 END EDIT

			//---------------------------------- 2p 1pi   ----------------------------------------------
			//Const int can be placed somewhere up after if for 2 protons F.H. 05.09.19
			const int N_2prot=2;
			//Variable might/could be placed in a more local context F.H. 05.09.19
			double Ecal_2p1pi_to2p0pi[N_2prot]={0};
			double p_miss_perp_2p1pi_to2p0pi[N_2prot]={0};

			if (num_pi_det == 1)
			{
				C2p1piALL++;

				TVector3 V3_1pi_corr;
				TLorentzVector V4_1pi_corr;
				double pion_acc_ratio = 1;

				if (choice == 1)
				{ //GENIE data
					//8.29.22 Edited to reflect better programming practices EDIT START
					//pion_acc_ratio = 0;//reset to 0 just to be safe for genie data
					//8.29.22 Corrected V3 assignment according to library for consistency - Also changed energy to Smeared NRG EDIT 
					//V4_1pi_corr.SetPxPyPzE(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]], TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[0]]*pf[ind_pi_phot[0]]));
					V4_1pi_corr.SetPxPyPzE(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]], Smeared_Epi[0]);
					V3_1pi_corr = V4_1pi_corr.Vect();
					double phi_pion = V3_1pi_corr.Phi();
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
					V3_1pi_corr.SetPhi(phi_pion); // Vec.Phi() is between (-180,180)
					V4_1pi_corr.SetPhi(phi_pion); //Set new Phi value to reflect V3_1pi_corr's phi

					double pion_theta = V3_1pi_corr.Theta();
					double pion_mom_corr = V3_1pi_corr.Mag();

					/*if (charge_pi[0] == 1)
					{ //acceptance for pi plus
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) 
						{
							CpiCut++; 
							continue; 
						}
					}
					else if (charge_pi[0] == -1)
					{    //acceptance for pi minus. using electron acceptance map
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) 
						{
							CpiCut++;
							continue; 
						}
					}
					else if (charge_pi[0] == 0)
					{    //acceptance for neutral, setting to 1 for now F.H. 09/24/19
						pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 2proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue; }
					*/				
				}

				double P_2p1pito2p0pi[2] = {0};
				double P_2p1pito1p1pi[2] = {0};
				double P_2p1pito1p0pi[2] = {0};
				double Ptot[2] = {0};
				double E_tot_2p[2] = {0};
				double p_perp_tot_2p[2] = {0};
				rotation->prot2_pi1_rot_func(V4_2prot_uncorr, V4_el_uncorr, V4_1pi_corr, V3_q, charge_pi[0], E_tot_2p,p_perp_tot_2p, P_2p1pito1p1pi);
				//double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec;
				//Is this correct in the following loop? F.H. 09/01/19
				if(charge_pi[0] == 1) {C2p1piPIPL++;}
				else if(charge_pi[0] == -1) {C2p1piPIMI++;}
				for(int z=0; z < N_2prot; z++)
				{ //looping over two protons
					if(charge_pi[0] == 1) //directly test the charge 8.30.22 EDIT
					{
						//11.3.21 EDIT: Added histogram fill for total events with at LEAST 1p1pi
						h1_E_cal_1p1pi_pipl_tot->Fill(E_tot_2p[z], histoweight);
						//---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------
						h1_E_tot_2p1pi_1p1pi_pipl->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h1_E_rec_2p1pi_1p1pi_pipl->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
					}//Minor bug fix: charge_pi for pimi
					else if(charge_pi[0] == -1) //Directly test the charge 8.30.22 EDIT
					{
						//11.3.21 EDIT: Added 1p1pi total event counter.
						h1_E_cal_1p1pi_pimi_tot->Fill(E_tot_2p[z],histoweight);
						//---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------
						h1_E_tot_2p1pi_1p1pi_pimi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h1_E_rec_2p1pi_1p1pi_pimi->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
					}

				}//filling the histograms for 2protons
			}//1pi requirement
		} //2prot requirement
	} //end of event loop (jentry)

	gStyle->SetOptFit(1);
	TH1F* h1_E_cal_pimi_sub = (TH1F*) h1_E_tot_pimi->Clone("h1_E_cal_pimi_sub");
	//Creates new histogram filled with floats and makes it a clone of h1_Ecal_pimi
	//Quotation marks are title of histogram when it is formed -> Should be same as variable name
	h1_E_cal_pimi_sub->Add(h1_E_tot_2p1pi_1p1pi_pimi, 1);
	h1_E_cal_pimi_sub->Write("process1");


	TH1F* h1_E_cal_pipl_sub = (TH1F*) h1_E_tot_pipl->Clone("h1_E_cal_pipl_sub");

	//2.24.22 -> If the weights are NEGATIVE -> These must be ADDED to reflect proper SUBTRACTION
	h1_E_cal_pipl_sub->Add(h1_E_tot_2p1pi_1p1pi_pipl, 1);

	//Weights of +1 to add positively

	gDirectory->Write("hist_Files", TObject::kOverwrite);
	// skim_tree->AutoSave();

	// --------------------------------------------------------------------------------------------------------

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Total Number of Events: " << nentries << std::endl;
	std::cout << std::endl << "One Pion Events = " << std::endl;
	/*std::cout << std::endl << "1p1pi ALL: " << C1p1piALL << std::endl;
	std::cout << std::endl << "1p1pi PIPL: " << C1p1piPIPL << std::endl;
	std::cout << std::endl << "1p1pi PIMI: " << C1p1piPIMI << std::endl;
	std::cout << std::endl << "2p1pi ALL: " <<C2p1piALL << std::endl;
	std::cout << std::endl << "2p1pi PIPL: " <<C2p1piPIPL << std::endl;
	std::cout << std::endl << "2p1pi PIMI: " <<C2p1piPIMI << std::endl;
	std::cout << std::endl << "3p1pi ALL: " <<C3p1piALL << std::endl;
	std::cout << std::endl << "3p1pi PIPL: " <<C3p1piPIPL << std::endl;
	std::cout << std::endl << "3p1pi PIMI: " <<C3p1piPIMI << std::endl;
	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Two Pion Events = " << std::endl;
	std::cout << std::endl << "1p2pi ALL: " <<C1p2piALL << std::endl;
	std::cout << std::endl << "1p2pi PIPL: " <<C1p2piPIPL << std::endl;
	std::cout << std::endl << "1p2pi PIMI: " <<C1p2piPIMI << std::endl;
	std::cout << std::endl << "2p2pi ALL: " << C2p2piALL << std::endl;
	std::cout << std::endl << "2p2pi PIPL: " << C2p2piPIPL << std::endl;
	std::cout << std::endl << "2p2pi PIMI: " <<C2p2piPIMI << std::endl;
	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Three Pion Events = " << std::endl;
	std::cout << std::endl << "1p3pi ALL: " <<C1p3piALL << std::endl;
	std::cout << std::endl << "1p3pi PIPL: " << C1p3piPIPL << std::endl;
	std::cout << std::endl << "1p3pi PIMI: " << C1p3piPIMI << std::endl;*/
	std::cout << std::endl << "Proton Acc Cut: " << CpCut << std::endl; 
	std::cout << std::endl << "Pion Acc Cut: " << CpiCut << std::endl; 
	std::cout << std::endl << "Electron Acc Cut: " << CeCut << std::endl;
	std::cout << std::endl << "Proton Fiducial Cut: " << CpFidCut << std::endl;
	std::cout << std::endl << "Electron Fiducial Cut: " << CeFidCut << std::endl;
	std::cout << std::endl << "Pion Fiducial Cut: " << CpiFidCut << std::endl; 
	std::cout << std::endl << "Proton Momentum Cut: " << CpMomCut << std::endl; 
	std::cout << std::endl << "Electron Momentum Cut: " << CeMomCut << std::endl; 
	std::cout << std::endl << "Pion Momentum Cut: " << CpiMomCut << std::endl; 
	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Initial # Events = " << fChain->GetEntries() << std::endl;
	std::cout << std::endl << "1e1p0pi Signal # Events = " << SignalEvents << std::endl;
	std::cout << std::endl << "Passing Rate = " << int(double(SignalEvents) / double(fChain->GetEntries())*100.) << " \%"<< std::endl << std::endl;

	/*std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "PMiss Fraction 1st Slice = " << int(double(PMiss_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "PMiss Fraction 2nd Slice = " << int(double(PMiss_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "PMiss Fraction 3rd Slice = " << int(double(PMiss_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "# Events With ECal Within 5\% of ETrue = " << ECalSignalEventsWithin5Perc << std::endl;
	std::cout << std::endl << "ECal 5% Fraction = " << int(double(ECalSignalEventsWithin5Perc) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 1st Slice = " << int(double(ECalSignalEventsWithin5Perc_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 2nd Slice = " << int(double(ECalSignalEventsWithin5Perc_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 3rd Slice = " << int(double(ECalSignalEventsWithin5Perc_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "# Events With EQE Within 5\% of ETrue = " << EQESignalEventsWithin5Perc << std::endl;
	std::cout << std::endl << "EQE 5% Fraction = " << int(double(EQESignalEventsWithin5Perc) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 1st Slice = " << int(double(EQESignalEventsWithin5Perc_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 2nd Slice = " << int(double(EQESignalEventsWithin5Perc_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 3rd Slice = " << int(double(EQESignalEventsWithin5Perc_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
*/
	if (choice == 1) {

		/*std::cout << std::endl << "QE Fractional Contribution = " << int(double(QESignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "MEC Fractional Contribution = " << int(double(MECSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "RES Fractional Contribution = " << int(double(RESSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "DIS Fractional Contribution = " << int(double(DISSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
*/
	}

}

//End Loop function

// -------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------

double genie_analysis::acceptance_c(double p, double cost, double phi, int particle_id,TFile* file_acceptance)
{

	//Redefinition of the phi angle
	// because the acceptance maps are defined between (-30,330)

	// Check that phi is between (0,360)

	//int redef = -30;
	int redef = 0;

	TH3D * acc;
	TH3D * gen;

	acc = (TH3D*)file_acceptance->Get("Accepted Particles");
	gen = (TH3D*)file_acceptance->Get("Generated Particles");

	//map 330 till 360 to [-30:0] for the acceptance map histogram
	if(phi > (2*TMath::Pi() - TMath::Pi()/6.) ) { phi -= 2*TMath::Pi(); }
	//Find number of generated events

	double pbin_gen = gen->GetXaxis()->FindBin(p);
	double tbin_gen = gen->GetYaxis()->FindBin(cost);
	double phibin_gen = gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
	double num_gen = gen->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

	//Find number of accepted events

	double pbin_acc = acc->GetXaxis()->FindBin(p);
	double tbin_acc = acc->GetYaxis()->FindBin(cost);
	double phibin_acc = acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
	double num_acc = acc->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

	double acc_ratio = (double)num_acc / (double)num_gen;
	double acc_err = (double)sqrt(acc_ratio*(1-acc_ratio)) / (double)num_gen;


	return acc_ratio;

}
