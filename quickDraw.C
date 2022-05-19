#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>

using namespace std;

double N_E_bins;
void divideByBinWidth(double NumBins, TH1F* h1);

int main(int argc, char* argv[])
{
  if(argc != 2)
  {
    cout << "Incorrect Use of file. Please use ./ProgramName filePath" << endl;
  }
  TFile *file_in;

  file_in = new TFile(argv[1]);


  TH1F* h1_E_cal_pimi_sub = (TH1F*) file_in-> Get("h1_E_cal_pimi_sub");
	TH1F* h1_E_cal_pipl_sub = (TH1F*) file_in-> Get("h1_E_cal_pipl_sub");
  TH1F* h1_E_tot_2p2pi_pipl_1step = (TH1F*) file_in-> Get("h1_E_tot_2p2pi_pipl_1step");
  TH1F* h1_E_tot_2p2pi_pipl_2step = (TH1F*) file_in-> Get("h1_E_tot_2p2pi_pipl_2step");
  TH1F* h1_E_tot_2p2pi_pimi_1step = (TH1F*) file_in-> Get("h1_E_tot_2p2pi_pimi_1step");
  TH1F* h1_E_tot_2p2pi_pimi_2step = (TH1F*) file_in-> Get("h1_E_tot_2p2pi_pimi_2step");

  gStyle->SetOptStat(0);

  divideByBinWidth(N_E_bins, h1_E_cal_pimi_sub);
  divideByBinWidth(N_E_bins, h1_E_cal_pipl_sub);
  divideByBinWidth(N_E_bins, h1_E_tot_2p2pi_pipl_1step);
  divideByBinWidth(N_E_bins, h1_E_tot_2p2pi_pipl_2step);
  divideByBinWidth(N_E_bins, h1_E_tot_2p2pi_pimi_1step);
  divideByBinWidth(N_E_bins, h1_E_tot_2p2pi_pimi_2step);

  TCanvas *c1 = new TCanvas("c1","",567,370);
  h1_E_cal_pimi_sub->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_pimi_sub->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_pimi_sub->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_pimi_sub->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_pimi_sub->SetLineColor(3);
  h1_E_cal_pimi_sub->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_pimi_sub->Draw("HIST");
  c1->Print("h1_E_cal_pimi_sub.png");

  TCanvas *c2 = new TCanvas("c2","",567,370);
  h1_E_cal_pipl_sub->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_pipl_sub->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_pipl_sub->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_pipl_sub->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_pipl_sub->SetLineColor(3);
  h1_E_cal_pipl_sub->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_pipl_sub->Draw("HIST");
  c2->Print("h1_E_cal_pipl_sub.png");

  TCanvas *c3 = new TCanvas("c2","",567,370);
  h1_E_tot_2p2pi_pipl_1step->GetXaxis()->SetTitle("E_{cal}");
  h1_E_tot_2p2pi_pipl_1step->GetXaxis()->SetTitleSize(0.05);
  h1_E_tot_2p2pi_pipl_1step->GetXaxis()->SetLabelSize(0.05);
  h1_E_tot_2p2pi_pipl_1step->GetYaxis()->SetLabelSize(0.05);
  h1_E_tot_2p2pi_pipl_1step->SetLineColor(3);
  h1_E_tot_2p2pi_pipl_1step->GetXaxis()->SetTitleOffset(0.7);
  h1_E_tot_2p2pi_pipl_1step->Draw("HIST");
  c2->Print("h1_E_tot_2p2pi_pipl_1step.png");

  TCanvas *c4 = new TCanvas("c2","",567,370);
  h1_E_tot_2p2pi_pipl_2step->GetXaxis()->SetTitle("E_{cal}");
  h1_E_tot_2p2pi_pipl_2step->GetXaxis()->SetTitleSize(0.05);
  h1_E_tot_2p2pi_pipl_2step->GetXaxis()->SetLabelSize(0.05);
  h1_E_tot_2p2pi_pipl_2step->GetYaxis()->SetLabelSize(0.05);
  h1_E_tot_2p2pi_pipl_2step->SetLineColor(3);
  h1_E_tot_2p2pi_pipl_2step->GetXaxis()->SetTitleOffset(0.7);
  h1_E_tot_2p2pi_pipl_2step->Draw("HIST");
  c2->Print("h1_E_tot_2p2pi_pipl_2step.png");

  TCanvas *c5 = new TCanvas("c2","",567,370);
  h1_E_tot_2p2pi_pimi_1step->GetXaxis()->SetTitle("E_{cal}");
  h1_E_tot_2p2pi_pimi_1step->GetXaxis()->SetTitleSize(0.05);
  h1_E_tot_2p2pi_pimi_1step->GetXaxis()->SetLabelSize(0.05);
  h1_E_tot_2p2pi_pimi_1step->GetYaxis()->SetLabelSize(0.05);
  h1_E_tot_2p2pi_pimi_1step->SetLineColor(3);
  h1_E_tot_2p2pi_pimi_1step->GetXaxis()->SetTitleOffset(0.7);
  h1_E_tot_2p2pi_pimi_1step->Draw("HIST");
  c2->Print("h1_E_tot_2p2pi_pimi_1step.png");

  TCanvas *c6 = new TCanvas("c2","",567,370);
  h1_E_tot_2p2pi_pimi_2step->GetXaxis()->SetTitle("E_{cal}");
  h1_E_tot_2p2pi_pimi_2step->GetXaxis()->SetTitleSize(0.05);
  h1_E_tot_2p2pi_pimi_2step->GetXaxis()->SetLabelSize(0.05);
  h1_E_tot_2p2pi_pimi_2step->GetYaxis()->SetLabelSize(0.05);
  h1_E_tot_2p2pi_pimi_2step->SetLineColor(3);
  h1_E_tot_2p2pi_pimi_2step->GetXaxis()->SetTitleOffset(0.7);
  h1_E_tot_2p2pi_pimi_2step->Draw("HIST");
  c2->Print("h1_E_tot_2p2pi_pimi_2step.png");
  return 0;
}

//Divides each bin by its bin width in order to area normalize the plot
//  this is helpful for comparing two different plots 
void divideByBinWidth(double NumBins, TH1F* h1)
{
  NumBins = h1->GetNbinsX();
  for(int i = 1; i<=NumBins; i++){
 	 h1->SetBinContent(i,h1->GetBinContent(i)/h1->GetBinWidth(i));
  }
}
