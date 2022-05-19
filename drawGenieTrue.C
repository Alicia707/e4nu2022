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


  TH1F* h1_E_cal_1p1pi_pimi_TRUE     = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pimi_TRUE");
	TH1F* h1_E_cal_1p1pi_pipl_TRUE     = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pipl_TRUE");
  TH1F* h1_E_cal_1p1pi_pimi_TRUE_qel = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pimi_TRUE_qel");
	TH1F* h1_E_cal_1p1pi_pipl_TRUE_qel = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pipl_TRUE_qel");
	TH1F* h1_E_cal_1p1pi_pimi_TRUE_mec = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pimi_TRUE_mec");
	TH1F* h1_E_cal_1p1pi_pipl_TRUE_mec = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pipl_TRUE_mec");
	TH1F* h1_E_cal_1p1pi_pimi_TRUE_res = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pimi_TRUE_res");
	TH1F* h1_E_cal_1p1pi_pipl_TRUE_res = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pipl_TRUE_res");
	TH1F* h1_E_cal_1p1pi_pimi_TRUE_dis = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pimi_TRUE_dis");
	TH1F* h1_E_cal_1p1pi_pipl_TRUE_dis = (TH1F*) file_in-> Get("h1_E_cal_1p1pi_pipl_TRUE_dis");

  gStyle->SetOptStat(0);

  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_TRUE);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_TRUE);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_TRUE_qel);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_TRUE_qel);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_TRUE_mec);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_TRUE_mec);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_TRUE_res);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_TRUE_res);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_TRUE_dis);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_TRUE_dis);

  TCanvas *c1 = new TCanvas("c1","",567,370);
  h1_E_cal_1p1pi_pimi_TRUE->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_TRUE->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pimi_TRUE->Draw("HIST");
  c1->Print("1p1pi_pimi_TRUE.png");

  TCanvas *c2 = new TCanvas("c2","",567,370);
  h1_E_cal_1p1pi_pipl_TRUE->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_TRUE->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pipl_TRUE->Draw("HIST");
  c2->Print("1p1pi_pipl_TRUE.png");

  TCanvas *c3 = new TCanvas("c3","",567,370);
  h1_E_cal_1p1pi_pipl_TRUE_qel->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_TRUE_qel->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_qel->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_qel->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_qel->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pipl_TRUE_qel->Draw("HIST");
  c3->Print("1p1pi_pipl_TRUE_qel.png");

  TCanvas *c4 = new TCanvas("c4","",567,370);
  h1_E_cal_1p1pi_pipl_TRUE_mec->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_TRUE_mec->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_mec->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_mec->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_mec->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pipl_TRUE_mec->Draw("HIST");
  c4->Print("1p1pi_pipl_TRUE_mec.png");

  TCanvas *c5 = new TCanvas("c5","",567,370);
  h1_E_cal_1p1pi_pipl_TRUE_res->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_TRUE_res->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_res->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_res->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_res->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pipl_TRUE_res->Draw("HIST");
  c5->Print("1p1pi_pipl_TRUE_res.png");

  TCanvas *c6 = new TCanvas("c6","",567,370);
  h1_E_cal_1p1pi_pipl_TRUE_dis->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_TRUE_dis->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_dis->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_dis->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_TRUE_dis->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pipl_TRUE_dis->Draw("HIST");
  c6->Print("1p1pi_pipl_TRUE_dis.png");

  TCanvas *c7 = new TCanvas("c7","",567,370);
  h1_E_cal_1p1pi_pimi_TRUE_qel->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_TRUE_qel->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_qel->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_qel->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_qel->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pimi_TRUE_qel->Draw("HIST");
  c7->Print("1p1pi_pimi_TRUE_qel.png");

  TCanvas *c8 = new TCanvas("c8","",567,370);
  h1_E_cal_1p1pi_pimi_TRUE_mec->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_TRUE_mec->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_mec->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_mec->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_mec->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pimi_TRUE_mec->Draw("HIST");
  c8->Print("1p1pi_pimi_TRUE_mec.png");

  TCanvas *c9 = new TCanvas("c9","",567,370);
  h1_E_cal_1p1pi_pimi_TRUE_res->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_TRUE_res->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_res->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_res->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_res->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pimi_TRUE_res->Draw("HIST");
  c9->Print("1p1pi_pimi_TRUE_res.png");

  TCanvas *c10 = new TCanvas("c10","",567,370);
  h1_E_cal_1p1pi_pimi_TRUE_dis->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_TRUE_dis->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_dis->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_dis->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_TRUE_dis->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pimi_TRUE_dis->Draw("HIST");
  c10->Print("1p1pi_pimi_TRUE_dis.png");
  return 0;
}

//Divides each bin by its bin width (dunno why its useful tbh)
void divideByBinWidth(double NumBins, TH1F* h1)
{
  NumBins = h1->GetNbinsX();
  for(int i = 1; i<=NumBins; i++){
 	 h1->SetBinContent(i,h1->GetBinContent(i)/h1->GetBinWidth(i));
  }
}
