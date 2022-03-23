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
void reflectOverX(double NumBins, TH1F* h1);
double findMax(TH1F* h1, TH1F* h2, TH1F* h3);
double findMin(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, TH1F* h5, TH1F* h6);

int main(int argc, char* argv[]){
  if(argc != 2)
  {
    cout << "Incorrect Use of file. Please use ./ProgramName filePath" << endl;
  }
  TFile *file_in;

  file_in = new TFile(argv[1]);

  //Subtraction histograms and all that jazz
  TH1F *h1_E_cal_pimi_sub         = (TH1F*)file_in->Get("h1_E_cal_pimi_sub");
  TH1F *h1_E_cal_pipl_sub         = (TH1F*)file_in->Get("h1_E_cal_pipl_sub");
  TH1F *h1_E_tot_2p1pi_1p1pi_pimi = (TH1F*)file_in->Get("h1_E_tot_2p1pi_1p1pi_pimi");
  TH1F *h1_E_tot_1p2pi_pimi       = (TH1F*)file_in->Get("h1_E_tot_1p2pi_pimi");
  TH1F *h1_E_tot_2p2pi_pimi       = (TH1F*)file_in->Get("h1_E_tot_2p2pi_pimi");
  TH1F *h1_E_tot_1p3pi_pimi       = (TH1F*)file_in->Get("h1_E_tot_1p3pi_pimi");
  TH1F *h1_E_tot_3p1pi_pimi       = (TH1F*)file_in->Get("h1_E_tot_3p1pi_pimi");
//pipl histos
  TH1F *h1_E_tot_2p1pi_1p1pi_pipl = (TH1F*)file_in->Get("h1_E_tot_2p1pi_1p1pi_pipl");
  TH1F *h1_E_tot_1p2pi_pipl       = (TH1F*)file_in->Get("h1_E_tot_1p2pi_pipl");
  TH1F *h1_E_tot_2p2pi_pipl       = (TH1F*)file_in->Get("h1_E_tot_2p2pi_pipl");
  TH1F *h1_E_tot_1p3pi_pipl       = (TH1F*)file_in->Get("h1_E_tot_1p3pi_pipl");
  TH1F *h1_E_tot_3p1pi_pipl       = (TH1F*)file_in->Get("h1_E_tot_3p1pi_pipl");

  //1p1pi Histogram stuff
  TH1F *h1_E_cal_1p1pi_pimi_tot   = (TH1F*)file_in->Get("h1_E_cal_1p1pi_pimi_tot");
  TH1F *h1_E_cal_1p1pi_pimi_ONLY  = (TH1F*)file_in->Get("h1_E_cal_1p1pi_pimi_ONLY");
  TH1F *h1_E_cal_1p1pi_pipl_tot   = (TH1F*)file_in->Get("h1_E_cal_1p1pi_pipl_tot");
  TH1F *h1_E_cal_1p1pi_pipl_ONLY  = (TH1F*)file_in->Get("h1_E_cal_1p1pi_pipl_ONLY");

  //variable to be used later ;)
  int min = 0;

  gStyle->SetOptStat(0);

  divideByBinWidth(N_E_bins, h1_E_cal_pimi_sub);
  divideByBinWidth(N_E_bins, h1_E_cal_pipl_sub);
  divideByBinWidth(N_E_bins, h1_E_tot_2p1pi_1p1pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_1p2pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_2p2pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_1p3pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_3p1pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_2p1pi_1p1pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_tot_1p2pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_tot_2p2pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_tot_1p3pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_tot_3p1pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_tot);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_ONLY);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_tot);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_ONLY);

  /*reflectOverX(N_E_bins, h1_E_tot_2p1pi_1p1pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_1p2pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_2p2pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_1p3pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_3p1pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_2p1pi_1p1pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_1p2pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_2p2pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_1p3pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_3p1pi_pipl);*/

  TCanvas *c1 = new TCanvas("c1","",567,370);
  h1_E_cal_1p1pi_pimi_tot->SetAxisRange(0,5.5,"X");
  double c1Min = h1_E_cal_pimi_sub->GetMinimum();
  double c1Max = findMax(h1_E_cal_1p1pi_pimi_tot, h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub);
  h1_E_cal_1p1pi_pimi_tot->SetMinimum(c1Min*1.25);
  h1_E_cal_1p1pi_pimi_tot->SetMaximum(c1Max*1.25);
  h1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
  h1_E_cal_pimi_sub->SetLineColor(3); //green
  h1_E_cal_pimi_sub->SetLineWidth(2);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pimi_tot->Draw("HIST");
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  c1->Print("1p1pi_pimiStuff_4GeV.png");

  TCanvas *c2 = new TCanvas("c2","",567,370);

  h1_E_cal_1p1pi_pipl_tot->SetAxisRange(0,5.5,"X");
  double c2Min = h1_E_cal_pipl_sub->GetMinimum();
  double c2Max = findMax(h1_E_cal_1p1pi_pipl_tot, h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub);
  h1_E_cal_1p1pi_pipl_tot->SetMinimum(c2Min*1.25);
  h1_E_cal_1p1pi_pipl_tot->SetMaximum(c2Max*1.25);
  h1_E_cal_1p1pi_pipl_tot->SetLineColor(1);
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2);
  h1_E_cal_pipl_sub->SetLineColor(3);
  h1_E_cal_pipl_sub->SetLineWidth(2);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pipl_tot->Draw("HIST");
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  c2->Print("1p1pi_piplStuff_4GeV.png");


  TCanvas *c3 = new TCanvas("c3","",567,370);

  h1_E_cal_1p1pi_pimi_tot->SetAxisRange(0,5.5,"X");
  //Make the minimum prettier :-)
  double c3Min = findMin(h1_E_cal_pimi_sub, h1_E_tot_2p1pi_1p1pi_pimi, h1_E_tot_1p2pi_pimi, h1_E_tot_2p2pi_pimi, h1_E_tot_1p3pi_pimi, h1_E_tot_3p1pi_pimi);
  double c3Max = findMax(h1_E_cal_1p1pi_pimi_tot, h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub);
  h1_E_cal_1p1pi_pimi_tot->SetMaximum(c3Max *1.25);
  h1_E_cal_1p1pi_pimi_tot->SetMinimum(c3Min *1.25);
  //Set colors and line styles for ALL 1p1pi, ONLY 1p1pi, and subtraction
  h1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
  h1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(2);
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
  h1_E_cal_pimi_sub->SetLineStyle(2);
  h1_E_cal_pimi_sub->SetLineColor(3); //green
  h1_E_cal_pimi_sub->SetLineWidth(2);
  //Plot members
  h1_E_tot_2p1pi_1p1pi_pimi->SetLineColor(4); //blue
  h1_E_tot_1p2pi_pimi->SetLineColor(40); //purple
  h1_E_tot_2p2pi_pimi->SetLineColor(6); //magenta
  h1_E_tot_1p3pi_pimi->SetLineColor(7); //cyan
  h1_E_tot_3p1pi_pimi->SetLineColor(8); //dirt green
  //Make the axes prettier
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleOffset(0.7);
  //Plot 1p1pi, 1p1pi AT LEAST, and sub
  h1_E_cal_1p1pi_pimi_tot->Draw("HIST");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
  //Plot members
  h1_E_tot_2p1pi_1p1pi_pimi->Draw("HIST SAME");
  h1_E_tot_1p2pi_pimi->Draw("HIST SAME");
  h1_E_tot_2p2pi_pimi->Draw("HIST SAME");
  h1_E_tot_1p3pi_pimi->Draw("HIST SAME");
  h1_E_tot_3p1pi_pimi->Draw("HIST SAME");
  //Create png
  c3->Print("pimiSubtraction_4GeV.png");

  TCanvas *c4 = new TCanvas("c4","",567,370);

  h1_E_cal_1p1pi_pipl_tot->SetAxisRange(0,5.5,"X");
  //Make the minimum prettier
  double c4Min = findMin(h1_E_cal_pipl_sub, h1_E_tot_2p1pi_1p1pi_pipl, h1_E_tot_1p2pi_pipl, h1_E_tot_2p2pi_pipl, h1_E_tot_1p3pi_pipl, h1_E_tot_3p1pi_pipl);
  double c4Max = findMax(h1_E_cal_1p1pi_pipl_tot, h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub);
  h1_E_cal_1p1pi_pipl_tot->SetMaximum(c4Max*1.25);
  h1_E_cal_1p1pi_pipl_tot->SetMinimum(c4Min*1.25);
  //Set colors and line styles for ALL 1p1pi, ONLY 1p1pi, and subtraction
  h1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //Black
  h1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(2);  //Dotted Line
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //Red
  h1_E_cal_pipl_sub->SetLineStyle(2);
  h1_E_cal_pipl_sub->SetLineColor(3); //green
  h1_E_cal_pipl_sub->SetLineWidth(2);
  //Plot members
  h1_E_tot_2p1pi_1p1pi_pipl->SetLineColor(4); //blue
  h1_E_tot_1p2pi_pipl->SetLineColor(40); //purple
  h1_E_tot_2p2pi_pipl->SetLineColor(6); //magenta
  h1_E_tot_1p3pi_pipl->SetLineColor(7); //cyan
  h1_E_tot_3p1pi_pipl->SetLineColor(8); //dirt green
  //Axes stuff
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleOffset(0.7);
  //Plot 1p1pi, 1p1pi AT LEAST, and sub
  h1_E_cal_1p1pi_pipl_tot->Draw("HIST");
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  //Plot members
  h1_E_tot_2p1pi_1p1pi_pipl->Draw("HIST SAME");
  h1_E_tot_1p2pi_pipl->Draw("HIST SAME");
  h1_E_tot_2p2pi_pipl->Draw("HIST SAME");
  h1_E_tot_1p3pi_pipl->Draw("HIST SAME");
  h1_E_tot_3p1pi_pipl->Draw("HIST SAME");
  //Create png
  c4->Print("piplSubtraction_4GeV.png");


  //Flip stuff over the x-axis such that we can see how much of etot the sub takes up
  reflectOverX(N_E_bins, h1_E_tot_2p1pi_1p1pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_1p2pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_2p2pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_1p3pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_3p1pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_2p1pi_1p1pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_1p2pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_2p2pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_1p3pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_3p1pi_pipl);

  // - - - 2p1pi Multiplicity
  TCanvas *c5 = new TCanvas("c5","",567,370);
  double c5Min = findMin(h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_2p1pi_1p1pi_pimi, h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_2p1pi_1p1pi_pimi);
  double c5Max = findMax(h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_2p1pi_1p1pi_pimi);
  h1_E_cal_1p1pi_pimi_ONLY->SetMinimum(c5Min *1.25);
  h1_E_cal_1p1pi_pimi_ONLY->SetMaximum(c5Max *1.1);
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2);
  h1_E_tot_2p1pi_1p1pi_pimi->SetLineColor(4);
  h1_E_cal_pimi_sub->SetLineColor(3);
  h1_E_cal_pimi_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST");
  h1_E_tot_2p1pi_1p1pi_pimi->Draw("HIST SAME");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  //create png
  c5->Print("2p1pi_pimi_Multiplicity_4GeV.png");

 //1p2pi pimi Multiplicity
  TCanvas *c6 = new TCanvas("c6","",567,370);
  double c6Min = findMin(h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_1p2pi_pimi, h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_1p2pi_pimi);
  h1_E_cal_1p1pi_pimi_ONLY->SetMinimum(c6Min *1.25);
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2);
  h1_E_tot_1p2pi_pimi->SetLineColor(40);
  h1_E_cal_pimi_sub->SetLineColor(3);
  h1_E_cal_pimi_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST");
  h1_E_tot_1p2pi_pimi->Draw("HIST SAME");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  //create png
  c6->Print("1p2pi_pimi_Multiplicity_4GeV.png");

  //2p2pi pimi Multiplicity
  TCanvas *c7 = new TCanvas("c7","",567,370);
  double c7Min = findMin(h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_2p2pi_pimi, h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_2p2pi_pimi);
  h1_E_cal_1p1pi_pimi_ONLY->SetMinimum(c7Min *1.25);
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2);
  h1_E_tot_2p2pi_pimi->SetLineColor(6);
  h1_E_cal_pimi_sub->SetLineColor(3);
  h1_E_cal_pimi_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST");
  h1_E_tot_2p2pi_pimi->Draw("HIST SAME");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  //create png
  c7->Print("2p2pi_pimi_Multiplicity_4GeV.png");

  TCanvas *c8 = new TCanvas("c8","",567,370);
  double c8Min = findMin(h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_1p3pi_pimi, h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_1p3pi_pimi);
  h1_E_cal_1p1pi_pimi_ONLY->SetMinimum(c8Min *1.25);
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2);
  h1_E_tot_1p3pi_pimi->SetLineColor(7);
  h1_E_cal_pimi_sub->SetLineColor(3);
  h1_E_cal_pimi_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST");
  h1_E_tot_1p3pi_pimi->Draw("HIST SAME");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  //create png
  c8->Print("1p3pi_pimi_Multiplicity_4GeV.png");

  TCanvas *c9 = new TCanvas("c9","",567,370);
  double c9Min = findMin(h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_3p1pi_pimi, h1_E_cal_1p1pi_pimi_ONLY, h1_E_cal_pimi_sub, h1_E_tot_3p1pi_pimi);
  h1_E_cal_1p1pi_pimi_ONLY->SetMinimum(c9Min *1.25);
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2);
  h1_E_tot_3p1pi_pimi->SetLineColor(8);
  h1_E_cal_pimi_sub->SetLineColor(3);
  h1_E_cal_pimi_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST");
  h1_E_tot_3p1pi_pimi->Draw("HIST SAME");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  //create png
  c9->Print("3p1pi_pimi_Multiplicity_4GeV.png");

  TCanvas *c10 = new TCanvas("c10","",567,370);
  double c10Min = findMin(h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_2p1pi_1p1pi_pipl, h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_2p1pi_1p1pi_pipl);
  double c10Max = findMax(h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_2p1pi_1p1pi_pipl);
  h1_E_cal_1p1pi_pipl_ONLY->SetMinimum(c10Min *1.25);
  h1_E_cal_1p1pi_pipl_ONLY->SetMaximum(c10Max *1.1);
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2);
  h1_E_tot_2p1pi_1p1pi_pipl->SetLineColor(4);
  h1_E_cal_pipl_sub->SetLineColor(3);
  h1_E_cal_pipl_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST");
  h1_E_tot_2p1pi_1p1pi_pipl->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  //create png
  c10->Print("2p1pi_pipl_Multiplicity_4GeV.png");

  TCanvas *c11 = new TCanvas("c11","",567,370);
  double c11Min = findMin(h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_1p2pi_pipl, h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_1p2pi_pipl);
  h1_E_cal_1p1pi_pipl_ONLY->SetMinimum(c11Min *1.25);
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2);
  h1_E_tot_1p2pi_pipl->SetLineColor(40);
  h1_E_cal_pipl_sub->SetLineColor(3);
  h1_E_cal_pipl_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST");
  h1_E_tot_1p2pi_pipl->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  //create png
  c11->Print("1p2pi_pipl_Multiplicity_4GeV.png");

  TCanvas *c12 = new TCanvas("c12","",567,370);
  double c12Min = findMin(h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_2p2pi_pipl, h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_2p2pi_pipl);
  h1_E_cal_1p1pi_pipl_ONLY->SetMinimum(c12Min *1.25);
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2);
  h1_E_tot_2p2pi_pipl->SetLineColor(6);
  h1_E_cal_pipl_sub->SetLineColor(3);
  h1_E_cal_pipl_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST");
  h1_E_tot_2p2pi_pipl->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  //create png
  c12->Print("2p2pi_pipl_Multiplicity_4GeV.png");

  TCanvas *c13 = new TCanvas("c13","",567,370);
  double c13Min = findMin(h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_1p3pi_pipl,h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_1p3pi_pipl);
  h1_E_cal_1p1pi_pipl_ONLY->SetMinimum(c13Min *1.25);
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2);
  h1_E_tot_1p3pi_pipl->SetLineColor(7);
  h1_E_cal_pipl_sub->SetLineColor(3);
  h1_E_cal_pipl_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST");
  h1_E_tot_1p3pi_pipl->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  //create png
  c13->Print("1p3pi_pipl_Multiplicity_4GeV.png");

  //3p1pi pipl Multiplicity
  TCanvas *c14 = new TCanvas("c14","",567,370);
  double c14Min = findMin(h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_3p1pi_pipl, h1_E_cal_1p1pi_pipl_ONLY, h1_E_cal_pipl_sub, h1_E_tot_3p1pi_pipl);
  h1_E_cal_1p1pi_pipl_ONLY->SetMinimum(c14Min *1.25);
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2);
  h1_E_tot_3p1pi_pipl->SetLineColor(8);
  h1_E_cal_pipl_sub->SetLineColor(3);
  h1_E_cal_pipl_sub->SetLineWidth(2);
  //Axis stuff
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleOffset(0.7);
  //Draw things :)
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST");
  h1_E_tot_3p1pi_pipl->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  //create png
  c14->Print("3p1pi_pipl_Multiplicity_4GeV.png");

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

//Reflects everything over the x axis to make it more visible
void reflectOverX(double NumBins, TH1F* h1)
{
  NumBins = h1->GetNbinsX();
  double binCont;
  for(int i = 1; i<NumBins; i++)
  {
    binCont = h1->GetBinContent(i);
    h1->SetBinContent(i,-1*binCont);
  }
}

//Function that finds the maximum of 3 histograms ;)
double findMax(TH1F* h1, TH1F* h2, TH1F* h3)
{
  double m1 = h1->GetMaximum();
  double m2 = h2->GetMaximum();
  double m3 = h3->GetMaximum();

  double ret = max(m1,m2);
  ret = max(ret,m3);

  return ret;
}

//Function that finds the minimum of all the histograms
double findMin(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, TH1F* h5, TH1F* h6)
{
  double data[] = {h1->GetMinimum(), h2->GetMinimum(),h3->GetMinimum(),h4->GetMinimum(),h5->GetMinimum(),h6->GetMinimum()};
  int arraySize = 6;
  double min = 0;
  for(int i = 0; i<arraySize; i++)
  {
    if(data[i] < min)
    {
      min = data[i];
    }
  }
  return min;

}

/*
 - Use subprograms to make plots with more lines on each plot
  - Another program that reads in the histogram and then makes pretty pictures
  - Make sure to label the signs the contributions come in as (added/subtracted)
    - If subtracted -> Plot it as negative for ease of understanding
      - Will be able to see graphically that blank + blank2 sums to subtraction histogram
      -
*/
