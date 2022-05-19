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
void fixPeak(TH1F* h1);
void reflectOverX(double NumBins, TH1F* h1);
void normalizeHisto(TH1F* h1);
double findMax(TH1F* h1, TH1F* h2, TH1F* h3);

int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    cout<< "Incorrect Use of file. Please use ./ProgramName target TargetEnergy first number" << endl;
  }

  TFile *data_file;
  TFile *genie_file;

  //arbitrary variable max to make the axes pretty
  double max = 10.;

  //char* dfileName = Form("/u/home/amand/RootWork/GenieAnalysis/data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  char* dfileName = Form("/mnt/c/Users/alici/Documents/Git/WorkingCode/PresentationStuff/drawFunctions/data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  data_file = new TFile(dfileName);

  //char* gfileName = Form("/u/home/amand/RootWork/GenieAnalysis/genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  char* gfileName = Form("/mnt/c/Users/alici/Documents/Git/WorkingCode/PresentationStuff/drawFunctions/genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  genie_file = new TFile(gfileName);

  //Pull in all of our histograms :)
  TH1F *dh1_E_cal_pimi_sub = (TH1F*)data_file->Get("h1_E_cal_pimi_sub");
  TH1F *gh1_E_cal_pimi_sub = (TH1F*)genie_file->Get("h1_E_cal_pimi_sub");
  TH1F *dh1_E_cal_pipl_sub = (TH1F*)data_file->Get("h1_E_cal_pipl_sub");
  TH1F *gh1_E_cal_pipl_sub = (TH1F*)genie_file->Get("h1_E_cal_pipl_sub");
  //1p1pi Histogram stuff (data)
  TH1F *dh1_E_cal_1p1pi_pimi_tot   = (TH1F*)data_file->Get("h1_E_cal_1p1pi_pimi_tot");
  TH1F *dh1_E_cal_1p1pi_pimi_ONLY  = (TH1F*)data_file->Get("h1_E_cal_1p1pi_pimi_ONLY");
  TH1F *dh1_E_cal_1p1pi_pipl_tot   = (TH1F*)data_file->Get("h1_E_cal_1p1pi_pipl_tot");
  TH1F *dh1_E_cal_1p1pi_pipl_ONLY  = (TH1F*)data_file->Get("h1_E_cal_1p1pi_pipl_ONLY");
  //genie
  TH1F *gh1_E_cal_1p1pi_pimi_tot   = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pimi_tot");
  TH1F *gh1_E_cal_1p1pi_pimi_ONLY  = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pimi_ONLY");
  TH1F *gh1_E_cal_1p1pi_pipl_tot   = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pipl_tot");
  TH1F *gh1_E_cal_1p1pi_pipl_ONLY  = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pipl_ONLY");
  //more genie
  TH1F *gh1_E_cal_1p1pi_pimi_TRUE = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pimi_TRUE");
  TH1F *gh1_E_cal_1p1pi_pipl_TRUE = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pipl_TRUE");

  //Divide by bin width for all histograms
  divideByBinWidth(N_E_bins, dh1_E_cal_pimi_sub);
  divideByBinWidth(N_E_bins, dh1_E_cal_pipl_sub);
  divideByBinWidth(N_E_bins, dh1_E_cal_1p1pi_pimi_tot);
  divideByBinWidth(N_E_bins, dh1_E_cal_1p1pi_pimi_ONLY);
  divideByBinWidth(N_E_bins, dh1_E_cal_1p1pi_pipl_tot);
  divideByBinWidth(N_E_bins, dh1_E_cal_1p1pi_pipl_ONLY);
  //genie stuff
  divideByBinWidth(N_E_bins, gh1_E_cal_pimi_sub);
  divideByBinWidth(N_E_bins, gh1_E_cal_pipl_sub);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pimi_tot);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pimi_ONLY);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pipl_tot);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pipl_ONLY);
  // True genie stuff
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pimi_TRUE);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pipl_TRUE);

  //Fix Bin Peak for All histograms
  fixPeak(dh1_E_cal_pimi_sub);
  fixPeak(dh1_E_cal_1p1pi_pimi_tot);
  fixPeak(dh1_E_cal_1p1pi_pimi_ONLY);
  //genie stuff
  fixPeak(gh1_E_cal_pimi_sub);
  fixPeak(gh1_E_cal_1p1pi_pimi_tot);
  fixPeak(gh1_E_cal_1p1pi_pimi_ONLY);
  fixPeak(gh1_E_cal_1p1pi_pimi_TRUE);
  fixPeak(gh1_E_cal_1p1pi_pipl_TRUE);

  //Normalize all of our histograms
  normalizeHisto(dh1_E_cal_pimi_sub);
  normalizeHisto(dh1_E_cal_pipl_sub);
  normalizeHisto(dh1_E_cal_1p1pi_pimi_tot);
  normalizeHisto(dh1_E_cal_1p1pi_pimi_ONLY);
  normalizeHisto(dh1_E_cal_1p1pi_pipl_tot);
  normalizeHisto(dh1_E_cal_1p1pi_pipl_ONLY);
  //genie stuff
  normalizeHisto(gh1_E_cal_pimi_sub);
  normalizeHisto(gh1_E_cal_pipl_sub);
  normalizeHisto(gh1_E_cal_1p1pi_pimi_tot);
  normalizeHisto(gh1_E_cal_1p1pi_pimi_ONLY);
  normalizeHisto(gh1_E_cal_1p1pi_pipl_tot);
  normalizeHisto(gh1_E_cal_1p1pi_pipl_ONLY);
  //true genie stuff
  normalizeHisto(gh1_E_cal_1p1pi_pimi_TRUE);
  normalizeHisto(gh1_E_cal_1p1pi_pipl_TRUE);

  gStyle->SetOptStat(0);

  if(strcmp(argv[3], "2")==0)
  {
    TCanvas *c1 = new TCanvas("c1","",567,370);
    //I will make the assumption that our maximum will be in genie plots since
    //there is more data to work with there
    max = findMax(gh1_E_cal_1p1pi_pimi_tot, gh1_E_cal_1p1pi_pimi_ONLY, gh1_E_cal_pimi_sub)*1.25;
    gh1_E_cal_1p1pi_pimi_tot->SetMaximum(max);

    gh1_E_cal_1p1pi_pimi_tot->SetAxisRange(0,3.5,"X");
    gh1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
    gh1_E_cal_1p1pi_pimi_tot->SetLineStyle(2);
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(2);
    gh1_E_cal_pimi_sub->SetLineColor(3); //green
    gh1_E_cal_pimi_sub->SetLineStyle(2);
    //data stuff
    dh1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
    dh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
    dh1_E_cal_pimi_sub->SetLineColor(3); //green
    //
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pimi_tot->Draw("HIST");
    gh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    gh1_E_cal_pimi_sub->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pimi_tot->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    dh1_E_cal_pimi_sub->Draw("HIST SAME");
    c1->Print("combined_1p1pi_pimiStuff_2GeV.png");

    TCanvas *c2 = new TCanvas("c2","",567,370);
    max = findMax(gh1_E_cal_1p1pi_pipl_tot, gh1_E_cal_1p1pi_pipl_ONLY, gh1_E_cal_pipl_sub)*1.25;
    gh1_E_cal_1p1pi_pipl_tot->SetMaximum(max);

    gh1_E_cal_1p1pi_pipl_tot->SetAxisRange(0,3.5,"X");
    gh1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //blue
    gh1_E_cal_1p1pi_pipl_tot->SetLineStyle(2);
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //pink
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(2);
    gh1_E_cal_pipl_sub->SetLineColor(3); //gross green
    gh1_E_cal_pipl_sub->SetLineStyle(2);
    //data stuff
    dh1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //black
    dh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //red
    //dh1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(2);
    dh1_E_cal_pipl_sub->SetLineColor(3); //green
    //dh1_E_cal_pipl_sub->SetLineStyle(3);
    //
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pipl_tot->Draw("HIST");
    gh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    gh1_E_cal_pipl_sub->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pipl_tot->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    dh1_E_cal_pipl_sub->Draw("HIST SAME");
    c2->Print("combined_1p1pi_piplStuff_2GeV.png");

    //True genie comparison
    TCanvas *c3 = new TCanvas("c3","",567,370);
    gh1_E_cal_1p1pi_pimi_ONLY->SetAxisRange(0,3.5,"X");
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(1);
    gh1_E_cal_1p1pi_pimi_TRUE->SetLineColor(2);
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(1);
    gh1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_ONLY->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pimi_TRUE->Draw("HIST");
    gh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    c3->Print("genie_TRUE_pimi.png");

    TCanvas *c4 = new TCanvas("c4","",567,370);
    gh1_E_cal_1p1pi_pipl_ONLY->SetAxisRange(0,3.5,"X");
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(1);
    gh1_E_cal_1p1pi_pipl_TRUE->SetLineColor(2);
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(1);
    gh1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_ONLY->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pipl_TRUE->Draw("HIST");
    gh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    c4->Print("genie_TRUE_pipl.png");

  }
  else if(strcmp(argv[3], "4")==0)
  {
    TCanvas *c1 = new TCanvas("c1","",567,370);
    max = findMax(gh1_E_cal_1p1pi_pimi_tot, gh1_E_cal_1p1pi_pimi_ONLY, gh1_E_cal_pimi_sub)*1.25;
    gh1_E_cal_1p1pi_pimi_tot->SetMaximum(max);

    gh1_E_cal_1p1pi_pimi_tot->SetAxisRange(0,5.5,"X");
    gh1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //blue
    gh1_E_cal_1p1pi_pimi_tot->SetLineStyle(2);
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //pink
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(2);
    gh1_E_cal_pimi_sub->SetLineColor(3); //gross green
    gh1_E_cal_pimi_sub->SetLineStyle(2);
    //data stuff
    dh1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
    dh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
    //dh1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(2);
    dh1_E_cal_pimi_sub->SetLineColor(3); //green
    //dh1_E_cal_pimi_sub->SetLineStyle(3);
    //
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pimi_tot->Draw("HIST");
    gh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    gh1_E_cal_pimi_sub->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pimi_tot->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    dh1_E_cal_pimi_sub->Draw("HIST SAME");
    c1->Print("combined_1p1pi_pimiStuff_4GeV.png");

    TCanvas *c2 = new TCanvas("c2","",567,370);
    max = findMax(gh1_E_cal_1p1pi_pipl_tot, gh1_E_cal_1p1pi_pipl_ONLY, gh1_E_cal_pipl_sub)*1.25;
    gh1_E_cal_1p1pi_pipl_tot->SetMaximum(max);

    gh1_E_cal_1p1pi_pipl_tot->SetAxisRange(0,5.5,"X");
    gh1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //blue
    gh1_E_cal_1p1pi_pipl_tot->SetLineStyle(2);
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //pink
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(2);
    gh1_E_cal_pipl_sub->SetLineColor(3); //gross green
    gh1_E_cal_pipl_sub->SetLineStyle(2);
    //data stuff
    dh1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //black
    dh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //red
    //dh1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(2);
    dh1_E_cal_pipl_sub->SetLineColor(3); //green
    //dh1_E_cal_pipl_sub->SetLineStyle(3);
    //
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pipl_tot->Draw("HIST");
    gh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    gh1_E_cal_pipl_sub->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pipl_tot->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    dh1_E_cal_pipl_sub->Draw("HIST SAME");
    c2->Print("combined_1p1pi_piplStuff_4GeV.png");

    //True genie comparison
    TCanvas *c3 = new TCanvas("c3","",567,370);
    gh1_E_cal_1p1pi_pimi_ONLY->SetAxisRange(0,5.5,"X");
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(1);
    gh1_E_cal_1p1pi_pimi_TRUE->SetLineColor(2);
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(1);
    gh1_E_cal_1p1pi_pimi_TRUE->SetLineColor(2);
    gh1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_ONLY->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_ONLY->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pimi_TRUE->Draw("HIST");
    gh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    c3->Print("genie_TRUE_pimi.png");

    TCanvas *c4 = new TCanvas("c4","",567,370);
    gh1_E_cal_1p1pi_pipl_ONLY->SetAxisRange(0,5.5,"X");
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(1);
    gh1_E_cal_1p1pi_pipl_TRUE->SetLineColor(2);
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(1);
    gh1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_ONLY->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_ONLY->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pipl_TRUE->Draw("HIST");
    gh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    c4->Print("genie_TRUE_pipl.png");

  }
  else{
    return 0;
  }
  return 0;
}

//Helper functions
//This divides by the bin width
void divideByBinWidth(double NumBins, TH1F* h1)
{
  NumBins = h1->GetNbinsX();
  for(int i = 1; i<=NumBins; i++){
 	 h1->SetBinContent(i,h1->GetBinContent(i)/h1->GetBinWidth(i));
  }
}
void fixPeak(TH1F*h1)
{
  int location = h1->GetMaximumBin();
  h1->SetBinContent(location, h1->GetBinContent(location)*.25);
}

//This function reflects over the x axis so it looks better
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
//Normalize Histogram function
void normalizeHisto(TH1F* h1)
{
  Double_t factor = 1.;
  h1->Scale(factor/h1->Integral());
}

double findMax(TH1F* h1, TH1F* h2, TH1F* h3)
{
  double m1 = h1->GetMaximum();
  double m2 = h2->GetMaximum();
  double m3 = h3->GetMaximum();

  double ret = max(m1,m2);
  ret = max(ret,m3);

  return ret;
}
