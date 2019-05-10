
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TROOT.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TCut.h"
#include <iostream>
#include <sstream>
#include "THStack.h"
#include "TMath.h"
#include <stdio.h>
#include "TAxis.h"


using std::cout;
using std::cin;
using std::endl;

void set_plot_style()
{
  /*
    const Int_t NRGBs = 3;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    */

  /*
  const Int_t Number = 3;
   Double_t Red[Number]    = { 1.00, 1.00, 0.00};
   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
   Double_t Blue[Number]   = { 0.00, 1.00, 1.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   Int_t nb=50;
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
  */
  
  Double_t stops[9] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000 };
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, stops, red, green, blue, 50);
  gStyle->SetNumberContours(50);
  
}

int Analyze_Sm147() {

  set_plot_style();

  int nbins = 1000;
  int K40start = 1000;
  int K40end = 2000;

  TFile * f1 = new TFile("/projects/cuore/data/TwoNu_DataRelease_Jan2019/Reduced_allDS_aggressive.root");

  TTree * t1 = (TTree*)f1->Get("outTree");
  
  TH1F* Sm147_detectors_inner = new TH1F("Sm147_detectors_inner", "Sm147_detectors_inner", 989, 0, 989);
  TH1F* Sm147_detectors_outer = new TH1F("Sm147_detectors_outer", "Sm147_detectors_outer", 989, 0, 989);
  TH1F* Sm147_detectors_all = new TH1F("Sm147_detectors_all", "Sm147_detectors_all", 100, 2325, 2310);
  
  TCut cutSm147 = "Energy < 2328 && Energy > 2305 && Multiplicity == 1 && Included == 1";
  TCut cut_inner = "Layer == 0";
  TCut cut_outer = "Layer == 1";
  
  t1->Draw("Channel >> Sm147_detectors_inner", cutSm147 && cut_inner, "goff");
  t1->Draw("Channel >> Sm147_detectors_outer", cutSm147 && cut_outer, "goff");

  t1->Draw("Energy >> Sm147_detectors_all", cutSm147, "goff");

  
  TCanvas * c2 = new TCanvas();
  Sm147_detectors_outer->Draw();
  Sm147_detectors_inner->Draw("SAME");
  Sm147_detectors_inner->SetLineColor(kRed);

  TCanvas * c3 = new TCanvas();
  Sm147_detectors_all->Draw();
  
}
