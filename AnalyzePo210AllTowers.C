
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

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TAxis.h"


using std::cout;
using std::cin;
using std::endl;
using namespace RooFit;

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

  const Int_t Number = 3;
   Double_t Red[Number]    = { 1.00, 1.00, 0.00};
   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
   Double_t Blue[Number]   = { 0.00, 1.00, 1.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   Int_t nb=50;
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);

   /*  
  Double_t stops[9] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000 };
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, stops, red, green, blue, 50);
  gStyle->SetNumberContours(50);
   */
}

int AnalyzePo210AllTowers() {

  set_plot_style();
  
  int nbins = 1000;
  int Po210start = 5000;
  int Po210end = 6000;

  //TFile * f1 = new TFile("/projects/cuore/data/ds3018_3021/Reduced_bkg_3018_3021.root");
  TFile * f1 = new TFile("/projects/cuore/data/TwoNu_DataRelease_Jan2019/Reduced_allDS_aggressive.root");
  TTree * t1 = (TTree*)f1->Get("outTree");
  
  TH1F* Po210 = new TH1F("Po210", "Po210", nbins, Po210start, Po210end);

  TH2F* Towers_Po210_hist = new TH2F("Towers_Po210_hist", "Towers_Po210_hist", 19, 1, 20, 13, 1, 14);
    
  TCut cutPo210 = "Energy < 5500 && Energy > 5250 && Included == 1";
  
  t1->Draw("Energy >> Po210", cutPo210, "goff");
  
  // Double Gaussian signal
RooRealVar x("x", "x", 0, 6000);
  RooRealVar mean("mean", "mean of gaussian", 5400, 0, 6000);
  RooRealVar sigma("sigma", "width of gaussian", 1, 0, 20);
  RooRealVar tail("tail", "tail of gaussian", 1, 0, 1000);
  RooRealVar norm("normal", "normalization", 0.5, 0 , 1);
  
  RooCBShape cball("cball", "crystal ball PDF", x, mean, sigma, tail, norm);
  
  RooRealVar mean2("mean2", "mean of secondary gaussian", 5500, 0, 6000);
  RooRealVar sigma2("sigma2", "sigma of secondary gaussian", 1, 0, 20);
  RooRealVar tail2("tail2", "tail of secondary gaussian", 1, 0, 1000);
  RooRealVar norm2("normal2", "secondary normalization", 0.5, 0 , 1);
  
  RooCBShape cball2("cball2", "crystal ball PDF", x, mean2, sigma2, tail2, norm2);

  RooRealVar mean3("mean3", "mean of third gaussian", 5310, 0, 6000);
  RooRealVar sigma3("sigma3", "sigma of third gaussian", 1, 0, 40);

  RooGaussian subgauss("subgauss", "subgaussian", x, mean3, sigma3);

    // Linear background
  RooRealVar a0("a0", "a0", 1, -10000, 10000);
  RooPolynomial p0("p0", "p0", x, RooArgList(a0));
  
  //Combined
  RooRealVar cball1frac("cball1frac", "fraction of main peak", 0.8, 0, 1);
  RooAddPdf doublepeak("doublepeak", "cball+cball2", RooArgList(cball, cball2), RooArgList(cball1frac));

  RooRealVar subgaussfrac("subgaussfrac", "fraction not in subgauss", 0.8, 0, 1);
  RooAddPdf triplepeak("triplepeak", "cball+cball2+gauss", RooArgList(doublepeak, subgauss), RooArgList(subgaussfrac));
  
  RooRealVar peakfrac("peakfrac", "fraction of gauss", 0.8, 0, 1);
  RooAddPdf peaklin("gausslin", "cball+p0", RooArgList(triplepeak, p0), RooArgList(peakfrac));


  x.setRange(5250, 5500);
  RooDataHist dataPo210("dataPo210", "Po-210", x, Po210);
  RooPlot * framePo210 = x.frame(Title("Po210"));

  mean.setVal(5340);
  mean.setRange(5300, 5400);
  sigma.setVal(7);

  mean2.setVal(5440);
  mean2.setRange(5430,5450);
  sigma2.setVal(7);

  mean3.setVal(5310);
  mean3.setRange(5300, 5335);
  sigma3.setVal(15);

  a0.setVal(0.1);
  a0.setRange(0,10);

  cball1frac.setVal(0.4);
  peakfrac.setVal(0.9);
  norm.setVal(0.8);
  norm2.setVal(0.9);
  subgaussfrac.setVal(0.7);
  tail.setVal(4);
  tail2.setVal(0.5);

  //a0.setVal(400);
  //a0.setRange(0,500);
  
  dataPo210.plotOn(framePo210);
  peaklin.fitTo(dataPo210);
  peaklin.plotOn(framePo210);
  peaklin.plotOn(framePo210, Components(p0), LineStyle(kDashed));
  peaklin.plotOn(framePo210, Components(subgauss), LineStyle(kDashed), LineColor(kRed));

  framePo210->Draw();
  framePo210->GetXaxis()->SetTitle("Energy");
  
  double flat = a0.getVal();
  double signal = peakfrac.getVal();
  
  //double signal = gaussfrac.getVal();
  
  //cout << flat << "\t" << linear << endl;

  //cout << Po210->Integral() <<endl;
  //cout << linear * 0.5 * (1550*1550 - 1350*1350) + flat * 200 << endl;
  cout << (1 - signal) * Po210->Integral() << endl;


  double Rate[13];
  double RateError[13];
  double Floor[13];
  double FloorError[13];
  double signal_floor[13];
  double signalerror_floor[13];
  double integral_floor[13];

  double FirstPeak_Rate[13];
  double SecondPeak_Rate[13];

  RooPlot * framePo210_floor[13];
  
  //TCanvas *c2 = new TCanvas();
  //c2->Divide(5,3);

  double livetimes[989]={0};
  // grab livetime info from txt file
  ifstream inFile;
  inFile.open("livetimes_2nu_2019.dat");
  std::string line;
  while (std::getline(inFile, line))
    {
      std::istringstream iss(line);
      int channel;
      double livetime;
      if (!(iss >> channel >> livetime)) break;
      else {
	cout << channel << "\t" << livetime << endl;
	livetimes[channel] += livetime;
      }
    }
  const int tower_floor = 13 * 19;
  double livetimes_towerfloor[tower_floor]={0};
  // turn into livetime by tower and floor
  for (int tower_i = 1; tower_i <=19; tower_i++)
    {
      for (int floor_i = 1; floor_i <=13; floor_i++)
	{
	// go over array by 13*tower + floor
	channel = 52 * (tower_i - 1) + 4 * (floor_i-1);
	cout << tower_i << "\t" << floor_i << "\t" << channel << endl;
	int towerfloor = 13 * (tower_i-1) + (floor_i-1);
	livetimes_towerfloor[towerfloor] = livetimes[channel] + livetimes[channel+1] + livetimes[channel+2] + livetimes[channel+3];
	cout << channel << "\t" << livetimes_towerfloor[towerfloor] << endl;
	}
    }
  // loop over all towers and floors
  
  for (int tower = 1; tower<=19; tower++)
    {
      std::ostringstream tower_ss;
      tower_ss << "Tower == " << tower;
      TCut cutTower = tower_ss.str().c_str();
  
      for (int floor = 1; floor <= 13; floor++)
	{
	  //c2->cd(floor);
	  std::ostringstream ss;
	  ss << "Floor == " << floor;
	  
	  cout << tower_ss.str() << endl;
	  cout  << ss.str() << endl;
      
	  if (Po210) delete Po210;
	  TH1F* Po210 = new TH1F("Po210", "Po210", nbins, Po210start, Po210end);
	  framePo210_floor[floor-1] = x.frame(Title(ss.str().c_str()));

	  TCut cutFloor = ss.str().c_str();
	  TCut fullCut = cutFloor && cutPo210 && cutTower;
	  t1->Draw("Energy >> Po210", fullCut, "goff");

	  RooDataHist dataPo210("dataPo210", "K-40", x, Po210);

	  mean.setVal(5340);
	  mean.setRange(5300, 5400);
	  sigma.setVal(7);

	  mean2.setVal(5440);
	  mean2.setRange(5430,5450);
	  sigma2.setVal(7);

	  mean3.setVal(5310);
	  mean3.setRange(5300, 5335);
	  sigma3.setVal(15);

	  a0.setVal(0.1);
	  a0.setRange(0,10);

	  cball1frac.setVal(0.4);
	  peakfrac.setVal(0.9);
	  norm.setVal(0.8);
	  norm2.setVal(0.9);
	  subgaussfrac.setVal(0.7);
	  tail.setVal(4);
	  tail2.setVal(0.5);
      
	  dataPo210.plotOn(framePo210_floor[floor-1], MarkerSize(0.2));
	  peaklin.fitTo(dataPo210);
	  peaklin.plotOn(framePo210_floor[floor-1]);
	  peaklin.plotOn(framePo210_floor[floor-1], Components(p0), LineStyle(kDashed));
	  peaklin.plotOn(framePo210_floor[floor-1], Components(subgauss), LineStyle(kDashed), LineColor(kRed));


	  framePo210_floor[floor-1]->Draw();
	  framePo210_floor[floor-1]->GetXaxis()->SetTitle("Energy");
	  //framePo210_floor[floor-1]->GetYaxis()->SetRangeUser(0,30);
      
	  double flat = a0.getVal();
	  signal_floor[floor-1] = peakfrac.getVal();
	  signalerror_floor[floor-1] = peakfrac.getError();

	  integral_floor[floor-1] = Po210->Integral();

	  towerfloor = 13 * (tower-1) + (floor-1);
	  Rate[floor-1] = Po210->Integral() / livetimes_towerfloor[towerfloor];
	  FirstPeak_Rate[floor-1] = (peakfrac.getVal() * (1 - subgaussfrac.getVal() * (1 - cball1frac.getVal())) * Po210->Integral()) / livetimes_towerfloor[towerfloor];
	  SecondPeak_Rate[floor-1] = (peakfrac.getVal() * (subgaussfrac.getVal() * (1 - cball1frac.getVal())) * Po210->Integral()) / livetimes_towerfloor[towerfloor];

	  Floor[floor-1] = floor;
	  RateError[floor-1] = sqrt(Po210->Integral()) / livetimes_towerfloor[towerfloor];
	  FloorError[floor-1] = 0;

	  Towers_Po210_hist->Fill(tower, floor, FirstPeak_Rate[floor-1] / SecondPeak_Rate[floor-1]);
	  Towers_Po210_hist->SetMinimum(0);
	  Towers_Po210_hist->SetMaximum(6);
	}
    }
  
    TCanvas * c4 = new TCanvas();
    c4->cd();
    framePo210_floor[0]->Draw();
    framePo210_floor[0]->GetXaxis()->SetTitle("Energy");
    framePo210_floor[0]->GetYaxis()->SetRangeUser(0,30);
    

    for (int j = 0; j < 13; j++)
      {
	cout << "floor: " << j+1 << " Events: " << integral_floor[j] << " Signal: " << signal_floor[j] << " Error: " << signalerror_floor[j] << endl; 
      }
  
    int n = 13;
    TGraphErrors * RatesByFloor = new TGraphErrors(n, Floor, Rate, FloorError, RateError);

    TCanvas * c3 = new TCanvas();
    c3->cd();
    RatesByFloor->Draw("APE");
    RatesByFloor->SetTitle("Po210 Rates by Floor");
    RatesByFloor->GetXaxis()->SetRangeUser(0,14);
    RatesByFloor->GetXaxis()->SetNdivisions(210,kTRUE);
    RatesByFloor->GetXaxis()->SetTitle("Floor");
    RatesByFloor->GetYaxis()->SetTitle("Events");
    RatesByFloor->GetYaxis()->SetTitleOffset(1.3);
    c3->SetGridy();

    TCanvas * c6 = new TCanvas();
    c6->cd();
    Towers_Po210_hist->Draw("COLZ");
    Towers_Po210_hist->SetTitle("Po210 Rates by Tower and Floor");
    Towers_Po210_hist->GetXaxis()->SetTitle("Tower [DAQ]");
    Towers_Po210_hist->GetYaxis()->SetTitle("Floor");

    Towers_Po210_hist->GetXaxis()->SetNdivisions(210, kTRUE);
    Towers_Po210_hist->GetYaxis()->SetNdivisions(210, kTRUE);

    /*
      c1->SaveAs("FloorAnalysis/Po210_sum.pdf");
      c1->SaveAs("FloorAnalysis/Po210_sum.C");
      c2->SaveAs("FloorAnalysis/Po210_floors.pdf");
      c2->SaveAs("FloorAnalysis/Po210_floors.C");
      c3->SaveAs("FloorAnalysis/Po210_rates.pdf");
      c3->SaveAs("FloorAnalysis/Po210_rates.C");
      c4->SaveAs("FloorAnalysis/Po210_floor0.pdf");
      c4->SaveAs("FloorAnalysis/Po210_floor0.C");
    */
}
