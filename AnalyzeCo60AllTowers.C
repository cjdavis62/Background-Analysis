
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

int AnalyzeCo60AllTowers() {

  set_plot_style();
  
  int nbins = 1000;
  int Co60start = 1000;
  int Co60end = 2000;

  TFile * f1 = new TFile("/projects/cuore/data/ds3018_3021/Reduced_bkg_3018_3021.root");
  TTree * t1 = (TTree*)f1->Get("outTree");
  
  TH1F* Co60 = new TH1F("Co60", "Co60", nbins, Co60start, Co60end);

  TH2F* Towers_Co60_hist = new TH2F("Towers_Co60_hist", "Towers_Co60_hist", 19, 1, 20, 13, 1, 14);
    
  TCut cutCo60 = "Energy < 1450 && Energy > 1250";
  
  t1->Draw("Energy >> Co60", cutCo60, "goff");
  
  // Double Gaussian signal
  RooRealVar x("x", "x", 0, 4000);
  RooRealVar mean("mean", "mean of gaussian", 1400, 0, 4000);
  RooRealVar sigma("sigma", "width of gaussian", 1, 0.3, 10);

  RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

  RooRealVar mean2("mean2", "mean of secondary gaussian", 1300, 0, 4000);
  RooRealVar sigma2("sigma2", "sigma of secondary gaussian", 1, 0.3, 10);
  
  RooGaussian gauss2("gauss2", "secondary gaussian PDF", x, mean2, sigma2);

  // Linear background
  RooRealVar a0("a0", "a0", 1, -10000, 10000);
  RooRealVar a1("a1", "a1", 1, -100, 100);
  RooPolynomial p1("p1", "p1", x, RooArgList(a0, a1), 0);

  //Combined
  //RooRealVar gaus1frac("g1frac", "fraction of main gaussian", 0.8, 0.6, 1);
  //RooAddPdf doublegauss("doublegaus", "gauss+gauss2", RooArgList(gauss, gauss2), RooArgList(gaus1frac));
  
  RooRealVar gaussfrac("gfrac", "fraction of gauss", 0.8, 0, 1);
  RooAddPdf gausslin("gausslin", "gauss+p1", RooArgList(gauss, p1), RooArgList(gaussfrac));


  x.setRange(1315, 1345);
  RooDataHist dataCo60("dataCo60", "K-40", x, Co60);
  RooPlot * frameCo60 = x.frame(Title("Co60"));

  mean.setVal(1330);
  mean.setRange(1320, 1340);
  sigma.setVal(3);

  //mean2.setVal(1120);
  //mean2.setRange(1110,1130);
  //sigma2.setVal(3);

  //a0.setVal(400);
  //a0.setRange(0,500);
  
  dataCo60.plotOn(frameCo60);
  gausslin.fitTo(dataCo60);
  gausslin.plotOn(frameCo60);
  gausslin.plotOn(frameCo60, Components(p1), LineStyle(kDashed));
  gausslin.plotOn(frameCo60, Components(gauss2), LineStyle(kDashed), LineColor(kRed));

  frameCo60->Draw();

  double linear = a1.getVal();
  double flat = a0.getVal();
  double signal = gaussfrac.getVal();
  
  //double signal = gaussfrac.getVal();
  
  //cout << flat << "\t" << linear << endl;

  //cout << Co60->Integral() <<endl;
  //cout << linear * 0.5 * (1550*1550 - 1350*1350) + flat * 200 << endl;
  cout << (1 - signal) * Co60->Integral() << endl;


  double Rate[13];
  double RateError[13];
  double Floor[13];
  double FloorError[13];
  double signal_floor[13];
  double signalerror_floor[13];
  double integral_floor[13];

  RooPlot * frameCo60_floor[13];
  
  //TCanvas *c2 = new TCanvas();
  //c2->Divide(5,3);

  double livetimes[988]={0};
  // grab livetime info from txt file
  ifstream inFile;
  inFile.open("exposures.txt");
  std::string line;
  while (std::getline(inFile, line))
    {
      std::istringstream iss(line);
      int channel, ds;
      double livetime;
      if (!(iss >> channel >> ds >> livetime)) break;
      else {
	cout << channel << "\t" << ds << "\t" << livetime << endl;
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
      
	  if (Co60) delete Co60;
	  TH1F* Co60 = new TH1F("Co60", "Co60", nbins, Co60start, Co60end);
	  frameCo60_floor[floor-1] = x.frame(Title(ss.str().c_str()));

	  TCut cutFloor = ss.str().c_str();
	  TCut fullCut = cutFloor && cutCo60 && cutTower;
	  t1->Draw("Energy >> Co60", fullCut, "goff");

	  RooDataHist dataCo60("dataCo60", "K-40", x, Co60);
      
	  mean.setVal(1330);
	  mean.setRange(1320, 1340);
	  sigma.setVal(3);

	  //mean2.setVal(1120);
	  //mean2.setRange(1110,1130);
	  //sigma2.setVal(3);

	  a0.setVal(1);
	  a1.setVal(1);

	  //gaus1frac.setVal(0.8);
	  gaussfrac.setVal(0.8);
      
	  dataCo60.plotOn(frameCo60_floor[floor-1], MarkerSize(0.2));
	  gausslin.fitTo(dataCo60);
	  gausslin.plotOn(frameCo60_floor[floor-1]);
	  gausslin.plotOn(frameCo60_floor[floor-1], Components(p1), LineStyle(kDashed));
	  //gausslin.plotOn(frameCo60_floor[floor-1], Components(gauss2), LineStyle(kDashed), LineColor(kRed));


	  frameCo60_floor[floor-1]->Draw();
	  frameCo60_floor[floor-1]->GetXaxis()->SetTitle("Energy");
	  frameCo60_floor[floor-1]->GetYaxis()->SetRangeUser(0,30);
      
	  double linear = a1.getVal();
	  double flat = a0.getVal();
	  signal_floor[floor-1] = gaussfrac.getVal();
	  signalerror_floor[floor-1] = gaussfrac.getError();

	  integral_floor[floor-1] = Co60->Integral();


	  towerfloor = 13 * (tower-1) + (floor-1);
	  Rate[floor-1] = Co60->Integral() / livetimes_towerfloor[towerfloor];
	  Floor[floor-1] = floor;
	  RateError[floor-1] = sqrt(Co60->Integral());
	  FloorError[floor-1] = 0;

	  Towers_Co60_hist->Fill(tower, floor, Rate[floor-1]);
	}
    }
  
    TCanvas * c4 = new TCanvas();
    c4->cd();
    frameCo60_floor[0]->Draw();
    frameCo60_floor[0]->GetXaxis()->SetTitle("Energy");
    frameCo60_floor[0]->GetYaxis()->SetRangeUser(0,30);
    

    for (int j = 0; j < 13; j++)
      {
	cout << "floor: " << j+1 << " Events: " << integral_floor[j] << " Signal: " << signal_floor[j] << " Error: " << signalerror_floor[j] << endl; 
      }
  
    int n = 13;
    TGraphErrors * RatesByFloor = new TGraphErrors(n, Floor, Rate, FloorError, RateError);

    TCanvas * c3 = new TCanvas();
    c3->cd();
    RatesByFloor->Draw("APE");
    RatesByFloor->SetTitle("Co60 Rates by Floor");
    RatesByFloor->GetXaxis()->SetRangeUser(0,14);
    RatesByFloor->GetXaxis()->SetNdivisions(210,kTRUE);
    RatesByFloor->GetXaxis()->SetTitle("Floor");
    RatesByFloor->GetYaxis()->SetTitle("Events");
    RatesByFloor->GetYaxis()->SetTitleOffset(1.3);
    c3->SetGridy();

    TCanvas * c6 = new TCanvas();
    c6->cd();
    Towers_Co60_hist->Draw("COLZ");
    Towers_Co60_hist->SetTitle("Co60 Rates by Tower and Floor");
    Towers_Co60_hist->GetXaxis()->SetTitle("Tower [DAQ]");
    Towers_Co60_hist->GetYaxis()->SetTitle("Floor");
  

    /*
      c1->SaveAs("FloorAnalysis/Co60_sum.pdf");
      c1->SaveAs("FloorAnalysis/Co60_sum.C");
      c2->SaveAs("FloorAnalysis/Co60_floors.pdf");
      c2->SaveAs("FloorAnalysis/Co60_floors.C");
      c3->SaveAs("FloorAnalysis/Co60_rates.pdf");
      c3->SaveAs("FloorAnalysis/Co60_rates.C");
      c4->SaveAs("FloorAnalysis/Co60_floor0.pdf");
      c4->SaveAs("FloorAnalysis/Co60_floor0.C");
    */
}
