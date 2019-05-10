
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


int AnalyzeCo60ByTower() {


 //get livetimes
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

  // get livetimes on each tower
  double livetime_tower[19] = {0};
  for (int channel = 0; channel < 988; channel++)
    {
      int tower_proxy = floor(channel/52.0);
      livetime_tower[tower_proxy] += livetimes[channel+1];
    }

  
  int nbins = 1000;
  int Co60start = 1000;
  int Co60end = 2000;

  //TFile * f1 = new TFile("/projects/cuore/data/ds3018_3021/Reduced_bkg_3018_3021.root");
  TFile * f1 = new TFile("/projects/cuore/data/TwoNu_DataRelease_Jan2019/Reduced_allDS_aggressive.root");

  TTree * t1 = (TTree*)f1->Get("outTree");
  
  TH1F* Co60 = new TH1F("Co60", "Co60", nbins, Co60start, Co60end);

    
  TCut cutCo60 = "Energy < 1450 && Energy > 1250 && Included == 1";
  
  t1->Draw("Energy >> Co60", cutCo60, "goff");

  /*
  TF1 * fit = new TF1("fit", "[0] * exp(-0.5 * ((x - [1]) / [2])^2) + [3] + x* [4]", -1, 3000);

  fit->SetParameter(0, 2000);
  fit->SetParameter(1, 1460);
  fit->SetParameter(2, 10);
  fit->SetParameter(3, 4000);
  fit->SetParameter(4, -2);
  
  Co60->Fit("fit", "", "", 1350, 1550);
  */
  
  // Double Gaussian signal
  RooRealVar x("x", "x", 0, 4000);
  RooRealVar mean("mean", "mean of gaussian", 1400, 0, 4000);
  RooRealVar sigma("sigma", "width of gaussian", 1, 0, 10);

  RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

  RooRealVar mean2("mean2", "mean of secondary gaussian", 1300, 0, 4000);
  RooRealVar sigma2("sigma2", "sigma of secondary gaussian", 1, 0, 10);
  
  RooGaussian gauss2("gauss2", "secondary gaussian PDF", x, mean2, sigma2);

  // Linear background
  RooRealVar a0("a0", "a0", 1, -10000, 10000);
  RooRealVar a1("a1", "a1", 1, -100, 100);
  RooPolynomial p1("p1", "p1", x, RooArgList(a0, a1), 0);

  //Combined
  RooRealVar gaus1frac("g1frac", "fraction of main gaussian", 0.8, 0, 1);
  RooAddPdf doublegauss("doublegaus", "gauss+gauss2", RooArgList(gauss, gauss2), RooArgList(gaus1frac));
  
  RooRealVar gaussfrac("gfrac", "fraction of gauss", 0.8, 0, 1);
  RooAddPdf gausslin("gausslin", "gauss+p1", RooArgList(doublegauss, p1), RooArgList(gaussfrac));


  x.setRange(1250, 1450);
  RooDataHist dataCo60("dataCo60", "Co-60", x, Co60);
  RooPlot * frameCo60 = x.frame(Title("Co60"));

  mean.setVal(1330);
  mean.setRange(1320, 1370);
  sigma.setVal(3);

  mean2.setVal(1320);
  mean2.setRange(1310,1330);
  sigma2.setVal(3);

  //a0.setVal(400);
  //a0.setRange(0,500);
  
  dataCo60.plotOn(frameCo60);
  gausslin.fitTo(dataCo60);
  gausslin.plotOn(frameCo60);
  gausslin.plotOn(frameCo60, Components(p1), LineStyle(kDashed));

  frameCo60->Draw();

  double linear = a1.getVal();
  double flat = a0.getVal();
  double signal = gaussfrac.getVal();
  
  //double signal = gaussfrac.getVal();
  
  //cout << flat << "\t" << linear << endl;

  //cout << Co60->Integral() <<endl;
  //cout << linear * 0.5 * (1550*1550 - 1350*1350) + flat * 200 << endl;
  cout << (1 - signal) * Co60->Integral() << endl;


  double Rate[19];
  double RateError[19];
  double Tower[19];
  double TowerError[19];
  double signal_tower[19];
  double signalerror_tower[19];
  double integral_tower[19];

  RooPlot * frameCo60_tower[19];

  TCanvas *c2 = new TCanvas();
  c2->Divide(5,4);
  
  for (int tower = 1; tower <= 19; tower++)
    {
      c2->cd(tower);

      std::ostringstream ss;
      ss << "Tower == " << tower;

      cout  << ss.str() << endl;
      
      if (Co60) delete Co60;
      TH1F* Co60 = new TH1F("Co60", "Co60", nbins, Co60start, Co60end);

      frameCo60_tower[tower-1] = x.frame(Title(ss.str().c_str()));


      mean.setVal(1330);
      mean.setRange(1320, 1370);
      sigma.setVal(3);

      mean2.setVal(1320);
      mean2.setRange(1310,1330);
      sigma2.setVal(3);

      a0.setVal(1);
      a1.setVal(1);

      gaus1frac.setVal(0.8);
      gaussfrac.setVal(0.8);
      
      TCut cutTower = ss.str().c_str();
      TCut fullCut = cutTower && cutCo60;
      t1->Draw("Energy >> Co60", fullCut, "goff");

      RooDataHist dataCo60("dataCo60", "Co-60", x, Co60);

      dataCo60.plotOn(frameCo60_tower[tower-1], MarkerSize(0.2));
      gausslin.fitTo(dataCo60);
      gausslin.plotOn(frameCo60_tower[tower-1]);
      gausslin.plotOn(frameCo60_tower[tower-1], Components(p1), LineStyle(kDashed));

      frameCo60_tower[tower-1]->Draw();
      frameCo60_tower[tower-1]->GetXaxis()->SetTitle("Energy");
      frameCo60_tower[tower-1]->GetYaxis()->SetRangeUser(0, 80);

      double linear = a1.getVal();
      double flat = a0.getVal();
      signal_tower[tower-1] = gaussfrac.getVal();
      signalerror_tower[tower-1] = gaussfrac.getError();

      integral_tower[tower-1] = Co60->Integral();
      
      Rate[tower-1] = signal_tower[tower-1] * Co60->Integral() / livetime_tower[tower-1];
      Tower[tower-1] = tower;
      RateError[tower-1] = signalerror_tower[tower-1] * Co60->Integral() / livetime_tower[tower-1];
      TowerError[tower-1] = 0;
    }
  
  TCanvas * c4 = new TCanvas();
  c4->cd();
  frameCo60_tower[0]->Draw();
  frameCo60_tower[0]->GetXaxis()->SetTitle("Energy");
  frameCo60_tower[0]->GetYaxis()->SetRangeUser(0,80);


  int n = 19;
  TGraphErrors * RatesByTower = new TGraphErrors(n, Tower, Rate, TowerError, RateError);

  TCanvas * c3 = new TCanvas();
  c3->cd();
  RatesByTower->Draw("APE");
  RatesByTower->SetTitle("Co60 Rates by Tower");
  RatesByTower->GetXaxis()->SetRangeUser(0,20);
  RatesByTower->GetXaxis()->SetNdivisions(210,kTRUE);
  RatesByTower->GetXaxis()->SetTitle("Tower");
  RatesByTower->GetYaxis()->SetTitle("Events");
  RatesByTower->GetYaxis()->SetTitleOffset(1.3);
  c3->SetGridy();

    //Save plots
  /*
  c1->SaveAs("TowerAnalysis/Co60_1330_sum.pdf");
  c1->SaveAs("TowerAnalysis/Co60_1330_sum.C");
  c2->SaveAs("TowerAnalysis/Co60_1330_towers.pdf");
  c2->SaveAs("TowerAnalysis/Co60_1330_towers.C");
  c3->SaveAs("TowerAnalysis/Co60_1330_rates.pdf");
  c3->SaveAs("TowerAnalysis/Co60_1330_rates.C");
  c4->SaveAs("TowerAnalysis/Co60_1330_tower0.pdf");
  c4->SaveAs("TowerAnalysis/Co60_1330_tower0.C");
  */

  
}
