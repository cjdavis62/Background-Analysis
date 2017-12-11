
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


int AnalyzePo210ByTower() {

  int nbins = 1000;
  int Po210start = 5000;
  int Po210end = 6000;

  TFile * f1 = new TFile("/projects/cuore/data/ds3018_3021/Reduced_bkg_3018_3021.root");
  TTree * t1 = (TTree*)f1->Get("outTree");
  
  TH1F* Po210 = new TH1F("Po210", "Po210", nbins, Po210start, Po210end);

    
  TCut cutPo210 = "Energy < 5500 && Energy > 5250";
  
  t1->Draw("Energy >> Po210", cutPo210, "goff");

  /*
  TF1 * fit = new TF1("fit", "[0] * exp(-0.5 * ((x - [1]) / [2])^2) + [3] + x* [4]", -1, 3000);

  fit->SetParameter(0, 2000);
  fit->SetParameter(1, 1460);
  fit->SetParameter(2, 10);
  fit->SetParameter(3, 4000);
  fit->SetParameter(4, -2);
  
  Po210->Fit("fit", "", "", 1350, 1550);
  */
  
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
  
  RooRealVar peakfrac("gfrac", "fraction of gauss", 0.8, 0, 1);
  RooAddPdf peaklin("gausslin", "cball+p0", RooArgList(triplepeak, p0), RooArgList(peakfrac));


  x.setRange(5250, 5500);
  RooDataHist dataPo210("dataPo210", "Co-60", x, Po210);
  RooPlot * framePo210 = x.frame(Title("Po210"));

  //subgaussfrac.setVal(0.99);
  //subgaussfrac.setRange(0.98,1.0);
  
  mean.setVal(5340);
  mean.setRange(5300, 5400);
  sigma.setVal(3);

  mean2.setVal(5440);
  mean2.setRange(5430,5450);
  sigma2.setVal(3);

  mean3.setVal(5310);
  mean3.setRange(5300, 5335);
  sigma3.setVal(1);
  
  a0.setVal(50);
  a0.setRange(0,100);

  TCanvas * c1 = new TCanvas();
  c1->cd();
  
  dataPo210.plotOn(framePo210);
  peaklin.fitTo(dataPo210);
  peaklin.plotOn(framePo210);
  peaklin.plotOn(framePo210, Components(p0), LineStyle(kDashed));
  peaklin.plotOn(framePo210, Components(subgauss), LineStyle(kDashed), LineColor(kRed));
  
  
  framePo210->Draw();
  framePo210->GetXaxis()->SetTitle("Energy");

  double flat = a0.getVal();
  double signal = peakfrac.getVal();
  
  cout << (1 - signal) * Po210->Integral() << endl;
  
  double Rate[19];
  double RateError[19];
  double Tower[19];
  double TowerError[19];
  double signal_tower[19];
  double integral_tower[19];

  //TH1F *Po210_tower[19];
  RooPlot * framePo210_tower[19];

  TCanvas *c2 = new TCanvas();
  c2->Divide(5,4);

  for (int tower = 1; tower <= 19; tower++)
    {
      std::ostringstream ss;
      ss << "Tower == " << tower;
      
      cout  << ss.str() << endl;
      
      c2->cd(tower);
      
      if (Po210) delete Po210;
      
      TH1F* Po210 = new TH1F("Po210", "Po210", nbins, Po210start, Po210end);
      framePo210_tower[tower-1] = x.frame(Title(ss.str().c_str()));
      
      TCut cutTower = ss.str().c_str();
      TCut fullCut = cutTower && cutPo210;
      t1->Draw("Energy >> Po210", fullCut, "goff");

      RooDataHist dataPo210("dataPo210", "Po-210", x, Po210);

      mean.setVal(5340);
      mean.setRange(5300, 5400);
      sigma.setVal(3);

      mean2.setVal(5440);
      mean2.setRange(5430,5450);
      sigma2.setVal(3);

      mean3.setVal(5310);
      sigma3.setVal(1);
  
      a0.setVal(50);
      a0.setRange(0,100);

      tail.setVal(1);
      norm.setVal(0.5);
      tail2.setVal(1);
      norm2.setVal(0.5);
      cball1frac.setVal(0.8);
      subgaussfrac.setVal(0.8);
      peakfrac.setVal(0.8);
      
      
      dataPo210.plotOn(framePo210_tower[tower-1], MarkerSize(0.2));
      peaklin.fitTo(dataPo210);
      peaklin.plotOn(framePo210_tower[tower-1]);
      peaklin.plotOn(framePo210_tower[tower-1], Components(p0), LineStyle(kDashed));
      peaklin.plotOn(framePo210_tower[tower-1], Components(p0), LineStyle(kDashed));
      peaklin.plotOn(framePo210_tower[tower-1], Components(subgauss), LineStyle(kDashed), LineColor(kRed));

      framePo210_tower[tower-1]->Draw();
      framePo210_tower[tower-1]->GetYaxis()->SetRangeUser(0,200);
      framePo210_tower[tower-1]->GetXaxis()->SetTitle("Energy");
      //      framePo210_tower[tower-1]->SetMarkerSize(0.3);
      double flat = a0.getVal();
      signal_tower[tower-1] = peakfrac.getVal();
      double signalerror = peakfrac.getError();

      integral_tower[tower-1] = Po210->Integral();
            
      Rate[tower-1] = signal_tower[tower-1] * Po210->Integral();
      Tower[tower-1] = tower;
      RateError[tower-1] = signalerror * Po210->Integral();
      TowerError[tower-1] = 0;
    }

  TCanvas * c4 = new TCanvas();
  c4->cd();
  framePo210_tower[0]->Draw();
  framePo210_tower[0]->GetYaxis()->SetRangeUser(0,200);
  framePo210_tower[0]->GetXaxis()->SetTitle("Energy");
  
  for (int j = 0; j < 19; j++)
    {
      cout << "tower: " << j+1 << " Events: " << integral_tower[j] << " Signal: " << signal_tower[j] << endl; 
    }

  int n = 19;
  TGraphErrors * RatesByTower = new TGraphErrors(n, Tower, Rate, TowerError, RateError);

  TCanvas * c3 = new TCanvas();
  c3->cd();
  RatesByTower->Draw("APE");
  RatesByTower->SetTitle("Po210 Rates by Tower");
  RatesByTower->GetXaxis()->SetRangeUser(0,20);
  RatesByTower->GetXaxis()->SetNdivisions(210,kTRUE);
  RatesByTower->GetXaxis()->SetTitle("Tower");
  RatesByTower->GetYaxis()->SetTitle("Events");
  RatesByTower->GetYaxis()->SetTitleOffset(1.3);
  c3->SetGridy();

  cout << signal << endl;


  //Save plots
  c1->SaveAs("TowerAnalysis/Po210_sum.pdf");
  c1->SaveAs("TowerAnalysis/Po210_sum.C");
  c2->SaveAs("TowerAnalysis/Po210_towers.pdf");
  c2->SaveAs("TowerAnalysis/Po210_towers.C");
  c3->SaveAs("TowerAnalysis/Po210_rates.pdf");
  c3->SaveAs("TowerAnalysis/Po210_rates.C");
  c4->SaveAs("TowerAnalysis/Po210_tower0.pdf");
  c4->SaveAs("TowerAnalysis/Po210_Tower0.C");

  
}
