
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


int AnalyzeCo60ByTower1170() {

  int nbins = 1000;
  int Co60start = 1000;
  int Co60end = 2000;

  TFile * f1 = new TFile("/projects/cuore/data/ds3018_3021/Reduced_bkg_3018_3021.root");
  TTree * t1 = (TTree*)f1->Get("outTree");
  
  TH1F* Co60 = new TH1F("Co60", "Co60", nbins, Co60start, Co60end);

    
  TCut cutCo60 = "Energy < 1250 && Energy > 1050";
  
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
  RooRealVar gauss1frac("g1frac", "fraction of main gaussian", 0.8, 0, 1);
  RooAddPdf doublegauss("doublegaus", "gauss+gauss2", RooArgList(gauss, gauss2), RooArgList(gauss1frac));
  
  RooRealVar gaussfrac("gfrac", "fraction of gauss", 0.8, 0, 1);
  RooAddPdf gausslin("gausslin", "gauss+p1", RooArgList(doublegauss, p1), RooArgList(gaussfrac));


  x.setRange(1050, 1250);
  RooDataHist dataCo60("dataCo60", "Co-60", x, Co60);
  RooPlot * frameCo60 = x.frame(Title("Co60"));

  mean.setVal(1170);
  mean.setRange(1160, 1180);
  sigma.setVal(3);

  mean2.setVal(1120);
  mean2.setRange(1110,1130);
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
  double signal = gaussfrac.getVal()*gauss1frac.getVal();
  
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

  TCanvas *c1 = new TCanvas();
  c1->Divide(5,4);
  
  for (int tower = 1; tower <= 19; tower++)
    {
      c1->cd(tower);
      
      if (Co60) delete Co60;
      TH1F* Co60 = new TH1F("Co60", "Co60", nbins, Co60start, Co60end);
      std::ostringstream ss;
      ss << "Tower == " << tower;

      cout  << ss.str() << endl;

      frameCo60_tower[tower-1] = x.frame(Title(ss.str().c_str()));
      
      TCut cutTower = ss.str().c_str();
      TCut fullCut = cutTower && cutCo60;
      t1->Draw("Energy >> Co60", fullCut, "goff");

      RooDataHist dataCo60("dataCo60", "Co-60", x, Co60);
      
      dataCo60.plotOn(frameCo60_tower[tower-1], MarkerSize(0.2));
      gausslin.fitTo(dataCo60);
      gausslin.plotOn(frameCo60_tower[tower-1]);
      gausslin.plotOn(frameCo60_tower[tower-1], Components(p1), LineStyle(kDashed));

      frameCo60_tower[tower-1]->Draw();

      double linear = a1.getVal();
      double flat = a0.getVal();
      signal_tower[tower-1] = gaussfrac.getVal();
      signalerror_tower[tower-1] = gaussfrac.getError();

      integral_tower[tower-1] = Co60->Integral();
      
      Rate[tower-1] = signal_tower[tower-1] * Co60->Integral();
      Tower[tower-1] = tower;
      RateError[tower-1] = signalerror_tower[tower-1] * Co60->Integral();
      TowerError[tower-1] = 0;
    }


  int n = 19;
  TGraphErrors * RatesByTower = new TGraphErrors(n, Tower, Rate, TowerError, RateError);

  TCanvas * c2 = new TCanvas();
  c2->cd();
  RatesByTower->Draw();

}
