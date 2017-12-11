
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


int AnalyzeK40Tower12() {

  int nbins = 1000;
  int K40start = 1000;
  int K40end = 2000;

  TFile * f1 = new TFile("/projects/cuore/data/ds3018_3021/Reduced_bkg_3018_3021.root");
  TTree * t1 = (TTree*)f1->Get("outTree");
  
  TH1F* K40 = new TH1F("K40", "K40", nbins, K40start, K40end);

    
  TCut cutK40 = "Energy < 1475 && Energy > 1445 && Tower==19";
  
  t1->Draw("Energy >> K40", cutK40, "goff");

  /*
  TF1 * fit = new TF1("fit", "[0] * exp(-0.5 * ((x - [1]) / [2])^2) + [3] + x* [4]", -1, 3000);

  fit->SetParameter(0, 2000);
  fit->SetParameter(1, 1460);
  fit->SetParameter(2, 10);
  fit->SetParameter(3, 4000);
  fit->SetParameter(4, -2);
  
  K40->Fit("fit", "", "", 1350, 1550);
  */
  
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
  RooRealVar gaus1frac("g1frac", "fraction of main gaussian", 0.95, 0.6, 1);
  RooAddPdf doublegauss("doublegaus", "gauss+gauss2", RooArgList(gauss, gauss2), RooArgList(gaus1frac));
  
  RooRealVar gaussfrac("gfrac", "fraction of gauss", 0.6, 0, 1);
  RooAddPdf gausslin("gausslin", "gauss+p1", RooArgList(doublegauss, p1), RooArgList(gaussfrac));


  x.setRange(1430, 1490);
  RooDataHist dataK40("dataK40", "K-40", x, K40);
  RooPlot * frameK40 = x.frame(Title("K40"));

  mean.setVal(1462);
  mean.setRange(1460, 1470);
  sigma.setVal(3);

  mean2.setVal(1452);
  mean2.setRange(1450,1455);
  sigma2.setVal(3);

  //a0.setVal(400);
  //a0.setRange(0,500);
  
  dataK40.plotOn(frameK40);
  gausslin.fitTo(dataK40);
  gausslin.plotOn(frameK40);
  gausslin.plotOn(frameK40, Components(p1), LineStyle(kDashed));
  gausslin.plotOn(frameK40, Components(gauss2), LineStyle(kDashed), LineColor(kRed));

  frameK40->Draw();

  double linear = a1.getVal();
  double flat = a0.getVal();
  double signal = gaussfrac.getVal();
  
  //double signal = gaussfrac.getVal();
  
  //cout << flat << "\t" << linear << endl;

  //cout << K40->Integral() <<endl;
  //cout << linear * 0.5 * (1550*1550 - 1350*1350) + flat * 200 << endl;
  cout << (1 - signal) * K40->Integral() << endl;


  double Rate[13];
  double RateError[13];
  double Floor[13];
  double FloorError[13];
  double signal_floor[13];
  double signalerror_floor[13];
  double integral_floor[13];

  RooPlot * frameK40_floor[13];
  
  TCanvas *c2 = new TCanvas();
  c2->Divide(5,3);
  
  for (int floor = 1; floor <= 13; floor++)
    {
      c2->cd(floor);
      
    
      std::ostringstream ss;
      ss << "Floor == " << floor;

      cout  << ss.str() << endl;

      
      if (K40) delete K40;
      TH1F* K40 = new TH1F("K40", "K40", nbins, K40start, K40end);
      frameK40_floor[floor-1] = x.frame(Title(ss.str().c_str()));

      TCut cutFloor = ss.str().c_str();
      TCut fullCut = cutFloor && cutK40;
      t1->Draw("Energy >> K40", fullCut, "goff");

      RooDataHist dataK40("dataK40", "K-40", x, K40);
      
      mean.setVal(1462);
      mean.setRange(1458, 1465);
      sigma.setVal(3);

      mean2.setVal(1452);
      mean2.setRange(1450,1460);
      sigma2.setVal(3);

      a0.setVal(1);
      a1.setVal(1);

      gaus1frac.setVal(0.97);
      gaussfrac.setVal(0.6);
      
      dataK40.plotOn(frameK40_floor[floor-1], MarkerSize(0.2));
      gausslin.fitTo(dataK40);
      gausslin.plotOn(frameK40_floor[floor-1]);
      gausslin.plotOn(frameK40_floor[floor-1], Components(p1), LineStyle(kDashed));
      gausslin.plotOn(frameK40_floor[floor-1], Components(gauss2), LineStyle(kDashed), LineColor(kRed));


      frameK40_floor[floor-1]->Draw();
      frameK40_floor[floor-1]->GetXaxis()->SetTitle("Energy");
      frameK40_floor[floor-1]->GetYaxis()->SetRangeUser(0,30);
      
      double linear = a1.getVal();
      double flat = a0.getVal();
      signal_floor[floor-1] = gaussfrac.getVal();
      signalerror_floor[floor-1] = gaussfrac.getError();

      integral_floor[floor-1] = K40->Integral();
      
      Rate[floor-1] = K40->Integral();
      Floor[floor-1] = floor;
      RateError[floor-1] = sqrt(K40->Integral());
      FloorError[floor-1] = 0;
    }

  TCanvas * c4 = new TCanvas();
  c4->cd();
  frameK40_floor[0]->Draw();
  frameK40_floor[0]->GetXaxis()->SetTitle("Energy");
  frameK40_floor[0]->GetYaxis()->SetRangeUser(0,30);


  for (int j = 0; j < 13; j++)
    {
      cout << "floor: " << j+1 << " Events: " << integral_floor[j] << " Signal: " << signal_floor[j] << " Error: " << signalerror_floor[j] << endl; 
    }
  
  int n = 13;
  TGraphErrors * RatesByFloor = new TGraphErrors(n, Floor, Rate, FloorError, RateError);

  TCanvas * c3 = new TCanvas();
  c3->cd();
  RatesByFloor->Draw("APE");
  RatesByFloor->SetTitle("K40 Rates by Floor");
  RatesByFloor->GetXaxis()->SetRangeUser(0,14);
  RatesByFloor->GetXaxis()->SetNdivisions(210,kTRUE);
  RatesByFloor->GetXaxis()->SetTitle("Floor");
  RatesByFloor->GetYaxis()->SetTitle("Events");
  RatesByFloor->GetYaxis()->SetTitleOffset(1.3);
  c3->SetGridy();

  /*
  c1->SaveAs("FloorAnalysis/K40_sum.pdf");
  c1->SaveAs("FloorAnalysis/K40_sum.C");
  c2->SaveAs("FloorAnalysis/K40_floors.pdf");
  c2->SaveAs("FloorAnalysis/K40_floors.C");
  c3->SaveAs("FloorAnalysis/K40_rates.pdf");
  c3->SaveAs("FloorAnalysis/K40_rates.C");
  c4->SaveAs("FloorAnalysis/K40_floor0.pdf");
  c4->SaveAs("FloorAnalysis/K40_floor0.C");
  */
}
