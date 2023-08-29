#include <iostream>
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
void SigmaGraphPlotter430()
{ 
  TCanvas * canvas_sigmas = new TCanvas("canvas_sigmas","canvas_sigmas",1200,1200);
  canvas_sigmas->Divide(2,2);

  TFile * file_2sigma = TFile::Open("output_run102_500k_2sigma.root");
  TFile * file_3sigma = TFile::Open("output_run102_500k_3sigma.root");
  TFile * file_4sigma = TFile::Open("output_run102_500k_4sigma.root");
  TFile * file_5sigma = TFile::Open("output_run102_500k_5sigma.root");

  TH2D * sigma2 = (TH2D *) file_2sigma->Get("h2_beam_XY");
  TH2D * sigma3 = (TH2D *) file_3sigma->Get("h2_beam_XY");
  TH2D * sigma4 = (TH2D *) file_4sigma->Get("h2_beam_XY");
  TH2D * sigma5 = (TH2D *) file_5sigma->Get("h2_beam_XY");

  canvas_sigmas->cd(1);
  gStyle->SetPalette(kRainBow);
  sigma2->SetTitle("2 Sigma Beamplot");
  sigma2->Draw("colz");

  TBox * box2 = new TBox(0.32-21/2,-3.66-21/2,0.32+21/2,-3.66+21/2);
  box2->SetFillStyle(0);
  box2->SetLineColor(kRed);
  box2->SetLineStyle(10);
  box2->SetLineWidth(2);
  box2->Draw();

  canvas_sigmas->cd(2);
  gStyle->SetPalette(kRainBow);
  sigma3->SetTitle("3 Sigma Beamplot");
  sigma3->Draw("colz");

  TBox * box3 = new TBox(0.32-21/2,-3.66-21/2,0.32+21/2,-3.66+21/2);
  box3->SetFillStyle(0);
  box3->SetLineColor(kRed);
  box3->SetLineStyle(10);
  box3->SetLineWidth(2);
  box3->Draw();
 
  canvas_sigmas->cd(3);
  gStyle->SetPalette(kRainBow);
  sigma4->SetTitle("4 Sigma Beamplot");
  sigma4->Draw("colz");

  TBox * box4 = new TBox(0.32-21/2,-3.66-21/2,0.32+21/2,-3.66+21/2);
  box4->SetFillStyle(0);
  box4->SetLineColor(kRed);
  box4->SetLineStyle(10);
  box4->SetLineWidth(2);
  box4->Draw();
 
  canvas_sigmas->cd(4);
  gStyle->SetPalette(kRainBow);
  sigma5->SetTitle("5 Sigma Beamplot");
  sigma5->Draw("colz");

  TBox * box5 = new TBox(0.32-21/2,-3.66-21/2,0.32+21/2,-3.66+21/2);
  box5->SetFillStyle(0);
  box5->SetLineColor(kRed);
  box5->SetLineStyle(10);
  box5->SetLineWidth(2);
  box5->Draw();

}
