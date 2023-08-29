#include <iostream>
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLatex.h"

void BeamspotGraphPlotter()
{
  TCanvas * canvas_beamspots = new TCanvas("canvas_beamplots","canvas_beamplots",1200,1200);
  canvas_beamspots->Divide(1,1);

  TFile * file_run52 = TFile::Open("output_run52_500k.root");
  TFile * file_run53 = TFile::Open("output_run53_500k.root");
  TFile * file_run54 = TFile::Open("output_run54_500k.root");
  TFile * file_run56 = TFile::Open("output_run56_500k.root");
  TFile * file_run57 = TFile::Open("output_run57_500k.root");
  TFile * file_run58 = TFile::Open("output_run58_500k.root");
  TFile * file_run59 = TFile::Open("output_run59_500k.root");
  TFile * file_run60 = TFile::Open("output_run60_500k.root");
  TFile * file_run61 = TFile::Open("output_run61_500k.root");

  TH2D * run52 = (TH2D *) file_run52->Get("h2_beam_XY");
  TH2D * run53 = (TH2D *) file_run53->Get("h2_beam_XY");
  TH2D * run54 = (TH2D *) file_run54->Get("h2_beam_XY");
  TH2D * run56 = (TH2D *) file_run56->Get("h2_beam_XY");
  TH2D * run57 = (TH2D *) file_run57->Get("h2_beam_XY");
  TH2D * run58 = (TH2D *) file_run58->Get("h2_beam_XY");
  TH2D * run59 = (TH2D *) file_run59->Get("h2_beam_XY");
  TH2D * run60 = (TH2D *) file_run60->Get("h2_beam_XY");
  TH2D * run61 = (TH2D *) file_run61->Get("h2_beam_XY");


  canvas_beamspots->cd(1);
  gStyle->SetPalette(kRainBow);
  run52->Draw("colz");
  run53->Draw("colz same");
  run54->Draw("colz same");
  run56->Draw("colz same");
  run57->Draw("colz same");
  run58->Draw("colz same");
  run59->Draw("colz same");
  run60->Draw("colz same");
  run61->Draw("colz same");

  TBox * box52 = new TBox(0.62-13/2,33.56-13/2,0.62+13/2,33.56+13/2);
  box52->SetFillStyle(0);
  box52->SetLineColor(kRed);
  box52->SetLineStyle(10);
  box52->SetLineWidth(2);
  box52->Draw();

  TLatex * text52 = new TLatex(0.62-13/2-5,33.56+13/2+0.5,"Run 52");
  text52->SetTextSize(0.025);
  text52->SetTextColor(kRed);
  text52->Draw();

  TBox * box53 = new TBox(39.95-13/2,33.32-13/2,39.95+13/2,33.32+13/2);
  box53->SetFillStyle(0);
  box53->SetLineColor(kRed);
  box53->SetLineStyle(10);
  box53->SetLineWidth(2);
  box53->Draw();

    TLatex * text53 = new TLatex(39.95-13/2-5,33.32+13/2+0.5,"Run 53");
  text53->SetTextSize(0.025);
  text53->SetTextColor(kRed);
  text53->Draw();

  TBox * box54 = new TBox(-39.76-13/2,32.94-13/2,-39.76+13/2,32.94+13/2);
  box54->SetFillStyle(0);
  box54->SetLineColor(kRed);
  box54->SetLineStyle(10);
  box54->SetLineWidth(2);
  box54->Draw();

  TLatex * text54 = new TLatex(-39.76-13/2-2,32.94+13/2+0.5,"Run 54");
  text54->SetTextSize(0.025);
  text54->SetTextColor(kRed);
  text54->Draw();

  TBox * box56 = new TBox(-29.32-13/2,-4.19-13/2,-29.32+13/2,-4.19+13/2);
  box56->SetFillStyle(0);
  box56->SetLineColor(kRed);
  box56->SetLineStyle(10);
  box56->SetLineWidth(2);
  box56->Draw();

  TLatex * text56 = new TLatex(-29.32-13/2-5,-4.19+13/2+0.5,"Run 56");
  text56->SetTextSize(0.025);
  text56->SetTextColor(kRed);
  text56->Draw();
  
  TBox * box57 = new TBox(-28.62-13/2,-35.43-13/2,-28.62+13/2,-35.43+13/2);
  box57->SetFillStyle(0);
  box57->SetLineColor(kRed);
  box57->SetLineStyle(10);
  box57->SetLineWidth(2);
  box57->Draw();

  TLatex * text57 = new TLatex(-28.62-13/2-5,-35.43+13/2+0.5,"Run 57");
  text57->SetTextSize(0.025);
  text57->SetTextColor(kRed);
  text57->Draw();
  
  TBox * box58 = new TBox(1.59-13/2,-36.14-13/2,1.59+13/2,-36.14+13/2);
  box58->SetFillStyle(0);
  box58->SetLineColor(kRed);
  box58->SetLineStyle(10);
  box58->SetLineWidth(2);
  box58->Draw();

  TLatex * text58 = new TLatex(1.59-13/2-5,-36.14+13/2+0.5,"Run 58");
  text58->SetTextSize(0.025);
  text58->SetTextColor(kRed);
  text58->Draw();
   
  TBox * box59 = new TBox(31.62-13/2,-36.20-13/2,31.63+13/2,-36.20+13/2);
  box59->SetFillStyle(0);
  box59->SetLineColor(kRed);
  box59->SetLineStyle(10);
  box59->SetLineWidth(2);
  box59->Draw();

  TLatex * text59 = new TLatex(31.62-13/2-5,-36.20+13/2+0.5,"Run 59");
  text59->SetTextSize(0.025);
  text59->SetTextColor(kRed);
  text59->Draw();

  TBox * box60 = new TBox(31.05-13/2,-4.59-13/2,31.05+13/2,-4.59+13/2);
  box60->SetFillStyle(0);
  box60->SetLineColor(kRed);
  box60->SetLineStyle(10);
  box60->SetLineWidth(2);
  box60->Draw();

  TLatex * text60 = new TLatex(31.05-13/2-5,-4.59+13/2+0.5,"Run 60");
  text60->SetTextSize(0.025);
  text60->SetTextColor(kRed);
  text60->Draw();

  TBox * box61 = new TBox(-2.02-13/2,-1.01-13/2,-2.02+13/2,-1.01+13/2);
  box61->SetFillStyle(0);
  box61->SetLineColor(kRed);
  box61->SetLineStyle(10);
  box61->SetLineWidth(2);
  box61->Draw();

  TLatex * text61 = new TLatex(-2.02-13/2-5,-1.01+13/2+0.5,"Run 61");
  text61->SetTextSize(0.025);
  text61->SetTextColor(kRed);
  text61->Draw();
  
}

