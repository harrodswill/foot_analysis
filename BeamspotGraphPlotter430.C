#include <iostream>
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLatex.h"

void BeamspotGraphPlotter430()
{
  TCanvas * canvas_beamspots = new TCanvas("canvas_beamplots","canvas_beamplots",1200,1200);
  canvas_beamspots->Divide(1,1);

  TFile * file_run93 = TFile::Open("output_run93_500k.root");
  TFile * file_run94 = TFile::Open("output_run94_500k.root");
  TFile * file_run95 = TFile::Open("output_run95_500k.root");
  TFile * file_run96 = TFile::Open("output_run96_500k.root");
  TFile * file_run97 = TFile::Open("output_run97_500k.root");
  TFile * file_run98 = TFile::Open("output_run98_500k.root");
  TFile * file_run100 = TFile::Open("output_run100_500k.root");
  TFile * file_run101 = TFile::Open("output_run101_500k.root");
  TFile * file_run102 = TFile::Open("output_run102_500k.root");

  TH2D * run93 = (TH2D *) file_run93->Get("h2_beam_XY");
  TH2D * run94 = (TH2D *) file_run94->Get("h2_beam_XY");
  TH2D * run95 = (TH2D *) file_run95->Get("h2_beam_XY");
  TH2D * run96 = (TH2D *) file_run96->Get("h2_beam_XY");
  TH2D * run97 = (TH2D *) file_run97->Get("h2_beam_XY");
  TH2D * run98 = (TH2D *) file_run98->Get("h2_beam_XY");
  TH2D * run100 = (TH2D *) file_run100->Get("h2_beam_XY");
  TH2D * run101 = (TH2D *) file_run101->Get("h2_beam_XY");
  TH2D * run102 = (TH2D *) file_run102->Get("h2_beam_XY");


  canvas_beamspots->cd(1);
  gStyle->SetPalette(kRainBow);
  run93->Draw("colz");
  run94->Draw("colz same");
  run95->Draw("colz same");
  run96->Draw("colz same");
  run97->Draw("colz same");
  run98->Draw("colz same");
  run100->Draw("colz same");
  run101->Draw("colz same");
  run102->Draw("colz same");

  TBox * box93 = new TBox(37.88-21/2,34.05-21/2,37.88+21/2,34.05+21/2);
  box93->SetFillStyle(0);
  box93->SetLineColor(kRed);
  box93->SetLineStyle(10);
  box93->SetLineWidth(2);
  box93->Draw();

  TLatex * text93 = new TLatex(37.88-21/2-5,34.05+21/2+0.5,"Run 93");
  text93->SetTextSize(0.025);
  text93->SetTextColor(kRed);
  text93->Draw();

  TBox * box94 = new TBox(-42.33-21/2,34.01-21/2,-42.33+21/2,34.01+21/2);
  box94->SetFillStyle(0);
  box94->SetLineColor(kRed);
  box94->SetLineStyle(10);
  box94->SetLineWidth(2);
  box94->Draw();

  TLatex * text94 = new TLatex(-42.33-21/2+4,34.01+21/2+0.5,"Run 94");
  text94->SetTextSize(0.025);
  text94->SetTextColor(kRed);
  text94->Draw();

  TBox * box95 = new TBox(-41.95-21/2,-3.85-21/2,-41.95+21/2,-3.85+21/2);
  box95->SetFillStyle(0);
  box95->SetLineColor(kRed);
  box95->SetLineStyle(10);
  box95->SetLineWidth(2);
  box95->Draw();

  TLatex * text95 = new TLatex(-41.95-21/2+4,-3.85+21/2+0.5,"Run 95");
  text95->SetTextSize(0.025);
  text95->SetTextColor(kRed);
  text95->Draw();

  TBox * box96 = new TBox(-35.44-21/2,-39.91-21/2,-35.44+21/2,-39.91+21/2);
  box96->SetFillStyle(0);
  box96->SetLineColor(kRed);
  box96->SetLineStyle(10);
  box96->SetLineWidth(2);
  box96->Draw();

  TLatex * text96 = new TLatex(-35.44-21/2-1,-39.91+21/2+0.5,"Run 96");
  text96->SetTextSize(0.025);
  text96->SetTextColor(kRed);
  text96->Draw();
  
  TBox * box97 = new TBox(1.58-21/2,-40.68-21/2,1.58+21/2,-40.68+21/2);
  box97->SetFillStyle(0);
  box97->SetLineColor(kRed);
  box97->SetLineStyle(10);
  box97->SetLineWidth(2);
  box97->Draw();

  TLatex * text97 = new TLatex(1.58-21/2-5,-40.68+21/2+0.5,"Run 97");
  text97->SetTextSize(0.025);
  text97->SetTextColor(kRed);
  text97->Draw();
  
  TBox * box98 = new TBox(39.46-21/2,-40.73-21/2,39.46+21/2,-40.73+21/2);
  box98->SetFillStyle(0);
  box98->SetLineColor(kRed);
  box98->SetLineStyle(10);
  box98->SetLineWidth(2);
  box98->Draw();

  TLatex * text98 = new TLatex(39.46-21/2-5,-40.73+21/2+0.5,"Run 98");
  text98->SetTextSize(0.025);
  text98->SetTextColor(kRed);
  text98->Draw();
   
  TBox * box100 = new TBox(38.97-21/2,-3.23-21/2,38.97+21/2,-3.23+21/2);
  box100->SetFillStyle(0);
  box100->SetLineColor(kRed);
  box100->SetLineStyle(10);
  box100->SetLineWidth(2);
  box100->Draw();

  TLatex * text100 = new TLatex(38.97-21/2-5,-3.23+21/2+0.5,"Run 100");
  text100->SetTextSize(0.025);
  text100->SetTextColor(kRed);
  text100->Draw();

  TBox * box101 = new TBox(3.32-21/2,34.00-21/2,3.32+21/2,34.00+21/2);
  box101->SetFillStyle(0);
  box101->SetLineColor(kRed);
  box101->SetLineStyle(10);
  box101->SetLineWidth(2);
  box101->Draw();

  TLatex * text101 = new TLatex(3.32-21/2-5,34.00+21/2+0.5,"Run 101");
  text101->SetTextSize(0.025);
  text101->SetTextColor(kRed);
  text101->Draw();

  TBox * box102 = new TBox(0.32-21/2,-3.66-21/2,0.32+21/2,-3.66+21/2);
  box102->SetFillStyle(0);
  box102->SetLineColor(kRed);
  box102->SetLineStyle(10);
  box102->SetLineWidth(2);
  box102->Draw();

  TLatex * text102 = new TLatex(0.32-21/2-5,-3.66+21/2+0.5,"Run 102");
  text102->SetTextSize(0.025);
  text102->SetTextColor(kRed);
  text102->Draw();
  
}

