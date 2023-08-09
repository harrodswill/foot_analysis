#include"libs.hh"
#include <TMathBase.h>
#define LENGTH(x) (sizeof x / sizeof *x)
#define NDETS LENGTH(foot_id)
using namespace std;
namespace{//don't change this!!!
  int foot_id[] = {15,16};
}
//Histogram binning and ranges
int Nbins = 500; double ymin = -50; double ymax = 400;
const double    FOOT_LENGTH   = 96.;//mm
const double    NSIGMA        = 6;

//Clustering data structures
typedef std::pair<int, double> strip_data;
typedef std::vector<strip_data> cluster; 
typedef std::vector<cluster> foot_data; 

Double_t pedestal[NDETS][640];
Double_t sigma[NDETS][640];
Double_t sigma_fine[NDETS][640];

bool is_good_strip(UInt_t det, UInt_t strip)
{
  switch(det){
    case 15:
      if(
	  strip==344
	) return false;
    case 16:
      if(
	  strip==16 || strip==114 || (strip>139 && strip<146) || strip==462 || strip==463 || strip==464 ||
	  strip==465 || strip==466 || strip==409 || strip==410 || strip==426 || strip==433 || strip==434 || strip==440 || strip==520
	) return false;
  }
  if((strip%64)>62 || (strip%64)<3) return false;//edge strips for every asic
  return true;
}

double get_cog(cluster c)//calculate center of gravity of a cluster
{
  double value = 0;  double esum = 0.;
  for(auto & s: c){
    value  += (s.first * s.second);
    esum += s.second;
  }
  value /= esum;//center of gravity
  return value;
}

double get_esum(cluster c)//calculate cluster sum
{
  double esum = 0.;
  for(auto & s: c){
    esum += s.second;
  }
  return esum;
}

void Check_Strip(UInt_t number, double energy, foot_data& fdata)
{
  cout << "\n\n---------------------------------------";
  cout << "\nChecking strip: " << number << "\t Energy: " << energy << endl;
  strip_data strip = std::make_pair(number,energy);
  cluster clust;
  if(fdata.size()==0)//no cluster yet, create new
  {
    clust.push_back(strip);
    fdata.push_back(clust);
    cout << "\n\t New cluster is created for this strip";
    return;
  }
  cluster    this_clust = fdata.back();
  strip_data this_strip = this_clust.back();
  if(abs(strip.first-this_strip.first)<2)//neighbour found 
  {
    cout << "\n\tStrip belong to exisitng cluster! Adding it...";
    fdata.back().push_back(strip);
    return;
  }
  else
  {
    cout << "\n\tStrip is a new cluster! Making it...";
    clust.clear();
    clust.push_back(strip);
    fdata.push_back(clust);
    return;
  }
}

void analyse(int firstEvent, int max_events, TChain * ch)
{
  TApplication* theApp = new TApplication("App", 0, 0);
  //---------- Definition of main  histograms ----------------
  TH2F * h2_peds_raw[NDETS];
  TH2F * h2_cal_coarse[NDETS];
  TH2F * h2_cal_fine[NDETS];
  TH2F * h2_baseline[NDETS];
  TH1D * h1_baseline[NDETS];
  TH1D * h1_peds[NDETS];
  TH1D * h1_sigma_raw[NDETS];
  TH1D * h1_sigma_fine[NDETS];
  TH1D * h1_cluster_e[NDETS];
  TH1I * h1_cluster_size[NDETS];
  TH2D * h2_cluster_e_vs_cog[NDETS];
  TH2D * h2_cluster_cog_vs_cog;
  TH2D * h2_beam_XY;
  TH1D * h1_beam_X;
  TH1D * h1_beam_Y;
  TH1I * h1_det_mul;
  for(int i=0; i<NDETS; i++)
  {
    h2_cal_coarse[i]   = new TH2F(Form("h2_cal_coarse_FOOT%d",foot_id[i]),Form("h2_cal_coarse_FOOT%d",foot_id[i]),640,1,641,Nbins,ymin,ymax);
    h2_cal_fine[i]     = new TH2F(Form("h2_cal_fine_FOOT%d",foot_id[i]),Form("h2_cal_fine_FOOT%d",foot_id[i]),640,1,641,Nbins,ymin,ymax);
    h2_peds_raw[i]     = new TH2F(Form("h%d",foot_id[i]),Form("h%d",foot_id[i]),640,1,641,1500,0,1500);
    h2_baseline[i]     = new TH2F(Form("h%d_baseline",foot_id[i]),Form("h%d_baseline",foot_id[i]),640,1,641,500,-100,100);
    h1_cluster_e[i]    = new TH1D(Form("h1_cluster_e_FOOT%d",foot_id[i]),Form("h1_cluster_e_FOOT%d",foot_id[i]),1000,-10,100);
    h1_cluster_size[i] = new TH1I(Form("h1_cluster_size_FOOT%d",foot_id[i]),Form("h1_cluster_size_FOOT%d",foot_id[i]),10,0,10);
    h2_cluster_e_vs_cog[i] = new TH2D(Form("h2_cluster_e_vs_cog_FOOT%d",foot_id[i]),Form("h2_cluster_e_vs_cog_FOOT%d",foot_id[i]),640,1,641,200,-10,100);

  }
  h1_det_mul = new TH1I("det_mul","det_mul",10,0,10);

  h2_beam_XY = new TH2D(Form("h2_beam_XY"),Form("h2_beam_XY"),100,-50,50,100,-50,50);
  //h2_beam_XY = new TH2D(Form("h2_beam_XY"),Form("h2_beam_XY"),640,0,640,640,0,640);
  h1_beam_X = new TH1D(Form("h1_beam_X"),Form("h1_beam_X"),640,0,640);
  h1_beam_Y = new TH1D(Form("h1_beam_Y"),Form("h1_beam_Y"),640,0,640);

  TCanvas * canvas_raw_peds = new TCanvas("canvas_raw_peds","canvas_raw_peds",1800,1200);
  canvas_raw_peds->Divide(2,1);

  TCanvas * canvas_baseline = new TCanvas("canvas_baseline","canvas_baseline",1800,1200);
  canvas_baseline->Divide(2,1);

  TCanvas * canvas_raw_sigma = new TCanvas("canvas_raw_sigma","canvas_raw_sigma",1800,1200);
  canvas_raw_sigma->Divide(2,1);

  TCanvas* canvas_coarse = new TCanvas("coarse","coarse",1800,1200);
  canvas_coarse->Divide(2,1);

  TCanvas* canvas_fine   = new TCanvas("fine","fine",1800,1200);
  canvas_fine->Divide(2,1);

  TCanvas* canvas_cluster_energy   = new TCanvas("cluster_energy","cluster_energy",1800,1200);
  canvas_cluster_energy->Divide(2,1);

  TCanvas* canvas_cluster_size   = new TCanvas("cluster_size","cluster_size",1800,1200);
  canvas_cluster_size->Divide(2,1);

  TCanvas* canvas_cluster_e_vs_cog   = new TCanvas("cluster_e_vs_cog","cluster_e_vs_cog",1800,1200);
  canvas_cluster_e_vs_cog->Divide(2,1);

  TCanvas* canvas_det_mul   = new TCanvas("det_mul","det_mul",1800,1200);

  TCanvas* canvas_XY   = new TCanvas("canvas_XY","canvas_XY",1200,1200);
  canvas_XY->Divide(2,2);
  //------ Define tree data for all foot detectors ----------
  UInt_t TRIGGER;
  UInt_t FOOT[NDETS];//multiplicity of foot detectors 
  UInt_t FOOTE[NDETS][640];
  UInt_t FOOTI[NDETS][640];
  UInt_t TPATv[1];
  ch->SetBranchAddress("TRIGGER",&TRIGGER);
  ch->SetBranchAddress("TPATv",TPATv);
  for(int i=0; i<NDETS; i++){
    cout << "\nTree branches being sorted out for FOOT: " << foot_id[i];
    TString bname   = Form("FOOT%d",  foot_id[i]);
    TString bname_E = Form("FOOT%dE", foot_id[i]);
    TString bname_I = Form("FOOT%dI", foot_id[i]);
    ch->SetBranchAddress(bname.Data(),&FOOT[i]);
    ch->SetBranchAddress(bname_E.Data(),FOOTE[i]);
    ch->SetBranchAddress(bname_I.Data(),FOOTI[i]);
  }
  //------ Define output tree  ----------
  TFile* _file = new TFile("output.root","Recreate","Write");
  UInt_t  F15_mul;
  Float_t F15_e[640];
  Float_t F15_cog[640];
  Float_t F15_size[640];
  UInt_t  F16_mul;
  Float_t F16_e[640];
  Float_t F16_cog[640];
  Float_t F16_size[640];
  TTree *tree = new TTree("tree","tree");
  tree->Branch("F15_mul"  , &F15_mul   , "F15_mul/i" );
  tree->Branch("F15_e"    ,  F15_e     , "F15_e[F15_mul]/F" );
  tree->Branch("F15_cog"  ,  F15_cog   , "F15_cog[F15_mul]/F" );
  tree->Branch("F15_size" ,  F15_size  , "F15_size[5_mul]/F" );
  tree->Branch("F16_mul"  , &F16_mul   , "F16_mul/i" );
  tree->Branch("F16_e"    ,  F16_e     , "F16_e[F16_mul]/F" );
  tree->Branch("F16_cog"  ,  F16_cog   , "F16_cog[F16_mul]/F" );
  tree->Branch("F16_size" ,  F16_size  , "F16_size[F16_mul]/F" );

  //----- Get starting set of pedestal data for all dets ------
  //
  int Nevents = ch->GetEntries();
  if(max_events>0) Nevents = max_events; 
  int     stat=0;
  for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
  {
    ch->GetEntry(ev);
    for(int f=0; f<NDETS; f++){
      for(int  j=0 ; j<640 ; j++){
	h2_peds_raw[f]->Fill(FOOTI[f][j],FOOTE[f][j]); 
      }
    }
    stat++;
    if(stat==1e4) break;
  }
  //-------  Slicing, fitting, saving raw pedestals and sigmas -------
  //
  TF1* foo = new TF1("foo","gaus",0,1000);
  for(int i=0; i<NDETS; i++)
  {
    cout << "\nPedestal and raw sigma data being prepared and plotted for FOOT: " << foot_id[i];
    h2_peds_raw[i]->FitSlicesY(foo,1,640,0,"QNR",0);
    h1_peds[i]      = (TH1D*)gDirectory->Get(Form("h%d_1",foot_id[i]))->Clone(Form("h1_peds_%d",  foot_id[i]));
    h1_sigma_raw[i] = (TH1D*)gDirectory->Get(Form("h%d_2",foot_id[i]))->Clone(Form("h1_sigma_raw_%d",foot_id[i]));

    canvas_raw_peds->cd(i+1);
    h2_peds_raw[i]->Draw("colz");
    h1_peds[i]->SetMarkerStyle(kFullCircle);
    h1_peds[i]->SetMarkerSize(1);
    h1_peds[i]->SetMarkerColor(kBlack);
    h1_peds[i]->SetLineColor(kBlack);
    h1_peds[i]->Draw("same");

    for(int c=0; c<10; c++)
    {
      TLine * l = new TLine(c*64,0,c*64,1500);
      l->Draw();
    } 

    canvas_raw_sigma->cd(i+1);
    h1_sigma_raw[i]->SetMarkerStyle(kFullCircle);
    h1_sigma_raw[i]->SetMarkerSize(0.2);
    h1_sigma_raw[i]->SetMarkerColor(kBlue);
    h1_sigma_raw[i]->SetLineColor(kBlue);
    h1_sigma_raw[i]->GetYaxis()->SetRangeUser(-2,10);
    h1_sigma_raw[i]->GetXaxis()->SetTitle("Strip No.");
    h1_sigma_raw[i]->GetYaxis()->SetTitle("ADC sigma");
    h1_sigma_raw[i]->SetTitle("Sigmas before baseline correction");
    h1_sigma_raw[i]->Draw();
    for(int c=0; c<10; c++)
    {
      TLine * l = new TLine(c*64,0,c*64,1500);
      l->Draw();
    } 


    for(int j=0; j<640; j++)
    {
      pedestal[i][j] = h1_peds[i]->GetBinContent(j+1);
      sigma[i][j]    = h1_sigma_raw[i]->GetBinContent(j+1);
      cout << "\nFOOT: " << foot_id[i] << "\tStrip: " << j << "\t Ped: " << pedestal[i][j] << "\tSig: " << sigma[i][j];
    } 
  } 

  //-------- Re-analyse with baseline correction to get fine sigmas -------
  double  mean_ssd = 0;
  double  asic_offset[10];
  double  signal = 0;
  double  signal_sum = 0;
  int     counter_asic =0;
  int     stat_baseline =0;
  for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
  {
    ch->GetEntry(ev);
    for(int f=0; f<NDETS; f++)//loop over all foots
    {
      //--------  Global base line correction in every FOOT in this event ---------
      mean_ssd=0; stat=0;
      for(int i=0; i<640; i++){
	if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
	signal = FOOTE[f][i] - pedestal[f][i];
	if(fabs(signal) > (10 * sigma[f][i])) continue; //possible hit candidate
	stat++;    mean_ssd += signal;
      }
      mean_ssd /= stat;
      if(fabs(mean_ssd)>10){ continue; }
      //------------ Calculating fine baseline correction for individual asics ---------
      stat=0; counter_asic=0;
      for(int i=0; i<10; i++){  asic_offset[i]=0.; }//reset asic baselines
      for(int i=0; i<640; i++){
	signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd;
	if(fabs(signal) < (4 * sigma[f][i]) && //ignore possible hit candidates
	    is_good_strip(foot_id[f],FOOTI[f][i]))//and bad strips
	{
	  stat++;
	  asic_offset[counter_asic] += signal;
	}
	if((FOOTI[f][i]%64)==0){//switch to next asic
	  //	  cout << "\nStrip: " << FOOTI[f][i] " for FOOT: " << foot_id[f] << "is corrected well for individual asics";
	  asic_offset[counter_asic] /= stat;
	  counter_asic++;  stat=0;
	}
      }
      //-------- Applying fine baseline correction and fill histogram ----------
      counter_asic=0;
      for(int i=0; i<640; i++)
      {
	if((FOOTI[f][i]%64) == 1 && FOOTI[f][i]>1) counter_asic++; 
	signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd - asic_offset[counter_asic];
	if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
	h2_baseline[f]->Fill(FOOTI[f][i],signal);
      }
    }//end detector loop
    stat_baseline++; if(stat_baseline==1e4) break;
  }//end evetloop

  //-------  Slicing, fitting, saving fine sigmas -------
  foo = new TF1("foo","gaus",-10,10);
  for(int i=0; i<NDETS; i++){
    cout << "\nFitting lines on FOOT: " << foot_id[i];
    h2_baseline[i]->FitSlicesY(foo,1,640,0,"QNR",0);
    h1_baseline[i]  = (TH1D*)gDirectory->Get(Form("h%d_baseline_1",foot_id[i]))->Clone(Form("h1_baseline_%d",foot_id[i]));
    h1_sigma_fine[i]  = (TH1D*)gDirectory->Get(Form("h%d_baseline_2",foot_id[i]))->Clone(Form("h1_sigma_fine_%d",foot_id[i]));

    canvas_baseline->cd(i+1);
    h2_baseline[i]->Draw("colz");
    h1_baseline[i]->SetMarkerStyle(kFullCircle);
    h1_baseline[i]->SetMarkerSize(0.1);
    h1_baseline[i]->SetMarkerColor(kRed);
    h1_baseline[i]->SetLineColor(kRed);
    h1_baseline[i]->Draw("same");

    canvas_raw_sigma->cd(i+1);
    h1_sigma_fine[i]->SetMarkerStyle(kFullCircle);
    h1_sigma_fine[i]->SetMarkerSize(0.2);
    h1_sigma_fine[i]->SetMarkerColor(kRed);
    h1_sigma_fine[i]->SetLineColor(kRed);
    h1_sigma_fine[i]->Draw("same");
    for(int j=0; j<640; j++){
      sigma_fine[i][j]    = h1_sigma_fine[i]->GetBinContent(j+1);
      cout << "\nFOOT: " << foot_id[i] << "\tStrip: " << j  << "\tFine sigma: " << sigma_fine[i][j];
    }
  }
  //--------- Analyzing  data ------------
  foot_data data[NDETS];//collection of clusters from all detectors 
  for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
  {
    cout << "\n-- Event # : " << ev;
    ch->GetEntry(ev);
    //--------  Global base line correction in every FOOT in this event ---------
    for(int f=0; f<NDETS; f++)//loop over all foots
    {
      //  cout << "\nPerforming global baseline correction for FOOT: " << foot_id[f];
      data[f].clear(); 

      mean_ssd=0; stat=0;
      for(int i=0; i<640; i++){
	if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
	signal = FOOTE[f][i] - pedestal[f][i];
	if(fabs(signal) > (10 * sigma[f][i])) continue; //possible hist candidates
	//	cout << "\nStrip: " << FOOTI[f][i] << "for FOOT: " << foot_id [f] <" is good for a global baseline correction, Energy: " << signal;	
	stat++;    mean_ssd += signal;
      }
      mean_ssd /= stat;
      if(fabs(mean_ssd)>10){
	cout << "\n--[WARNING]: In FOOT " << foot_id[f] << "Mean ssd is " << mean_ssd << "which shows the results need calibration." << endl;
      }

      //------- Fine baseline correction for individual asics ---------
      stat=0; counter_asic=0;
      for(int i=0; i<10; i++){  asic_offset[i]=0.; }//reset asic baselines
      for(int i=0; i<640; i++)
      {
	signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd;
	if(fabs(signal) < (4 * sigma[f][i]) && //ignore possible hit candidates
	    is_good_strip(foot_id[f],FOOTI[f][i]))//and bad strips
	{
	  stat++;
	  asic_offset[counter_asic] += signal;
	}
	if((FOOTI[f][i]%64)==0){//switch to next asic
	  asic_offset[counter_asic] /= stat;
	  counter_asic++;  stat=0;
	}
      }

      //-------- baseline correction and ASIC sum------------
      counter_asic=0; signal_sum=0;
      for(int i=0; i<640; i++)
      {
	if((FOOTI[f][i]%64) == 1 && FOOTI[f][i]>1)
	{
	  signal_sum=0;
	  counter_asic++;
	}
	if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
	signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd - asic_offset[counter_asic];
	//      cout << "\nStrip: " << FOOTI[f][i] << " for FOOT: " << foot_id[f] < " has been baseline corrected including the asic sum, Energy: " << signal;
	signal_sum += signal;
      }

      //-------- Apply baseline correction, fill cluster data and histos ------------
      counter_asic=0; signal_sum=0;
      bool is_good_asic=false;
      for(int i=0; i<640; i++)
      {
	if((FOOTI[f][i]%64) == 1 && FOOTI[f][i]>1) counter_asic++;
	if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
	signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd - asic_offset[counter_asic];

	h2_cal_coarse[f]->Fill(FOOTI[f][i], (FOOTE[f][i] - pedestal[f][i]));
	h2_cal_fine[f]->Fill(FOOTI[f][i], signal);

	if(signal>(NSIGMA * sigma_fine[f][i]))
	{
	  cout << "\nStrip: " << FOOTI[f][i] << " for FOOT: " << foot_id[f] << " has passed the threshold condition, this has Energy: " << signal ;
	  Check_Strip(FOOTI[f][i], signal, data[f]);
	}
      }
    }//end loop detectors

    //-------- Inspect how many detectors observed hits
    int mul=0;
    for(int f=0; f<NDETS; f++)
    {
      if(FOOT[f]<640) continue;
      if(data[f].size()>0) mul++;
    }
    h1_det_mul->Fill(mul);

    //---------- Plot beam profile
    double X, Y, xp;
    for(auto & c0: data[0])
    {
      for(auto & c1: data[1])
      {   
	xp = get_cog(c1);
	X = (-1)*(xp*FOOT_LENGTH/640.- FOOT_LENGTH/2.);
	Y = (get_cog(c0)*FOOT_LENGTH/640.- FOOT_LENGTH/2.);
	h2_beam_XY->Fill(X,Y);
	//h2_beam_XY->Fill(get_cog(c1),get_cog(c0));
	h1_beam_X->Fill(X);
	h1_beam_Y->Fill(Y);

      }
    }
    //--------- Fill output tree 
    F15_mul = 0;
    F16_mul = 0; 
    for(int f=0; f<NDETS; f++)
    {
      for(auto & c0: data[f])
      {
	switch(f)//filling up output tree
	{
	  case 0: 
	    F15_e[F15_mul]   = get_esum(c0);
	    F15_cog[F15_mul] = get_cog(c0);
	    F15_size[F15_mul]= c0.size();
	    F15_mul++;
	    break;

	  case 1: 
	    F16_e[F16_mul]   = get_esum(c0);
	    F16_cog[F16_mul] = get_cog(c0);
	    F16_size[F16_mul]= c0.size();
	    F16_mul++;
	    break;
	}                    
	h1_cluster_e[f]->Fill( get_esum(c0) );
	h1_cluster_size[f]->Fill( c0.size() );
	h2_cluster_e_vs_cog[f]->Fill( get_cog(c0), get_esum(c0));
      }
    }
    //tree->Fill();
  }//end of eventloop

  //------- Plotting everything -------
  for(int f=0; f<NDETS; f++)
  {
    canvas_coarse->cd(f+1);
    gPad->SetLogz();
    h2_cal_coarse[f]->Draw("colz");

    canvas_fine->cd(f+1);

    h2_cal_fine[f]->Draw("colz");

    canvas_cluster_energy->cd(f+1);
    gPad->SetLogy();
    h1_cluster_e[f]->Draw();

    canvas_cluster_size->cd(f+1);
    h1_cluster_size[f]->Draw();

    canvas_cluster_e_vs_cog->cd(f+1);
    h2_cluster_e_vs_cog[f]->Draw("colz");
  }
  canvas_XY->cd(1);

  gPad->SetLogz();
  h2_beam_XY->Draw("colz");

  canvas_XY->cd(2);
  h1_beam_X->Draw();

  canvas_XY->cd(3);
  h1_beam_Y->Draw();

  canvas_det_mul->cd();
  h1_det_mul->Draw();

  tree->AutoSave();
  _file->Close();

  theApp->Run();
  return;
} 

int main(Int_t argc, Char_t* argv[])
{
  gRandom = new TRandom3();
  gRandom->SetSeed(0);
  gROOT->Macro("rootlogon.C");
  gStyle->SetPalette(kRainBow);

  TChain * ch = new TChain("h101");
  ch->Add("/u/lndgst02/william/roots/run93_unpacked.root");
  //ch->Add("test.root");
  analyse(0,-1,ch);
  //analyse(0,5e4,ch);
  return 0;
}
