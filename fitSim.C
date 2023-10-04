#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "fitAll_mod_original.h"

FILE *fFitResult;
void fitSim( int iniRing = 0, int finRing = 225, Double_t fitminNa = 250, Double_t fitmaxNa = 2500, string anelN="All", string simLife="55Ni_54Ni_6+_cascade_970.root"){

  char *simFileName[4]={(char*)simLife.c_str(),(char*)"55Ni_54Ni_4+_cascade.root",(char*)"55Ni_54Ni_2+.root",(char*)"55Ni_54Ni_511keV.root"}; // you can add as many (char*)"PATH/nome_..." as the number of peaks you want to fit
  char *expFileName[1]={(char*)Form("result.root")};
  
  //Define the binning
  
  std::string tempName;
  std::string input;
  int minBin = 0;
  int maxBin = 6000; //SP
  int binning = 10;
  int numBin = (maxBin-minBin)/binning;
  const int numsims=4; //number of peaks to be fitted
  
  //-----------------------------

  const int rebin=2; //<-------------------------------REBIN
  
  //Id starts with one!!!
  int daliIDMinNa = iniRing;
  int daliIDMaxNa = finRing; //backward crystals are from 0 to 105, forward are from 106 to 225
 
  //-----------------------------
  //Some style
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetOptLogy(0);
  TGaxis::SetMaxDigits(5);
  

  //-----------------------------------
  //Read the input files
  cout<<"Opening simulation file ..."<<endl;
  TFile *sim[numsims];
  for(int i=0;i<numsims;i++){
   
    sim[i] = TFile::Open(simFileName[i]);
    if(sim[i]==nullptr){
      std::cerr << "Error: Cannot open Simulation File: " << tempName.c_str() << ", please check it!" << std::endl;
      return;
    } else{
      std::cout<<"File "<<simFileName[i]<<" opened"<<std::endl;
    }
  }
  

  TFile *expe = TFile::Open(expFileName[0]);
  if(expe==nullptr) {
    std::cerr << "Error: Cannot open experimental data file" <<expFileName[0] <<  std::endl;
    return;
    } else{
      std::cout<<"File "<<expFileName[0]<<" opened"<<std::endl;
  }
  
  //----------------------------------------
  //Create the histograms to be read
  //Simulated peaks
  
  TH1F **hsim = new TH1F*[numsims];
  TH2F **hsim_id_doppler = new TH2F*[numsims];
  
  Double_t norm;
  for(int i=0;i<numsims;i++){
    std::cout<<"Simulation number= "<<i<<std::endl;
  
    hsim_id_doppler[i] = (TH2F*)sim[i]->Get("crystal_fired_doppler_addback");
    tempName = Form("crystal_fired_doppler_addback_py[%i]",i);
    hsim[i] = (TH1F*)hsim_id_doppler[i]->ProjectionY(tempName.c_str(),daliIDMinNa,daliIDMaxNa);
    

    norm = hsim[i]->Integral(1,600);
    
    for(int j=1;j<=600;j++){
        hsim[i]->SetBinContent(j, hsim[i]->GetBinContent(j)/norm );
        hsim[i]->SetBinError(j, hsim[i]->GetBinError(j)/norm );
        gSystem->ProcessEvents();
    }
    
    hsim[i]->Scale(1000);
    
  }

//----------------------------------------
//The experimental spectra:
  TH1F **hexp = new TH1F*[numsims];
  TH2F **hexp_id_doppler = new TH2F*[numsims];

  for(int i=0;i<numsims;i++){
 
    hexp_id_doppler[i] = (TH2F*)expe->Get("h_cut55i_cut54_AB_mle4");
    tempName = Form("h_exp_py[%i]",i);
    hexp[i] = (TH1F*)hexp_id_doppler[i]->ProjectionY(tempName.c_str(),daliIDMinNa+1,daliIDMaxNa+1);
    
    norm = hexp[i]->Integral(1,600);
    
    for(int j=1;j<=600;j++){
        hexp[i]->SetBinContent(j, hexp[i]->GetBinContent(j)/norm );
        hexp[i]->SetBinError(j, hexp[i]->GetBinError(j)/norm );
        gSystem->ProcessEvents();
    }
    
    hexp[i]->Scale(1000);
    
  }


//------------------------------------------------
//Give some nice, fancy style to the histos
for(int i=0;i<numsims;i++) 
    {
      hsim[i]->Rebin(rebin); //<-----------------------
    }
  //Experimental spectra
    
  for(int i=0;i<numsims;i++) {
    hexp[i]->Rebin(rebin);   //<-----------------------
    hexp[i]->SetStats(0);
    hexp[i]->SetFillColor(0);
    hexp[i]->SetLineColor(kBlack);
    hexp[i]->SetLineWidth(2);
    hexp[i]->GetXaxis()->SetRangeUser(0,2500);
    //hexp[i]->GetYaxis()->SetRangeUser(0,150);
    hexp[i]->GetXaxis()->SetNdivisions(7);
    hexp[i]->GetYaxis()->SetNdivisions(7);
    hexp[i]->GetYaxis()->SetTitle(Form("Counts (%3.1f keV/bin)",hexp[i]->GetBinWidth(2)));
    hexp[i]->GetXaxis()->SetTitle("Energy [keV]");
    hexp[i]->GetXaxis()->SetTitleOffset(1.2);
    hexp[i]->GetYaxis()->SetTitleOffset(1.2);
    hexp[i]->GetXaxis()->SetTitleFont(132);
    hexp[i]->GetYaxis()->SetTitleFont(132);
    hexp[i]->GetXaxis()->SetTitleSize(0.065);
    hexp[i]->GetYaxis()->SetTitleSize(0.065);
    hexp[i]->GetXaxis()->SetLabelSize(0.065);
    hexp[i]->GetYaxis()->SetLabelSize(0.065);
    hexp[i]->GetYaxis()->SetDecimals();
    hexp[i]->GetXaxis()->SetDecimals();
    //How to get error bar for each bin
    hexp[i]->SetDefaultSumw2(kTRUE);
    hexp[i]->SetTitle("");
  }
  
  //-------------------------------------------------
  //Create the graphs to fit. One per simulated peak
  peakg = new TGraph*[numsims];
  for(int j = 0 ; j < numsims ; j++)
    peakg[j] = new TGraph(hsim[j]);
  
  //-------------------------------------------
  //Define the function to fit and the parameter range
  TF1 *whole1;
  switch(numsims){
    case 1:
      whole1 = new TF1("whole1",ex_respf1,fitminNa,fitmaxNa,4+2*numsims);
      break;
    case 2:
      whole1 = new TF1("whole1",ex_respf2,fitminNa,fitmaxNa,4+2*numsims);
      break;
    case 3:
      whole1 = new TF1("whole1",ex_respf3,fitminNa,fitmaxNa,4+2*numsims);
      break;
    case 4:
      whole1 = new TF1("whole1",ex_respf4,fitminNa,fitmaxNa,4+2*numsims);
      break;
    case 5:
      whole1 = new TF1("whole1",ex_respf5,fitminNa,fitmaxNa,4+2*numsims);
      break;
    case 6:
      whole1 = new TF1("whole1",ex_respf6,fitminNa,fitmaxNa,4+2*numsims);
      break;
    default:
      std::cerr << "Number of sims " << numsims << "not foreseen (1-6)" << std::endl;
      return;
  }

  
  whole1->SetParameter(0,-0.231347918537524);
  whole1->SetParameter(2,3.37233e+00);
  //parameter 1,3 expo slopes
  whole1->SetParameter(1,2.64E-03);
  whole1->SetParameter(3,-8.18725e-04);
 
  whole1->SetParLimits(0,-2.24051e+01,0);
  whole1->SetParLimits(1,1.64614e-03,10.20576e-03);
  whole1->SetParLimits(2,1,5);
  whole1->SetParLimits(3,-2.17835e-03,0);
 
 
  for(int j = 0 ; j < numsims ; j++){
    whole1->SetParLimits(5+2*j,0,1);
    whole1->SetParLimits(4+2*j,-10,10);
  }
  whole1->SetParLimits(10,-15,15);
  
  whole1->FixParameter(4,0);
  whole1->FixParameter(6,0);
  whole1->FixParameter(8,0);

  whole1->SetLineColor(kRed);
  whole1->SetLineWidth(4);

  //-------------------------------------------
  //Fit!
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);
  TVirtualFitter::SetMaxIterations(10000);
  whole1->SetNpx(6000);
  int res = hexp[0]->Fit(whole1,"RL");
    
  cout << res << "\t" << whole1->GetNDF() << "\t"<< whole1->GetChisquare() << "\t" << whole1->GetChisquare()/whole1->GetNDF() << "\t" << whole1->GetProb() << endl;
  //---------------------------------------
  
  //Drawing the result
  TF1 **indPeaks = new TF1*[numsims];
  Int_t ColorPalette[] = {kBlue,kGreen+2, kPink, kTeal+1, kCopper, kWaterMelon};
  for(Int_t j = 0 ; j < numsims ; j++){
    indPeaks[j] = new TF1(Form("peak%if",j), "respf",150,fitmaxNa,3);
    indPeaks[j]->SetParameter(2,j );
    indPeaks[j]->SetParameter(0,whole1->GetParameter(nbExponentialPar+2*j)   );
    indPeaks[j]->SetParameter(1,whole1->GetParameter(nbExponentialPar+2*j+1) );
    indPeaks[j]->SetLineColor(ColorPalette[j]);
    indPeaks[j]->SetLineWidth(3);
    indPeaks[j]->SetLineStyle(2);
    indPeaks[j]->Draw("hist same");
    indPeaks[j]->GetYaxis()->SetDecimals();
    indPeaks[j]->GetXaxis()->SetDecimals();
    indPeaks[j]->SetNpx(6000);
    
  }

  TF1 *expon1= new TF1( "expon1",expf1,150,fitmaxNa,4);
  expon1->SetParameters(whole1->GetParameter(0),
                        whole1->GetParameter(1),
                        whole1->GetParameter(2),
                        whole1->GetParameter(3));
  expon1->SetLineColor(kBlack);
  expon1->SetLineWidth(3);
  expon1->SetLineStyle(7);
  //hexp[0]->Draw("hist E");
  hexp[0]->Draw("hist E same");
  expon1->Draw("same");
  whole1->Draw("same");
  
  
  
  auto legend = new TLegend(0.38,0.6,0.88,0.9);
  legend->AddEntry((TObject*)0, Form("%s",simFileName[0]));
  legend->AddEntry((TObject*)0, Form("%s",anelN.c_str()));
  legend->AddEntry((TObject*)0, Form("%d  %d  %.1f",res,whole1->GetNDF(),whole1->GetChisquare()) , "");
  legend->Draw();
  
}

void rings(string simLife="55Ni_54Ni_6+_cascade_970.root"){
    TCanvas *c1 = new TCanvas("c1","bgLeft, peak, bgRight",1400,1400);
  c1->Divide(3,3);
    c1->cd(1);
    fitSim( 208, 227, 370, 2800, "114-122deg",simLife.c_str());
    c1->Update();
    gSystem->ProcessEvents();
    c1->cd(2);
    fitSim(  163, 207, 370, 2700, "93-102deg",simLife.c_str());
    c1->cd(3);
    fitSim( 135, 162, 370, 2700, "85deg",simLife.c_str());
    c1->cd(4);
    fitSim( 107, 134, 370, 2700, "77deg",simLife.c_str());
    c1->cd(5);
    fitSim( 87, 106, 370, 2700, "69deg",simLife.c_str());
    c1->cd(6);
    fitSim( 67, 86, 370, 2700, "59deg",simLife.c_str());
    c1->cd(7);
    fitSim( 47, 66, 370, 2700, "48deg",simLife.c_str());
    c1->cd(8);
    fitSim( 20, 46, 370, 3000, "28-35deg",simLife.c_str());
    c1->cd(9);
    fitSim( 0, 19, 370, 2700, "16-22deg",simLife.c_str());
    
}


