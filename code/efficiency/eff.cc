//
// ROOT's Macro to plot toghether two graphs with different y-scales
//
// Input file format:
//  PMTnr
//  Nnew
//  WorkHV
//  HV	PMT Tnear	T&near  Tfar    T&far
//
// Luigi Pertoldi, 13/04/2016
//
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

void eff( std::string file_name ) {
    
    int time = 500;                         // acquisition time (s)
    int n = 0;                              // total number of points
    int PMTnr;                              // PMT number
    int Nnew = 0;                           // number of additional points (far)
    float WorkHV;
    ifstream file(file_name.c_str());       // input file
	
    float HV      [40];                     // voltage
    float PMT     [40];                     // counts
    float TN      [40];                     // telescope near PMT
    float TandN   [40];                     // telescope near PMT & PMT (coincidence)
    float effN    [40];                     // efficiency near PMT
    float TF      [40];                     // telescope far from PMT
    float TandF   [40];                     // telescope far from PMT & PMT (coincidence)
    float effF    [40];                     // efficiency far from PMT
    
    float HVN     [40];
    float HVF     [40];
    float HVc     [40];
    
    // errors
    float PMTErr  [40];
    float TErrN   [40];
    float TandErrN[40];
    float effErrN [40];
    float TErrF   [40];
    float TandErrF[40];
    float effErrF [40];
    
    float nullErr [40];
    
    // retrieve PMT number for labels
    file >> PMTnr;
    ostringstream convert;
    convert << PMTnr;
    string PMTname = convert.str();
    
    file >> Nnew;
    file >> WorkHV;
    
    // initialization
    for ( int i = 0 ; i < 40 ; i++ ) PMT[i] = TN[i]= TandN[i] = TF[i] = TandF[i] = HVN[i] = HVF[i] = HVc[i] = 1;
    
    // retrieve data from file
    int i_N = 0;
    int i_F = 0;
    int i_c = 0;
    while (  file >> HV[n] >> PMT[i_c] >> TN[i_N] >> TandN[i_N] >> TF[i_F] >> TandF[i_F]  ) {
	
	if ( PMT[i_c] != 1 )                    { HVc[i_c] = HV[n]; i_c++; }
        if ( TN[i_N]  != 1 || TandN[i_N] != 1 ) { HVN[i_N] = HV[n]; i_N++; }
        if ( TF[i_F]  != 1 || TandF[i_F] != 1 ) { HVF[i_F] = HV[n]; i_F++; }
        n++;	
    }
    
    // find max counts for PMT (for axis ranges)
    int max = 0;
    int maxHV = 0, minHV = 20000;
    for ( int i = 0 ; i < n ; i++ ) {
        
        if ( HV[i] > maxHV ) {
            maxHV = HV[i];
            do maxHV++; while ( maxHV % 50 != 0 );  // find a nice range
        }
        if ( HV[i] < minHV ) {
            minHV = HV[i];
            do minHV--; while ( minHV % 50 != 0 );  // find a nice range
        }
        if ( i < i_c && PMT[i] > max ) max = PMT[i];
    }
    
    // compute errors & fill arrays
    for ( int i = 0 ; i < n ; i++ ) {
        
        effN[i] = TandN[i] *1./ TN[i];
        effF[i] = TandF[i] *1./ TF[i];
        
        nullErr[i] = 0;
        
        TErrN[i]    = TMath::Sqrt(TN[i]);
        TErrF[i]    = TMath::Sqrt(TF[i]);
        TandErrN[i] = TMath::Sqrt(TandN[i]);
        TandErrF[i] = TMath::Sqrt(TandF[i]);
        //effErrN[i]  = TMath::Sqrt(effN[i]*(1-effN[i])/TN[i]);
            effErrN[i] = TMath::Sqrt(effN[i]/TN[i]);
        //effErrF[i]  = TMath::Sqrt(effF[i]*(1-effF[i])/TF[i]); 
            effErrF[i] = TMath::Sqrt(effF[i]/TF[i]);
        
        PMTErr[i] = TMath::Sqrt(PMT[i]) *1./ time;  // -> frequencies
        PMT[i] = PMT[i]*1./ time;
    }
    
	// create objects for graphs
    TGraphErrors* countsPMT  = new TGraphErrors( i_c      , HVc , PMT  , nullErr , PMTErr  );
    TGraphErrors* effNGr     = new TGraphErrors( i_N      , HVN , effN , nullErr , effErrN );
    TGraphErrors* effFGr     = new TGraphErrors( i_F-Nnew , HVF , effF , nullErr , effErrF );
    
    TGraphErrors* effFGrnew  = new TGraphErrors( Nnew );
    
    for ( int i = i_F - Nnew ; i < i_F ; i++ ) {
        
        effFGrnew->SetPoint( i - i_F + Nnew , HVF[i]     , effF[i]    );
        effFGrnew->SetPointError( i - i_F + Nnew , nullErr[i] , effErrF[i] );
    }
    
    countsPMT->Sort();
    effNGr   ->Sort();
    effFGr   ->Sort();
    
    // set style
    countsPMT->SetName("countsPMT");
    countsPMT->SetTitle("Counts w.r.t. HV");
    countsPMT->SetMarkerColor(kBlack);
    countsPMT->SetMarkerStyle(20);
    countsPMT->SetMarkerSize(1);
    countsPMT->SetLineStyle(2);
    countsPMT->SetLineColor(kBlack);
    
    effNGr->SetName("effNGr");
    effNGr->SetTitle("Efficiency Near w.r.t. HV");
    effNGr->SetMarkerColor(kRed);
    effNGr->SetMarkerStyle(23);
    effNGr->SetMarkerSize(1.3);
    effNGr->SetLineColor(kRed);
    
    effFGr->SetName("effFGr");
    effFGr->SetTitle("Efficiency Far w.r.t. HV");
    effFGr->SetMarkerColor(kBlue);
    effFGr->SetMarkerStyle(22);
    effFGr->SetMarkerSize(1.3);
    effFGr->SetLineColor(kBlue);
    
    effFGrnew->SetMarkerStyle(22);
    effFGrnew->SetMarkerSize(1.3);
    effFGrnew->SetMarkerColor(kGreen+2);
    effFGrnew->SetLineColor(kGreen+2);
    
    ///////////////// DRAW //////////////////
    
    // draw efficiency graph
    TCanvas* can = new TCanvas( "can" , "PMT Efficiency and Counts" , 0 , 500 , 1050 , 725 );
    //TCanvas* can = new TCanvas( "can" , "PMT Efficiency and Counts" , 0 , 500 , 650 , 725 );

    //TCanvas* can = new TCanvas( "can" , "PMT Efficiency and Counts" ,1000 , 500 , 2050 , 725 );	// AppleTV extended screen

    TPad* pad = new TPad( "pad" , "" , 0 , 0 , 1 , 1 );
    pad->SetTitle("PMT Efficiency and Counts");
    pad->SetGrid();
    pad->Draw();
    pad->cd();
    
    TH1F* hr = pad->DrawFrame( minHV , 0 , maxHV , 1.2 ); // X,Y ranges for efficiency
    hr->SetXTitle("HV [V]");
    hr->SetYTitle("Efficiency");
    pad->GetFrame()->SetBorderSize(12);
    
    effNGr->Draw("P");
    effFGr->Draw("PSAME");
    effFGrnew->Draw("PSAME");

    // draw counts pad
    can->cd();
    TPad* overlay = new TPad( "overlay" , "" , 0 , 0 , 1 , 1 );
    overlay->SetTitle("PMT Efficiency and Counts");
    overlay->SetFillStyle(4000);
    overlay->SetFillColor(0);
    overlay->SetFrameFillStyle(4000);
    overlay->Draw();
    overlay->cd();
    double xmin = pad->GetUxmin();
    double xmax = pad->GetUxmax();
    double ymin = 0;
    double ymax = (int)(max*1./500)+100;
    TH1F* hframe = overlay->DrawFrame( xmin , ymin , xmax , ymax );
    hframe->GetXaxis()->SetLabelOffset(99);
    hframe->GetYaxis()->SetLabelOffset(99);
    hframe->GetYaxis()->SetTickSize(0.0);
    
    countsPMT->Draw("PXC");
    
    // draw counts axis
    TGaxis* axis = new TGaxis( xmax , ymin , xmax , ymax , ymin , ymax , 510 , "+L" );
    axis->SetLabelSize(0.035);
    axis->SetTitleSize(0.035);
    axis->SetLabelFont(42);
    axis->SetTitleFont(42);
    axis->SetTitle("Counts [Hz]");
    axis->Draw("P");
    
    // set legend
    TLegend* leg = new TLegend( 0.66 , 0.15 , 0.86 , 0.27 );
    leg->AddEntry( effNGr , "Efficiency Near", "P");
    leg->AddEntry( effFGr , "Efficiency Far", "P");
    leg->AddEntry( countsPMT , "Counts", "P");
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0,0.7);
    leg->SetTextSize(0.035);
    leg->Draw();
    
    // title
    string title = "Efficiency and Counts for PMT" + PMTname;
    TText* t = new TText( 0.5 , 0.93 , title.c_str());
    t->SetNDC();
    t->SetTextAlign(21);
    t->Draw();

    // line with working voltage
    TLine * line = new TLine( WorkHV , 0 , WorkHV , ymax );
	line->SetLineColor(kGreen+3);
	line->SetLineWidth(3);
	line->SetLineStyle(10);
	line->Draw();
    
    // save .pdf
    string fileTitle = "eff" + PMTname + ".pdf";
    can->SaveAs(fileTitle.c_str());
    
    return; 
}
