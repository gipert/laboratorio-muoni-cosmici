// MonteCarlo per la misura della vita media dei muoni cosmici in alluminio
//
// Authors: Mattia Faggin, Davide Piras, Luigi Pertoldi
//
// Prova:
// ./montecarlo --help
//

// nel codice utilizziamo 3 metodi per stimare la baseline: 1-media, 2-fit likelihood, 3-fit chi2

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <chrono>

#include "TApplication.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGraphErrors.h"

#include "../../ProgressBar/progressbar.h"

// questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
#define	Begin     18	// inizio istogramma
#define StartBase 2600	// punto dell'istogramma in cui comincia la baseline
#define End       3904	// fine istogramma
#define Nsim      500   // numero simulazioni
#define beginFit  18	// inizio fit esponenziale

TRandom3 r;
int counts;

struct distRange{
    double minimo;
    double massimo;
};
    
void fitGaus(TH1D& h);
distRange getRange(std::vector<double> v);

//--------------- main ---------------------
int main( int argc, char* argv[] ) {

    // salvo gli argomenti prima di darli in pasto al TApplication
    // perchè lui li modifica
    std::vector<std::string> args;
    args.reserve(argc);
    for ( int i = 0; i < argc; i++ ) args.push_back(argv[i]);

    if ( argc == 2 && args[1] == "--help" ) {
        std::cout << std::endl
                  << "MonteCarlo per la misura della vita media dei muoni cosmici in alluminio." << std::endl
                  << "Autori: Mattia Faggin, Davide Piras, Luigi Pertoldi" << std::endl << std::endl
                  << "Utilizzo:" << std::endl 
                  << "    $ ./montecarlo [Baseline] [tau] [Integrale] [rebinFactor]" << std::endl << std::endl;
        return 0;
    }

    if ( argc < 3 ) {
        std::cout << "Pochi argomenti! Se non ti ricordi c'è l'opzione '--help'" << std::endl
                  << "Termino l'esecuzione..." << std::endl;
        return 0;
    }

    float B       = std::stof(args[1]);
    double tau	  = std::stoi(args[2]);
    double integrale  = std::stoi(args[3]);
    int RebFactor = std::stoi(args[4]);

    double A = (-integrale/tau)*1/(exp(-End/tau)-exp(Begin/tau));

    TApplication Root("App",&argc,argv); 
    gStyle->SetOptStat(0);

    // METODO 1
    std::vector<double> vFitTauL;
    vFitTauL.reserve(Nsim);
    std::vector<double> vFitAL;
    vFitAL.reserve(Nsim);
    std::vector<double> vFitBL;
    vFitBL.reserve(Nsim);
    std::vector<double> vFitErrTauL;
    vFitErrTauL.reserve(Nsim);
    std::vector<double> vFitErrAL;
    vFitErrAL.reserve(Nsim);
    std::vector<double> vFitErrBL;
    vFitErrBL.reserve(Nsim);

    // METODO 2
    std::vector<double> vFitTauL2;
    vFitTauL2.reserve(Nsim);
    std::vector<double> vFitAL2;
    vFitAL2.reserve(Nsim);
    std::vector<double> vFitBL2;
    vFitBL2.reserve(Nsim);
    std::vector<double> vFitErrTauL2;
    vFitErrTauL2.reserve(Nsim);
    std::vector<double> vFitErrAL2;
    vFitErrAL2.reserve(Nsim);
    std::vector<double> vFitErrBL2;
    vFitErrBL2.reserve(Nsim);

    // simulazione baseline+exp
    TH1D total("total","total",4096,0,4096);

    TF1 fitFunc("fitFunc","[0]*TMath::Exp(-x/[1])+[2]",beginFit,End);
    fitFunc.SetParName(0,"A");
    fitFunc.SetParName(1,"tau");
    fitFunc.SetParName(2,"B");
    TF1 fitFunc2("fitFunc2","[0]*TMath::Exp(-x/[1])+[2]",beginFit,End);
    fitFunc2.SetParName(0,"A");
    fitFunc2.SetParName(1,"tau");
    fitFunc2.SetParName(2,"B");

    std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)";
    TCanvas can( "METODO1", (canName + " METODO 1").c_str(), 1200 , 700 );
    can.Divide(3,2);
    can.cd(1);

    // ciclo delle simulazioni
    ProgressBar bar(Nsim);
    bar.Init();
    for(int k=0; k<Nsim; k++)
    {
	    bar.Update(k);
        r.SetSeed(0);
    	// simulazione della baseline
    	TH1D baseline("baseline","baseline",4096,0,4096);
    	for ( int k = 0; k < B*(End-Begin); k++ )
    	{
        	baseline.Fill(r.Uniform(Begin,End));
    	}
    	// rebin 
    	baseline.Rebin(RebFactor);
    	//simulazione esponenziale
       	TH1D exponential("exponential","exponential",4096,0,4096);
       	double val = 0;
        for ( int k = 0; k < A*tau*exp(-Begin/tau); k++ )
    	{
            val = r.Exp(tau);
            while (val<Begin || val>End) val = r.Exp(tau);
		    exponential.Fill(val);
    	}
    	// rebin 
    	exponential.Rebin(RebFactor);
 
    	total.Add(&baseline,&exponential);
	//total.Draw();
        // METODO 1: funzione completa
        // parameter setting
	//fitFunc.SetParameter(0,A);
	fitFunc.SetParameter(1,tau);
	//fitFunc.SetParameter(2,B);

///////////////// FIT UNICO //////////////////// 
      total.Fit("fitFunc","RLQN");
/////////////////////////////////////////////////        
        
        vFitAL.push_back(fitFunc.GetParameter("A"));
        vFitBL.push_back(fitFunc.GetParameter("B"));
        vFitTauL.push_back(fitFunc.GetParameter("tau"));
        vFitErrAL.push_back(fitFunc.GetParError(0));
        vFitErrTauL.push_back(fitFunc.GetParError(1));
        vFitErrBL.push_back(fitFunc.GetParError(2));
        //std::cout<<"\nB = " << fitFunc.GetParameter("B") << " +- " << fitFunc.GetParError(2);

/////////// METODO 2: pol0 per B, quindi B-fixing e fit con fitFunc
	TFitResultPtr basePtr = total.Fit("pol0","LSNQ","",StartBase,End);
        fitFunc2.SetParameter(2,basePtr->Parameter(0));
        fitFunc2.SetParameter(1,tau);

/////////// FIT RESIDUO ////////////////
        total.Fit("fitFunc2","RLQN","",beginFit,End);
        vFitBL2.push_back(fitFunc2.GetParameter(2));
        vFitErrBL2.push_back(fitFunc2.GetParError(2));
        vFitAL2.push_back(fitFunc2.GetParameter("A"));
        vFitTauL2.push_back(fitFunc2.GetParameter("tau"));
        vFitErrAL2.push_back(fitFunc2.GetParError(0));
        vFitErrTauL2.push_back(fitFunc2.GetParError(1));
    }

    // range distribuzioni METODO 1
    distRange rFitTauL    = getRange(vFitTauL);
    distRange rFitAL      = getRange(vFitAL);
    distRange rFitBL      = getRange(vFitBL);
    distRange rFitErrTauL = getRange(vFitErrTauL);
    distRange rFitErrAL   = getRange(vFitErrAL);
    distRange rFitErrBL   = getRange(vFitErrBL);

    // range distribuzioni METODO 2
    distRange rFitTauL2    = getRange(vFitTauL2);
    distRange rFitAL2      = getRange(vFitAL2);
    distRange rFitBL2      = getRange(vFitBL2);
    distRange rFitErrTauL2 = getRange(vFitErrTauL2);
    distRange rFitErrAL2   = getRange(vFitErrAL2);
    distRange rFitErrBL2   = getRange(vFitErrBL2);

    // creo gli istogrammi delle distribuzioni METODO 1
    TH1D hFitTauL 	    ( "hFitTauL" 	, "#tau", 100, rFitTauL.minimo, rFitTauL.massimo );
    TH1D hFitAL   	    ( "hFitAL"   	, "A"   , 100, rFitAL.minimo, rFitAL.massimo );
    TH1D hFitBL   	    ( "hFitBL"   	, "B"   , 100, rFitBL.minimo, rFitBL.massimo );
    TH1D hFitErrBL 	    ( "hFitErrBL"   , "Errori #sigma_{B}"   , 100, rFitErrBL.minimo, rFitErrBL.massimo); 
    TH1D hFitErrTauL 	( "hFitErrTauL" , "Errori #sigma_{#tau}", 100, rFitErrTauL.minimo, rFitErrTauL.massimo); 
    TH1D hFitErrAL   	( "hFitErrAL"   , "Errori #sigma_{A}"   , 100, rFitErrAL.minimo, rFitErrAL.massimo); 

    // creo gli istogrammi delle distribuzioni METODO 2
    TH1D hFitTauL2 	    ( "hFitTauL2" 	 , "#tau", 100, rFitTauL2.minimo, rFitTauL2.massimo );
    TH1D hFitAL2   	    ( "hFitAL2"   	 , "A"   , 100, rFitAL2.minimo, rFitAL2.massimo );
    TH1D hFitBL2   	    ( "hFitBL2"   	 , "B"   , 100, rFitBL2.minimo, rFitBL2.massimo );
    TH1D hFitErrBL2 	( "hFitErrBL2"   , "Errori #sigma_{B}"   , 100, rFitErrBL2.minimo, rFitErrBL2.massimo); 
    TH1D hFitErrTauL2 	( "hFitErrTauL2" , "Errori #sigma_{#tau}", 100, rFitErrTauL2.minimo, rFitErrTauL2.massimo); 
    TH1D hFitErrAL2   	( "hFitErrAL2"   , "Errori #sigma_{A}"   , 100, rFitErrAL2.minimo, rFitErrAL2.massimo); 

    // riempio gli istogrammi METODO 1
    for ( int i = 0; i < Nsim; i++ ) {        
        hFitTauL.Fill(vFitTauL.at(i));
        hFitAL.Fill(vFitAL.at(i));
        hFitBL.Fill(vFitBL.at(i));
	    hFitErrTauL.Fill(vFitErrTauL.at(i));
	    hFitErrAL.Fill(vFitErrAL.at(i));
        hFitErrBL.Fill(vFitErrBL.at(i));
    }
    // riempio gli istogrammi METODO 2
    for ( int i = 0; i < Nsim; i++ ) {        
        hFitTauL2.Fill(vFitTauL2.at(i));
        hFitAL2.Fill(vFitAL2.at(i));
        hFitBL2.Fill(vFitBL2.at(i));
	    hFitErrTauL2.Fill(vFitErrTauL2.at(i));
	    hFitErrAL2.Fill(vFitErrAL2.at(i));
        hFitErrBL2.Fill(vFitErrBL2.at(i));
    }


    std::cout << std::endl 
	      << "---------------------------------------------------------------------------------" 
	      << "\n./montecarlo.out " << B << " " << tau << " " << integrale << " " << RebFactor <<std::endl 
	      << "\nEventi totali   " << integrale
	      << "\nEventi baseline " << B*(End-Begin)
	      << "\nEventi exp      " << integrale - B*(End-Begin)
	      << "\nValori in INPUT:"
	      << "\n	Baseline  " << B
	      << "\n	tau       " << tau
	      << "\n	Integrale " << integrale
	      << "\n	RebFactor " << RebFactor << std::endl
          << "\nSimulazioni MC effettuate:                         "  << Nsim
	      << "\nIntervallo di generazione:                         [" << Begin     << ", " << End << "]"
          << "\nIntervallo di fit METODO 1:                        [" << beginFit  << ", " << End << "]" 
	      << "\nIntervallo di fit baseline METODO 2:               [" << StartBase << ", " << End << "]"
	      << "\nIntervallo di fit complessivo (B fixed) METODO 2:  [" << beginFit  << ", " << End << "]" << std::endl;

    // fit gaussiani delle distribuzioni METODO 1
    std::cout<<"\n METODO 1" << std::endl;
    fitGaus(hFitTauL);
    fitGaus(hFitAL);
    fitGaus(hFitBL);
    fitGaus(hFitErrTauL);
    fitGaus(hFitErrAL);
    fitGaus(hFitErrBL);

        hFitTauL.Draw();
        TLine lTau(tau,0,tau,hFitTauL.GetMaximum());
	lTau.SetLineColor(kRed);
        lTau.Draw();
    can.cd(2);
        hFitAL.Draw();
        TLine lA(A,0,A,hFitAL.GetMaximum());
	lA.SetLineColor(kRed);
        lA.Draw();
    can.cd(3);
        hFitBL.Draw();
        TLine lB(B,0,B,hFitBL.GetMaximum());
	lB.SetLineColor(kRed);
        lB.Draw();
    can.cd(4);
        hFitErrTauL.Draw();
    can.cd(5);
        hFitErrAL.Draw();
    can.cd(6);
        hFitErrBL.Draw();

    // fit gaussiani delle distribuzioni METODO 2
    std::cout<<"\n METODO 2" << std::endl;
    fitGaus(hFitTauL2);
    fitGaus(hFitAL2);
    fitGaus(hFitBL2);
    fitGaus(hFitErrTauL2);
    fitGaus(hFitErrAL2);
    fitGaus(hFitErrBL2);
    TCanvas can2( "METODO2", (canName + " METODO 2").c_str(), 1200 , 700 );
    can2.Divide(3,2);
    can2.cd(1);
        hFitTauL2.Draw();
        TLine lTau2(tau,0,tau,hFitTauL2.GetMaximum());
	    lTau2.SetLineColor(kGreen+3);
        lTau2.SetLineWidth(3);
        lTau2.Draw();
    can2.cd(2);
        hFitAL2.Draw();
        TLine lA2(A,0,A,hFitAL2.GetMaximum());
	    lA2.SetLineColor(kGreen+3);
        lA2.SetLineWidth(3);
        lA2.Draw();
    can2.cd(3);
        hFitBL2.Draw();
        TLine lB2(B,0,B,hFitBL2.GetMaximum());
	    lB2.SetLineColor(kGreen+3);
        lB2.SetLineWidth(3);
        lB2.Draw();
    can2.cd(4);
        hFitErrTauL2.Draw();
    can2.cd(5);
        hFitErrAL2.Draw();
    can2.cd(6);
        hFitErrBL2.Draw();

    Root.Run();

    return 0;
}
//--------------- fine main ----------------------


void fitGaus(TH1D& h){
   int dim = h.GetNbinsX();

   std::string name = h.GetName();
   std::cout << std::endl << name <<std::endl;
   if (h.GetBinContent(0)>10 || h.GetBinContent(dim+1)>10)
   {
	std::cout <<"Underflow " << h.GetBinContent(0)
		  <<"\nOverflow  " << h.GetBinContent(dim+1)
		  <<"\nTroppi dati fuori range: fit non fatto" << std::endl;;
	return;
   }
   double min = h.GetBinCenter(1);
   double max = h.GetBinCenter(dim);
   TFitResultPtr p = h.Fit("gaus","SQ","",min,max);
   double mean = p->Parameter(1);
   double err  = p->ParError(1);
   double sigma = p->Parameter(2);
   double errsigma = p->ParError(2);

   std::cout << "mean = "       << mean  << " +- " << err << std::endl
	         << "sigma = "      << sigma << " +- "<< errsigma << std::endl
             << "mean/sigma = " << "("   << sigma*100/mean << "%)" << std::endl;
   return;
}


distRange getRange(std::vector<double> v){
   double min = 100000;
   double max = 0;
   double entry;
   int N = v.size();
   for(int k=0; k<N; k++)
   {
        entry = v.at(k);
	if(min==0 || min>entry )	min = entry;
	if(max==0 || max<entry )	max = entry;
   }
   //auto min = std::min_element(v.begin(),v.end());
   //auto max = std::max_element(v.begin(),v.end());
   distRange d {min,max};
   return d;
}
