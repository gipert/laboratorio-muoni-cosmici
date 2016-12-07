// MonteCarlo per la misura della vita media dei muoni cosmici in alluminio
//
// Authors: Mattia Faggin, Davide Piras, Luigi Pertoldi
//
// Prova:
// ./montecarlo --help
//

// nel codice generiamo due esponenziali e la baseline, cercando quindi il metodo migliore per stimare i loro parametri

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
#include "TGraphErrors.h"

#include "../../ProgressBar/progressbar.h"

// questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
#define	Begin     0	// inizio istogramma
#define StartBase 2600	// punto dell'istogramma in cui comincia la baseline --> FISSATO dai test su fit pol0 della baseline (baselineStart.cc)
#define End       3904	// fine istogramma
#define Nsim      100   // numero simulazioni
#define beginFit  860	// inizio fit esponenziale dei muoni lunghi


TRandom3 r;
int counts;

struct distRange{
    double minimo;
    double massimo;
};
    
struct dataBase {
    double media;
    double fitL;
    double fitC;

    double errMedia;
    double errFitL;
    double errFitC;

    TH1D hBase;
};

struct dataExp{
    double tauL;
    double tauC;
    double AL;
    double AC;

    double errTauL;
    double errTauC;
    double errAL;
    double errAC;

    TH1D* hExp;
};

dataBase simulateBase( float B, int rebin );
dataExp  simulateExp( double tau, double A, int rebin );
void fitGaus(TH1D h);
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
                  << "    $ ./montecarlo [Baseline] [tau] [Integrale] [taucorto] [R] [rebinFactor] " << std::endl << std::endl;
        return 0;
    }

    if ( argc < 7 ) {
        std::cout << "Pochi argomenti! Se non ti ricordi c'è l'opzione '--help'" << std::endl
                  << "Termino l'esecuzione..." << std::endl;
        return 0;
    }

    float B          = std::stof(args[1]); // valore tipico per 1 settimana: 2 (circa)
    double tau	     = std::stof(args[2]); // valore vero = 429 canali
	double integrale = std::stof(args[3]); // valore tipico per 1 settimana: 70000 (circa)
	double taucorto  = std::stof(args[4]); // valore vero = 172 canali (circa)
	double R         = std::stof(args[5]); // valore vero: 1.261
	int RebFactor    = std::stoi(args[6]);
	double A         = (integrale - B*(End - Begin)) / (tau + taucorto / R);
	double Aminus    = A/R;	

    TApplication Root("App",&argc,argv); 
    gErrorIgnoreLevel = kError; // toglie i warning

    // METODO 2
    std::vector<double> vFitTauL2;
    vFitTauL2.reserve(Nsim);
    std::vector<double> vFitTauShortL2;
    vFitTauShortL2.reserve(Nsim);
    std::vector<double> vFitBL2;
    vFitBL2.reserve(Nsim);
    std::vector<double> vFitRL2;
    vFitRL2.reserve(Nsim);
    
    std::vector<double> vFitErrTauL2;
    vFitErrTauL2.reserve(Nsim);
    std::vector<double> vFitErrTauShortL2;
    vFitErrTauShortL2.reserve(Nsim);
    std::vector<double> vFitErrBL2;
    vFitErrBL2.reserve(Nsim);
    std::vector<double> vFitErrRL2;
    vFitErrRL2.reserve(Nsim);

    // simulazione baseline+exp
    TH1D total("total","total",4096,0,4096);
	TH1D total2("total2", "total2", 4096, 0, 4096);

    TF1 fitFunc("fitFunc","[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])+[4]",Begin,End);
    fitFunc.SetParName(0,"Aminus");
    fitFunc.SetParName(1,"tauShort");
    fitFunc.SetParName(2,"Aplus");
    fitFunc.SetParName(3,"tauLong");
    fitFunc.SetParName(4,"B");
    TF1 fitFunc2("fitFunc2","[0]*TMath::Exp(-x/[1])+[2]",beginFit,End);
    fitFunc2.SetParName(0,"A");
    fitFunc2.SetParName(1,"tau");
    fitFunc2.SetParName(2,"B");

    // ProgressBar
    ProgressBar bar(Nsim);
    bar.Init();
    // ciclo delle simulazioni
    double val      = 0;
    double entina   = 0;    // gigi, digita "Valentina Nappi" su google e vedi cosa esce...
    double ratio    = 0;
    double errRatio = 0;
    for(int k=0; k<Nsim; k++)
    {

        bar.Update(k);
	    r.SetSeed(k+1);
    	// simulazione della baseline
    	TH1D baseline("baseline","baseline",4096,0,4096);
    	for ( int k = 0; k < B*(End-Begin); k++ )
    	{
        	baseline.Fill(r.Uniform(Begin,End));
    	}
    	// rebin 
    	baseline.Rebin(RebFactor);
    	//simulazione primo esponenziale
   	    TH1D exponential("exponential","exponential",4096,0,4096);
    	for ( int k = 0; k < A*tau; k++ )
    	{
    	val = r.Exp(tau);
    	while (val<Begin || val>End) val = r.Exp(tau);
		exponential.Fill(val);
    	}
    	// rebin 
    	exponential.Rebin(RebFactor);
 
        // somma istogrammi exp+ e baseline
		total.Add(&baseline, &exponential);

		//simulazione secondo esponenziale
		TH1D exponential2("exponential2", "exponential2", 4096, 0, 4096);
		for (int k = 0; k < A*taucorto/R; k++)
		{
		    entina = r.Exp(taucorto);
		    while (entina<Begin || entina>End) entina = r.Exp(taucorto);
			exponential2.Fill(entina);
		}
		// rebin 
		exponential2.Rebin(RebFactor);
		
		// somma istogrammi exp+ e baseline (già in total) e exp-
		total2.Add(&total, &exponential2);

        // METODO 2: pol0 per B, quindi B-setting e fit con fitFunc2, quindi fit complessivo
	    TFitResultPtr basePtr = total2.Fit("pol0","LSNQ","",StartBase,End);
        fitFunc2.SetParameter(2,basePtr->Parameter(0));
        fitFunc2.SetParameter(1,tau);

        // FIT (exp+)+Base direttamente su total2 con fitfunc2 (solo 1 exp) 
        total2.Fit("fitFunc2","LQN","",beginFit,End);
        fitFunc.SetParameter("tauShort",taucorto);
        //fitFunc.SetParameter("Aplus",fitFunc2.GetParameter("A"));
        fitFunc.SetParameter("tauLong",fitFunc2.GetParameter("tau"));
        fitFunc.SetParameter("B",fitFunc2.GetParameter("B"));
        // FIT totale su total2 con fitfunc (ambedue le exp)
        total2.Fit("fitFunc","LQN","",Begin+1,End);
        
        /*std::cout << "\nfitFunc2.GetParameter(\"A\") (1o fit)     = " << fitFunc2.GetParameter("A")
                  << "\nfitFunc.GetParameter(\"Aplus\")           = " << fitFunc.GetParameter("Aplus")
                  << "\nfitFunc.GetParameter(\"Aminus\")          = " << fitFunc.GetParameter("Aminus") << std::endl;*/ 
        
        ratio    = fitFunc.GetParameter("Aplus")/fitFunc.GetParameter("Aminus");
        errRatio = sqrt( pow(fitFunc.GetParError(2),2)+pow(ratio*fitFunc.GetParError(0),2) )/fitFunc.GetParameter("Aminus");
        
        vFitTauShortL2.push_back(fitFunc.GetParameter("tauShort"));
        vFitErrTauShortL2.push_back(fitFunc.GetParError(1)); 
        vFitTauL2.push_back(fitFunc.GetParameter("tauLong"));
        vFitErrTauL2.push_back(fitFunc.GetParError(3));       
        vFitBL2.push_back(fitFunc.GetParameter("B"));
        vFitErrBL2.push_back(fitFunc.GetParError(4));
        vFitRL2.push_back(ratio);
        vFitErrRL2.push_back(errRatio);
        //std::cout << "\ntau primo fit =" << fitFunc2.GetParameter("tau") << std::endl;
    }
    
    //histo drawing
	total2.Draw();
	fitFunc.Draw("same");
	//std::cout << "\ntau primo fit =" << fitFunc2.GetParameter("tau") << std::endl;
	/*std::cout << "Aminus =" << fitFunc.GetParameter("Aminus")
	          << "tau-   =" << fitFunc.GetParameter("tauShort")
	          << "Aplus  =" << fitFunc.GetParameter("Aplus")
	          << "tau+   =" << fitFunc.GetParameter("tauLong")
	          << "B      =" << fitFunc.GetParameter("B") << std::endl;*/

    // range distribuzioni METODO 2
    distRange rFitTauL2         = getRange(vFitTauL2);
    distRange rFitTauShortL2    = getRange(vFitTauShortL2);
    distRange rFitBL2           = getRange(vFitBL2);
    distRange rFitRL2           = getRange(vFitRL2);
    
    distRange rFitErrTauL2          = getRange(vFitErrTauL2);
    distRange rFitErrTauShortL2     = getRange(vFitErrTauShortL2);
    distRange rFitErrBL2            = getRange(vFitErrBL2);
    distRange rFitErrRL2            = getRange(vFitErrRL2); 

    // creo gli istogrammi delle distribuzioni METODO 2
    TH1D hFitTauL2 	        ( "hFitTauL2" 	        , "Likelihood (#tau_{+})"   , 100, rFitTauL2.minimo, rFitTauL2.massimo );
    TH1D hFitTauShortL2     ( "hFitTauShortL2"      , "Likelihood (#tau_{-})"   , 100, rFitTauShortL2.minimo, rFitTauShortL2.massimo );
    TH1D hFitBL2   	        ( "hFitBL2"   	        , "Likelihood (B)"          , 100, rFitBL2.minimo, rFitBL2.massimo );
    TH1D hFitRL2            ( "hFitRL2"             , "Likelihood (R)"          , 100, rFitRL2.minimo, rFitRL2.massimo );   
 
    TH1D hFitErrTauL2 	    ( "hFitErrTauL2"        , "Likelihood (#tau-{+})"   , 100, rFitErrTauL2.minimo, rFitErrTauL2.massimo); 
    TH1D hFitErrTauShortL2  ( "hFitErrTauShortL2"   , "Likelihood (#tau_{-})"   , 100, rFitErrTauShortL2.minimo, rFitErrTauShortL2.massimo); 
    TH1D hFitErrBL2 	    ( "hFitErrBL2"          , "Likelihood (B)"          , 100, rFitErrBL2.minimo, rFitErrBL2.massimo);
    TH1D hFitErrRL2 	    ( "hFitErrRL2"          , "Likelihood (R)"          , 100, rFitErrRL2.minimo, rFitErrRL2.massimo);

    // riempio gli istogrammi METODO 2
    for ( int i = 0; i < Nsim; i++ ) {        
        hFitTauL2.Fill(vFitTauL2.at(i));
        hFitTauShortL2.Fill(vFitTauShortL2.at(i));
        hFitBL2.Fill(vFitBL2.at(i));
        hFitRL2.Fill(vFitRL2.at(i));
        
	    hFitErrTauL2.Fill(vFitErrTauL2.at(i));
	    hFitErrTauShortL2.Fill(vFitErrTauShortL2.at(i));
        hFitErrBL2.Fill(vFitErrBL2.at(i));
        hFitErrRL2.Fill(vFitErrRL2.at(i));
    }

    std::cout << std::endl 
	      << "---------------------------------------------------------------------------------" 
	      << "\n./montecarlo.out " << B << " " << tau << " " << integrale << " " << taucorto 
		  << " " << R << " " << RebFactor <<std::endl 
	        << "\nEventi totali    " << integrale
	      << "\nEventi baseline  " << B*(End-Begin)
	      << "\nEventi exp lungo " << A*tau
		  << "\nEventi exp corto " << A*taucorto/R
	      << "\nValori in INPUT:"
	      << "\n	Baseline  " << B
	      << "\n	tau       " << tau
	      << "\n	integrale " << integrale
		  << "\n	taucorto  " << taucorto
		  << "\n	R         " << R
		  << "\n	RebFactor " << RebFactor
		  << "\nA+        " << A 
		  << "\nA-        " << Aminus << std::endl
          << "\nSimulazioni MC effettuate:                         "  << Nsim
	      << "\nIntervallo di generazione:                         [" << Begin     << ", " << End << "]" 
	      << "\nIntervallo di fit baseline METODO 2:               [" << StartBase << ", " << End << "]"
	      << "\nIntervallo di fit (exp+)+base METODO 2:            [" << beginFit  << ", " << End << "]"
	      << "\nIntervallo di fit (exp-)+(exp+)+base METODO 2:     [" << Begin     << ", " << End << "]" << std::endl;
    std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)";

    // fit gaussiani delle distribuzioni METODO 2
    std::cout<<"\n METODO 2" << std::endl;
    fitGaus(hFitTauL2);
    fitGaus(hFitTauShortL2);
    fitGaus(hFitBL2);
    fitGaus(hFitRL2);
    
    fitGaus(hFitErrTauL2);
    fitGaus(hFitErrTauShortL2);
    fitGaus(hFitErrBL2);
    fitGaus(hFitErrRL2);
    
    TCanvas can2( "METODO2", (canName + " METODO 2").c_str(), 1200 , 700 );
    can2.Divide(4,2);
    can2.cd(1);
        hFitTauL2.Draw();
        TLine lTau2(tau,0,tau,hFitTauL2.GetMaximum());
	    lTau2.SetLineColor(kRed);
        lTau2.Draw();
    can2.cd(2);
        hFitTauShortL2.Draw();
        TLine lTauShort2(taucorto,0,taucorto,hFitTauShortL2.GetMaximum());
	    lTauShort2.SetLineColor(kRed);
        lTauShort2.Draw();
    can2.cd(3);
        hFitBL2.Draw();
        TLine lB2(B,0,B,hFitBL2.GetMaximum());
	    lB2.SetLineColor(kRed);
        lB2.Draw();
    can2.cd(4);
        hFitRL2.Draw();
        TLine lR2(R,0,R,hFitRL2.GetMaximum());
	    lR2.SetLineColor(kRed);
        lR2.Draw();
    can2.cd(5);
        hFitErrTauL2.Draw();
    can2.cd(6);
        hFitErrTauShortL2.Draw();
    can2.cd(7);
        hFitErrBL2.Draw();
    can2.cd(8);
        hFitErrRL2.Draw();
        
    Root.Run();

    return 0;
}
//--------------- fine main ----------------------
 
void fitGaus(TH1D h){
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
   TFitResultPtr p = h.Fit("gaus","SQN","",min,max);
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
   double min = 0;
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

