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
#include "math.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TFile.h"

#include "../../ProgressBar/progressbar.h"

// questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
#define	Begin     0	    // inizio istogramma
#define StartBase 2600	// punto dell'istogramma in cui comincia la baseline --> FISSATO dai test su fit pol0 della baseline (baselineStart.cc)
#define End       3904	// fine istogramma
#define Nsim      500   // numero simulazioni
#define beginFit  860	// inizio fit esponenziale dei muoni lunghi

struct distRange{
    double minimo;
    double massimo;
};

void fitGaus(TH1D& h);
distRange getRange(std::vector<double>& v);

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
                  << "    $ ./montecarlo [Baseline] [tau] [Integrale] [taucorto] [R] [rebinFactor] [OPZIONI]" << std::endl << std::endl
                  << "OPZIONI:" << std::endl
                  << "    --save-fit      : Salva i risultati grafici dei singoli Nsim fit in un file savedfits.root" << std::endl 
                  << "    --no-comp-check : Non effettua i controlli di compatibilità" << std::endl << std::endl;
        return 0;
    }

    if ( argc < 7 ) {
        std::cout << "Pochi argomenti! Se non ti ricordi c'è l'opzione '--help'" << std::endl
                  << "Termino l'esecuzione..." << std::endl;
        return 0;
    }
    
    bool saveFits = false;
    if ( argc == 8 && args[7] == "--save-fit" ) {
        saveFits = true;
        std::cout << "Salvo i grafici dei singoli fit in savedfits.root." << std::endl;
    }
    
    bool compCheck = true;
    if ( argc == 8 && args[7] == "--no-comp-check" ) {
        compCheck = false;
        std::cout << "Non faccio i controlli di compatibilità." << std::endl;
    }

    if ( argc == 9 && (( args[7] == "--save-fit" && args[8] == "--no-comp-check" ) || 
                       ( args[8] == "--save-fit" && args[7] == "--no-comp-check" )) ) {
        saveFits = true;
        compCheck = false;
        std::cout << "Salvo i grafici dei singoli fit in savedfits.root." << std::endl;
        std::cout << "Non faccio i controlli di compatibilità." << std::endl;
    }

    float B          = std::stof(args[1]); // valore tipico per 1 settimana: 1 (circa)
    double tau	     = std::stof(args[2]); // valore vero = 429 canali
	double integrale = std::stof(args[3]); // valore tipico per 1 settimana: 80000 (circa)
	double taucorto  = std::stof(args[4]); // valore vero = 172 canali (circa)
	double R         = std::stof(args[5]); // valore vero: 1.261
	int RebFactor    = std::stoi(args[6]);
	double A         = (integrale - B*(End - Begin)) / (tau + taucorto / R);
    double Aminus    = A/R;
	
	TApplication Root("App",&argc,argv);
	
    //gErrorIgnoreLevel = kError; // toglie i warning
    //gStyle->SetOptStat(0);

    // METODO 2
    std::vector<double> vFitTauL2;
    std::vector<double> vFitTauShortL2;
    std::vector<double> vFitBL2;
    std::vector<double> vFitRL2;
    
    std::vector<double> vFitErrTauL2;
    std::vector<double> vFitErrTauShortL2;
    std::vector<double> vFitErrBL2;
    std::vector<double> vFitErrRL2;

    // simulazione baseline+exp
    TH1D total("total","total",4096,0,4096);
		total.Rebin(RebFactor);
	TH1D total2("total2", "total2", 4096, 0, 4096);
		total2.Rebin(RebFactor);
        total2.SetStats(kFALSE);

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
    //TF1 prelFit ("prelFit","[0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x)",Begin,End); // funzione per fit preliminare       
    
    double val      = 0;
    double entina   = 0;
    double ratio    = 0;
    double errRatio = 0;

    TRandom3 r;
    int countTauPlusErr = 0;   // contatori per i fit con errore grande
    int countTauMinErr = 0;
    int countRErr = 0;
    
    int countTauPlusComp = 0;   // contatori per i fit incompatibili
    int countTauMinComp = 0;
    int countRComp = 0;
    
    int countTot = 0;
       
    // limiti di rigetto
    int compLimit = 3;
    double errLimit = 0.2;

    // ProgressBar
    ProgressBar bar(Nsim);
    bar.Init();
    
    TFile * file;
    TCanvas * c;
    if (saveFits) {
        file = new TFile("savedfits.root","RECREATE");
        c = new TCanvas("c2","c2",1);
        c->cd();
    }

////////////// CICLO MENSILE ////////////////
    for( int k = 0; k < Nsim; k++ ) {

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

        // fit preliminare per settare tau+ e tau- la 1a volta
        TFitResultPtr ptr = total2.Fit("expo","SNQ");
        
        double tauSet = -1/(ptr->Parameter(1));
        
        // METODO 2: pol0 per B, quindi B-setting e fit con fitFunc2, quindi fit complessivo
	    TFitResultPtr basePtr = total2.Fit("pol0","LSNQ","",StartBase,End);
        fitFunc2.SetParameter(2,basePtr->Parameter(0));
        fitFunc2.SetParameter(1,tauSet);
        //fitFunc2.SetParameter(1,tauPlusSet);

        // FIT (exp+)+Base direttamente su total2 con fitfunc2 (solo 1 exp) 
        total2.Fit("fitFunc2","LQN","",beginFit,End);
        fitFunc.SetParameter("tauShort",tauSet);
        //fitFunc.SetParameter("tauShort",tauMinSet);
        fitFunc.SetParameter("tauLong",fitFunc2.GetParameter("tau"));
        fitFunc.SetParameter("B",fitFunc2.GetParameter("B"));

        // FIT totale su total2 con fitfunc (ambedue le exp)
        total2.Fit("fitFunc","LQN","",Begin,End);
        
        if (saveFits) {
            total2.Draw();
            fitFunc.Draw("SAME");
            c->Write();
        }

        ratio    = fitFunc.GetParameter("Aplus")/fitFunc.GetParameter("Aminus");
        errRatio = sqrt( pow(fitFunc.GetParError(2),2)+pow(ratio*fitFunc.GetParError(0),2) )/fitFunc.GetParameter("Aminus");
    
        if (compCheck) {       
            // calcolo errori relativi e compatibilità       
            double compTauPlus = TMath::Abs(fitFunc.GetParameter("tauLong") - tau)/fitFunc.GetParError(3);
            double compTauMin  = TMath::Abs(fitFunc.GetParameter("tauShort") - taucorto)/fitFunc.GetParError(1);
            double compR       = TMath::Abs(ratio - R)/errRatio;

            double errRelTauPlus = fitFunc.GetParError(3)/fitFunc.GetParameter("tauLong");
            double errRelTauMin  = fitFunc.GetParError(1)/fitFunc.GetParameter("tauShort");
            double errRelR       = errRatio/ratio;

            // pushback se compatibili e con errori piccoli
            if ( compTauPlus   > compLimit ) countTauPlusComp++;
            if ( errRelTauPlus > errLimit  ) countTauPlusErr++; 
            if ( compTauMin    > compLimit ) countTauMinComp++;
            if ( errRelTauMin  > errLimit  ) countTauMinErr++; 
            if ( compR         > compLimit ) countRComp++;
            if ( errRelR       > errLimit  ) countRErr++;

            if ( compTauPlus > compLimit || errRelTauPlus > errLimit ||
                 compTauMin  > compLimit || errRelTauMin  > errLimit ||
                 compR       > compLimit || errRelR       > errLimit ) countTot++;
            
            else {
                vFitTauShortL2.push_back(fitFunc.GetParameter("tauShort"));
                vFitErrTauShortL2.push_back(fitFunc.GetParError(1)); 
                vFitTauL2.push_back(fitFunc.GetParameter("tauLong"));
                vFitErrTauL2.push_back(fitFunc.GetParError(3));       
                vFitBL2.push_back(fitFunc.GetParameter("B")/RebFactor);
                vFitErrBL2.push_back(fitFunc.GetParError(4)/RebFactor);
                vFitRL2.push_back(ratio);
                vFitErrRL2.push_back(errRatio);
            }
        }
        else {
            vFitTauShortL2.push_back(fitFunc.GetParameter("tauShort"));
            vFitErrTauShortL2.push_back(fitFunc.GetParError(1)); 
            vFitTauL2.push_back(fitFunc.GetParameter("tauLong"));
            vFitErrTauL2.push_back(fitFunc.GetParError(3));       
            vFitBL2.push_back(fitFunc.GetParameter("B")/RebFactor);
            vFitErrBL2.push_back(fitFunc.GetParError(4)/RebFactor);
            vFitRL2.push_back(ratio);
            vFitErrRL2.push_back(errRatio);
        }
    }
///////////// FINE CICLO MENSILE //////////////////
    if (saveFits) { 
        file->Close();
        delete file;
        delete c;
    }
	
    // range distribuzioni METODO 2
    distRange rFitTauL2         = getRange(vFitTauL2);
        std::cout << "\nhFitTauL2:         [" << rFitTauL2.minimo << ", " << rFitTauL2.massimo << "]" << std::endl;
    distRange rFitTauShortL2    = getRange(vFitTauShortL2);
        std::cout << "hFitTauShortL2:    [" << rFitTauShortL2.minimo << ", " << rFitTauShortL2.massimo << "]" << std::endl;
    distRange rFitBL2           = getRange(vFitBL2);
        std::cout << "hFitBL2:           [" << rFitBL2.minimo << ", " << rFitBL2.massimo << "]" << std::endl;
    distRange rFitRL2           = getRange(vFitRL2);
        std::cout << "hFitRL2:           [" << rFitRL2.minimo << ", " << rFitRL2.massimo << "]" << std::endl;
    
    distRange rFitErrTauL2      = getRange(vFitErrTauL2);
        std::cout << "hFitErrTauL2:      [" << rFitErrTauL2.minimo << ", " << rFitErrTauL2.massimo << "]" << std::endl;
    distRange rFitErrTauShortL2 = getRange(vFitErrTauShortL2);
        std::cout << "hFitErrTauShortL2: [" << rFitErrTauShortL2.minimo << ", " << rFitErrTauShortL2.massimo << "]" << std::endl;
    distRange rFitErrBL2        = getRange(vFitErrBL2);
        std::cout << "hFitErrBL2:        [" << rFitErrBL2.minimo << ", " << rFitErrBL2.massimo << "]" << std::endl;
    distRange rFitErrRL2        = getRange(vFitErrRL2); 
        std::cout << "hFitErrRL2:        [" << rFitErrRL2.minimo << ", " << rFitErrRL2.massimo << "]" << std::endl;

    // creo gli istogrammi delle distribuzioni METODO 2
    TH1D hFitTauL2 	       ( "hFitTauL2" 	     , "#tau_{+}"          , 100, rFitTauL2.minimo        , rFitTauL2.massimo*1.01 );
    TH1D hFitTauShortL2    ( "hFitTauShortL2"    , "#tau_{-}"          , 100, rFitTauShortL2.minimo   , rFitTauShortL2.massimo*1.01 );
    TH1D hFitBL2   	       ( "hFitBL2"   	     , "B"                 , 100, rFitBL2.minimo          , rFitBL2.massimo*1.01 );
    TH1D hFitRL2           ( "hFitRL2"           , "R"                 , 100, rFitRL2.minimo          , rFitRL2.massimo*1.01 );   
 
    TH1D hFitErrTauL2 	   ( "hFitErrTauL2"      , "#tau_{+} - errori" , 100, rFitErrTauL2.minimo     , rFitErrTauL2.massimo*1.01); 
    TH1D hFitErrTauShortL2 ( "hFitErrTauShortL2" , "#tau_{-} - errori" , 100, rFitErrTauShortL2.minimo, rFitErrTauShortL2.massimo*1.01); 
    TH1D hFitErrBL2 	   ( "hFitErrBL2"        , "B - errori"        , 100, rFitErrBL2.minimo       , rFitErrBL2.massimo*1.01);
    TH1D hFitErrRL2 	   ( "hFitErrRL2"        , "R - errori"        , 100, rFitErrRL2.minimo       , rFitErrRL2.massimo*1.01);

    // riempio gli istogrammi METODO 2
    for ( int i = 0; i < vFitTauL2.size(); i++ )       hFitTauL2.Fill(vFitTauL2[i]);
    for ( int i = 0; i < vFitTauShortL2.size(); i++ )  hFitTauShortL2.Fill(vFitTauShortL2[i]);
    for ( int i = 0; i < vFitBL2.size(); i++ )         hFitBL2.Fill(vFitBL2[i]);
    for ( int i = 0; i < vFitRL2.size(); i++ )         hFitRL2.Fill(vFitRL2[i]);
       
    for ( int i = 0; i < vFitErrTauL2.size(); i++ )      hFitErrTauL2.Fill(vFitErrTauL2[i]);
    for ( int i = 0; i < vFitErrTauShortL2.size(); i++ ) hFitErrTauShortL2.Fill(vFitErrTauShortL2[i]);
    for ( int i = 0; i < vFitErrBL2.size(); i++ )        hFitErrBL2.Fill(vFitErrBL2[i]);
    for ( int i = 0; i < vFitErrRL2.size(); i++ )        hFitErrRL2.Fill(vFitErrRL2[i]);
    
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
		  << "\n	A+        " << A 
		  << "\n	A-        " << Aminus << std::endl
          << "\nSimulazioni MC effettuate:                         "  << Nsim
	      << "\nIntervallo di generazione:                         [" << Begin     << ", " << End << "]" 
	      << "\nIntervallo di fit baseline METODO 2:               [" << StartBase << ", " << End << "]"
	      << "\nIntervallo di fit (exp+)+base METODO 2:            [" << beginFit  << ", " << End << "]"
	      << "\nIntervallo di fit (exp-)+(exp+)+base METODO 2:     [" << Begin     << ", " << End << "]" << std::endl << std::endl;

    // fit gaussiano delle distribuzioni
    fitGaus(hFitTauL2);
    fitGaus(hFitErrTauL2);
	
    fitGaus(hFitTauShortL2);
    fitGaus(hFitErrTauShortL2);
    
    fitGaus(hFitRL2);
    fitGaus(hFitErrRL2);
    
    fitGaus(hFitBL2); 
    fitGaus(hFitErrBL2);
	
    std::cout << "Valori con errore relativo < " << errLimit*100 << "% :" << std::endl
              << "tauLong:  " << Nsim-countTauPlusErr << "/" << Nsim << " [" << (Nsim-countTauPlusErr)*100./Nsim << "%]" << std::endl
              << "tauShort: " << Nsim-countTauMinErr  << "/" << Nsim << " [" << (Nsim-countTauMinErr)*100./Nsim  << "%]" << std::endl
              << "R:        " << Nsim-countRErr       << "/" << Nsim << " [" << (Nsim-countRErr)*100./Nsim       << "%]" << std::endl << std::endl
              << "Valori con compatibilità < " << compLimit << " :" << std::endl
              << "tauLong:  " << Nsim-countTauPlusComp << "/" << Nsim << " [" << (Nsim-countTauPlusComp)*100./Nsim << "%]" << std::endl
              << "tauShort: " << Nsim-countTauMinComp  << "/" << Nsim << " [" << (Nsim-countTauMinComp)*100./Nsim  << "%]" << std::endl
              << "R:        " << Nsim-countRComp       << "/" << Nsim << " [" << (Nsim-countRComp)*100./Nsim       << "%]" << std::endl << std::endl
              << "Efficienza complessiva:  " << Nsim-countTot << "/" << Nsim << " [" << (Nsim-countTot)*100./Nsim  << "%]" << std::endl;

    // DRAW SECTION   
   std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)";
   TCanvas can2( "METODO2", (canName + " - Likelihood").c_str(), 1200 , 700 );
    can2.Divide(4,2);
    can2.cd(1);
        hFitTauL2.Draw();
        TLine lTau2(tau,0,tau,hFitTauL2.GetMaximum());
	    lTau2.SetLineColor(kGreen+3);
	    lTau2.SetLineWidth(3);
        lTau2.Draw();
    can2.cd(2);
        hFitTauShortL2.Draw();
        TLine lTauShort2(taucorto,0,taucorto,hFitTauShortL2.GetMaximum());
	    lTauShort2.SetLineColor(kGreen+3);
	    lTauShort2.SetLineWidth(3);
        lTauShort2.Draw();
    can2.cd(3);
        hFitBL2.Draw();
        TLine lB2(B,0,B,hFitBL2.GetMaximum());
	    lB2.SetLineColor(kGreen+3);
	    lB2.SetLineWidth(3);
        lB2.Draw();
    can2.cd(4);
        hFitRL2.Draw();
        TLine lR2(R,0,R,hFitRL2.GetMaximum());
	    lR2.SetLineColor(kGreen+3);
	    lR2.SetLineWidth(3);
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
//----------- fine main ---------------

void fitGaus(TH1D& h){
   
    int dim = h.GetNbinsX();
    double min = h.GetBinCenter(1);
    double max = h.GetBinCenter(dim);
    //TF1 f( "gaus", "[0]*exp(-0.5*((TMath::Abs(x)-[1])/[2])**2)",0,1000);
    //    f.SetParameter(1,(max+min)/2);
    //    f.SetParameter(2,(max-min)/2);

    TFitResultPtr p = h.Fit("gaus","SQ","",min,max);

    std::cout << h.GetName() << std::endl;
    if ( p!= -1 ) {
        double mean     = p->Parameter(1);
        double err      = p->ParError(1);
        double sigma    = p->Parameter(2);
        double errsigma = p->ParError(2);

        std::cout << "mean = "       << mean  << " +- " << err << std::endl
	              << "sigma = "      << sigma << " +- "<< errsigma << std::endl
                  << "mean/sigma = " << "("   << sigma*100/mean << "%)" << std::endl << std::endl;
    }
    
    else std::cout << "Fit fallito." << std::endl << std::endl;

    return;
}

distRange getRange(std::vector<double>& v){
   double min = 10000;
   double max = 0;
   double entry;
   int N = v.size();
   for ( int k = 0; k < N; k++ ) {
        entry = v[k];
	    if( min > entry ) min = entry;
	    if( max < entry ) max = entry;
   }
   
   if ( (min == 10000 && max == 0) || min >= max ) {
       min = 0;
       max = 10;
   }
   
   distRange d {min,max};
   return d;
}
