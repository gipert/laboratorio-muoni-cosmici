/* 
 * Studio della bontà della ricostruzione del fit exponenziale lungo in exp+exp+cost al variare del punto
 * di partenza del fit exp. Modalità di fit: METODO 2: fitto prima solo la baseline e successivamente
 * la funzione intera settando il valore della baseline, e successivamente la funzione intera settando i 
 * valori dei parametri precedenti.
 *
 * Created: 09/12/2016
 * Authors: Mattia Faggin, Davide Piras, Luigi Pertoldi
 * 
 */

#include <iostream>
#include <vector>
#include <string>

#include "TApplication.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TLine.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "../../ProgressBar/progressbar.h"

double m = 5.12162E-03;
double q = -7.6994E-01;

// ----> VALORI CON CUI GIOCARE <----
#define StartBase  2600// inizio del fit baseline, già scelto
#define	Begin     18	// inizio istogramma
#define End       3904	// fine istogramma
#define Nsim      30 // numero punti nel plot
#define beginFit  860	// inizio fit esponenziale -> qui si gioca
#define range     300*2// intervallo di valori attorno a beginFit (canali)


TRandom3 r;
int counts;

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
                  << "$ ./baselineStart [Baseline] [tau] [Integrale] [taucorto] [R] [rebinFactor] " << std::endl << std::endl;
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
	int RebFactor = std::stoi(args[6]);
	double Aminus = (integrale - B*(End - Begin)) / (R*tau*exp(-Begin / tau) + taucorto*exp(-Begin / taucorto));
	double A = Aminus*R;

    TApplication Root("App",&argc,argv); 

    // simulazione baseline+exp
    TH1D total("total","total",4096,0,4096);
	TH1D total2("total2", "total2", 4096, 0, 4096);

	TF1 fitFunc2("fitFunc2", "[0]*TMath::Exp(-x/[1])+[2]", beginFit, End);
	fitFunc2.SetParName(0, "A");
	fitFunc2.SetParName(1, "tau");
	fitFunc2.SetParName(2, "B");

    int Startfit = beginFit;
    Startfit    += range/2;
    
    TProfile profileB("profileB","TProfile baseline",range,beginFit-range/2,beginFit+range/2,"s");
    TProfile profileA("profileA","TProfile Aplus",range,beginFit-range/2,beginFit+range/2,"s");
	TProfile profiletau("profiletau", "Ricostruzione #tau_{+} con TProfile", range, (beginFit - range / 2)*m+q, (beginFit + range / 2)*m+q, "s");
    profileB.SetMinimum(B-1);
    profileB.SetMaximum(B+1);
	profileA.SetMinimum(A - 50);
	profileA.SetMaximum(A + 50);
	profiletau.SetMinimum((tau - 22)*m);
	profiletau.SetMaximum((tau + 5)*m);

    //if ( beginFit-range/2 < 0 || beginFit+range/2 > 4096 ) { std::cout << "Hai esagerato. Termino..." << std::endl; return 0; }

	// ciclo delle simulazioni
	double val = 0;
	double entina = 0;
	double ratio = 0;
	double errRatio = 0;

    ProgressBar bar(Nsim);
    bar.Init();
    for(int j=0; j<Nsim; j++)
    {
        bar.Update(j);
        for( int i = 0; i < 100; i++ ) {
            
            r.SetSeed(i+1);
			// simulazione della baseline
			TH1D baseline("baseline", "baseline", 4096, 0, 4096);
			for (int k = 0; k < B*(End - Begin); k++)
			{
				baseline.Fill(r.Uniform(Begin, End));
			}
			// rebin 
			baseline.Rebin(RebFactor);
			//simulazione primo esponenziale
			TH1D exponential("exponential", "exponential", 4096, 0, 4096);
			for (int k = 0; k < A*tau*exp(-Begin / tau); k++)
			{
				val = r.Exp(tau);
				while (val<Begin || val>End) val = r.Exp(tau);
				exponential.Fill(val);
			}
			// rebin 
			exponential.Rebin(RebFactor);

			// rebin 
			exponential.Rebin(RebFactor);
        
            // sommo le due distribuzioni
    	    total.Add( &baseline , &exponential );

			//simulazione secondo esponenziale
			TH1D exponential2("exponential2", "exponential2", 4096, 0, 4096);
			for (int k = 0; k < Aminus*taucorto*exp(-Begin / taucorto); k++)
			{
				entina = r.Exp(taucorto);
				while (entina<Begin || entina>End) entina = r.Exp(taucorto);
				exponential2.Fill(entina);
			}
			// rebin 
			exponential2.Rebin(RebFactor);

			// somma istogrammi exp+ e baseline (già in total) e exp-
			total2.Add(&total, &exponential2);
        
            // prendo i risultati del fit, METODO 2: 'pol0' per fittare la baseline
	        TFitResultPtr basePtr = total2.Fit("pol0","LSNQ","",StartBase,End);
			fitFunc2.SetParameter(2, basePtr->Parameter(0));
			fitFunc2.SetParameter(1, tau);

			// FIT (exp+)+Base direttamente su total2 con fitfunc2 (solo 1 exp) 
			total2.Fit("fitFunc2", "LQN", "", Startfit, End);
	        profileB.Fill(Startfit,fitFunc2.GetParameter(2));
			profileA.Fill(Startfit, fitFunc2.GetParameter(0));
			profiletau.Fill(Startfit*m+q, fitFunc2.GetParameter(1)*m);

        }
              Startfit -= range/Nsim;
    }
    
    //total2.Draw();
    
    TCanvas c("c", "Simulations",1200,900);
    /*c.Divide(1,3);
    c.cd(1);
    profileB.Draw();
    TLine lB( profileB.GetXaxis()->GetXmin(), B, profileB.GetXaxis()->GetXmax(), B );
    lB.SetLineColor(kRed);
    lB.Draw();
    c.cd(2);
	profileA.Draw();
	TLine lA(profileA.GetXaxis()->GetXmin(), A, profileA.GetXaxis()->GetXmax(), A);
	lA.SetLineColor(kRed);
	lA.Draw();
	c.cd(3);*/
	profiletau.GetXaxis()->SetTitle("Ritardo [#mu s]");
	profiletau.GetYaxis()->SetTitle("#tau_{+} [#mu s]");
	profiletau.SetStats(0);
	gPad->SetGrid();
	profiletau.Draw("E1");
	TLine ltau(profiletau.GetXaxis()->GetXmin(), tau*m, profiletau.GetXaxis()->GetXmax(), tau*m);
	ltau.SetLineColor(kRed);
	ltau.Draw();
    
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
		  << "\n	A         " << A << std::endl;

    Root.Run();
    return 0;
}
