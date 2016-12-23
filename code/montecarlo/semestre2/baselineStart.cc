/* 
 * Studio della bontà della ricostruzione della baseline in exp+cost al variare del punto
 * di partenza del fit. Modalità di fit: METODO 2: fitto prima solo la baseline e successivamente
 * la funzione intera settando il valore della baseline.
 *
 * Created: 29/11/2016
 * Authors: Mattia Faggin, Davide Piras, Luigi Pertoldi
 * 
 * To Do: 1. Il range è simmetrico, forse è più comodo che non lo sia
 *        2. Va sistemato per farlo funzionare con rebin > 1
 *        3. Adesso stiamo facendo simulazioni diverse per ogni midValue, ha senso?
 */

#include <iostream>
#include <vector>
#include <string>

#include "TApplication.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TLine.h"
#include "TStyle.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "../../ProgressBar/progressbar.h"

// ----> VALORI CON CUI GIOCARE <----
#define midValue  2538  // =5tau
#define range     1000*2// intervallo di valori attorno a midValue (canali)

#define	Begin     0	// inizio istogramma
#define End       3904	// fine istogramma
#define Nsim      30 // numero punti nel plot
#define beginFit  860	// inizio fit esponenziale

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
	int RebFactor    = std::stoi(args[6]);
	double A         = (integrale - B*(End - Begin)) / (tau + taucorto / R);

    TApplication Root("App",&argc,argv); 
    gStyle->SetOptStat(0);

    std::vector<double> vFitBL2;
    std::vector<double> vFitErrBL2;
    std::vector<double> vStartBase;

    // simulazione baseline+exp
    TH1D total("total","total",4096,0,4096);

    int StartBase = midValue;
    StartBase    -= range/2;
    
    TCanvas c("c", "Simulations",1200,900);
    TProfile profileB("profileB","Ricostruzione baseline b (TProfile)",range,midValue-range/2,midValue+range/2,"s");
        profileB.SetMinimum(B-0.2);
        //profileB.SetMaximum(B+0.8);
        profileB.GetXaxis()->SetTitle("canale");
        profileB.GetYaxis()->SetTitle("baseline");
    TProfile profileErrB("profileErrB","Incertezza #sigma_{b} (TProfile)",range,midValue-range/2,midValue+range/2,"s");
        profileErrB.SetMinimum(0.02);
        //profileErrB.SetMaximum(0.045);
        profileErrB.GetXaxis()->SetTitle("canale");
        profileErrB.GetYaxis()->SetTitle("errore");

    if ( midValue-range/2 < 0 || midValue+range/2 > 4096 ) { std::cout << "Hai esagerato. Termino..." << std::endl; return 0; }

    // ciclo delle simulazioni
    std::cout << "Run: ";
    ProgressBar bar(Nsim);
    bar.Init();
    for(int j=0; j<Nsim; j++)
    {
        bar.Update(j);
        for( int i = 0; i < Nsim; i++ ) {
            
            r.SetSeed(i+1);
    	
         // simulazione della baseline
    	    TH1D baseline("baseline","baseline",4096,0,4096);
    	    for ( int k = 0; k < B*(End-Begin); k++ ) baseline.Fill(r.Uniform(Begin,End));
            // rebin 
    	    //baseline.Rebin(RebFactor);
    	
         //simulazione esponenziale
   	        TH1D exponential("exponential","exponential",4096,0,4096);
    	    double val;
            for ( int k = 0; k < A*tau; k++ ) {
                val = r.Exp(tau);
                while (val<Begin || val>End)    val = r.Exp(tau);
                /*if ( val > Begin && val < End )*/ exponential.Fill(val);
            }
    	    // rebin eventuale
    	    //exponential.Rebin(RebFactor);
        
            // sommo le due distribuzioni
    	    total.Add( &baseline , &exponential );
        
            // prendo i risultati del fit, METODO 2: 'pol0' per fittare la baseline
	        TFitResultPtr basePtr = total.Fit("pol0","LSNQ","",StartBase,End);
	        profileB.Fill(StartBase,basePtr->Parameter(0));
	        profileErrB.Fill(StartBase,basePtr->ParError(0));
        }
            /*vFitBL2.push_back(basePtr->Parameter(0));
            vFitErrBL2.push_back(basePtr->ParError(0));
            vStartBase.push_back(StartBase);*/
        
            StartBase += range/Nsim;
    }
    
    c.Divide(1,2);
    c.cd(1);
    gPad->SetGrid();
    profileB.Draw("E1");
    TLine l10( profileB.GetXaxis()->GetXmin(), B, profileB.GetXaxis()->GetXmax(), B );
    l10.SetLineColor(kRed);
    l10.Draw();
    c.cd(2);
    gPad->SetGrid();
    profileErrB.Draw("E1");

    c.SaveAs("baselineStart.cc.pdf");
    
    /*TCanvas can("can", "StartBase Analisys",1200,900);
    can.Divide(1,2);

    // definisco il grafico
    TGraphErrors gr( Nsim, &vStartBase[0] , &vFitBL2[0] , 0 , &vFitErrBL2[0] );
    gr.SetTitle("Baseline ricostruita al variare di StartBase");
    gr.GetXaxis()->SetTitle("channel");
    gr.GetYaxis()->SetTitle("baseline ricostruita");
    gr.SetMarkerStyle(23);
    
    // linea sul valore vero della baseline B
    TLine l10( gr.GetXaxis()->GetXmin(), B, gr.GetXaxis()->GetXmax(), B );
    l10.SetLineColor(kRed);
    TLine lTick( 2538, gr.GetYaxis()->GetXmin(), 2538, gr.GetYaxis()->GetXmax());
    lTick.SetLineColor(kRed);
    
    can.cd(1);
    gr.Draw("AP");
    l10.Draw();
    lTick.Draw();

    // istogramma per evidenzare il range sulla distribuzione totale
    TH1F hSel("hSel", "Selezionati", 4096 , 0, 4096);
    hSel.SetFillStyle(3001);
    hSel.SetFillColor(kBlue);
    for ( int i = midValue-range/2+1; i < midValue+range/2+1; i++ ) hSel.SetBinContent( i, total.GetBinContent(i));
    
    // sistemo total per essere disegnato
    total.SetStats(kFALSE);
    total.SetTitle("Spettro simulato");
    total.SetXTitle("channel");
    total.SetYTitle("counts");
    
    TLine lTick2( 2538, 0, 2538, total.GetMaximum());
    lTick2.SetLineColor(kRed);

    can.cd(2);
    total.Draw();
    hSel.Draw("SAME");
    lTick2.Draw();*/

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
