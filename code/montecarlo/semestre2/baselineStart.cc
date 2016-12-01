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
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TLine.h"
#include "TText.h"
#include "TGraphErrors.h"

// ----> VALORI CON CUI GIOCARE <----
#define midValue  2538  // =5tau
#define range     1000*2// intervallo di valori attorno a midValue (canali)

#define	Begin     160	// inizio istogramma
#define End       3904	// fine istogramma
#define Nsim      100   // numero punti nel plot
#define beginFit  160	// inizio fit esponenziale

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
                  << "MonteCarlo per la misura della vita media dei muoni cosmici in alluminio. Analisi del punto startBase." << std::endl
                  << "Autori: Mattia Faggin, Davide Piras, Luigi Pertoldi" << std::endl << std::endl
                  << "Utilizzo:" << std::endl 
                  << "    $ ./montecarlo [Baseline] [tau] [A]" << std::endl << std::endl;
        return 0;
    }

    if ( argc < 3 ) {
        std::cout << "Pochi argomenti! Se non ti ricordi c'è l'opzione '--help'" << std::endl
                  << "Termino l'esecuzione..." << std::endl;
        return 0;
    }

    float B       = std::stof(args[1]);
    double tau	  = std::stoi(args[2]);
    double A	  = std::stoi(args[3]);
    //int RebFactor = std::stoi(args[4]);

    TApplication Root("App",&argc,argv); 

    std::vector<double> vFitBL2;
    std::vector<double> vFitErrBL2;
    std::vector<double> vStartBase;

    // simulazione baseline+exp
    TH1D total("total","total",4096,0,4096);

    int StartBase = midValue;
    StartBase    -= range/2;

    if ( midValue-range/2 < 0 || midValue+range/2 > 4096 ) { std::cout << "Hai esagerato. Termino..." << std::endl; return 0; }

    // ciclo delle simulazioni
    std::cout << "Run: ";
    for( int k = 0; k < Nsim; k++ ) {

        r.SetSeed(1);
    	
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
            if ( val > Begin && val < End ) exponential.Fill(val);
        }
    	// rebin eventuale
    	//exponential.Rebin(RebFactor);
        
        // sommo le due distribuzioni
    	total.Add( &baseline , &exponential );
        
        // prendo i risultati del fit, METODO 2: 'pol0' per fittare la baseline
	    TFitResultPtr basePtr = total.Fit("pol0","LSNQ","",StartBase,End);
        
        vFitBL2.push_back(basePtr->Parameter(0));
        vFitErrBL2.push_back(basePtr->ParError(0));
        vStartBase.push_back(StartBase);
        
        StartBase += range/Nsim;
    }
    
    TCanvas can("can", "StartBase Analisys",1200,900);
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
    lTick2.Draw();

    Root.Run();
    return 0;
}
