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
#include "omp.h"

#include "TApplication.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"

// questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
#define	Begin     160	// inizio istogramma
#define StartBase 2138	// punto dell'istogramma in cui comincia la baseline
#define End       3904	// fine istogramma
#define Nsim      100  // numero simulazioni

TRandom3 r;
int counts;
 
struct dataBase {
    double media;
    double fitL;
    double fitC;

    double errMedia;
    double errFitL;
    double errFitC;
};

dataBase simulateBase( float B, int rebin );

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
                  << "    $ ./montecarlo [Baseline] [rebinFactor]" << std::endl << std::endl;
        return 0;
    }

    if ( argc < 3 ) {
        std::cout << "Pochi argomenti! Se non ti ricordi c'è l'opzione '--help'" << std::endl
                  << "Termino l'esecuzione..." << std::endl;
        return 0;
    }

    float B       = std::stof(args[1]);
    int RebFactor = std::stoi(args[2]);
    counts = B*(End-Begin);
    std::cout << "Eventi totali: " << counts << std::endl;
    
    TApplication Root("App",&argc,argv); 

    float sig = 2;
    TH1D hMedia( "hMedia", "Media"     , 100, B-sig, B+sig );
    TH1D hFitL ( "hFitL" , "Likelihood", 100, B-sig, B+sig );
    TH1D hFitC ( "hFitC" , "Chi2"      , 100, B-sig, B+sig );

    auto start = std::chrono::high_resolution_clock::now();
    dataBase data;

#pragma omp parallel for private(data,r,counts) shared(hMedia, hFitL, hFitC)
    for ( int i = 0; i < Nsim; i++ ) {        
        r.SetSeed(i);
        data = simulateBase( B, RebFactor );
#pragma omp critical
{
        hMedia.Fill(data.media);
        hFitL.Fill(data.fitL);
        hFitC.Fill(data.fitC);
}
    }
    auto time = std::chrono::high_resolution_clock::now() - start;
    std::cout << "CPU time: " << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << " msec." << std::endl;

    std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)"; 
    TCanvas can( "can", canName.c_str(), 1 );
    can.Divide(2,2);

    // TLine ............

    can.cd(1);
        hMedia.Draw();
        gPad->PaintLine(B,0,B,100);
    can.cd(2);
        hFitL.Draw();
    can.cd(3);
        hFitC.Draw();
    can.cd(4);
        //baseline.Draw();

    Root.Run();

    return 0;
}

dataBase simulateBase( float B, int rebin ) {
  
    // simulazione della baseline
    TH1F baseline("baseline","baseline",4096,0,4096);
    for ( int k = 0; k < counts; k++ )
    {
        baseline.Fill(r.Uniform(Begin,End));
	    //if (c>StartBase) k++;	
	    // 'counts' è l'integrale della baseline solo nella parte finale dello spettro
	    // il contatore k è incrementato solo se il numero generato è oltre StartBase
    }
    // rebin dell'istogramma
    baseline.Rebin(rebin); 

    // 1 - metodo della MEDIA
    //vediamo quali bin compongono la baseline nello spettro vero
    int N = 0;
    double Mean = 0;
    double ErrMean = 0;
    // calcolo della media
    for(int j=1; j<=(baseline.GetNbinsX()); j++)
    {
        double BinCenter = baseline.GetBinCenter(j);
	    if( BinCenter>StartBase && BinCenter<End) 
	    {
		    N++;
		    Mean += baseline.GetBinContent(j);
	    }
    }
    Mean = Mean/N;
    // calcolo l'errore della media
    for(int i=1; i<=(baseline.GetNbinsX()); i++)
    {
	    if( (baseline.GetBinCenter(i))>StartBase ) 
	    {
		    ErrMean += pow(baseline.GetBinContent(i)-Mean,2);
	    }
    }
    ErrMean = ErrMean/(N-1);
    ErrMean = sqrt(ErrMean/N); // la singola misura come errore ha lo scarto quadratico medio

    // 2 - metodo del FIT LIKELIHOOD
    TFitResultPtr BaseLikePtr = baseline.Fit("pol0","LSQ","",StartBase,End);
    //TF1 retta("retta","[0]+[1]*x",0,4096);
    //retta.FixParameter(1,0);
    //TFitResultPtr BaseLikePtr = baseline.Fit("retta","LSQ","",StartBase,End);
    // 3 - metodo del FIT CHI2
    TFitResultPtr BaseChi2Ptr = baseline.Fit("pol0","SQ","",StartBase,End);
    //TFitResultPtr BaseChi2Ptr = baseline.Fit("retta","SQ","",StartBase,End);
/*
    std::cout <<"\n=============== RISULTATI ================" //45
              <<"\nMedia baseline          = " << Mean                      << " +- " << ErrMean
              <<"\nFit Likelihood baseline = " << BaseLikePtr->Parameter(0) << " +- " << BaseLikePtr->ParError(0)
              <<"\nFit Chi2 baseline       = " << BaseChi2Ptr->Parameter(0) << " +- " << BaseChi2Ptr->ParError(0) << std::endl;
*/    
    dataBase data { 
                    Mean, 
                    BaseLikePtr->Parameter(0), 
                    BaseChi2Ptr->Parameter(0), 
                    ErrMean, 
                    BaseLikePtr->ParError(0), 
                    BaseChi2Ptr->ParError(0) 
                    };
    //TCanvas canvas("canvas","canvas",1);
    //canvas.cd();
    //baseline.Draw();
    //canvas.SaveAs("can.C");                
    return data;
}
