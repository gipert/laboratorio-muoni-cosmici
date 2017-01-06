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

#include "../../ProgressBar/progressbar.h"

// questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
#define	Begin     160	// inizio istogramma
#define StartBase 2138	// punto dell'istogramma in cui comincia la baseline
#define End       3904	// fine istogramma
#define Nsim      2000   // numero simulazioni

#define beginFit  1500	// inizio fit esponenziale


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
                  << "    $ ./montecarlo [Baseline] [tau] [A] [rebinFactor]" << std::endl << std::endl;
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
    int RebFactor = std::stoi(args[4]);

    TApplication Root("App",&argc,argv);
    auto oldErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);

    counts = B*(End - Begin);

// nuova versione
// I vari vector sono i numeri generati dalla funzione simulateExp. Successivamente, di questi numeri individuo il minimo e il massimo, in modo tale da poter poi definire correttamente il range degli istogrammi
    std::vector<double> vFitBMean;
    vFitBMean.reserve(Nsim);
	std::vector<double> vFitBL;
	vFitBL.reserve(Nsim);
	std::vector<double> vFitBC;
	vFitBC.reserve(Nsim);
	std::vector<double> vFitErrBMean;
	vFitErrBMean.reserve(Nsim);
	std::vector<double> vFitErrBL;
	vFitErrBL.reserve(Nsim);
	std::vector<double> vFitErrBC;
	vFitErrBC.reserve(Nsim);

    std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)";
    TCanvas can( "can", canName.c_str(), 1 );
    // simulo le esponenziali, inserendo le stime dei parametri nei vector
    dataBase data;

    ProgressBar bar(Nsim);
    bar.Init();

    for ( int i = 0; i < Nsim; i++ ) {
        bar.Update(i);
        r.SetSeed(i);
        data = simulateBase( B, RebFactor );
		vFitBMean.push_back(data.media);
		vFitBL.push_back(data.fitL);
		vFitBC.push_back(data.fitC);
		vFitErrBMean.push_back(data.errMedia);
		vFitErrBL.push_back(data.errFitL);
		vFitErrBC.push_back(data.errFitC);
    }

    // recupero minimo e massimo dei range 
	distRange rFitBMean = getRange(vFitBMean);
	distRange rFitBL = getRange(vFitBL);
	distRange rFitBC = getRange(vFitBC);
	distRange rFitErrBMean = getRange(vFitErrBMean);
	distRange rFitErrBL = getRange(vFitErrBL);
	distRange rFitErrBC = getRange(vFitErrBC);

    // creo gli istogrammi
	TH1D hFitBMean("hFitBMean", "Media", 50, rFitBMean.minimo, rFitBMean.massimo);
	TH1D hFitBL("hFitBL", "Likelihood", 50, rFitBL.minimo, rFitBL.massimo);
	TH1D hFitBC("hFitBC", "#chi^{2}", 50, rFitBC.minimo, rFitBC.massimo);
	TH1D hFitErrBMean("hFitErrBMean", "Errori Media", 50, rFitErrBMean.minimo, rFitErrBMean.massimo);
	TH1D hFitErrBL("hFitErrBL", "Errori Likelihood", 50, rFitErrBL.minimo, rFitErrBL.massimo);
	TH1D hFitErrBC("hFitErrBC", "Errori #chi^{2}", 50, rFitErrBC.minimo, rFitErrBC.massimo);

    // riempio gli istogrammi
    for ( int i = 0; i < Nsim; i++ ) {        
		hFitBMean.Fill(vFitBMean.at(i));
		hFitBL.Fill(vFitBL.at(i));
		hFitBC.Fill(vFitBC.at(i));
		hFitErrBMean.Fill(vFitErrBMean.at(i));
		hFitErrBL.Fill(vFitErrBL.at(i));
		hFitErrBC.Fill(vFitErrBC.at(i));		
    }

    std::cout << std::endl 
	      << "---------------------------------------------------------------------------------" 
	      << "\n./montecarlo.out " << B << " " << tau << " " << A << " " << RebFactor <<std::endl 
	      << "\nEventi totali" << counts
	      << "\nValori in INPUT:"
	      << "\n	Baseline  " << B
	      << "\n	tau       " << tau
	      << "\n	A	  " << A
	      << "\n	RebFactor " << RebFactor
              << "\nSimulazioni MC effettuate " << Nsim << std::endl;
	      //<< "\nIntervallo di fit dell'esponenziale [" << beginFit << "," << StartBase << "]" << std::endl;

    // fit gaussiani delle distribuzioni

    TLine line(B,0,B,120);
    line.SetLineColor(kGreen+3);
    line.SetLineWidth(3);
    can.Divide(3,2);
    can.cd(1);
	fitGaus(hFitBMean);
        hFitBMean.GetXaxis()->SetTitle("B");
        hFitBMean.GetYaxis()->SetTitle("counts");
        hFitBMean.Draw();
        line.Draw("same");
    can.cd(2);
	fitGaus(hFitBL);
        hFitBL.GetXaxis()->SetTitle("B");
        hFitBL.GetYaxis()->SetTitle("counts");
        hFitBL.Draw();
        line.Draw("same");
    can.cd(3);
	fitGaus(hFitBC);
        hFitBC.GetXaxis()->SetTitle("B");
        hFitBC.GetYaxis()->SetTitle("counts");
        hFitBC.Draw();
        line.Draw("same");
   //   TCanvas canErr( "canErr", (canName + "ERRORI").c_str() , 1 );
    //canErr.Divide(2,2);
    can.cd(4);
	fitGaus(hFitErrBMean);
        hFitErrBMean.GetXaxis()->SetTitle("#sigma_{B}");
        hFitErrBMean.GetYaxis()->SetTitle("counts");
        hFitErrBMean.Draw();
    can.cd(5);
	fitGaus(hFitErrBL);
        hFitErrBL.GetXaxis()->SetTitle("#sigma_{B}");
        hFitErrBL.GetYaxis()->SetTitle("counts");
        hFitErrBL.Draw();
    can.cd(6);
	fitGaus(hFitErrBC);
        hFitErrBC.GetXaxis()->SetTitle("#sigma_{B}");
        hFitErrBC.GetYaxis()->SetTitle("counts");
        hFitErrBC.Draw();
  
    Root.Run();

    return 0;
}
//--------------- fine main ----------------------

dataBase simulateBase( float B, int rebin ) {
  
    // simulazione della baseline
    TH1D baseline("baseline","baseline",4096,0,4096);
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
    TFitResultPtr BaseLikePtr = baseline.Fit("pol0","LSQN","",StartBase,End);
    // 3 - metodo del FIT CHI2
    TFitResultPtr BaseChi2Ptr = baseline.Fit("pol0","SQN","",StartBase,End);

   /*std::cout <<"\n=============== RISULTATI ================" //45
              <<"\nMedia baseline          = " << Mean                      << " +- " << ErrMean
              <<"\nFit Likelihood baseline = " << BaseLikePtr->Parameter(0) << " +- " << BaseLikePtr->ParError(0)
              <<"\nFit Chi2 baseline       = " << BaseChi2Ptr->Parameter(0) << " +- " << BaseChi2Ptr->ParError(0) << std::endl;*/
    
    dataBase data { 
                    Mean, 
                    BaseLikePtr->Parameter(0), 
                    BaseChi2Ptr->Parameter(0), 
                    ErrMean, 
                    BaseLikePtr->ParError(0), 
                    BaseChi2Ptr->ParError(0),
		    baseline 
                    };
              
    return data;
}

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

   std::cout << "mean = " << mean << "+-" << err<<std::endl 
             << "sigma = "<< sigma << "(" << sigma*100/mean << "%)" << std::endl;
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
   distRange d {min,max};
   return d;
}
