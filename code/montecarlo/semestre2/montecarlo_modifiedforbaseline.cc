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

// questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
#define	Begin     160	// inizio istogramma
#define StartBase 2138	// punto dell'istogramma in cui comincia la baseline
#define End       3904	// fine istogramma
#define Nsim      500   // numero simulazioni

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

    counts = B*(End - Begin);
	//counts = A*tau;
    //std::cout << "Eventi totali: " << counts << std::endl;




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
	/*std::vector<double> vFitTauL;
	vFitTauL.reserve(Nsim);
    std::vector<double> vFitTauC;
    vFitTauC.reserve(Nsim);
    std::vector<double> vFitAL;
    vFitAL.reserve(Nsim);
    std::vector<double> vFitAC;
    vFitAC.reserve(Nsim);
    std::vector<double> vFitErrTauL;
    vFitErrTauL.reserve(Nsim);
    std::vector<double> vFitErrTauC;
    vFitErrTauC.reserve(Nsim);
    std::vector<double> vFitErrAL;
    vFitErrAL.reserve(Nsim);
    std::vector<double> vFitErrAC;
    vFitErrAC.reserve(Nsim);*/

    // simulo le esponenziali, inserendo le stime dei parametri nei vector
    dataBase data;
    for ( int i = 0; i < Nsim; i++ ) {        
        r.SetSeed(i);
        data = simulateBase( B, RebFactor );
		vFitBMean.push_back(data.media);
		vFitBL.push_back(data.fitL);
		vFitBC.push_back(data.fitC);
		vFitErrBMean.push_back(data.errMedia);
		vFitErrBL.push_back(data.errFitL);
		vFitErrBC.push_back(data.errFitC);
		/*vFitTauL.push_back(data.tauL);
        vFitTauC.push_back(data.tauC);
        vFitAL.push_back(data.AL);
 	vFitAC.push_back(data.AC);
	vFitErrTauL.push_back(data.errTauL);
	vFitErrTauC.push_back(data.errTauC);
	vFitErrAL.push_back(data.errAL);
	vFitErrAC.push_back(data.errAC);*/
    }

    // recupero minimo e massimo dei range 
	distRange rFitBMean = getRange(vFitBMean);
	distRange rFitBL = getRange(vFitBL);
	distRange rFitBC = getRange(vFitBC);
	distRange rFitErrBMean = getRange(vFitErrBMean);
	distRange rFitErrBL = getRange(vFitErrBL);
	distRange rFitErrBC = getRange(vFitErrBC);
	/*distRange rFitTauL = getRange(vFitTauL);
    distRange rFitTauC = getRange(vFitTauC);
    distRange rFitAL = getRange(vFitAL);
    distRange rFitAC = getRange(vFitAC);
    distRange rFitErrTauL = getRange(vFitErrTauL);
    distRange rFitErrTauC = getRange(vFitErrTauC);
    distRange rFitErrAL = getRange(vFitErrAL);
    distRange rFitErrAC = getRange(vFitErrAC);*/

    // creo gli istogrammi
    
	TH1D hFitBMean("hFitBMean", "Media (B)", 100, rFitBMean.minimo, rFitBMean.massimo);
	TH1D hFitBL("hFitBL", "Likelihood (B)", 100, rFitBL.minimo, rFitBL.massimo);
	TH1D hFitBC("hFitBC", "Chi (B)", 100, rFitBC.minimo, rFitBC.massimo);
	TH1D hFitErrBMean("hFitErrBMean", "ErrMedia (B)", 100, rFitErrBMean.minimo, rFitErrBMean.massimo);
	TH1D hFitErrBL("hFitErrBL", "ErrLikelihood (B)", 100, rFitErrBL.minimo, rFitErrBL.massimo);
	TH1D hFitErrBC("hFitErrBC", "ErrChi2 (B)", 100, rFitErrBC.minimo, rFitErrBC.massimo);
	/*TH1D hFitTauL 	( "hFitTauL" 	, "Likelihood (#tau)", 100, rFitTauL.minimo, rFitTauL.massimo );
    TH1D hFitTauC 	( "hFitTauC" 	, "Chi2 (#tau)"      , 100, rFitTauC.minimo, rFitTauC.massimo );
    TH1D hFitAL   	( "hFitAL"   	, "Likelihood (A)"   , 100, rFitAL.minimo, rFitAL.massimo );
    TH1D hFitAC   	( "hFitAC"   	, "Chi2 (A)"         , 100, rFitAC.minimo, rFitAC.massimo );
    TH1D hFitErrTauL 	( "hFitErrTauL" , "Likelihood (#tau)", 100, rFitErrTauL.minimo, rFitErrTauL.massimo); 
    TH1D hFitErrTauC 	( "hFitErrTauC" , "Chi2 (#tau)"      , 100, rFitErrTauC.minimo, rFitErrTauC.massimo); //2.27-sigerr, 2.27+sigerr );
    TH1D hFitErrAL   	( "hFitErrAL"   , "Likelihood (A)"   , 100, rFitErrAL.minimo, rFitErrAL.massimo); //0.707-sigerr,0.707+sigerr );
    TH1D hFitErrAC   	( "hFitErrAC"   , "Chi2 (A)"         , 100, rFitErrAC.minimo, rFitErrAC.massimo); //0.73-sigerr,0.73+sigerr );*/
    // riempio gli istogrammi
    for ( int i = 0; i < Nsim; i++ ) {        
		hFitBMean.Fill(vFitBMean.at(i));
		hFitBL.Fill(vFitBL.at(i));
		hFitBC.Fill(vFitBC.at(i));
		hFitErrBMean.Fill(vFitErrBMean.at(i));
		hFitErrBL.Fill(vFitErrBL.at(i));
		hFitErrBC.Fill(vFitErrBC.at(i));		
		/*hFitTauL.Fill(vFitTauL.at(i));
        hFitTauC.Fill(vFitTauC.at(i));
        hFitAL.Fill(vFitAL.at(i));
 	hFitAC.Fill(vFitAC.at(i));
	hFitErrTauL.Fill(vFitErrTauL.at(i));
	hFitErrTauC.Fill(vFitErrTauC.at(i));
	hFitErrAL.Fill(vFitErrAL.at(i));
	hFitErrAC.Fill(vFitErrAC.at(i));*/
    }


// vecchia versione
/*
    // da sistemare al variare dei parametri
    // se come estremi metto 0,0 l'istogramma è centrato in automatico
    float sig = 4;
    float sigerr = 0.05;
    TH1D hFitTauL 	( "hFitTauL" 	, "Likelihood (#tau)", 100, 426, 432 );
    TH1D hFitTauC 	( "hFitTauC" 	, "Chi2 (#tau)"      , 100, 422, 427 );
    TH1D hFitAL   	( "hFitAL"   	, "Likelihood (A)"   , 100, 995, 1004 );
    TH1D hFitAC   	( "hFitAC"   	, "Chi2 (A)"         , 100, 1000, 1010 );
    TH1D hFitErrTauL 	( "hFitErrTauL" , "Likelihood (#tau)", 100, 0.719,0.727); //2.29-sigerr, 2.29+sigerr );
    TH1D hFitErrTauC 	( "hFitErrTauC" , "Chi2 (#tau)"      , 100, 0.701,0.714); //2.27-sigerr, 2.27+sigerr );
    TH1D hFitErrAL   	( "hFitErrAL"   , "Likelihood (A)"   , 100, 2.224,2.242); //0.707-sigerr,0.707+sigerr );
    TH1D hFitErrAC   	( "hFitErrAC"   , "Chi2 (A)"         , 100, 2.226,2.244); //0.73-sigerr,0.73+sigerr );

    auto start = std::chrono::high_resolution_clock::now();
    dataExp data;

    for ( int i = 0; i < Nsim; i++ ) {        
        r.SetSeed(i);
        data = simulateExp( tau, A, RebFactor );
        hFitTauL.Fill(data.tauL);
        hFitTauC.Fill(data.tauC);
        hFitAL.Fill(data.AL);
 	hFitAC.Fill(data.AC);
	hFitErrTauL.Fill(data.errTauL);
	hFitErrTauC.Fill(data.errTauC);
	hFitErrAL.Fill(data.errAL);
	hFitErrAC.Fill(data.errAC);
    }

    auto time = std::chrono::high_resolution_clock::now() - start;
    std::cout << "CPU time: " << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << " msec." << std::endl;
*/




    //(data.hExp)->Draw();
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
    std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)";

    // fit gaussiani delle distribuzioni
	fitGaus(hFitBMean);
	fitGaus(hFitBL);
	fitGaus(hFitBC);
	fitGaus(hFitErrBMean);
	fitGaus(hFitErrBL);
	fitGaus(hFitErrBC);
	/*fitGaus(hFitTauL);
    fitGaus(hFitTauC);
    fitGaus(hFitAL);
    fitGaus(hFitAC);
    fitGaus(hFitErrTauL);
    fitGaus(hFitErrTauC);
    fitGaus(hFitErrAL);
    fitGaus(hFitErrAC); */

    TLine line(B,0,B,20);
    line.SetLineColor(kRed);
    TCanvas can( "can", canName.c_str(), 1 );
    can.Divide(3,2);
    can.cd(1);
        hFitBMean.GetXaxis()->SetTitle("B");
        hFitBMean.GetYaxis()->SetTitle("counts");
        hFitBMean.Draw();
        line.Draw("same");
    can.cd(2);
        hFitBL.GetXaxis()->SetTitle("B");
        hFitBL.GetYaxis()->SetTitle("counts");
        hFitBL.Draw();
        line.Draw("same");
    can.cd(3);
        hFitBC.GetXaxis()->SetTitle("B");
        hFitBC.GetYaxis()->SetTitle("counts");
        hFitBC.Draw();
        line.Draw("same");
   //   TCanvas canErr( "canErr", (canName + "ERRORI").c_str() , 1 );
    //canErr.Divide(2,2);
    can.cd(4);
        hFitErrBMean.GetXaxis()->SetTitle("#sigma_{B}");
        hFitErrBMean.GetYaxis()->SetTitle("counts");
        hFitErrBMean.Draw();
    can.cd(5);
        hFitErrBL.GetXaxis()->SetTitle("#sigma_{B}");
        hFitErrBL.GetYaxis()->SetTitle("counts");
        hFitErrBL.Draw();
    can.cd(6);
        hFitErrBC.GetXaxis()->SetTitle("#sigma_{B}");
        hFitErrBC.GetYaxis()->SetTitle("counts");
        hFitErrBC.Draw();
  
		/*can.Divide(2, 2);
		can.cd(1);
		hFitTauL.Draw();
		can.cd(2);
		hFitTauC.Draw();
		can.cd(3);
		hFitAL.Draw();
		can.cd(4);
		hFitAC.Draw();
		TCanvas canErr("canErr", (canName + "ERRORI").c_str(), 1);
		canErr.Divide(2, 2);
		canErr.cd(1);
		hFitErrTauL.Draw();
		canErr.cd(2);
		hFitErrTauC.Draw();
		canErr.cd(3);
		hFitErrAL.Draw();
		canErr.cd(4);
		hFitErrAC.Draw();*/
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


dataExp  simulateExp( double tau, double A, int rebin ){

    TH1D* exponential = new TH1D ("exponential","exponential",4096,0,4096);
    for ( int k = 0; k < counts; k++ )
    {
	exponential->Fill(r.Exp(tau));
    }
    exponential->Rebin(rebin);
    
    TF1 fExp("fExp","[0]*TMath::Exp(-x/[1])",0,4096);
    fExp.SetParameters(A,tau);

    // fit Likelihood
    double tauLike;
    double errTauLike;
    double ALike;
    double errALike;
    exponential->Fit("fExp","LQN","",beginFit,StartBase);
    tauLike     = fExp.GetParameter(1);
    errTauLike  = fExp.GetParError(1);
    ALike       = fExp.GetParameter(0);
    errALike	= fExp.GetParError(0);

    // fit Chi2
    double tauChi;
    double errTauChi;
    double AChi;
    double errAChi;
    exponential->Fit("fExp","QN","",beginFit,StartBase);
    tauChi      = fExp.GetParameter(1);
    errTauChi   = fExp.GetParError(1);
    AChi        = fExp.GetParameter(0);
    errAChi 	= fExp.GetParError(0);

    dataExp data {
		tauLike,
		tauChi,
		ALike,
		AChi,
		errTauLike,
		errTauChi,
		errALike,
    		errAChi,
		exponential
                 };
    return data; 
}

 
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

// main baseline
/*
    counts = B*(End-Begin);
    std::cout << "Eventi totali: " << counts << std::endl;
    float sig = 2;
    TH1D hMedia( "hMedia", "Media"     , 100, B-sig, B+sig );
    TH1D hFitL ( "hFitL" , "Likelihood", 100, B-sig, B+sig );
    TH1D hFitC ( "hFitC" , "Chi2"      , 100, B-sig, B+sig );

    auto start = std::chrono::high_resolution_clock::now();
    dataBase data;
    for ( int i = 0; i < Nsim; i++ ) {        
        r.SetSeed(i);
        data = simulateBase( B, RebFactor );
        hMedia.Fill(data.media);
        hFitL.Fill(data.fitL);
        hFitC.Fill(data.fitC);
    }
    auto time = std::chrono::high_resolution_clock::now() - start;
    std::cout << "CPU time: " << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << " msec." << std::endl;

    std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)"; 
    TCanvas can( "can", canName.c_str(), 1 );
    can.Divide(3);

    // TLine ............

    can.cd(1);
        hMedia.Draw();
        gPad->PaintLine(B,0,B,100);
    can.cd(2);
        hFitL.Draw();
    can.cd(3);
        hFitC.Draw();*/



