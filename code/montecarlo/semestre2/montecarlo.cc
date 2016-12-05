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
#define	Begin     160	// inizio istogramma
#define StartBase 2700	// punto dell'istogramma in cui comincia la baseline
#define End       3904	// fine istogramma
#define Nsim      500  // numero simulazioni
#define beginFit  860	// inizio fit esponenziale


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
	

    TApplication Root("App",&argc,argv); 
    gErrorIgnoreLevel = kError; // toglie i warning

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
	TH1D total2("total2", "total2", 4096, 0, 4096);

    TF1 fitFunc("fitFunc","[0]*TMath::Exp(-x/[1])+[2]",beginFit,End);
    fitFunc.SetParName(0,"A");
    fitFunc.SetParName(1,"tau");
    fitFunc.SetParName(2,"B");
    TF1 fitFunc2("fitFunc2","[0]*TMath::Exp(-x/[1])+[2]",beginFit,End);
    fitFunc2.SetParName(0,"A");
    fitFunc2.SetParName(1,"tau");
    fitFunc2.SetParName(2,"B");

    // ProgressBar
    ProgressBar bar(Nsim);
    bar.Init();
    // ciclo delle simulazioni
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
		exponential.Fill(r.Exp(tau));
    	}
    	// rebin 
    	exponential.Rebin(RebFactor);
 
		total.Add(&baseline, &exponential);

		//simulazione secondo esponenziale
		TH1D exponential2("exponential2", "exponential2", 4096, 0, 4096);
		for (int k = 0; k < A*taucorto/R; k++)
		{
			exponential2.Fill(r.Exp(taucorto));
		}
		// rebin 
		exponential2.Rebin(RebFactor);
		
		total2.Add(&total, &exponential2);
		
        
	
	// METODO 1: funzione completa
        // parameter setting
	//fitFunc.SetParameter(0,A);
	fitFunc.SetParLimits(1, tau-20, tau+20);
	//fitFunc.SetParameter(2,B);
	////////////////// FIT UNICO //////////////////// 
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
	total2.Draw();
	
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
    TH1D hFitTauL 	    ( "hFitTauL" 	, "Likelihood (#tau)", 100, rFitTauL.minimo, rFitTauL.massimo );
    TH1D hFitAL   	    ( "hFitAL"   	, "Likelihood (A)"   , 100, rFitAL.minimo, rFitAL.massimo );
    TH1D hFitBL   	    ( "hFitBL"   	, "Likelihood (B)"   , 100, rFitBL.minimo, rFitBL.massimo );
    TH1D hFitErrBL 	    ( "hFitErrBL"   , "Likelihood (B)"   , 100, rFitErrBL.minimo, rFitErrBL.massimo); 
    TH1D hFitErrTauL 	( "hFitErrTauL" , "Likelihood (#tau)", 100, rFitErrTauL.minimo, rFitErrTauL.massimo); 
    TH1D hFitErrAL   	( "hFitErrAL"   , "Likelihood (A)"   , 100, rFitErrAL.minimo, rFitErrAL.massimo); 

    // creo gli istogrammi delle distribuzioni METODO 2
    TH1D hFitTauL2 	    ( "hFitTauL2" 	 , "Likelihood (#tau)", 100, rFitTauL2.minimo, rFitTauL2.massimo );
    TH1D hFitAL2   	    ( "hFitAL2"   	 , "Likelihood (A)"   , 100, rFitAL2.minimo, rFitAL2.massimo );
    TH1D hFitBL2   	    ( "hFitBL2"   	 , "Likelihood (B)"   , 100, rFitBL2.minimo, rFitBL2.massimo );
    TH1D hFitErrBL2 	( "hFitErrBL2"   , "Likelihood (B)"   , 100, rFitErrBL2.minimo, rFitErrBL2.massimo); 
    TH1D hFitErrTauL2 	( "hFitErrTauL2" , "Likelihood (#tau)", 100, rFitErrTauL2.minimo, rFitErrTauL2.massimo); 
    TH1D hFitErrAL2   	( "hFitErrAL2"   , "Likelihood (A)"   , 100, rFitErrAL2.minimo, rFitErrAL2.massimo); 

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
		  << "\n	A         " << A << std::endl
          << "\nSimulazioni MC effettuate:                         "  << Nsim
	      << "\nIntervallo di generazione:                         [" << Begin     << ", " << End << "]"
          << "\nIntervallo di fit METODO 1:                        [" << beginFit  << ", " << End << "]" 
	      << "\nIntervallo di fit baseline METODO 2:               [" << StartBase << ", " << End << "]"
	      << "\nIntervallo di fit complessivo (B fixed) METODO 2:  [" << beginFit  << ", " << End << "]" << std::endl;
    std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)";

    // fit gaussiani delle distribuzioni METODO 1
    std::cout<<"\n METODO 1" << std::endl;
    fitGaus(hFitTauL);
    fitGaus(hFitAL);
    fitGaus(hFitBL);
    fitGaus(hFitErrTauL);
    fitGaus(hFitErrAL);
    fitGaus(hFitErrBL);
    TCanvas can( "METODO1", (canName + " METODO 1").c_str(), 1200 , 700 );
    can.Divide(3,2);
    can.cd(1);
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
	lTau2.SetLineColor(kRed);
        lTau2.Draw();
    can2.cd(2);
        hFitAL2.Draw();
        TLine lA2(A,0,A,hFitAL2.GetMaximum());
	lA2.SetLineColor(kRed);
        lA2.Draw();
    can2.cd(3);
        hFitBL2.Draw();
        TLine lB2(B,0,B,hFitBL2.GetMaximum());
	lB2.SetLineColor(kRed);
        lB2.Draw();
    can2.cd(4);
        hFitErrTauL2.Draw();
    can2.cd(5);
        hFitErrAL2.Draw();
    can2.cd(6);
        hFitErrBL2.Draw();


/*    fitGaus(hFitTauL);
    fitGaus(hFitTauC);
    fitGaus(hFitAL);
    fitGaus(hFitAC);
    fitGaus(hFitErrTauL);
    fitGaus(hFitErrTauC);
    fitGaus(hFitErrAL);
    fitGaus(hFitErrAC); 

    TCanvas can( "can", canName.c_str(), 1 );
    can.Divide(2,2);
    can.cd(1);
        hFitTauL.Draw();
    can.cd(2);
        hFitTauC.Draw();
    can.cd(3);
        hFitAL.Draw();
    can.cd(4);
        hFitAC.Draw();
    TCanvas canErr( "canErr", (canName + "ERRORI").c_str() , 1 );
    canErr.Divide(2,2);
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

/*
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
*/
 
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



