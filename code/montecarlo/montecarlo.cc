// nel codice utilizziamo 3 metodi per stimare la baseline: 1-media, 2-fit likelihood, 3-fit chi2

#include <iostream>
#include <cmath>
#include "TApplication.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

// questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
#define	Begin 160	// inizio istogramma
#define StartBase 2138	// punto dell'istogramma in cui comincia la baseline
#define End 3904	// fine istogramma

int main(int argc, char* argv[]) {

    TApplication Root("App",&argc,argv);   
    TCanvas can("can","Simulazione MC",800,500);

    // simulazione della baseline
    TH1F baseline("baseline","baseline",4096,0,4096);
    int counts;
    std::cout<<"\nConteggi = ";
    std::cin>>counts;
    int seed;
    std::cout<<"\nSeme = ";
    std::cin>>seed;		// inserendo da tastiera semi diversi non abbiamo sempre lo stesso istogramma
    TRandom3 r;
    r.SetSeed(seed);
    int k=0;
    int c=0;
    while(k<counts)
    {
        c=r.Uniform(Begin,End);
        baseline.Fill(c);
	if (c>StartBase) k++;	
	// 'counts' è l'integrale della baseline solo nella parte finale dello spettro
	// il contatore k è incrementato solo se il numero generato è oltre StartBase
    }
    // rebin dell'istogramma
    int RebFactor;
    std::cout<<"\nRebin = ",
    std::cin>>RebFactor;
    baseline.Rebin(RebFactor); 

    // 1 - metodo della MEDIA
    //vediamo quali bin compongono la baseline nello spettro vero
    int N = 0;
    float Mean = 0;
    float ErrMean = 0;
    for(int j=1; j<=(baseline.GetNbinsX()); j++)
    {
        float BinCenter = baseline.GetBinCenter(j);
	if( BinCenter>StartBase ) 
	{
		N++;
		Mean += baseline.GetBinContent(j);
	}
    }
    Mean = Mean/N;
    for(int i=1; i<=(baseline.GetNbinsX()); i++)
    {
	if( (baseline.GetBinCenter(i))>StartBase ) 
	{
		ErrMean += pow(baseline.GetBinContent(i)-Mean,2);
	}
    }
    ErrMean =ErrMean/(N-1);
    ErrMean =sqrt(ErrMean/N); // la singola misura come errore ha lo scarto quadratico medio

    // 2 - metodo del FIT LIKELIHOOD
    TFitResultPtr BaseLikePtr = baseline.Fit("pol0","LS","",StartBase,End);
    // 3 - metodo del FIT CHI2
    TFitResultPtr BaseChi2Ptr = baseline.Fit("pol0","S","",StartBase,End);

    std::cout<<"\n***RISULTATI***";
    std::cout<<"\nMedia baseline = "<< Mean << " +- "<< ErrMean;
    std::cout<<"\nFit Likelihood baseline = "<< BaseLikePtr->Parameter(0) << " +- "<< BaseLikePtr->ParError(0);
    std::cout<<"\nFit Chi2 baseline = "<< BaseChi2Ptr->Parameter(0) << " +- "<< BaseChi2Ptr->ParError(0)<<std::endl;

    baseline.Draw();
    Root.Run();
    return 0;
}


