#include <iostream>
#include <cmath>
#include "TApplication.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"

#define	Begin 160
#define StartBase 2138
#define End 3904

int main(int argc, char* argv[]) {

    TApplication Root("App",&argc,argv);   

    TCanvas can("can","Simulazione MC",800,500);

    TH1F baseline("baseline","baseline",4096,0,4096);
    int counts;

    // questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
    /*int Begin	  = 160;
    int StartBase = 2138;
    int End	  = 3904;*/

    std::cout<<"\nConteggi = ";
    std::cin>>counts;
    TRandom3 r;
    int k=0;
    int c=0;
    while(k<counts)
    {
        c=r.Uniform(Begin,End);
        baseline.Fill(c);
	if (c>StartBase) k++;	// 'counts' è l'integrale della baseline solo nella parte finale dello spettro
    }
    
    int RebFactor;
    std::cout<<"\nRebin = ",
    std::cin>>RebFactor;
    baseline.Rebin(RebFactor); 

    //vediamo quali bin compongono la baseline nello spettro vero
    int N = 0;
    float Mean = 0;
    float ErrMean = 0;
    for(int j=1; j<=(baseline.GetNbinsX()); j++)
    {
        float BinCenter = baseline.GetBinCenter(j);
	//std::cout<<"\n"<<j<<"   BinCenter = "<<BinCenter;
	if( BinCenter>StartBase ) 
	{
		N++;
		Mean += baseline.GetBinContent(j);
	}
        //std::cout<<"   N = "<<N<<"   Mean = "<<Mean;
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

    std::cout<<"\nMedia baseline = "<< Mean << " +- "<< ErrMean << std::endl;

    baseline.Draw();

    Root.Run();
    return 0;
}


