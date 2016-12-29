#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TText.h"

#define	m	1//5.121E-03
#define q	0//-7.70E-01//-1.635	// dati del fit di calibrazione canale - tempo (in microsecondi!)

void openADC( std::string fileName , int rebin = 1 );

int main( int argc, char** argv) {

    // salvo gli argomenti prima di darli in pasto al TApplication
    // perchè lui li modifica
    std::vector<std::string> args;
    args.reserve(argc);
    for ( int i = 0; i < argc; i++ ) args.push_back(argv[i]);

    if ( argc == 2 && args[1] == "--help" ) {
        std::cout << std::endl
                  << "Visualizzatore di dati per l'ADC (misura della vita media dei muoni cosmici in alluminio)." << std::endl
                  << "Autori: Mattia Faggin, Davide Piras, Luigi Pertoldi" << std::endl << std::endl
                  << "Utilizzo:" << std::endl 
                  << "    $ ./lifetimeAnalysis [filelist]      [rebinFactor]" << std::endl
                  << "    $ ./lifetimeAnalysis [singletxtfile] [rebinFactor]" << std::endl << std::endl;
        return 0;
    }

    if ( argc < 3 ) {
        std::cout << "Pochi argomenti! Se non ti ricordi c'è l'opzione '--help'" << std::endl
                  << "Termino l'esecuzione..." << std::endl;
        return 0;
    }

    TApplication Root("app" , &argc , argv );
    openADC( args[1] , std::stoi(args[2]));
    Root.Run();

    return 0;
}

void openADC( std::string fileName , int rebin) {

	int rebinFactor = rebin;
    int totevents = 0;
    TH1D * data;

    if ( fileName == "filelist" ) {

        std::ifstream file(fileName.c_str());
        std::vector<std::ifstream> files;
        std::string name;
        while ( file >> name ) files.emplace_back(name.c_str());
	
	    Double_t min = 0*m + q;
	    Double_t max = 4096*m + q;
	    data = new TH1D( "data" , "data" , 4096 , min , max);
		    data->SetTitle("");
		    data->SetYTitle("Counts");
		    data->SetXTitle("Delay (#mus)");
		    data->SetLineColor(4);
		    data->SetStats(kFALSE);

	    int height;
	    std::string line;

    	// data retrieving
        for ( int j = 0 ; j < files.size() ; j++ ) {
    	    for ( int i = 0 ; i < 8 ; i++ ) { 
                std::getline(files[j],line); 
                } // skip first lines
        }
    
        for ( int j = 0 ; j < files.size() ; j++ ) {
    	    for ( int i = 0 ; i < 4096 ; i++ ) {
    
    		    files[j] >> height;
    		    data->SetBinContent( i , data->GetBinContent(i)+height );
    		    totevents += height; // tot events
    	    }
        files[j].close();
        }
    }


    else {
    
        std::ifstream file1(fileName.c_str());

	    Double_t min = 0*m + q;
	    Double_t max = 4096*m + q;
	    data = new TH1D( "data" , "data" , 4096 , min , max);
		    data->SetTitle("");
		    data->SetYTitle("Counts");
		    data->SetXTitle("Delay (#mus)");
		    data->SetLineColor(4);
		    data->SetStats(kFALSE);

	    int height;
	    std::string line;

	    // data retrieving
	    for ( int i = 0 ; i < 8 ; i++ ) std::getline(file1,line); // skip first lines
	    for ( int i = 0 ; i < 4096 ; i++ ) {
	
		    file1 >> height;
		    data->SetBinContent( i , height );
		    totevents += height; // tot events
	    }
    }
	
	data->Rebin(rebinFactor); // rebin histograms

	TCanvas * can = new TCanvas( "ADCdataCan" , "Adc data" , 1500 , 600);
	can->cd();
	//can->SetLogy();
	data->Draw();

	// dump number of events
	std::ostringstream convert;
	convert << totevents;
	std::string eventsstring = convert.str();

	TText * text = new TText( 0.5 , 0.94 , (fileName + " -- " + eventsstring + " events").c_str());
	text->SetTextSize(0.04);
	text->SetNDC();
	text->SetTextAlign(22);
	text->Draw();

	return;
}
