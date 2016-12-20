/* Programma per l'analisi dati finale (2° semestre)
 * Mattia Faggin, Davide Piras, Luigi Pertoldi
 * Created: 20 dic 2016
 *
 * Usage: ./analisiFinale [filelist] [rebin] 
 *
 * inputFile format: elenco dei file .txt con i dati, attenzione a non 
 *                   aggiungere spazi oppure newline characters
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <string>

#include "TH1D.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TApplication.h"

int main ( int argc , char** argv ) {

    // salvo gli argomenti prima di darli in pasto al TApplication perchè lui li modifica
    std::vector<std::string> args;
    args.reserve(argc);
    for ( int i = 0; i < argc; i++ ) args.push_back(argv[i]);

    if ( argc == 2 && args[1] == "--help" ) {
        std::cout << std::endl
                  << "Analisi Dati Finale per la misura della vita media dei muoni cosmici in alluminio." << std::endl
                  << "Autori: Mattia Faggin, Davide Piras, Luigi Pertoldi" << std::endl << std::endl
                  << "Utilizzo:" << std::endl 
                  << "    $ ./lifetimeAnalysis [filelist] [rebin]" << std::endl << std::endl;
        return 0;
    }

    if ( argc < 3 ) {
        std::cout << "Pochi argomenti! Se non ti ricordi c'è l'opzione '--help'" << std::endl
                  << "Termino l'esecuzione..." << std::endl;
        return 0;
    }

    std::ifstream file(args[1].c_str());
    std::vector<std::ifstream> files;
    std::string name;
    while ( file >> name ) files.push_back(std::ifstream(name.c_str()));

    int rebin = std::stoi(args[2]);

    // definiamo gli estremi dell'istogramma
    TH1D data( "data" , "data" , 4096 , 0 , 4096);
        //data.SetTitle("Spettro globale");
        data.SetYTitle("dN/dx");
        data.SetXTitle("Delay (canali)");
        data.SetLineColor(kBlack);
        data.SetStats(kFALSE);
        //data.GetYaxis()->SetTitleOffset(1.5);
        data.SetLineColor(kBlue);

    std::string line;
    int totevents = 0;
    int height = 0;

    // data retrieving
    for ( int j = 0 ; j < files.size() ; j++ ) {
        for ( int i = 0 ; i < 8 ; i++ ) std::getline(files[j],line);    // skip first lines
    }
    
    // initialize histogram (maybe not needed...)
    for ( int i = 0 ; i <= 4096 ; i++ ) data.SetBinContent( i, 0 );
    
    // fill histograms
    for ( int j = 0 ; j < files.size() ; j++ ) {
        for ( int i = 1 ; i <= 4096 ; i++ ) {
    
            files[j] >> height;
            data.SetBinContent( i , data.GetBinContent(i)+height );
            totevents += height; // tot events
        }
        files[j].close();   // close the file
    }

    // shift data to left
    int shift = 150;
    for ( int i = 1; i <= 4096; i++ ) {
        if ( i <= 4096-shift ) data.SetBinContent(i, data.GetBinContent(i+shift));
        if ( i >  4096-shift ) data.SetBinContent(i, 0);
    }

    // dump number of events in title
    std::ostringstream convert;
    convert << totevents;
    std::string eventsstring =  convert.str() + " eventi";
    data.SetTitle(eventsstring.c_str());

    // REBIN
    data.Rebin(rebin);

    //TApplication app("app", &argc, argv);
    TCanvas c( "c", "Analisi Dati", 1200 , 700);
    c.cd();
    c.SetGrid();

/* =========== FIT ============ */
    
    // define variables
    int begin = 168-shift;
    int end = 4096-shift;
    
    int midBase = 2600;
    int midExp = 860;
    
    TF1 fitFunc1("fitFunc1", "[0]*TMath::Exp(-x/[1])+[2]", midExp, end);
        fitFunc1.SetParName(0,"A");
        fitFunc1.SetParName(1,"tau");
        fitFunc1.SetParName(2,"B");

    TF1 fitFunc2("fitFunc2", "[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])+[4]", begin, end);
        fitFunc2.SetParName(0,"Aminus");
        fitFunc2.SetParName(1,"tauShort");
        fitFunc2.SetParName(2,"Aplus");
        fitFunc2.SetParName(3,"tauLong");
        fitFunc2.SetParName(4,"B");

    // fit preliminare per settare tau+ e tau
    TFitResultPtr ptrL = data.Fit("expo","SNQ", "", midExp, midBase);
    double tauSetL = -1/(ptrL->Parameter(1));
    
    TFitResultPtr ptrS = data.Fit("expo","SNQ", "", begin, midExp);
    double tauSetS = -1/(ptrS->Parameter(1));
    
    // pol0 per B
    TFitResultPtr basePtr = data.Fit("pol0", "LSNQ", "", midBase, end);
    // set par for next fit
    fitFunc1.SetParameter(2,basePtr->Parameter(0));
    fitFunc1.SetParameter(1,tauSetL);

    // exp + B
    data.Fit("fitFunc1", "LQN", "", midExp, end);
    // set par for next fit
    fitFunc2.SetParameter("tauShort", tauSetS);
    fitFunc2.SetParameter("tauLong",  fitFunc1.GetParameter("tau"));
    fitFunc2.SetParameter("B",        fitFunc1.GetParameter("B"));
 
    // exp + exp + B
    data.Fit("fitFunc2", "LQ", "", begin, end);

    // retrieve fit results
    double tauLong  = fitFunc2.GetParameter("tauLong");     // tau+
    double tauShort = fitFunc2.GetParameter("tauShort");    // tau-
    double Aminus   = fitFunc2.GetParameter("Aminus")/rebin;
    double Aplus    = fitFunc2.GetParameter("Aplus")/rebin;
    double B        = fitFunc2.GetParameter("B")/rebin;
    
    double tauLongErr  = fitFunc2.GetParError(3);
    double tauShortErr = fitFunc2.GetParError(1);
    double AminusErr   = fitFunc2.GetParError(0)/rebin;
    double AplusErr    = fitFunc2.GetParError(2)/rebin;
    double BErr        = fitFunc2.GetParError(4)/rebin;

    double R = Aplus/Aminus;
    double RErr = R * sqrt( AplusErr*AplusErr/(Aplus*Aplus) + AminusErr*AminusErr/(Aminus*Aminus));

    // calibration
    double m = 5.12162E-03;

    tauLong *= m;
    tauShort *= m;
 
    tauLongErr *= m;
    tauShortErr *= m;

    // dump results:
    std::cout << std::endl << "=========== Risultati Fit ===========" << std::endl << std::endl
              << "Tau+ = ( " << tauLong  << " +- " << tauLongErr  << " ) us" << std::endl
              << "Tau- = ( " << tauShort << " +- " << tauShortErr << " ) us" << std::endl
              << "A+   = "   << Aplus    << " +- " << AplusErr               << std::endl
              << "A-   = "   << Aminus   << " +- " << AminusErr              << std::endl
              << "R    = "   << R        << " +- " << RErr                   << std::endl
              << "B    = "   << B        << " +- " << BErr                   << std::endl
              << std::endl << "=====================================" << std::endl << std::endl;

    data.Draw();
    c.SaveAs("spectrumFit.pdf");
    //app.Run();

    return 0;
}
