#include <iostream>
#include <vector>
#include <string>

#include "../../ProgressBar/progressbar.h"

std::vector<bool> montecarlo( float B, double tau, double integrale, double taucorto, double R, int RebFactor );

int main( int argc, char* argv[] ) {
/*
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

    float B          = std::stof(args[1]); // valore tipico per 1 settimana: 1 (circa)
    double tau	     = std::stof(args[2]); // valore vero = 429 canali
	double integrale = std::stof(args[3]); // valore tipico per 1 settimana: 80000 (circa)
	double taucorto  = std::stof(args[4]); // valore vero = 172 canali (circa)
	double R         = std::stof(args[5]); // valore vero: 1.261
	int RebFactor    = std::stoi(args[6]);
	double A         = (integrale - B*(End - Begin)) / (tau + taucorto / R);
	double Aminus    = A/R;	
*/
    
   std::vector<bool> vOut;
   std::vector<int>  vCounter(4,0);
   //std::vector<std::vector<std::vector<std::vector<bool>>>> matrixBool;
   // matrixBool[i]         taulungo
   //             [j]       taucorto
   //               [l]     R
   //                 [k]   Integrale
    
    int progress = 1;

    for ( int tauL = 250; tauL < 650; tauL+=180 ) {
        for ( int tauS = 100; tauS < 250; tauS+=72 ) {
            for ( double R = 0.621; R < 2; R+=0.6 ) {
                for ( int I = 80000; I < 410000; I+=160000 ) {
                    int B = I/80000;
                    std::cout << "\nB=" << B << "   tauL=" << tauL << "   tauS=" << tauS << "   I=" << I << "   R=" << R << std::endl;
                    std::cout << "Simulazione " << progress << "/81 " << std::flush;
                    
                    
                    vOut = montecarlo( B, tauL, I, tauS, R, 1 );
                    for ( int i = 0; i < 4; i++ ) vCounter[i] += vOut[i]; 
                    std::cout << "\nvOut={" << vOut[0] << "," << vOut[1] << "," << vOut[2] << "," << vOut[3] << std::endl;
                    progress++;
                    std::cout << std::endl;
                }
            }
        }
    }
    
    std::cout << "Efficienza tau+: " << 100 - vCounter[0]*100./progress << "%" << std::endl
              << "Efficienza tau-: " << 100 - vCounter[1]*100./progress << "%" << std::endl
              << "Efficienza R   : " << 100 - vCounter[2]*100./progress << "%" << std::endl
              << "Efficienza Tot : " << vCounter[3]*100./progress << "%" << std::endl;
    
    return 0;
}
