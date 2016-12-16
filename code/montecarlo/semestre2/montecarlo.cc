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

#include <sstream>

#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TStyle.h"

#include "../../ProgressBar/progressbar.h"

// questi valori andranno poi aggiustati (inizio istogramma, inizio baseline, fine istogramma)
#define	Begin     0	    // inizio istogramma
#define StartBase 2600	// punto dell'istogramma in cui comincia la baseline --> FISSATO dai test su fit pol0 della baseline (baselineStart.cc)
#define End       3904	// fine istogramma
#define Nsim      500   // numero simulazioni
#define beginFit  860	// inizio fit esponenziale dei muoni lunghi

struct distRange{
    double minimo;
    double massimo;
};

TFitResultPtr fitGaus(TH1D h);
distRange getRange(std::vector<double> v);

//--------------- main ---------------------
std::vector<bool> montecarlo( float B, double tau, double integrale, double taucorto, double R, int RebFactor  ) {

    double A = (integrale - B*(End - Begin)) / (tau + taucorto / R);

    gErrorIgnoreLevel = kError; // toglie i warning
    //gStyle->SetOptStat(0);

    // METODO 2
    std::vector<double> vFitTauL2;
    std::vector<double> vFitTauShortL2;
    std::vector<double> vFitBL2;
    std::vector<double> vFitRL2;
    
    std::vector<double> vFitErrTauL2;
    std::vector<double> vFitErrTauShortL2;
    std::vector<double> vFitErrBL2;
    std::vector<double> vFitErrRL2;

    // simulazione baseline+exp
    TH1D total("total","total",4096,0,4096);
	TH1D total2("total2", "total2", 4096, 0, 4096);

    TF1 fitFunc("fitFunc","[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])+[4]",Begin,End);
        fitFunc.SetParName(0,"Aminus");
        fitFunc.SetParName(1,"tauShort");
        fitFunc.SetParName(2,"Aplus");
        fitFunc.SetParName(3,"tauLong");
        fitFunc.SetParName(4,"B");
    TF1 fitFunc2("fitFunc2","[0]*TMath::Exp(-x/[1])+[2]",beginFit,End);
        fitFunc2.SetParName(0,"A");
        fitFunc2.SetParName(1,"tau");
        fitFunc2.SetParName(2,"B");
    
    double val      = 0;
    double entina   = 0;
    double ratio    = 0;
    double errRatio = 0;

    TRandom3 r;
    int countTauPlus = 0;   // contatori per i successi dell'intero codice
    int countTauMin = 0;
    int countR = 0;
    int countTot = 0;

    // ProgressBar
    ProgressBar bar(Nsim);
    bar.Init();

////////////// CICLO MENSILE ////////////////
    for( int k = 0; k < Nsim; k++ ) {

        bar.Update(k);
	    r.SetSeed(0);

    	// simulazione della baseline
    	TH1D baseline("baseline","baseline",4096,0,4096);
    	for ( int k = 0; k < B*(End-Begin); k++ )
    	{
        	baseline.Fill(r.Uniform(Begin,End));
    	}
    	// rebin 
    	//baseline.Rebin(RebFactor);
    	//simulazione primo esponenziale
   	    TH1D exponential("exponential","exponential",4096,0,4096);
    	for ( int k = 0; k < A*tau; k++ )
    	{
    	    val = r.Exp(tau);
    	    while (val<Begin || val>End) val = r.Exp(tau);
		    exponential.Fill(val);
    	}
    	// rebin 
    	//exponential.Rebin(RebFactor);
 
        // somma istogrammi exp+ e baseline
		total.Add(&baseline, &exponential);

		//simulazione secondo esponenziale
		TH1D exponential2("exponential2", "exponential2", 4096, 0, 4096);
		for (int k = 0; k < A*taucorto/R; k++)
		{
		    entina = r.Exp(taucorto);
		    while (entina<Begin || entina>End) entina = r.Exp(taucorto);
			exponential2.Fill(entina);
		}
		// rebin 
		exponential2.Rebin(RebFactor);
		
		// somma istogrammi exp+ e baseline (già in total) e exp-
		total2.Add(&total, &exponential2);

        // fit preliminare per settare tau+ e tau- la 1a volta
        //TFitResultPtr ptr = total2.Fit("expo","SNQ");
        //double tauSet = -1/(ptr->Parameter(1));

        // METODO 2: pol0 per B, quindi B-setting e fit con fitFunc2, quindi fit complessivo
	    TFitResultPtr basePtr = total2.Fit("pol0","LSNQ","",StartBase,End);
        fitFunc2.SetParameter(2,basePtr->Parameter(0));
        //fitFunc2.SetParameter(1,tauSet);
        fitFunc2.SetParameter(1,tau);

        // FIT (exp+)+Base direttamente su total2 con fitfunc2 (solo 1 exp) 
        total2.Fit("fitFunc2","LQN","",beginFit,End);
        //fitFunc.SetParameter("tauShort",tauSet);
        fitFunc.SetParameter("tauShort",taucorto);
        fitFunc.SetParameter("tauLong",fitFunc2.GetParameter("tau"));
        fitFunc.SetParameter("B",fitFunc2.GetParameter("B"));

        // FIT totale su total2 con fitfunc (ambedue le exp)
        total2.Fit("fitFunc","LQN","",Begin,End);
        
        ratio    = fitFunc.GetParameter("Aplus")/fitFunc.GetParameter("Aminus");
        errRatio = sqrt( pow(fitFunc.GetParError(2),2)+pow(ratio*fitFunc.GetParError(0),2) )/fitFunc.GetParameter("Aminus");
       
        // calcolo errori relativi e compatibilità       
        double compTauPlus = TMath::Abs(fitFunc.GetParameter("tauLong") - tau)/fitFunc.GetParError(3);
        double compTauMin  = TMath::Abs(fitFunc.GetParameter("tauShort") - taucorto)/fitFunc.GetParError(1);
        double compR       = TMath::Abs(ratio - R)/errRatio;

        double errRelTauPlus = fitFunc.GetParError(3)/fitFunc.GetParameter("tauLong");
        double errRelTauMin  = fitFunc.GetParError(1)/fitFunc.GetParameter("tauShort");
        double errRelR       = errRatio/ratio;
        
        // limiti di rigetto
        int compLimit = 3;
        double errLimit = 0.2;

        // pushback se compatibili e con errori piccoli
        if ( compTauPlus > compLimit || errRelTauPlus > errLimit ) countTauPlus++; 
        if ( compTauMin  > compLimit || errRelTauMin  > errLimit ) countTauMin++; 
        if ( compR       > compLimit || errRelR       > errLimit ) countR++;

        if ( compTauPlus > compLimit || errRelTauPlus > errLimit ||
             compTauMin  > compLimit || errRelTauMin  > errLimit ||
             compR       > compLimit || errRelR       > errLimit ) countTot++;
        else {
            vFitTauShortL2.push_back(fitFunc.GetParameter("tauShort"));
            vFitErrTauShortL2.push_back(fitFunc.GetParError(1)); 
            vFitTauL2.push_back(fitFunc.GetParameter("tauLong"));
            vFitErrTauL2.push_back(fitFunc.GetParError(3));       
            vFitBL2.push_back(fitFunc.GetParameter("B"));
            vFitErrBL2.push_back(fitFunc.GetParError(4));
            vFitRL2.push_back(ratio);
            vFitErrRL2.push_back(errRatio);
        }
    }
    
    // range distribuzioni METODO 2
    distRange rFitTauL2         = getRange(vFitTauL2);
    distRange rFitTauShortL2    = getRange(vFitTauShortL2);
    distRange rFitBL2           = getRange(vFitBL2);
    distRange rFitRL2           = getRange(vFitRL2);
    
    distRange rFitErrTauL2      = getRange(vFitErrTauL2);
    distRange rFitErrTauShortL2 = getRange(vFitErrTauShortL2);
    distRange rFitErrBL2        = getRange(vFitErrBL2);
    distRange rFitErrRL2        = getRange(vFitErrRL2); 

    // creo gli istogrammi delle distribuzioni METODO 2
    TH1D hFitTauL2 	       ( "hFitTauL2" 	     , "Likelihood (#tau_{+})" , 100, rFitTauL2.minimo     , rFitTauL2.massimo );
    TH1D hFitTauShortL2    ( "hFitTauShortL2"    , "Likelihood (#tau_{-})" , 100, rFitTauShortL2.minimo, rFitTauShortL2.massimo );
    TH1D hFitBL2   	       ( "hFitBL2"   	     , "Likelihood (B)"        , 100, rFitBL2.minimo       , rFitBL2.massimo );
    TH1D hFitRL2           ( "hFitRL2"           , "Likelihood (R)"        , 100, rFitRL2.minimo       , rFitRL2.massimo );   
 
    TH1D hFitErrTauL2 	   ( "hFitErrTauL2"      , "Likelihood (#tau_{+})" , 100, rFitErrTauL2.minimo     , rFitErrTauL2.massimo); 
    TH1D hFitErrTauShortL2 ( "hFitErrTauShortL2" , "Likelihood (#tau_{-})" , 100, rFitErrTauShortL2.minimo, rFitErrTauShortL2.massimo); 
    TH1D hFitErrBL2 	   ( "hFitErrBL2"        , "Likelihood (B)"        , 100, rFitErrBL2.minimo       , rFitErrBL2.massimo);
    TH1D hFitErrRL2 	   ( "hFitErrRL2"        , "Likelihood (R)"        , 100, rFitErrRL2.minimo       , rFitErrRL2.massimo);

    // riempio gli istogrammi METODO 2
    for ( int i = 0; i < vFitTauL2.size(); i++ ) {        
        hFitTauL2.Fill(vFitTauL2.at(i));
        hFitTauShortL2.Fill(vFitTauShortL2.at(i));
        hFitBL2.Fill(vFitBL2.at(i));
        hFitRL2.Fill(vFitRL2.at(i));
        
	    hFitErrTauL2.Fill(vFitErrTauL2.at(i));
	    hFitErrTauShortL2.Fill(vFitErrTauShortL2.at(i));
        hFitErrBL2.Fill(vFitErrBL2.at(i));
        hFitErrRL2.Fill(vFitErrRL2.at(i));
    }
/*
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
		  << "\n	A+        " << A 
		  << "\n	A-        " << Aminus << std::endl
          << "\nSimulazioni MC effettuate:                         "  << Nsim
	      << "\nIntervallo di generazione:                         [" << Begin     << ", " << End << "]" 
	      << "\nIntervallo di fit baseline METODO 2:               [" << StartBase << ", " << End << "]"
	      << "\nIntervallo di fit (exp+)+base METODO 2:            [" << beginFit  << ", " << End << "]"
	      << "\nIntervallo di fit (exp-)+(exp+)+base METODO 2:     [" << Begin     << ", " << End << "]"
          << "\nEfficienza tau: "         << (Nsim-countTauPlus)*100./Nsim << "%"
          << "\nEfficienza taucorto: "    << (Nsim-countTauMin)*100./Nsim  << "%"
          << "\nEfficienza R: "           << (Nsim-countR)*100./Nsim       << "%"
          << "\nEfficienza complessiva: " << (Nsim-countTot)*100./Nsim     << "%" << std::endl;

*/
    // fit gaussiani delle distribuzioni METODO 2
    double errLimit = 0.2;
    int compLimit = 3;
    int factorLimit = 2;

    // recupero risultati fit e controllo errori e campatibilità
    std::vector<bool> vBool(4, true);
    TFitResultPtr p, pErr;

    p    = fitGaus(hFitTauL2);
    pErr = fitGaus(hFitErrTauL2);
    
    std::cout << "\np " << p << "   pErr " << pErr << std::endl;
    
    if ( p == -1 || pErr == -1 ) // il valore -1 emerge quando il fit non converge 
    {
        std::cout << "Tau+ fit failed. Bool set to 'false'." << std::endl;
        vBool[0] = false;
    }
/*
    if ( !(p == -1 || pErr == -1) &&
          (pErr->Parameter(1)/p->Parameter(1) > errLimit || 
           TMath::Abs(p->Parameter(1)-tau)/pErr->Parameter(1) > compLimit ||
           p->Parameter(2) > factorLimit*pErr->Parameter(1)                 ) ) vBool[0] = false;*/
    if ( !(p == -1 || pErr == -1) && pErr->Parameter(1)/p->Parameter(1) > errLimit)
    {
        std::cout << "Errore relativo troppo grande (" << pErr->Parameter(1)/p->Parameter(1) << ">" << errLimit << ")" << std::endl;
        vBool[0] = false;
    } 
    if ( !(p == -1 || pErr == -1) && TMath::Abs(p->Parameter(1)-tau)/pErr->Parameter(1) > compLimit)
    {
        std::cout << "Compatibilità troppo grande (" << TMath::Abs(p->Parameter(1)-tau)/pErr->Parameter(1) << ">" << compLimit << ")" << std::endl;
        vBool[0] = false;
    }  
    if ( !(p == -1 || pErr == -1) && p->Parameter(2) > factorLimit*pErr->Parameter(1))
    {
        std::cout << "Sigma troppo grande rispetto all'errore (" << p->Parameter(2) << ">" << factorLimit*pErr->Parameter(1) << ")" << std::endl;
        vBool[0] = false;
    }   
 
    p    = fitGaus(hFitTauShortL2);
    pErr = fitGaus(hFitErrTauShortL2);
    
    std::cout << "p " << p << "   pErr " << pErr << std::endl;    
    
    if ( p == -1 || pErr == -1 ) 
    {
        std::cout << "Tau- fit failed. Bool set to 'false'." << std::endl;
        vBool[1] = false;
    }
/*    
    if ( !(p == -1 || pErr == -1) &&
          (pErr->Parameter(1)/p->Parameter(1) > errLimit || 
           TMath::Abs(p->Parameter(1)-tau)/pErr->Parameter(1) > compLimit ||
           p->Parameter(2) > factorLimit*pErr->Parameter(1)                ) ) vBool[1] = false;
*/
    if ( !(p == -1 || pErr == -1) && pErr->Parameter(1)/p->Parameter(1) > errLimit)
    {
        std::cout << "Errore relativo troppo grande (" << pErr->Parameter(1)/p->Parameter(1) << ">" << errLimit << ")" << std::endl;
        vBool[1] = false;
    } 
    if ( !(p == -1 || pErr == -1) && TMath::Abs(p->Parameter(1)-taucorto)/pErr->Parameter(1) > compLimit)
    {
        std::cout << "Compatibilità troppo grande (" << TMath::Abs(p->Parameter(1)-taucorto)/pErr->Parameter(1) << ">" << compLimit << ")" << std::endl;
        vBool[1] = false;
    }  
    if ( !(p == -1 || pErr == -1) && p->Parameter(2) > factorLimit*pErr->Parameter(1))
    {
        std::cout << "Sigma troppo grande rispetto all'errore (" << p->Parameter(2) << ">" << factorLimit*pErr->Parameter(1) << ")" << std::endl;
        vBool[1] = false;
    }

    p    = fitGaus(hFitRL2);
    pErr = fitGaus(hFitErrRL2);
    
    std::cout << "p " << p << "   pErr " << pErr << std::endl;
    
    if ( p == -1 || pErr == -1 ) 
    {
        std::cout << "R fit failed. Bool set to 'false." << std::endl;
        vBool[2] = false;
    }
/*
    if ( !(p == -1 || pErr == -1) &&
          (pErr->Parameter(1)/p->Parameter(1) > errLimit || 
           TMath::Abs(p->Parameter(1)-tau)/pErr->Parameter(1) > compLimit ||
           p->Parameter(2) > factorLimit*pErr->Parameter(1)                ) ) vBool[2] = false;
*/
    if ( !(p == -1 || pErr == -1) && pErr->Parameter(1)/p->Parameter(1) > errLimit)
    {
        std::cout << "Errore relativo troppo grande (" << pErr->Parameter(1)/p->Parameter(1) << ">" << errLimit << ")" << std::endl;
        vBool[2] = false;
    } 
    if ( !(p == -1 || pErr == -1) && TMath::Abs(p->Parameter(1)-R)/pErr->Parameter(1) > compLimit)
    {
        std::cout << "Compatibilità troppo grande (" << TMath::Abs(p->Parameter(1)-R)/pErr->Parameter(1) << ">" << compLimit << ")" << std::endl;
        vBool[2] = false;
    }  
    if ( !(p == -1 || pErr == -1) && p->Parameter(2) > factorLimit*pErr->Parameter(1))
    {
        std::cout << "Sigma troppo grande rispetto all'errore (" << p->Parameter(2) << ">" << factorLimit*pErr->Parameter(1) << ")" << std::endl;
        vBool[2] = false;
    }

    if ( vBool[0] == false || vBool[1] == false || vBool[2] == false ) vBool[3] = false;
 
    p = fitGaus(hFitBL2); 
    pErr = fitGaus(hFitErrBL2);
    
    // DRAW
    TFile f("matrix_canvas.root","UPDATE");
    
    std::ostringstream strs;
    strs << B << tau << integrale << taucorto << R << RebFactor;
    std::string str = strs.str();
    
    std::string set = "set" + str;
    std::string canName = "Simulazione MC (" + std::to_string(Nsim) + " simulazioni)";
    TCanvas can2( ("METODO2_"+set).c_str(), (canName + " METODO 2").c_str(), 1200 , 700 );
    can2.Divide(4,2);
    can2.cd(1);
        hFitTauL2.Draw();
        TLine lTau2(tau,0,tau,hFitTauL2.GetMaximum());
	    lTau2.SetLineColor(kRed);
        lTau2.Draw();
    can2.cd(2);
        hFitTauShortL2.Draw();
        TLine lTauShort2(taucorto,0,taucorto,hFitTauShortL2.GetMaximum());
	    lTauShort2.SetLineColor(kRed);
        lTauShort2.Draw();
    can2.cd(3);
        hFitBL2.Draw();
        TLine lB2(B,0,B,hFitBL2.GetMaximum());
	    lB2.SetLineColor(kRed);
        lB2.Draw();
    can2.cd(4);
        hFitRL2.Draw();
        TLine lR2(R,0,R,hFitRL2.GetMaximum());
	    lR2.SetLineColor(kRed);
        lR2.Draw();
    can2.cd(5);
        hFitErrTauL2.Draw();
    can2.cd(6);
        hFitErrTauShortL2.Draw();
    can2.cd(7);
        hFitErrBL2.Draw();
    can2.cd(8);
        hFitErrRL2.Draw();
        
    // save can2 into ROOT file
    f.WriteTObject(&can2);
    f.Close();
    // save .pdf
    can2.SaveAs("plot.pdf");
   
    return vBool;
}
//--------------- fine main ----------------------
 
TFitResultPtr fitGaus(TH1D h){
   int dim = h.GetNbinsX();

   std::string name = h.GetName();

   double min = h.GetBinCenter(1);
   double max = h.GetBinCenter(dim);
   TFitResultPtr p = h.Fit("gaus","SQN","",min,max);
   //double mean = p->Parameter(1);
   //double err  = p->ParError(1);
   //double sigma = p->Parameter(2);
   //double errsigma = p->ParError(2);

//   std::cout << "mean = "       << mean  << " +- " << err << std::endl
//	         << "sigma = "      << sigma << " +- "<< errsigma << std::endl
//             << "mean/sigma = " << "("   << sigma*100/mean << "%)" << std::endl;
   return p;
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
