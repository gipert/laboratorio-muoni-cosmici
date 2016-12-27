#define randomV  0.06
#define scaleV  0.02
#define randomT  0.06

#include <iostream>
#include <math.h>
#include <fstream>

/*#include "Tension.h"
#include "Time.h"*/
#include <TGraph.h>

struct Tension{
	double value;
	double error;
	double Vdiv;

};
struct Time{
	double value;
	double error;
	double secdiv;
};

using namespace std;
void cal(string name)
{
	ifstream in(name.c_str());
	char l;							// indicatore del tipo di dato (tensione, tempo, canale)

	const Int_t n = 50; 					// dimensione array
	Double_t tens[n];					// array tensioni
	Double_t errtens[n];					// array errori tensioni
	Double_t time[n];					// array tempi
	Double_t errtime[n];					// array errori tempi
	Double_t channel[n];					// array canali
	//Double_t* errchan = new Double_t[n];			// array errori canali	--> non hanno errore...

	Int_t v = 0;						// contatori caselle array (v = tensioni) (t = tempi) (c = canali)
	Int_t t = 0;
	Int_t c = 0;
    float fwhm;
	while(in>> l)
	{
		if (l=='V' || l=='v') 				// tensioni
		{
			Tension* V = new Tension;
			in >> V->value >> V->Vdiv;
			//V->error = sqrt((randomV)*(V->Vdiv)*(randomV)*(V->Vdiv) + (scaleV)*(V->value)*(scaleV)*(V->value));
			V->error = (V->Vdiv/10)/(sqrt(12));
			tens[v] = V->value;
			errtens[v] = V->error;
			cout<<"V = " << tens[v] << " +- " << errtens[v] << endl;
			v++;	
		}
		if (l=='T' || l=='t')				// tempi
		{
			Time*    T = new Time;
			in >> T->value >> T->secdiv;
			//T->error = sqrt((randomT)*(T->secdiv)*(randomT)*(T->secdiv));
			T->error = (T->secdiv/10)/(sqrt(12));
			time[t] = T->value;
			errtime[t] = T->error;
			cout<<"T = " << time[t] << " +- " << errtime[t] << endl;
			t++;	
		}
		if (l=='C' || l=='c')				// canali
		{
			in >> channel[c] >> fwhm;
			cout<<"C = " << channel[c] << endl;
			c++;
		}
	}

        TF1 *linear    = new TF1("linear","pol1",0,100); 
	linear->SetParNames("Intercept","Slope");
	Double_t par[2];
	TF1* base = new TF1("base", "0",0,5000);
	if (v == 0)
	{
		TCanvas *c1 = new TCanvas("c1","Calibration",1100,400); 
        	c1->Divide(2,1,0.01,0.01);
		c1->cd(1);
		TGraphErrors* g = new TGraphErrors(t,channel,time,0,errtime);
		g->Draw("A*");
		g->Fit("linear");
		linear->GetParameters(&par[0]);
		Double_t ERRq = linear ->GetParError(0);
		Double_t ERRm = linear ->GetParError(1);
		Double_t* res    = new Double_t [t];
		Double_t* ERRres = new Double_t [t];
		for(int k=0; k<t; k++)										// residui
		{
			res[k] = time[k] - (par[0]+channel[k]*par[1]);
			ERRres[k] = sqrt(errtime[k]*errtime[k]+channel[k]*channel[k]*ERRm*ERRm+ERRq*ERRq);	
		}

    		g->SetTitle("Channel vs. Time");
   	        g->SetMarkerColor(kBlue);
   	        g->SetMarkerStyle(22);
   		g->SetMarkerSize(0.8);
   	        g->SetLineColor(kBlue);
		g->GetXaxis()->SetTitle("Channel");
		g->GetYaxis()->SetTitle("Time (s)");

		c1->cd(2);
		TGraphErrors* residui = new TGraphErrors(t,channel,res,0,ERRres);
		residui->Draw("A*");
		base->Draw("same");
    		residui->SetTitle("Residuals");
   	        residui->SetMarkerColor(kBlue);
   	        residui->SetMarkerStyle(22);
   		residui->SetMarkerSize(0.8);
   	        residui->SetLineColor(kBlue);
		residui->GetXaxis()->SetTitle("Channel");
		residui->GetYaxis()->SetTitle("Residuals (s)");
	}
	if (t == 0)
	{
		TCanvas *c1 = new TCanvas("c1","ADC Calibration",1100,400); 
        	c1->Divide(2,1,0.01,0.01);
		c1->cd(1);
		TGraphErrors* g = new TGraphErrors(v,channel,tens,0,errtens);
		g->Draw("A*");
		g->Fit("linear");
		linear->GetParameters(&par[0]);
		Double_t ERRq = linear ->GetParError(0);
		Double_t ERRm = linear ->GetParError(1);
		Double_t* res    = new Double_t [v];
		Double_t* ERRres = new Double_t [v];
		for(int k=0; k<c; k++)										// residui
		{
			res[k] = tens[k] - (par[0]+channel[k]*par[1]);
			ERRres[k] = sqrt(errtens[k]*errtens[k]+channel[k]*channel[k]*ERRm*ERRm+ERRq*ERRq);	
			cout << "res[" << k << "]= " <<res[k]<< endl;
			cout << "ERRres[" << k << "]= "<<ERRres[k] << endl;
		}

    		g->SetTitle("Channel vs. Tension");
   	        g->SetMarkerColor(kBlue);
   	        g->SetMarkerStyle(22);
   		g->SetMarkerSize(0.8);
   	        g->SetLineColor(kBlue);
		g->GetXaxis()->SetTitle("Channel");
		g->GetYaxis()->SetTitle("Tension (V)");

		c1->cd(2);
		TGraphErrors* residui = new TGraphErrors(v,channel,res,0,ERRres);
		residui->Draw("A*");
		base->Draw("same");
    		residui->SetTitle("Residuals");
   	        residui->SetMarkerColor(kBlue);
   	        residui->SetMarkerStyle(22);
   		residui->SetMarkerSize(0.8);
   	        residui->SetLineColor(kBlue);
		residui->GetXaxis()->SetTitle("Channel");
		residui->GetYaxis()->SetTitle("Residuals (V)");
	}
	if (c == 0)
	{
		TCanvas *c1 = new TCanvas("c1","TAC Calibration",1100,400); 
        	c1->Divide(2,1,0.01,0.01);
		c1->cd(1);
		TGraphErrors* g = new TGraphErrors(t,time,tens,errtime,errtens);
		g->Draw("A*");
		g->Fit("linear");
		linear->GetParameters(&par[0]);
		Double_t ERRq = linear ->GetParError(0);
		Double_t ERRm = linear ->GetParError(1);
		Double_t* res    = new Double_t [t];
		Double_t* ERRres = new Double_t [t];
		for(int k=0; k<t; k++)										// residui
		{
			res[k] = tens[k] - (par[0]+time[k]*par[1]);
			ERRres[k] = sqrt(errtens[k]*errtens[k]+par[1]*par[1]*errtime[k]*errtime[k]+time[k]*time[k]*ERRm*ERRm+ERRq*ERRq);	
		}

    		g->SetTitle("Time vs. Tension");
   	        g->SetMarkerColor(kBlue);
   	        g->SetMarkerStyle(22);
   		g->SetMarkerSize(0.8);
   	        g->SetLineColor(kBlue);
		g->GetXaxis()->SetTitle("Time (s)");
		g->GetYaxis()->SetTitle("Tension (V)");

		c1->cd(2);
		TGraphErrors* residui = new TGraphErrors(t,time,res,errtime,ERRres);
		residui->Draw("A*");
		base->Draw("same");
    		residui->SetTitle("Residuals");
   	        residui->SetMarkerColor(kBlue);
   	        residui->SetMarkerStyle(22);
   		residui->SetMarkerSize(0.8);
   	        residui->SetLineColor(kBlue);
		residui->GetXaxis()->SetTitle("Time (s)");
		residui->GetYaxis()->SetTitle("Residuals (V)");
		
	}
	return;
}



