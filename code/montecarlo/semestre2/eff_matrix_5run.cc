void eff_matrix ()
{
	const int n = 5;
	float effTauPlus [n] = {97.5309,93.8272,86.4198,84.1463,67.9012};
	float effTauMinus[n] = {97.5309,93.8272,87.6543,84.1463,67.9012};
	float effR       [n] = {97.5309,93.8272,82.716,78.0488,66.6667};
	float effTot     [n] = {97.5309,92.5926,82.716,78.0488,65.4321};
	
	TH1F* Tplus  = new TH1F ("Tplus",  "Efficienza", 5, 0.5, 5.5);
	TH1F* Tminus = new TH1F ("Tminus", "Efficienza", 5, 0.5, 5.5);
	TH1F* TR     = new TH1F ("TR",     "Efficienza"       , 5, 0.5, 5.5);
	TH1F* Ttot   = new TH1F ("Ttot",   "Efficienza"  , 5, 0.5, 5.5);
	Tplus   ->SetStats(0);
	Tminus  ->SetStats(0);
	TR      ->SetStats(0);
	Ttot    ->SetStats(0);
	Tplus   ->GetXaxis()->SetTitle("Ampiezza del range di variazione dei parametri (%)");
	Tplus   ->GetYaxis()->SetTitle("Efficienza (%)");
	Tplus	->GetYaxis()->SetTitleOffset(1);
	Tminus  ->GetXaxis()->SetTitle("Ampiezza del range di variazione dei parametri (%)");
	Tminus  ->GetYaxis()->SetTitle("Efficienza (%)");
	Tminus  ->GetYaxis()->SetTitleOffset(1);
	TR      ->GetXaxis()->SetTitle("Ampiezza del range di variazione dei parametri (%)");
	TR      ->GetYaxis()->SetTitle("Efficienza (%)");
	TR		->GetYaxis()->SetTitleOffset(1);
	Ttot    ->GetXaxis()->SetTitle("Ampiezza del range di variazione dei parametri (%)");
	Ttot    ->GetYaxis()->SetTitle("Efficienza (%)");
	Ttot	->GetYaxis()->SetTitleOffset(1);
	Tplus   ->GetYaxis()->SetRangeUser(50,100);
	Tminus  ->GetYaxis()->SetRangeUser(50,100);
	TR      ->GetYaxis()->SetRangeUser(50,100);
	Ttot    ->GetYaxis()->SetRangeUser(50,100);
       Tplus   ->GetXaxis()->SetNdivisions(5);
	   Tminus  ->GetXaxis()->SetNdivisions(5);
        TR      ->GetXaxis()->SetNdivisions(5);
        Ttot    ->GetXaxis()->SetNdivisions(5);
		char str[10];
	for (int k=0; k<5; k++)
	{
			int m = (k + 1) * 10;
			sprintf (str, "%d", m);
	        Tplus   ->SetBinContent(k+1,effTauPlus[k]);
	        Tminus  ->SetBinContent(k+1,effTauMinus[k]);
	        TR      ->SetBinContent(k+1,effR[k]);
	        Ttot    ->SetBinContent(k+1,effTot[k]);	     
			Tplus->GetXaxis()->SetBinLabel(k+1, str);
			Tminus->GetXaxis()->SetBinLabel(k + 1, str);
			TR->GetXaxis()->SetBinLabel(k + 1, str);
			Ttot->GetXaxis()->SetBinLabel(k + 1, str);
	}
	
	TCanvas *c1 = new TCanvas("c1","Efficienze",700,800);
	//c1      ->SetFillColor(390);
	gPad    ->SetGrid();
	
	Tplus   ->SetBarWidth(0.2);
	Tplus   ->SetFillColor(433);
	Tplus   ->SetBarOffset(0.1);
	Tminus  ->SetFillColor(417);
	Tminus  ->SetBarWidth(0.2);
	Tminus  ->SetBarOffset(0.3);
	TR      ->SetFillColor(401);
	TR      ->SetBarWidth(0.2);
	TR      ->SetBarOffset(0.5);
	Ttot    ->SetFillColor(633);
	Ttot    ->SetBarWidth(0.2);
	Ttot    ->SetBarOffset(0.7);
	
        Tplus   ->Draw("bar1");
	Tminus  ->Draw("bar1,same");
	TR      ->Draw("bar1,same");
	Ttot    ->Draw("bar1,same");
	
	TLegend *legend = new TLegend(0.60,0.67,0.81,0.85);
        legend->AddEntry(Tplus  ,"#tau^{+}","f");
        legend->AddEntry(Tminus ,"#tau^{-}","f");
        legend->AddEntry(TR     ,"R"       ,"f");
        legend->AddEntry(Ttot   ,"Totale"  ,"f");
        legend->SetFillColor(390);
        legend->Draw();
	
        return;
}
