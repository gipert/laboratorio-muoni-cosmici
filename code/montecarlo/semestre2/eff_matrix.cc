void eff_matrix ()
{
	float effTauPlus [5]  {84.1463,82.716,93.8272,87.6543,88.8889};
	float effTauMinus[5]  {84.1463,85.1852,93.8272,90.1235,90.1235};
	float effR       [5]  {78.0488,77.7778,88.8889,86.4198,86.4198};
	float effTot     [5]  {78.0488,77.7778,87.6543,86.4198,83.9506};
	
	TH1F* Tplus  = new TH1F ("Tplus",  "Efficienza", 5, 0.5, 5.5);
	TH1F* Tminus = new TH1F ("Tminus", "Efficienza", 5, 0.5, 5.5);
	TH1F* TR     = new TH1F ("TR",     "Efficienza"       , 5, 0.5, 5.5);
	TH1F* Ttot   = new TH1F ("Ttot",   "Efficienza"  , 5, 0.5, 5.5);
	Tplus   ->SetStats(0);
	Tminus  ->SetStats(0);
	TR      ->SetStats(0);
	Ttot    ->SetStats(0);
	Tplus   ->GetXaxis()->SetTitle("Run");
	Tplus   ->GetYaxis()->SetTitle("Efficienza (%)");
	Tminus  ->GetXaxis()->SetTitle("Run");
	Tminus  ->GetYaxis()->SetTitle("Efficienza (%)");
	TR      ->GetXaxis()->SetTitle("Run");
	TR      ->GetYaxis()->SetTitle("Efficienza (%)");
	Ttot    ->GetXaxis()->SetTitle("Run");
	Ttot    ->GetYaxis()->SetTitle("Efficienza (%)");
	Tplus   ->GetYaxis()->SetRangeUser(75,95);
	Tminus  ->GetYaxis()->SetRangeUser(75,95);
	TR      ->GetYaxis()->SetRangeUser(75,95);
	Ttot    ->GetYaxis()->SetRangeUser(75,95);
        Tplus   ->GetXaxis()->SetNdivisions(5);
        Tminus  ->GetXaxis()->SetNdivisions(5);
        TR      ->GetXaxis()->SetNdivisions(5);
        Ttot    ->GetXaxis()->SetNdivisions(5);
	
	for (int k=0; k<5; k++)
	{
	        Tplus   ->SetBinContent(k+1,effTauPlus[k]);
	        Tminus  ->SetBinContent(k+1,effTauMinus[k]);
	        TR      ->SetBinContent(k+1,effR[k]);
	        Ttot    ->SetBinContent(k+1,effTot[k]);	        
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
	
	TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
        legend->AddEntry(Tplus  ,"#tau^{+}","f");
        legend->AddEntry(Tminus ,"#tau^{-}","f");
        legend->AddEntry(TR     ,"R"       ,"f");
        legend->AddEntry(Ttot   ,"Totale"  ,"f");
        legend->SetFillColor(390);
        legend->Draw();
	
        return;
}
