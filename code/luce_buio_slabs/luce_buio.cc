// ROOTMacro per il rate di conteggio dei PMT con luce e buio
// Usage: .x luce_buio.cc("nome_file_con_estensione")
// Formattazione file: 	1col: nr. PMT
// 			2col: conteggi luce
// 			3col: conteggi buio

void luce_buio( string name ) {

	TCanvas * can = new TCanvas("luce_buio","LuceBuio",1200,600);
		can->SetGrid();
	TGraph * grluce = new TGraph( name.c_str() ,"%lg %lg %*lg");
	  	grluce->SetLineStyle(2);
		grluce->SetLineColor(2);
	        grluce->SetMarkerStyle(22);
    grluce->SetMarkerSize(1.3);
	        grluce->SetMarkerColor(2);
	TGraph * grbuio = new TGraph( name.c_str() ,"%lg %*lg %lg");
		grbuio->SetLineStyle(2);
		grbuio->SetLineColor(4);
		grbuio->SetMarkerStyle(22);
    grbuio->SetMarkerSize(1.3);
		grbuio->SetMarkerColor(4);
	
		double x,x_,y,y_;
	for ( int i = 0 ; i < 14 ; i++ ) {
		grluce->GetPoint(i,x,y);
		grbuio->GetPoint(i,x_,y_);
		grluce->SetPoint(i,x,y/100);
		grbuio->SetPoint(i,x_,y_/100);
	}	
		grluce->SetTitle("");
		grluce->GetHistogram()->SetXTitle("PMT number");
		grluce->GetHistogram()->SetYTitle("Counts [Hz]");	

	can->cd();
	grluce->Draw("APL");
	grbuio->Draw("PLSAME");

	TLegend * leg = new TLegend(0.66,0.15,0.86,0.27);
	leg->AddEntry(grluce,"Luce","P");
	leg->AddEntry(grbuio,"Buio","P");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.035);
	leg->Draw();

return;
}

