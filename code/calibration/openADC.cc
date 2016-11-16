#define	m	1//5.144E-03
#define q	0//7.734E-01-1.635	// dati del fit di calibrazione canale - tempo (in microsecondi!)

TH1F * data; 

void openADC( string fileName1 , int rebin = 1) {

	int rebinFactor = rebin;

	ifstream file1(fileName1.c_str());
	
	Double_t min = 0*m + q;
	Double_t max = 4096*m + q;
	data = new TH1F( "data" , "data" , 4096 , min , max);
		data->SetTitle("");
		data->SetYTitle("Counts");
		data->SetXTitle("Delay (#mus)");
		data->SetLineColor(4);
		data->SetStats(kFALSE);

	int height1;
	string line;
	int totevents = 0;

	// data retrieving
	for ( int i = 0 ; i < 8 ; i++ ) std::getline(file1,line); // skip first lines
	for ( int i = 0 ; i < 4096 ; i++ ) {
	
		file1 >> height1;
		data->SetBinContent( i , height1 );
		totevents += height1; // tot events
	}
	
	data->Rebin(rebinFactor); // rebin histograms

	TCanvas * can = new TCanvas( "ADCdataCan" , "Adc data" , 1500 , 600);
	can->cd();
	//can->SetLogy();
	data->Draw();

	// dump number of events
	ostringstream convert;
	convert << totevents;
	string eventsstring = convert.str();

	TText * text = new TText( 0.5 , 0.94 , (eventsstring + " events").c_str());
	text->SetTextSize(0.04);
	text->SetNDC();
	text->SetTextAlign(22);
	text->Draw();

	return;
}

void openADC( char* fileName1 , char*fileName2 , int rebin = 1) {

	int rebinFactor = rebin;

	ifstream file1(fileName1);
	ifstream file2(fileName2);
	
	Double_t min = 0*m + q;
	Double_t max = 4096*m + q;
	TH1F * data = new TH1F( "data" , "data" , 4096 , min , max);
		data->SetTitle("");
		data->SetYTitle("Counts");
		data->SetXTitle("Delay (#mus)");
		data->SetLineColor(4);
		data->SetStats(kFALSE);

	int height1, height2;
	string line;
	int totevents = 0;

	// data retrieving
	for ( int i = 0 ; i < 8 ; i++ ) { std::getline(file1,line); std::getline(file2,line); } // skip first lines
	for ( int i = 0 ; i < 4096 ; i++ ) {
	
		file1 >> height1;
		file2 >> height2;
		data->SetBinContent( i , height1+height2 );
		totevents += height1+height2; // tot events
	}

	data   ->Rebin(rebinFactor); // rebin histograms

	TCanvas * can = new TCanvas( "ADCdataCan" , "Adc data" , 1500 , 600);
	can->cd();
	can->SetLogy();
	data->Draw();

	// dump number of events
	ostringstream convert;
	convert << totevents;
	string eventsstring = convert.str();

	TText * text = new TText( 0.5 , 0.94 , (eventsstring + " events").c_str());
	text->SetTextSize(0.04);
	text->SetTextAlign(22);
	text->SetNDC();
	text->Draw();
}
