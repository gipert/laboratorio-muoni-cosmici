#define	m	1 //5.144E-03
#define calib 5.144E-03
#define sigmacalib 6.015E-6		// dati del fit di calibrazione canale - tempo (in microsecondi!)
#define q	0-150 //7.734E-01-1.635	 
#define fitMin  30 //30
#define fitMax  2138
#define Tauptrue 2.1969
#define sTauptrue 0.0000021
#define Taumtrue 0.88
#define sTaumtrue 0.01
#define RatioT 1.261
#define sRatioT 0.009
#define sigmaxstar 0 //incertezza sullo 0

// midValue ottimale = 1000
// rebin ottimale = 12
												
void lifetime2nocalib_Bkgr_sub( char* fileName1 , char* fileName2 , char* fileName3 ,/*char* filename4 ,*/ int rebin = 1 , double midValue = 1000 ) {		
	ifstream file1(fileName1);
	ifstream file2(fileName2);
	ifstream file3(fileName3);
	//ifstream file4(filename4);
	
	// definiamo gli estremi dell'istogramma
	Double_t min = 0*m + q;
	Double_t max = 4096*m + q;
	TH1F * data = new TH1F( "data" , "data" , 4096 , min , max);
		data->SetTitle("Spettro globale");
		data->SetYTitle("dN/dx");
		data->SetXTitle("Delay (canali)");
		data->SetLineColor(kBlack);
		data->SetStats(kFALSE);
		data->GetYaxis()->SetTitleOffset(1.5);
		data->SetLineColor(kBlue);

	TCanvas * can = new TCanvas( "ADCdataCan" , "ADC data" , 1500 , 600 );
	can->Divide(2);
	can->cd(1);
	gPad->SetGrid();

	int height1, height2, height3/*, height4*/;
	string line;
	int totevents = 0;

	// data retrieving
	for ( int i = 0 ; i < 8 ; i++ ) 
	{ 
		std::getline(file1,line); 
		std::getline(file2,line); 
		std::getline(file3,line); 
		//std::getline(file4,line); 
	} // skip first lines
	for ( int i = 0 ; i < 4096 ; i++ ) {
	
		file1 >> height1;
		file2 >> height2;
		file3 >> height3;
		//file4 >> height4;
		data->SetBinContent( i , height1+ height2+ height3/*+ height4*/);
		//totevents += height1+height2+height3+height4; // tot events
	}
/*
	// compute mean background
	int binThr = 2500-q; // lower level
	double mean = 0;
	for ( int i = binThr ; i < 4096 ; i++ ) {
		
		mean += data->GetBinContent(i);
		

		}

	mean      /= 4096 - binThr; // bkg subtraction for a single bin
	double errmean = mean/ (sqrt(4096-binThr));

	totevents -= (int)(mean*data->GetSize());// bkg events subtraction from total number

	// background subtraction
	int newBinContent = 0;
	for ( int i = 0 ; i < data->GetSize() - 2 ; i++ ) {

		newBinContent = data->GetBinContent(i) - mean;
		data->SetBinContent( i , newBinContent );
	}
*/

	// dump number of events
	ostringstream convert;
	convert << totevents;
	string eventsstring = convert.str();

	TText * text = new TText( 0.5 , 0.94 , (eventsstring + " events").c_str());
	text->SetTextSize(0.04);
	text->SetTextAlign(22);
	text->SetNDC();
	//text->Draw();
	
	//rebin
	data->Rebin(rebin);
	//data->Scale(1./rebinFactor);

	
	// compute mean background after rebinning
	int binThr = data->GetSize()-2;
	int k = 0;
	double center = 3001;
	while (center>3000)
	{
		center = data->GetBinCenter(binThr);
		binThr--;
		k++;
	}

	double mean = 0;
	int j;
	for(j = binThr; j < (data->GetSize()-2) ; j++)
	{

		mean += data->GetBinContent(j);	
		//cout << "data->GetBinContent(j) = " << data->GetBinContent(j) << endl;
		//cout << "mean = " << mean << endl;
	}

	cout << "\ndata->GetBinCenter("<<binThr<<") = " << data->GetBinCenter(binThr) << endl;
	cout << "Passaggi = "<< k << endl;
	cout << "binThr = " << binThr << endl;
	cout << "Canali = " << data->GetSize()-2 << endl;

	cout << "Sum = " << mean << endl;
	mean = mean/k;
	cout << "Bkgr = " << mean << endl;
	// background subtraction
	int newBinContent = 0;
	for ( int i = 0 ; i < data->GetSize() - 2 ; i++ ) {

		newBinContent = data->GetBinContent(i) - mean;
		data->SetBinContent( i , newBinContent );
				

	}

	// fit
	//TF1 * fminu = new TF1( "fminu" , "[0]+[1]*x", 0 , 20);
	//TF1 * fplus = new TF1( "fplus" , "[0]*exp(-x/[1])+[2]", q , 4096+q);
	TF1 * fplus = new TF1( "fplus" , "[0]*exp(-x/[1])", q , 4096+q);
	//fplus->FixParameter(2,mean);
	fplus->SetParameter(0,300);
	fplus->SetParameter(1,417);
	//fminu->SetLineColor(kBlue);

	//TFitResultPtr resultminu = dataLog->Fit( "fminu" , "SQR" , "" , fitMin , midValue );
	TFitResultPtr resultplus = data->Fit( "fplus" , "SQR" , "" , midValue , fitMax );
	
	//double tauM = -1./resultminu->Value(1);
	double tauP = resultplus->Value(1);
	double Aplus   = resultplus->Value(0);
	double sigmaAplus = sqrt(pow(resultplus->Error(0),2)+ pow(Aplus*sigmaxstar/tauP,2)) / rebin;
	double Nplus = Aplus * tauP / rebin; 
	double sigmaNplus = sqrt(pow(tauP, 2)*pow(resultplus->Error(0), 2) + pow(Aplus, 2)*pow(resultplus->Error(1), 2)+pow(Aplus,2)*pow(sigmaxstar,2) +2 * tauP*Aplus*resultplus->CovMatrix(0, 1))/rebin; 
	double sigmataupcalib = sqrt(pow(calib, 2)*pow(resultplus->Error(1), 2) + pow(resultplus->Value(1), 2)*pow(sigmacalib, 2));


	cout << endl;
	cout << "Primo fit a spanne (DX):\n";
	//cout << "tau_- = " << tauM << " +- " << resultminu->Error(1)/(tauM*tauM) << endl;
	cout << "tau+ = " << tauP << " +- " << resultplus->Error(1) << endl;
	cout << "tau+ calibrato = " << tauP*calib << " +- " << sigmataupcalib << endl;
	cout << "A+ = " << Aplus / rebin << " +- " << sigmaAplus << endl;
	//cout << "N+ = " << Nplus << " +- " << sigmaNplus << endl;
	cout << "Covariance check = " << resultplus->CovMatrix(0, 1) << endl<< endl;

	// draw
	//TCanvas * can = new TCanvas( "ADCdataCan" , "ADC data" , 1500 , 600 );
	//TCanvas * canLog = new TCanvas( "ADCdataCan (log)" , "ADC data (log)" , 1000 , 500 , 2050 , 725 ); // AppleTV
	
//	can->Divide(2);
	//can->cd();
	//gPad->SetGrid();
	//data->Draw("");
	//fplus->DrawCopy("SAME");	
	//canLog->SaveAs("datalog.png");
	
	// mu+ subtraction:
	// second fit to find A, t0
	//TF1 * f = new TF1( "f" , "log([0])+[1]/[2]-x/[2]" , 0 , 20); 
	//f->FixParameter( 2 , tauP );
	//f->SetParNames( "A" , "t0" );

	//TFitResultPtr result = dataLog->Fit ( "f" , "SQN" , "" , midValue , fitMax );
	
	//cout << "Secondo fit della seconda parte dello spettro (tau_-) per la stima dell'esponenziale da sottrarre:\n";
	//cout << "A = " << A << " +- " << result->Error(0) << endl;
	//cout << "t0 = " << t0 << " +- " << result->Error(1) << endl << endl;

	can->cd(2);
	gPad->SetGrid();

	TH1F *  dataSub = new TH1F( "dataSub" , "DataSub" , data->GetSize()-2 , min , max );
		dataSub->SetTitle("Spettro dei #mu^{-}");
		dataSub->SetXTitle("Delay (canali)");
		dataSub->SetYTitle("dN/dt");
		dataSub->SetStats(kFALSE);
		dataSub->GetYaxis()->SetTitleOffset(1.5);

	// subtraction
	for ( int i = 0 ; i < dataSub->GetSize()-2 ; i++ ) 
	{
		double content = data->GetBinContent(i) - fplus->Eval(data->GetBinCenter(i));
		dataSub->SetBinContent( i , content );
	}

	// rebinning & rescaling
	//dataSub->Rebin(rebin2);
	//dataSub->Scale(1./rebin2);
	
//	can->cd(2);

	// fitting the mu-
	//TF1 * fminu2 = new TF1( "fminu2" , "[0]*exp(-x/[1])+[2]*exp(-x/[3])+[4]", q , 4096+q);
	TF1 * fminu2 = new TF1( "fminu2" , "[0]*exp(-x/[1])", q , 4096+q);
	fminu2->SetParameter(0,140);
	fminu2->SetParameter(1,194);
	TFitResultPtr resultminu2 = dataSub->Fit( "fminu2" , "SQR+" , "" , fitMin , midValue );

	double tauM2 = resultminu2->Value(1);
	double Aminus = resultminu2->Value(0);
	double sigmaAminus = sqrt(pow(resultminu2->Error(0), 2) + pow(Aminus*sigmaxstar / tauM2, 2)) / rebin;
	double Nminus = Aminus * tauM2/ rebin;
	double sigmaNminus = sqrt(pow(tauM2, 2)*pow(resultminu2->Error(0), 2) + pow(Aminus, 2)*pow(resultminu2->Error(1), 2) + pow(Aminus, 2)*pow(sigmaxstar, 2) + 2 * tauM2*Aminus*resultminu2->CovMatrix(0, 1))/rebin; //qui non si è tenuto conto dell'xstar
	double sigmataumcalib = sqrt(pow(calib, 2)*pow(resultminu2->Error(1), 2) + pow(resultminu2->Value(1), 2)*pow(sigmacalib, 2));

	cout << "Stima finale di tau- (SX):\n";
	cout << "tau- = " << tauM2 << " +- " << resultminu2->Error(1) << endl;
	cout << "tau- calibrato = " << tauM2*calib << " +- " << sigmataumcalib << endl;
	cout << "A- = " << Aminus / rebin << " +- " <<sigmaAminus << endl; 
	//cout << "N- = " << Nminus << " +- " << sigmaNminus << endl;
	cout << "Covariance check = " << resultminu2->CovMatrix(0, 1) << endl<< endl;
	
	double ratio = Nplus / Nminus;
	double sigmaratio = sqrt(pow(sigmaNplus,2) / pow(Nminus,2) + pow(Nplus,2)*pow(sigmaNminus,2) / pow(Nminus,4));
	double ratioA = Aplus / Aminus;
	double sigmaratioA = sqrt(pow(sigmaAplus, 2) / pow(Aminus/rebin, 2) + pow(Aplus/rebin, 2)*pow(sigmaAminus, 2) / pow(Aminus/rebin, 4));
	//cout << "Population ratio = " << ratio << " +- " << sigmaratio << endl;
	cout << "Ratio of A factors: " << ratioA << " +- " << sigmaratioA << endl;

	cout << endl;
	cout << "Compatibilita'" << endl;
	cout << "tau+: " << -(tauP*calib - Tauptrue) / sqrt(pow(sTauptrue, 2)+pow(sigmataupcalib, 2)) << endl;
	cout << "tau-: " << (tauM2*calib - Taumtrue) / sqrt(pow(sTaumtrue, 2)+pow(sigmataumcalib, 2)) << endl;
	cout << "Ratio: " << (ratioA - RatioT) / sqrt(pow(sigmaratioA, 2) + pow(sRatioT, 2)) << endl;

	double deltap, deltam, deltaA;
	deltap = tauP*calib-Tauptrue;
	deltam = tauM2*calib-Taumtrue;
	deltaA = ratioA - RatioT;
	cout<< "Delta+ = " << deltap << endl;
	cout<< "Delta- = " << deltam << endl;
	cout<< "DeltaA = " << deltaA << endl;

	//cout << "Population ratio: " << (ratio - RatioT) / sqrt(pow(sigmaratio, 2) + pow(sRatioT, 2)) << endl;

	
	// draw
	//TCanvas * can = new TCanvas( "canNew" , "canNew" , 1 ) ;
	//data->DrawCopy();

	//cout << "Il numero totale di eventi escluso il fondo e': " << data->Integral(-q/rebin, 4096/rebin) - mean*(4096+2*q) << endl; // - mean*(4096+q) 

	return;
}
