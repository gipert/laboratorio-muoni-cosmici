int T = 100;		// defining some parameters: TIME
double L = 1.83;		//LENGTH
double D = 0.20;		//WIDTH
double H = 0.08;		//DISTANCE
double sp = 0.027;		//THICKNESS
int FLUSSO = 130;		//FLUX
const int N_events = 4758;
//1000;		//NUMBER OF EVENTS
const double PI = atan(1) * 4;
int M = 4 / PI;		// this is the constant in front of the distribution; see later

void Sim1D(int s)
{
	double x[N_events], z[N_events], zfin[N_events], xfin[N_events], theta[N_events]; // all the necessary arrays
	TRandom *rand = new TRandom();
	rand->SetSeed(s); // to set the seed
	for (int i = 0; i < N_events; i++) {
		//rand->SetSeed(0);		//forget about this, it's wrong; I kept it only to remember not to use it
		//if (s == 1) x[i] = rand->Uniform(-L*s/H, L+ L*s / H);		//I tried this, but it didn't work; I rather used the uniform distribution even fot the z axis
		x[i] = rand->Uniform(0, L);		//generates N_events values along X
		//y[i] = rand->Uniform(0, D);
		z[i] = 0;
		//phi[i] = rand->Uniform(0, 2 * PI);
	}
	int k = 0;		// now we try and find the cosine squared distribution
	double uniform, u;
	while (k < N_events) {
		uniform = 0;
		u = 2;
		while (M*u > M*cos(uniform)*cos(uniform)) {		//rejection condition
														//rand->SetSeed(0);
			u = rand->Uniform(0, 1);
			uniform = rand->Uniform(-PI/2, PI / 2);
		};
		theta[k] = uniform;
		k++;
	}
	double counter = 0;		// counter for coincidence events
	for (int i = 0; i < N_events; i++) {
		xfin[i] = x[i] - (H + z[i]) * tan(theta[i]);		// these should be correct
		zfin[i] = -H;
		//yfin[i] = y[i] - (H + z[i]) * tan(theta[i]) * sin(phi[i]);
		if (0 <= xfin[i] && xfin[i] <= L) counter++;
	}

	cout << N_events << " " << counter << " " << counter / N_events * 100 << "%" << endl;		// this ends the program



	TCanvas *c1 = new TCanvas("c1", "", 200, 10, 600, 400);		//finally the plot
	c1->SetFillColor(10);
	c1->SetGrid();
	TGraph *grafico = new TGraph(N_events, x, z);
	TGraph *grafico2 = new TGraph(N_events, xfin, zfin);
	grafico2->SetMarkerColor(9);
	//grafico2->SetMarkerStyle(23);
	grafico2->SetFillStyle(0);
	grafico->SetMarkerColor(2);
	//grafico->SetMarkerStyle(19);	
	grafico->SetFillStyle(0);
	grafico->SetTitle("S1");
	grafico2->SetTitle("S2");
	TMultiGraph *mg = new TMultiGraph("mg", "Slabs coincidence, 1D, H = 0.08 m");
	mg->Add(grafico, "*");
	mg->Add(grafico2, "*");
	mg->Draw("A");
	mg->GetXaxis()->SetTitle("X-Axis [m]");
	mg->GetYaxis()->SetTitle("Z-Axis [m]");
	mg->SetMinimum(-0.1.);
	mg->SetMaximum(0.02);
	mg->Draw("A*");
	c1->BuildLegend();
	return;
	}

