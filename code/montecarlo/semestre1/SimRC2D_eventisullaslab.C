int T = 100;		// defining some parameters: TIME
double L = 1.83;		//LENGTH
double D = 0.20;		//WIDTH
double H = 0.08;		//DISTANCE
double sp = 0.027;		//THICKNESS
int FLUSSO = 130;		//FLUX
const int N_events = 4758;
//1000;		//NUMBER OF EVENTS
const double PI = atan(1) * 4;
int M = 4/PI;		// this is the constant in front of the distribution; see later

void SimRC2D_eventisullaslab (int s)
{
	double x[N_events], y[N_events], z[N_events], zfin[N_events], xfin[N_events], yfin[N_events], phi[N_events], theta[N_events]; // all the necessary arrays
	TRandom *rand = new TRandom();
	rand->SetSeed(s); // to set the seed
	for (int i = 0; i < N_events; i++) {
		//rand->SetSeed(0);		//forget about this, it's wrong; I kept it only to remember not to use it
		//if (s == 1) x[i] = rand->Uniform(-L*s/H, L+ L*s / H);		//I tried this, but it didn't work; I rather used the uniform distribution even fot the z axis
		x[i] = rand->Uniform(0, L);		//generates N_events values along X
		y[i] = rand->Uniform(0, D);
		z[i] = 0;
		phi[i] = rand->Uniform(0, 2 * PI);
	}
	int k = 0;		// now we try and find the cosine squared distribution
	double uniform, u;
	while (k < N_events) {
		uniform = 0;
		u = 2;
		while (M*u > M*cos(uniform)*cos(uniform)*sin(uniform)) {		//rejection condition
			//rand->SetSeed(0);
			u = rand->Uniform(0, 1);
			uniform = rand->Uniform(0, PI / 2);
		};
		theta[k] = uniform;
		k++;
	}
	double counter = 0;		// counter for coincidence events
	for (int i = 0; i < N_events; i++) {
		xfin[i] = x[i] - (H + z[i]) * tan(theta[i]) * cos(phi[i]);		// these should be correct
		yfin[i] = y[i] - (H + z[i]) * tan(theta[i]) * sin(phi[i]);
		if (0 <= xfin[i] && xfin[i] <= L && yfin[i] >= 0 && yfin[i] <= D) counter++;
	}
	
	cout << N_events << " " << counter << " " << counter/N_events * 100  << "%" << endl;		// this ends the program
	
	const int n = 2 * N_events;		//I write the files all in one file
	double xfinal[n];
	double yfinal[n];
	double zfinal[n];
	for (int i = 0; i < N_events; i++) {
		xfinal[i] = x[i];
		yfinal[i] = y[i];
		zfinal[i] = z[i];
		xfinal[i + N_events] = xfin[i];
		yfinal[i + N_events] = yfin[i];
		zfinal[i + N_events] = -H;
	}
		
	TCanvas *c1 = new TCanvas("c1", "", 200, 10, 600, 400);		//finally the plot
	c1->SetFillColor(10);
	c1->SetGrid();
	TGraph2D *grafico = new TGraph2D(n, xfinal, yfinal, zfinal);
	grafico->SetTitle("Slabs coincidence, 2D, H = 0.08 m");
	grafico->SetMarkerColor(2);
	grafico->GetXaxis()->SetTitle("X-Axis [m]");
	grafico->GetYaxis()->SetTitle("Y-Axis [m]");
	grafico->GetZaxis()->SetTitle("Z-Axis [m]");
	grafico->Draw("PCOL");
	return;
}
	

