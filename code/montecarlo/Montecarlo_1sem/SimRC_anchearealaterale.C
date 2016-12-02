int T = 100;		// defining some parameters: TIME
double L = 1.83;		//LENGTH
double D = 0.20;		//WIDTH
double H = 0.08;		//DISTANCE
double sp = 0.027;		//THICKNESS
int FLUSSO = 130;		//FLUX
double square = 0.5;	//SIZE OF THE SQUARE AROUND THE SLAB
int numero = int((L * D + 2 * sp*D + 2 * sp*L) * T * FLUSSO);		//CONSIDERING THE WHOLE AREA
const int N_events =  numero;		//NUMBER OF EVENTS 
const double PI = atan(1) * 4;
int M = 4/PI;		// this is the constant in front of the distribution; see later

void SimRC_anchearealaterale (int seed)
{
	double x[N_events], y[N_events], z[N_events], zfin[N_events], xfin[N_events], yfin[N_events], phi[N_events], theta[N_events]; // all the necessary arrays
	double uniform, u;
	TRandom3 *rand = new TRandom();
	rand->SetSeed(seed); // to set the seed
	const int n = 2 * N_events;		//I write the files all in one file to plot in the end
	double xfinal[n];
	double yfinal[n];
	double zfinal[n];
	double counterge = 0;		//to count good events
	double countentries = 0;		//to count the total entries
	int k = 0;		//dynamic counter
	int i = 0;
	while (i > -1) {
		//rand->SetSeed(0);		//forget about this, it's wrong; I kept it only to remember not to use it
		countentries++;
		if (counterge >= N_events) break;
		x[i] = rand->Uniform(-square,L+square);		//generates N_events values along X
		y[i] = rand->Uniform(-square, D+square);
		phi[i] = rand->Uniform(0, 2*PI);
		uniform = 0;
		u = 2;
		while (M*u > M*cos(uniform)*cos(uniform)) {		//rejection condition
			u = rand->Uniform(0, 1);
			uniform = rand->Uniform(0, PI / 2);
		};
		theta[i] = uniform;
		if (x[i] >= 0 && x[i] <= L && y[i] >= 0 && y[i] <= D) {
			counterge++;
			z[i] = 0;	
			xfinal[i] = x[i];
			yfinal[i] = y[i];
			k++;
			cout << k << endl;
			i++;
		}
		else if (x[i] < 0 && phi[i] >PI / 2 && phi[i] < PI && (y[i] - x[i] * tan(phi[i])) > 0 && (y[i] - x[i] * tan(phi[i])) < D && theta[i] > atan(x[i] / (cos(phi[i]) * sp))) {
			counterge++;
			z[i] = -x[i] / (cos(phi[i]) * tan(theta[i]));
			xfinal[i] = 0;
			yfinal[i] = (y[i] - x[i] * tan(phi[i])); 
			k++;
			cout << k << endl;
			i++;
		}
		else if (x[i] < 0 && phi[i] >PI / 2 && phi[i] < PI && (y[i] + x[i] * tan(PI-phi[i])) > D && x[i]-(D-y[i])/ tan(PI-phi[i]) > 0 && x[i] - (D - y[i]) / tan(PI-phi[i]) <L && theta[i] > atan((y[i]-D) / (sin(PI-phi[i]) * sp))) {
			counterge++;
			z[i] = -(y[i] - D) / (sin(PI - phi[i]) * tan(theta[i]));
			yfinal[i] = D; 
			xfinal[i] = x[i] - (D - y[i]) / tan(PI-phi[i]); 
			k++;
			cout << k << endl;
			i++;
		}
		else if (x[i] < 0 && phi[i] >PI && phi[i] < 3*PI/2 && (y[i] - x[i] * tan(phi[i])) > 0 && (y[i] - x[i] * tan(phi[i])) < D && theta[i] > atan(-x[i] / (cos(PI-phi[i]) * sp))) {
			counterge++;
			z[i] = -x[i] / (cos(phi[i]) * tan(theta[i]));
			xfinal[i] = 0; 
			yfinal[i] = (y[i] - x[i] * tan(phi[i]));
			k++;
			cout << k << endl;
			i++;
		}
		else if (x[i] < 0 && phi[i] >PI && phi[i] < 3 * PI / 2 && x[i] - y[i] / tan(phi[i]) > 0 && x[i] - y[i] / tan(phi[i]) < L&& (y[i] - x[i] * tan(phi[i])) < 0 && theta[i] > atan(-y[i] / (sin(phi[i]-PI) * sp))) {
			counterge++;
			z[i] = y[i] / (sin(phi[i] - PI) * tan(theta[i]));
			yfinal[i] = 0;
			xfinal[i] = x[i] - y[i] / tan(phi[i]); 
			k++;
			cout << k << endl;
			i++;
		}
		else if (x[i] > L && phi[i] > 0 && phi[i] < PI / 2 && (y[i] - (x[i] - L) * tan(phi[i])) > 0 && (y[i] - (x[i] - L) * tan(phi[i])) < D && theta[i] > atan((x[i] - L) / (cos(phi[i]) * sp))) {
			counterge++;
			z[i] = -(x[i] - L) / (cos(phi[i]) * tan(theta[i]));
			xfinal[i] = L; 
			yfinal[i] = y[i] - (x[i] - L) * tan(phi[i]); 
			k++;
			cout << k << endl;
			i++;
		}
		else if (x[i] > L && phi[i] > 0 && phi[i] < PI / 2 && (y[i] - (x[i] - L) * tan(phi[i])) > D && x[i] - (y[i] - D) / tan(phi[i]) > 0 && x[i] - (y[i] - D) / tan(phi[i]) < L && theta[i] > atan((y[i]-D) / (sin(phi[i]) * sp))) {
			counterge++;
			z[i] = -(y[i] - D) / (sin(phi[i]) * tan(theta[i]));
			yfinal[i] = D; 
			xfinal[i] = x[i] - (y[i] - D) / tan(phi[i]); 
			k++;
			cout << k << endl;
			i++;
		}
		else if (x[i] > L && phi[i] > 3 * PI / 2 && (y[i] + (x[i]-L) / tan(phi[i] - 3 * PI / 2)) > 0  && (y[i] + (x[i]-L) / tan(phi[i] - 3 * PI / 2)) < D && theta[i] > atan((x[i]-L) / (sin(phi[i] - 3 * PI / 2) * sp))) {
			counterge++;
			z[i]  = -(x[i] - L) / (sin(phi[i] - 3 * PI / 2) * tan(theta[i]));
			xfinal[i] = L; 
			yfinal[i] = y[i] + (x[i] - L) / tan(phi[i] - 3 * PI / 2); 
			k++;
			cout << k << endl;
			i++;
		}
		else if (x[i] > L && phi[i] > 3 * PI / 2 && x[i] + y[i] * tan(phi[i] - 3 * PI / 2) > 0 && x[i] + y[i] * tan(phi[i] - 3 * PI / 2) < L && (y[i] + (x[i] - L) / tan(phi[i] - 3 * PI / 2)) < 0 && theta[i] > atan(-y[i] / (sin(2*PI-phi[i]) * sp))) {
			counterge++;
			z[i] = y[i] / (sin(2 * PI - phi[i]) * tan(theta[i]));
			yfinal[i] = 0; 
			xfinal[i] = x[i] + y[i] * tan(phi[i]-3*PI/2); 
			k++;
			cout << k << endl;
			i++;
		}
		else if (y[i] < 0 && phi[i] > PI && phi[i] < 3 * PI/2 && (x[i] - y[i] / tan(phi[i]-PI)) > 0 && (x[i] - y[i] / tan(phi[i]-PI)) < L && theta[i] > atan(-y[i] / (sin(phi[i]-PI) * sp))) {
			counterge++;
			z[i] = y[i] / (sin(phi[i]-PI) * tan(theta[i]));
			xfinal[i] = x[i] - y[i] / tan(phi[i] - PI);
			k++;
			cout << k << endl;
			yfinal[i] = 0;
			i++;
		}
		//else if (y[i] < 0 && phi[i] > PI && phi[i] < 3 * PI / 2 && (x[i] - y[i] / tan(phi[i] - PI)) > 0 - D / (tan(phi[i] - PI)) && (x[i] - y[i] / tan(phi[i] - PI)) < 0 && theta[i] > atan(-x[i] / (cos(phi[i]-PI) * sp))) {
		//counterge++;
		//z[i] = -x[i] / (cos(phi[i] - PI) * tan(theta[i]));
		//xfinal[i] = 0;
		//yfinal[i] = -x[i]*tan(phi[i]-PI)+y[i];//correggere
		//k++;
		//cout << k << endl;
		//i++;
		//}
		else if (y[i] < 0 && phi[i] > 3 * PI / 2 && phi[i] < 2 * PI && (x[i] + y[i] / tan(2 * PI - phi[i])) > 0 && (x[i] + y[i] / tan(2*PI-phi[i])) < L && theta[i] > atan(-y[i] / (sin(2*PI-phi[i]) * sp))) {
			counterge++;
			z[i] = y[i] / (sin(2*PI-phi[i]) * tan(theta[i]));
			xfinal[i] = x[i] + y[i] / tan(2 * PI - phi[i]); 
			yfinal[i] = 0; 
			k++;
			cout << k << endl;
			i++;
		}
		//else if (y[i] < 0 && phi[i] > 3 * PI / 2 && (x[i] + y[i] / tan(2 * PI - phi[i])) > L && (x[i] + y[i] / tan(phi[i])) < L + D / (tan(2 * PI - phi[i])) && theta[i] > atan((x[i] - L) / (cos(2 * PI - phi[i]) * sp))) {
		//counterge++;
		//z[i] = (x[i] - L) / (cos(2 * PI - phi[i]) * tan(theta[i]));
		//xfinal[i] = L;		
		//yfinal[i] = (x[i] - L) * tan(2 * PI - phi[i]) + y[i];//correggere
		//k++;
		//cout << k << endl;
		//i++;
		//}
		else if (y[i] > D && phi[i] > 0 && phi[i]  < PI / 2 && (x[i] - (y[i] - D) / tan(phi[i])) > 0 && (x[i] - (y[i] - D) / tan(phi[i])) < L && theta[i] > atan((y[i] - D) / (sin(phi[i]) * sp))) {
			counterge++;
			z[i] = -(y[i] - D) / (sin(phi[i]) * tan(theta[i]));
			xfinal[i] = x[i] - (y[i] - D) / tan(phi[i]);
			yfinal[i] = D;
			k++;
			cout << k << endl;
			i++;
		}
		//else if (y[i] > D && phi[i] > 0 && phi[i]  < PI / 2 && (x[i] - (y[i] - D) / tan(phi[i])) > L && (x[i] - (y[i] - D) / tan(phi[i])) < L + D / (tan(phi[i])) && theta[i] > atan((x[i]-L) / (cos(phi[i]) * sp))) {
		//counterge++;
		//z[i] = (x[i] - L) / (cos(phi[i]) * tan(theta[i]));
		//xfinal[i] = L;		//correggere
		//yfinal[i] = y[i] - (x[i]-L)/tan(phi[i]);
		//k++;
		//cout << k << endl;
		//i++;
		//}
		else if (y[i] > D && phi[i] > PI / 2 && phi[i] < PI && (x[i] + (y[i] - D) / tan(PI - phi[i])) > 0 && (x[i] + (y[i] - D) / tan(PI - phi[i])) < L  && theta[i] > atan ( (y[i] - D) / (sin(PI-phi[i]) * sp))) {
			counterge++;
			z[i] = -(y[i] - D) / (sin(PI-phi[i]) * tan(theta[i])); 
			xfinal[i] = x[i] + (y[i] - D) / tan(PI - phi[i]);
			yfinal[i] = D;
			k++;
			cout << k << endl;
			i++;
		}
		//else if (y[i] > D && phi[i] > PI / 2 && phi[i] < PI && (x[i] + (y[i] - D) / tan(PI - phi[i])) > 0 - D / (tan(PI - phi[i])) && (x[i] + (y[i] - D) / tan(PI - phi[i])) < 0  && theta[i] > atan(-x[i] / (cos(PI-phi[i]) * sp))) {
		//counterge++;
		//z[i] = -x[i] / (cos(PI - phi[i]) * tan(theta[i]));
		//xfinal[i] = 0;
		//yfinal[i] = y[i]+x[i]/tan(PI-phi[i]);//correggere
		//k++;
		//cout << k << endl;
		//i++;
		//}
		else continue;		//calculations need to be checked

	}
	
	double counter = 0;		// counter for coincidence events
	double counterz = 0;
	double counterzcoin = 0;
	for (int i = 0; i < N_events; i++) {
		xfin[i] = x[i] - (H + z[i]) * tan(theta[i]) * cos(phi[i]);		// these should be correct
		yfin[i] = y[i] - (H + z[i]) * tan(theta[i]) * sin(phi[i]);
		if (0 <= xfin[i] && xfin[i] <= L && yfin[i] >= 0 && yfin[i] <= D) counter++;
		if (z[i] != 0) counterz++;		//counts the event fallen on the sp thick side
		if (0 <= xfin[i] && xfin[i] <= L && yfin[i] >= 0 && yfin[i] <= D && z[i] != 0) counterzcoin++;
	}
	
	cout << "Il numero di eventi e': " << N_events << endl << "Il numero di eventi in coincidenza e': " << counter << endl << "La percentuale in coincidenza e': " << counter / N_events * 100 << "%" << endl << "Il numero di eventi totali controllato e': " << countentries << endl << "Il numero totale buono di eventi e': " << counterge << endl << "Il numero totale di eventi passato lateralmente e': " << counterz << endl << "Il numero totale di eventi passati lateralmente e in coincidenza e' : " << counterzcoin << endl;		// this ends the program
	
	/*for (int i = 0; i < N_events; i++) {
		zfinal[i] = z[i];
		xfinal[i + N_events] = xfin[i];
		yfinal[i + N_events] = yfin[i];
		zfinal[i + N_events] = -H;
	}
		
	TCanvas *c1 = new TCanvas("c1", "", 200, 10, 600, 400);		//finally the plot
	c1->SetFillColor(10);
	c1->SetGrid();
	TGraph2D *grafico = new TGraph2D(n, xfinal, yfinal, zfinal);
	grafico->SetTitle("Slabs coincidence TER, H = 0.08, sp = 0.027, square = 10");
	grafico->SetMarkerColor(2);
	grafico->GetXaxis()->SetTitle("X-Axis");
	grafico->GetYaxis()->SetTitle("Y-Axis");
	grafico->GetZaxis()->SetTitle("Z-Axis");
	grafico->Draw("PCOL");*/
	return;
} 
	
