
int T = 100;		// defining some parameters: TIME
double L = 1.83;		//LENGTH
double D = 0.20;		//WIDTH
double H = 0.08;		//DISTANCE
double sp = 0.027;		//THICKNESS
int FLUSSO = 130;		//FLUX
double square = 1;	//SIZE OF THE SQUARE AROUND THE SLAB
const int N_events = 109538;
//int ((2*square+L) * (2*square+D) * T * FLUSSO);		//NUMBER OF EVENTS
const double PI = atan(1) * 4;
int M = 4/PI;		// this is the constant in front of the distribution; see later

void SimRC_eventisuareatotale(int seed)
{
	double x[N_events], y[N_events], z[N_events], zfin[N_events], x1final[N_events], y1final[N_events], xfin[N_events], yfin[N_events], phi[N_events], theta[N_events], x1[N_events], y1[N_events], z1[N_events], phi1[N_events], theta1[N_events]; // all the necessary array
	double uniform, u;
	TRandom3 *rand = new TRandom();
	rand->SetSeed(seed); // to set the seed
	int counterge = 0;		//to count good events
	int k = 0;		//dynamic counter 
	for (int i = 0; i < N_events; i++) {
		//rand->SetSeed(0);		//forget about this, it's wrong; I kept it only to remember not to use it
		x[i] = rand->Uniform(-square,L+square);		//generates N_events values along X
		y[i] = rand->Uniform(-square, D+square);
		phi[i] = rand->Uniform(0, 2*PI);
		uniform = 0;
		u = 2;
		while (M*u > M*cos(uniform)*cos(uniform)*sin(uniform)) {		//rejection condition
			u = rand->Uniform(0, 1);
			uniform = rand->Uniform(0, PI / 2);
		};
		theta[i] = uniform;
		if (x[i] >= 0 && x[i] <= L && y[i] >= 0 && y[i] <= D) {
			z[i] = 0;
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = x[i];
			y1final[counterge] = y[i];
			counterge++; 
			k++;
			//cout << k << endl;
		}
		else if (x[i] < 0 && phi[i] >PI / 2 && phi[i] < PI && (y[i] - x[i] * tan(phi[i])) > 0 && (y[i] - x[i] * tan(phi[i])) < D && theta[i] > atan(x[i] / (cos(phi[i]) * sp))) {
			z[i] = -x[i] / (cos(phi[i]) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = 0;
			y1final[counterge] = (y[i] - x[i] * tan(phi[i]));
			counterge++;
			k++;
			//cout << k << endl;			
		}
		else if (x[i] < 0 && phi[i] >PI / 2 && phi[i] < PI && (y[i] + x[i] * tan(PI-phi[i])) > D && x[i]-(D-y[i])/ tan(PI-phi[i]) > 0 && x[i] - (D - y[i]) / tan(PI-phi[i]) <L && theta[i] > atan((y[i]-D) / (sin(PI-phi[i]) * sp))) {
			z[i] = -(y[i] - D) / (sin(PI - phi[i]) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			y1final[counterge] = D;
			x1final[counterge] = x[i] - (D - y[i]) / tan(PI - phi[i]);
			counterge++;
			k++;
			//cout << k << endl;
		}
		else if (x[i] < 0 && phi[i] >PI && phi[i] < 3*PI/2 && (y[i] - x[i] * tan(phi[i])) > 0 && (y[i] - x[i] * tan(phi[i])) < D && theta[i] > atan(-x[i] / (cos(PI-phi[i]) * sp))) {
			z[i] = -x[i] / (cos(phi[i]) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = 0;
			y1final[counterge] = (y[i] - x[i] * tan(phi[i]));
			counterge++;
			k++;
			//cout << k << endl;
		}
		else if (x[i] < 0 && phi[i] >PI && phi[i] < 3 * PI / 2 && x[i] - y[i] / tan(phi[i]) > 0 && x[i] - y[i] / tan(phi[i]) < L&& (y[i] - x[i] * tan(phi[i])) < 0 && theta[i] > atan(-y[i] / (sin(phi[i]-PI) * sp))) {
			z[i] = y[i] / (sin(phi[i] - PI) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			y1final[counterge] = 0;
			x1final[counterge] = x[i] - y[i] / tan(phi[i]);
			counterge++;
			k++;
			//cout << k << endl;
		}
		else if (x[i] > L && phi[i] > 0 && phi[i] < PI / 2 && (y[i] - (x[i] - L) * tan(phi[i])) > 0 && (y[i] - (x[i] - L) * tan(phi[i])) < D && theta[i] > atan((x[i] - L) / (cos(phi[i]) * sp))) {
			z[i] = -(x[i] - L) / (cos(phi[i]) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = L;
			y1final[counterge] = y[i] - (x[i] - L) * tan(phi[i]);
			counterge++;
			k++;
			//cout << k << endl;
		}
		else if (x[i] > L && phi[i] > 0 && phi[i] < PI / 2 && (y[i] - (x[i] - L) * tan(phi[i])) > D && x[i] - (y[i] - D) / tan(phi[i]) > 0 && x[i] - (y[i] - D) / tan(phi[i]) < L && theta[i] > atan((y[i]-D) / (sin(phi[i]) * sp))) {
			z[i] = -(y[i] - D) / (sin(phi[i]) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			y1final[counterge] = D;
			x1final[counterge] = x[i] - (y[i] - D) / tan(phi[i]);
			counterge++;
			k++;
			//cout << k << endl;
		}
		else if (x[i] > L && phi[i] > 3 * PI / 2 && (y[i] + (x[i]-L) / tan(phi[i] - 3 * PI / 2)) > 0  && (y[i] + (x[i]-L) / tan(phi[i] - 3 * PI / 2)) < D && theta[i] > atan((x[i]-L) / (sin(phi[i] - 3 * PI / 2) * sp))) {
			z[i]  = -(x[i] - L) / (sin(phi[i] - 3 * PI / 2) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = L;
			y1final[counterge] = y[i] + (x[i] - L) / tan(phi[i] - 3 * PI / 2);
			counterge++;
			k++;
			//cout << k << endl;
		}
		else if (x[i] > L && phi[i] > 3 * PI / 2 && x[i] + y[i] * tan(phi[i] - 3 * PI / 2) > 0 && x[i] + y[i] * tan(phi[i] - 3 * PI / 2) < L && (y[i] + (x[i] - L) / tan(phi[i] - 3 * PI / 2)) < 0 && theta[i] > atan(-y[i] / (sin(2*PI-phi[i]) * sp))) {
			z[i] = y[i] / (sin(2 * PI - phi[i]) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			y1final[counterge] = 0;
			x1final[counterge] = x[i] + y[i] * tan(phi[i] - 3 * PI / 2);
			counterge++;
			k++;
			//cout << k << endl;
		}
		else if (y[i] < 0 && phi[i] > PI && phi[i] < 3 * PI/2 && (x[i] - y[i] / tan(phi[i]-PI)) > 0 && (x[i] - y[i] / tan(phi[i]-PI)) < L && theta[i] > atan(-y[i] / (sin(phi[i]-PI) * sp))) {
			z[i] = y[i] / (sin(phi[i]-PI) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = x [i] - y[i] / tan(phi[i] - PI);
			y1final[counterge] = 0;
			counterge++;
			k++;
			//cout << k << endl;
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
			z[i] = y[i] / (sin(2*PI-phi[i]) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = x[i] + y[i] / tan(2 * PI - phi[i]);
			y1final[counterge] = 0;
			counterge++;
			k++;
			//cout << k << endl;
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
			z[i] = -(y[i] - D) / (sin(phi[i]) * tan(theta[i]));
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = x[i] - (y[i] - D) / tan(phi[i]);
			y1final[counterge] = D;
			counterge++;
			k++;
			//cout << k << endl;
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
			z[i] = -(y[i] - D) / (sin(PI-phi[i]) * tan(theta[i])); 
			x1[counterge] = x[i];
			y1[counterge] = y[i];
			z1[counterge] = z[i];
			theta1[counterge] = theta[i];
			phi1[counterge] = phi[i];
			x1final[counterge] = x[i] + (y[i] - D) / tan(PI - phi[i]);
			y1final[counterge] = D;
			counterge++;
			k++;
			//cout << k << endl;
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
	for (int i = 0; i < counterge; i++) {
		xfin[i] = x1[i] - (H + z1[i]) * tan(theta1[i]) * cos(phi1[i]);		// these should be correct
		yfin[i] = y1[i] - (H + z1[i]) * tan(theta1[i]) * sin(phi1[i]);
		if (0 <= xfin[i] && xfin[i] <= L && yfin[i] >= 0 && yfin[i] <= D) counter++;
		if (z1[i] != 0) counterz++;		//counts the event fallen on the sp thick side
		if (0 <= xfin[i] && xfin[i] <= L && yfin[i] >= 0 && yfin[i] <= D && z[i] != 0) counterzcoin++;
	}
	
	cout << "Il numero di eventi e': " << counterge << endl << "Il numero di eventi in coincidenza e': " << counter << endl << "La percentuale in coincidenza e': " << counter/counterge * 100  << "%" << endl << "Il numero di eventi totali controllato e': " << N_events << endl << "Il numero totale buono di eventi e': " << counterge << endl << "Il numero totale di eventi passato lateralmente e': " << counterz << endl << "Il numero totale di eventi passati lateralmente e in coincidenza e' : " << counterzcoin << endl;		// this ends the program
	
	/*const int n1 = 2*counterge;		//I write the files all in one file to plot in the end
	double xfinal1[n1];
	double yfinal1[n1];
	double zfinal1[n1];
	for (int i = 0; i < counterge; i++) {
		xfinal1[i] = x1final[i];
		yfinal1[i] = y1final[i];
		zfinal1[i] = z1[i];
		xfinal1[i + counterge] = xfin[i];
		yfinal1[i + counterge] = yfin[i];
		zfinal1[i + counterge] = -H;
	}
		
	TCanvas *c1 = new TCanvas("c1", "", 200, 10, 600, 400);		//finally the plot
	c1->SetFillColor(10);
	c1->SetGrid();
	TGraph2D *grafico = new TGraph2D(n1, xfinal1, yfinal1, zfinal1);
	grafico->SetTitle("Slabs coincidence 3D, H = 0.08 m, square = 10 m");
	grafico->SetMarkerColor(2);
	grafico->GetXaxis()->SetTitle("X-Axis [m]");
	grafico->GetYaxis()->SetTitle("Y-Axis [m]");
	grafico->GetZaxis()->SetTitle("Z-Axis [m]");
	grafico->Draw("PCOL");*/
	return;
} 
	
