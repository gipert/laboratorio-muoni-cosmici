#include<iostream>
#include<math.h>

#include<TRandom.h>
#include<TCanvas.h>
#include<TH1F.h>
#include<TApplication.h>

using namespace std;

TRandom randRoot;
double cos2( double low , double up );
void prog( double height );

int main() {
	
	//for ( int i = 1 ; i <= 10 ; i++ ) {
		
	//	prog(i*0.10);
	//}
	
	prog(0.08);
	return 0;
}

void prog( double heigth ) {
	
	double xMin = 0;
	double xMax = 1.83;
	double yMin = 0;
	double yMax = 0.20;
	double h = heigth;

	double x, x_;
	double y, y_;
	double theta;
	double phi;

	randRoot.SetSeed(2);

	//TH1F hist( "cos2" , "cos2" , 100 , 0 , M_PI/2 );

	int N = 1000000;
	int n = 0;	
	for ( int i = 0 ; i < N ; i++ ) {
		
		x     = randRoot.Uniform(xMin,xMax);
		y     = randRoot.Uniform(yMin,yMax);
		phi   = randRoot.Uniform(0,2*M_PI);
		theta = cos2(0,M_PI/2);
		//hist.Fill(theta);
		
		x_ = -h*tan(theta)*cos(phi) + x;
		y_ = -h*tan(theta)*sin(phi) + y;

		if ( x_ > xMin && x_ < xMax && y_ > yMin && y_ < yMax ) n++;
	}

	//TApplication app( "app" , &argc , argv );
	//TCanvas can("cos2","cos2",1);
	//can.cd();
	//hist.Draw();
	//can.SaveAs("cos2.png");
	//app.Run(kFALSE);
	
	//cout << "Eventi generati: " << N << endl;
	//cout << "Eventi in coincidenza: " << n << endl;
	cout << flush << "Distanza slabs: " << heigth << " m, percentuale di coincidenze: " << (n*1./N)*100 << "%" << endl;

	return;
}

double cos2( double low , double up ) {
	
	double u1, u2, x;

	do { 	u1 = randRoot.Uniform(low,up);
		u2 = randRoot.Uniform(0,1);

		x = cos(u1)*cos(u1)*sin(u1);
	}
	while ( x < u2 );

	return u1;
}
