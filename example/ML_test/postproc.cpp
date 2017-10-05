#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>
#include<stdlib.h>

using namespace std;

int main(){

	int stat;
	string line;
	stat = system("cp nozzle.dat solution_flow.dat");
	stat = system("SU2_SOL config_CFD.cfg");
	double sr, tw, p1, p2, p3, beta;
	ifstream i1, i2, i3, i4;
	i1.open("verification.dat");
	i2.open("sample_features.dat");
	i3.open("beta_for_su2.dat");
	i4.open("flow.dat");
	unsigned long index;
	ofstream of("final_flow.dat");
	getline(i4, line);
	of<<line<<endl;
	getline(i4, line);
	of<<line<<"\"strain_rate\"\"TauWallNearest\"\"WallDist\"\"Beta\""<<endl;
	getline(i4, line);
	of<<line<<endl;
	while(i3>>index>>beta){

		i1>>index>>sr>>tw;
		i2>>index>>p1>>p2>>p3;
		getline(i4, line);
		of<<line<<scientific<<setprecision(10)<<'\t'<<sr<<'\t'<<tw<<'\t'<<p3<<'\t'<<beta<<'\t'<<endl;		

	}
	while(getline(i4, line)){

		of<<line<<endl;

	}

}