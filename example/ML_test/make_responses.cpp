#include<cstdio>
#include<fstream>
#include<cmath>

using namespace std;

int main(){

	double ML, BL;
	ifstream infile("results.out");
	infile>>ML;
	infile.close();
	infile.open("Baseline/results.out");
	infile>>BL;
	infile.close();
	ofstream outfile("responses.out");
	outfile<<abs(BL-ML)<<endl;
	outfile<<ML<<endl;
	outfile.close();

}