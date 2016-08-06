#define _USE_MATH_DEFINES

#include<iostream>
#include<vector>
#include<string>
#include<math.h>
#include<limits>
#include<fstream>
#include <time.h>

#include<boost/property_tree/ptree.hpp>
#include<boost/property_tree/ini_parser.hpp>

#include<omp.h>

using namespace std;

class granule {
public:
//--------------- These parameters are given in .inp1 files--------------------------------
	double a; //Granular radius
//--------------- End parameters given by .inp1 files-------------------------------------

//----------------These parameters are givin in .inp2 files -----------------------------
	double pulsedur; //pulse duration in seconds
	double alpha; // Thermal Diffusivity of medium in cm**2/sec
	double cond; // Thermal conductivity of medium in j/cm-Csec
//-------------- End Parameters from .inp2 files ----------------------------------------
	double A0; //used in (orig,in,surf,out,vanilla)temp functions, calculated in bigthrem from other input parameters
	double OutTemp(double r, double t);
	double SurfTemp(double t);
	double InTemp(double r, double t);
	double OrigTemp(double t);
	double Temp(double r, double t);
	double MelTemp(double r, double t);
	double tmin, PulseSep;
	int PulseNum;
	vector<double> PulseVec;
};

class BigTherm {
	public:
		BigTherm(char* configfilename);
		vector<int> elapsed, tarray, guess;
		int k, s, improf, n,TempFlag, SeedFlag;
		int xnum, ynum, znum, tnum, melnum, rnum;
		vector<int> test;
		double BoxLength, BoxDepth, MelDens;
		double xmin, xmax, xsize;
		double ymin, ymax, ysize;
		double zmin, zmax, zsize;
		double tmin, tmax, tsize;
		vector< vector< double > > Pos;
		vector< double >  weight, PulseVec;
		vector< vector< double > > VecTemp; 
		vector< double > rsize;
		double  dx;
		double ran1, ran2, ran3;
		double irr,Tbody, absorb;
		double Temp, x, y, z, t, r;
		double sigma, corspot, Focus;
		double spotsize,MelTemp,trans,dobs;
		double Qtot, Qlost, Qkept;
		granule granny;
		boost::property_tree::ptree configfile;
		string MelPlacementFile;
		string TempOutputFile;
		string DamageOutputFile;
		double BigThermTempAt(double x, double y, double z, double t);
		// this function is to build the function again if you change any variables that affect weight or vectemp
		void ReConstruct(void);
		void BigTherm2(void);
		void BigTherm2(double outputxmin, double outputxmax,int outxnum, double outputymin, double outputymax, int outynum, double outputzmin, double outputzmax, int outznum);
		int PulseNum;
		double PulseSep;
		void TakataDamage(void);
		int TakataDamageAt(double x, double y, double z);
		int DamageInfo;
		int numthreads;
};
void takatadamage(char* configfilename);

