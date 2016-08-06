#include"TGlib.h"

double Gmma(double aa, double bb, double cc){
if(bb == 0.0) {
//-------------gigamma function ---------------------------
if(cc > 30.0)
	cc = 30.0;
//-------------gammaser function -------------------
//-------------modified(un-normalized) version of GSER from numerical recipies --------------
//modification of GSER is ignoring(removing)the "gln" variable, and the function that calculates it gammln(a). this is used in the sum evaluating the exponential for del < sum *EPS, as subraction by the natural log of the complete gamma function inside the exponential acts the same as dividing by the same value.
	const int ITMAX = 100;
	const double EPS = numeric_limits<double>::epsilon();
	double sum, del,ap;
	if(cc <= 0.0){
		if(cc<0.0){
			return -1;
			//fail in a mysterious way
		}	
		else
			return 0.0;
	}
	else{
		ap = aa;
		del = sum = 1.0/aa;
		for(int n = 0; n <ITMAX; n++){
			ap += 1.0;
			del *= cc/ap;
			sum += del;
			if(fabs(del) < fabs(sum)*EPS){
				return sum*exp(-cc+aa*log(cc));
			}
		}
	}
}
}

double LagL(double a, double b, double c){
double d, gm1, gm2, Gamma, answer;
d = -c;
gm1 = -0.5;
gm2 = 0.0;
Gamma = Gmma(gm1,gm2,d);
answer = - (Gamma * sqrt(d) * (1.0 - 2*c/3) + (2*exp(c))/3)/M_PI;
return answer;
}
//calculates the temperature outside the granual
double granule::OutTemp(double r, double t){

double Outterm[4], gm[2], lg[2];

double OutTheta;

gm[0] = 0.5;

gm[1] = 0.0;

lg[0] = 1.5;

lg[1] = -0.5;

double gamma = Gmma(gm[0],gm[1],(((-a + r)*(-a + r))/(4.0*alpha*t)));
	Outterm[0] = ((a*t*((sqrt(M_PI)*((-a + r)*(-a + r) + 2.0*alpha*t))/(2.0*alpha*t) - (2.0*(-a + r)*(0.5*exp((-((-a + r)*(-a +r))/(4.0*alpha*t))) + (sqrt(((-a + r)*(-a +r))/(alpha*t))*(((-a + r)*(-a +r)) + 2.0*alpha*t)* gamma)/(4.0*((-a + r)*(-a+r)))))/(sqrt(alpha)*sqrt(t))))/(sqrt(M_PI)));

	Outterm[1] = (a*a*a - 3.0*a*a*r + 3.0*a*r*r - (r*r*r) + 6.0*a*alpha*t - 6.0*alpha*r*t + 6.0*(alpha*sqrt(alpha))*sqrt(M_PI)*(t*sqrt(t))*LagL(lg[0],lg[1],-((-a + r)*(-a + r))/(4.0*alpha*t)))/(6.0*alpha);

	Outterm[2] = (-a*a*a - 3.0*a*a*r - 3.0*a*r*r - r*r*r - 6.0*a*alpha*t - 6.0*alpha*r*t + 6.0*(alpha*sqrt(alpha))*sqrt(M_PI)*(t*sqrt(t))*LagL(lg[0],lg[1],(-(a + r)*(a +r))/(4.0*alpha*t)))/(6.0*alpha);

	Outterm[3] = (a*t*((sqrt(M_PI) * ((a + r)*(a + r) + 2.0*alpha*t))/(2.0*alpha*t) - (2.0*(a + r)* (0.5*exp(-((a + r)*(a + r))/(4.0*alpha*t)) + (sqrt(((a + r)*(a + r))/(alpha*t))*((a + r)*(a + r) + 2.0*alpha*t)*Gmma(gm[0], gm[1], ((a + r)*(a + r))/(4.0*alpha*t)))/(4.0*(a + r)*(a + r))))/(sqrt(alpha)*sqrt(t))))/(sqrt(M_PI));

OutTheta = alpha*A0*(Outterm[0] - Outterm[1] + Outterm[2] + Outterm[3])/(2.0*cond);

return OutTheta/r;
}
// calculates temperature on the surface of the granual
double granule::SurfTemp(double t){

double Surffnc[4];

double SurfTheta, gm[2], lg[2];

gm[0] = 0.5;

gm[1] = 0.0;

lg[0] = 1.5;

lg[1] = -0.5;

Surffnc[0] = t;

Surffnc[1] = t*(sqrt(M_PI)*(2.0*a*a + alpha*t)/(alpha*t) - 4.0*a*(0.5*exp(-(a*a)/(alpha*t)) + sqrt((a*a)/(alpha*t))*(2.0*a*a + alpha*t)*Gmma(gm[0],gm[1],(a*a)/(alpha*t))/(4.0*a*a))/(sqrt(alpha)*sqrt(t)))/(sqrt(M_PI));

Surffnc[2] = 4.0*(t*sqrt(t))/(3.0*sqrt(M_PI));

Surffnc[3] = -4.0*(a*a*a)/(3.0*alpha*sqrt(alpha)) - 2.0*a*t/sqrt(alpha) + sqrt(M_PI)*t*sqrt(t)*LagL(lg[0],lg[1],-((a*a)/(alpha*t)));

Surffnc[4] = a*alpha*A0*(Surffnc[0] + Surffnc[1])/(2.0*cond) + alpha*sqrt(alpha)*A0*(Surffnc[3] - Surffnc[2])/(2.0*cond);

return Surffnc[4]/a;
}
// Calculates temperature inside the granual
double granule::InTemp(double r, double t){

double Interm[5], gm[2], lg[2];

double Intheta;

gm[0] = 0.5;

gm[1] = 0.0;

lg[0] = 1.5;

lg[1] = -0.5;

Interm[0] = r*t;

Interm[1] = (-a*a*a + 3.0*a*a*r - 3.0*a*r*r + r*r*r - 6.0*a*alpha*t + 6.0*alpha*r*t + 6.0*alpha*sqrt(alpha)*sqrt(M_PI)*t*sqrt(t)*LagL(lg[0],lg[1],-((a - r)*(a - r))/(4.0*alpha*t)))/(12.0*alpha);

Interm[2] = (-a*a*a - 3.0*a*a*r - 3.0*a*r*r - r*r*r - 6.0*a*alpha*t - 6.0*alpha*r*t + 6.0*alpha*sqrt(alpha)*sqrt(M_PI)*t*sqrt(t)*LagL(lg[0],lg[1],-((a + r)*(a + r))/(4.0*alpha*t)))/(12.0*alpha);

Interm[3] = a*t*(sqrt(M_PI)*(((a - r)*(a - r)) + 2.0*alpha*t)/(2.0*alpha*t) - 2.0*(a - r)*(0.5*exp(-((a - r)*(a - r))/(4.0*alpha*t)) + (a - r)/sqrt(alpha*t)*((a - r)*(a - r) + 2.0*alpha*t)*Gmma(gm[0],gm[1],((a - r)*(a - r))/(4.0*alpha*t))/(4.0*((a - r)*(a - r))))/sqrt(alpha*t))/(2.0*sqrt(M_PI));

Interm[4] = a*t*(sqrt(M_PI)*(((a + r)*(a + r)) + 2.0*alpha*t)/(2.0*alpha*t) - 2.0*(a + r)*(0.5*exp(-((a + r)*(a + r))/(4.0*alpha*t)) + (a + r)/sqrt(alpha*t)*(((a + r)*(a + r)) + 2.0*alpha*t)*Gmma(gm[0],gm[1],((a + r)*(a + r))/(4.0*alpha*t))/(4.0*((a + r)*(a + r))))/sqrt(alpha*t))/(2.0*sqrt(M_PI));

Intheta = alpha*A0*(Interm[0] - Interm[1] + Interm[2] - Interm[3] + Interm[4])/cond;

return Intheta/r;
}

double granule::OrigTemp(double t){
double Origfnc[3], gm[2], lg[2];

gm[0] = 0.5;

gm[1] = 0.0;

lg[0] = 1.5;

lg[1] = -0.5;

Origfnc[0] = t;

Origfnc[1] = t*(sqrt(M_PI)*(a*a + 2.0*alpha*t)/(2.0*alpha*t) - 2.0*a*(0.5*exp(-a*a/(4.0*alpha*t)) + sqrt(a*a/(alpha*t))*(a*a + 2.0*alpha*t)*Gmma(gm[0],gm[1],a*a/(4.0*alpha*t))/(4.0*a*a))/(sqrt(alpha)*sqrt(t)))/sqrt(M_PI);

Origfnc[2] = sqrt(t)*(-(a*sqrt(M_PI)/(sqrt(alpha)*sqrt(t))) - sqrt(a*a/(alpha*t))*Gmma(-gm[0],gm[1],a*a/(4.0*alpha*t))/2.0)/sqrt(M_PI);

return alpha*A0*(Origfnc[0] - Origfnc[1] - (a*Origfnc[2])/sqrt(alpha))/cond;
}
//
double granule::Temp(double r, double t){

if(t == 0.0)
	return 0.0;
else{
	if(r == 0.0)
		return OrigTemp(t);
	else{
		if((r>0.0)&&(r<a))
			return InTemp(r,t);
		else{
			if(r == a)
				return SurfTemp(t);
			else{
				if(r>a)
					return OutTemp(r,t);
			}
		}
	}
}

}

double granule::MelTemp(double r, double t){
	double times;
	double ReturnTemp;
	for(int i =0; i<PulseVec.size(); i++){
		times = t - PulseVec[i];
//This part here is basically the original MelTemp Function from FORTRAN, here it is looped over the number of pulses that are given
		if( times < 0.0){
			ReturnTemp += 0.0;
		}
		else{
			if(times > pulsedur)
				ReturnTemp += Temp(r,times) - Temp(r,times-pulsedur);
			else
				ReturnTemp += Temp(r,times);
		}
	}
	return ReturnTemp;
/*
if( t < 0.0)
	return 0.0;
if(t > pulsedur)
	return Temp(r,t) - Temp(r,t-pulsedur);
else
	return Temp(r,t);
*/
}

BigTherm::BigTherm(char* configfilename){
boost::property_tree::ptree configfile;
boost::property_tree::ini_parser::read_ini(configfilename, configfile);
BoxLength=configfile.get<double>("Bigtherm.BoxLength");
BoxDepth=configfile.get<double>("Bigtherm.BoxDepth");
granny.a=configfile.get<double>("Bigtherm.GranularRadius");
xnum=configfile.get<int>("Bigtherm.xnum");
xmin=configfile.get<double>("Bigtherm.xmin");
xmax=configfile.get<double>("Bigtherm.xmax");
ynum=configfile.get<int>("Bigtherm.ynum");
ymin=configfile.get<double>("Bigtherm.ymin");
ymax=configfile.get<double>("Bigtherm.ymax");
znum=configfile.get<int>("Bigtherm.znum");
zmin=configfile.get<double>("Bigtherm.zmin");
zmax=configfile.get<double>("Bigtherm.zmax");
tnum=configfile.get<int>("Bigtherm.tnum");
tmin=configfile.get<double>("Bigtherm.tmin");
tmax=configfile.get<double>("Bigtherm.tmax");
MelDens=configfile.get<double>("Bigtherm.MelDens");
rnum=configfile.get<int>("Bigtherm.rnum");
absorb=configfile.get<double>("Bigtherm.absorb");
irr=configfile.get<double>("Bigtherm.irr");
granny.pulsedur=configfile.get<double>("Bigtherm.pulsedur");
granny.alpha=configfile.get<double>("Bigtherm.alpha");
granny.cond=configfile.get<double>("Bigtherm.cond");
corspot=configfile.get<double>("Bigtherm.corspot");
spotsize=configfile.get<double>("Bigtherm.spotsize");
trans=configfile.get<double>("Bigtherm.trans");
improf=configfile.get<int>("Bigtherm.improf");
SeedFlag=configfile.get<int>("Bigtherm.SeedFlag");
PulseNum=configfile.get<int>("Bigtherm.PulseNum");
PulseSep=configfile.get<double>("Bigtherm.PulseSep");
PulseVec.push_back(tmin);
for(int i =1; i<PulseNum; i++){
	PulseVec.push_back(PulseVec[i-1] + granny.pulsedur + PulseSep);
}
granny.PulseVec = PulseVec;
//---------------------------------------------------------------
//This bit will need some options, if it should be here at all
TempOutputFile=configfile.get<string>("Bigtherm.TempOutputFile");
DamageOutputFile=configfile.get<string>("Bigtherm.DamageOutputFile");
DamageInfo=configfile.get<int>("Bigtherm.DamageInfo");
numthreads=configfile.get<int>("BigTherm.numthreads",1);
if(numthreads<1)
	numthreads = 1;

//////////////////////////////////////////////////////////////////
//end loading
//---------------------------------------------------------------
//////////////////////////////////////////////////
if((spotsize>xmax) || (spotsize>ymax))
//	cout<<"box size too small"<<endl;
xsize = (xmax - xmin)/xnum;

      ysize = (ymax - ymin)/ynum;

      zsize = (zmax - zmin)/znum;

      tsize = (tmax - tmin)/tnum;

melnum = 0;
melnum = configfile.get<int>("Bigtherm.Melnum", 0);
if(melnum == 0){
	MelDens=configfile.get<double>("Bigtherm.MelDens");
	melnum = MelDens * BoxLength*BoxLength * BoxDepth;
}
else
	MelDens = 0;
//      melnum = MelDens * BoxLength*BoxLength * BoxDepth;
//cout<<xsize<<" "<<ysize<<" "<<zsize<<endl;
//	cout<<xsize<<" "<<ysize<<" "<<zsize<<" "<<tsize<<" "<<melnum<<endl;
//look into adding option for inputting this number directly instead of using this ass backwards method
/*
!***********************************************************************
!
!     Determine the image profile on the retina.  
!     Use image profile technique determined by
!     the value of improf:
!
!                  improf=0  <->  Gaussian Image
!                  improf=3  <->  Top Hat Source
!                  improf=4  <->  Annular Beam 
!
!     The image diameter is then determined by the value of
!     spotsize.  The amount of focusing
!     that has taken place is a function of the retinal image diameter
!     as well. For annular beams, which assume a 37 per cent central
!     obscuration, the diameter of the central obscuration, dobs, is
!     spotsize*DSQRT(0.37).
!
!***********************************************************************
*/
	//cout<<melnum<<endl;
	dobs = spotsize * sqrt(0.37);
	//cout<<dobs<<endl;

	if(improf ==  0) 
		Focus = 2 * (corspot*corspot) / (spotsize*spotsize);

        if (improf == 3)
            Focus = (corspot*corspot) / (spotsize*spotsize);

	if(improf == 4)
            Focus = (corspot*corspot) / (spotsize*spotsize);

	//cout<<Focus<<endl;
         granny.A0 = 3 * irr * trans * Focus * (1 - (1 - exp(-2 * granny.a * absorb)*   (1 + 2 * granny.a * absorb))/(2 * granny.a*granny.a * absorb*absorb) )/ (4 * granny.a * granny.pulsedur);

	//cout<<granny.A0<<endl;

         sigma = 8.0 / (spotsize*spotsize);
	//cout<<sigma<<endl;
      if (tmax < granny.pulsedur)
          Qtot = granny.A0 * 4.0 * M_PI * granny.a*granny.a*granny.a * tmax / 3.0;
      else
          Qtot = granny.A0 * 4.0 * M_PI * granny.a*granny.a*granny.a * granny.pulsedur / 3.0;
	
	//cout<<Qtot<<endl;
//----------------------------------
//need to have pos vector filled before here
vector<double> dummy;
// computing weighting of the fluence for the placed granule
if(SeedFlag == 0){
	MelPlacementFile=configfile.get<std::string>("Bigtherm.MelPlacementFile");
	fstream posfile(MelPlacementFile.c_str());
	if(posfile.good()){
		double in;
		string throwaway;
		while(!posfile.eof()){
			posfile>>in;
			dummy.push_back(in);
			posfile>>in;
			dummy.push_back(in);
				posfile>>in;
			dummy.push_back(in);
			Pos.push_back(dummy);
			dummy.clear();
			getline(posfile,throwaway);
		}
	}
	else{
		cout<<"Placement File not found, using random seed"<<endl;
		SeedFlag = 1;

	}
}
double thermcon;
if(SeedFlag == 1){

	srand(time(NULL));
	for( int j = 0; j<melnum; j++) {
		if(j== 0)
			Pos.clear();
		ran1 = (rand() % 101)/100.0;
		ran2 = (rand() % 101)/100.0; 
		ran3 = (rand() % 101)/100.0; 
		thermcon = (ran1 * BoxLength) - (BoxLength/2.0);
		dummy.push_back(thermcon);
		thermcon = (ran2 * BoxLength) - (BoxLength/2.0);
		dummy.push_back(thermcon);
		thermcon = (ran3 * BoxDepth) - (BoxDepth/2.0);
		dummy.push_back(thermcon);
		Pos.push_back(dummy);
		dummy.clear();
	}
}

for( int j = 0; j<melnum; j++) {

	r = sqrt(Pos[j][0]*Pos[j][0] + Pos[j][1]*Pos[j][1]);

	if (improf == 0)
               weight.push_back(exp(-sigma * r*r));
	else if (improf == 3){
		if (r <=  spotsize/2.0){
			weight.push_back(1.0);
		}
		else{
			weight.push_back(0.0);
		}
	}
	else if (improf == 4){
		if( r > spotsize/2.0)
			weight.push_back(0.0);
		else { 
			if(r < dobs/2.0)
				weight.push_back(0.0);
			else
				weight.push_back(1.0);
		}

	}
}
for(int i = 0; i<=tnum; i++){

	t = tmin + i * tsize;
	rsize.push_back(0.0);

        if (t <= 0.00001) 
            rsize[i] = (0.0005/rnum);
         else
            rsize[i] = (0.0005*(7.0 + log10(t))/rnum);
	
	for(int j = 0; j<=rnum; j++){

		r = j * rsize[i];

       		dummy.push_back(granny.MelTemp(r,t)); 

       		if(dummy[j] + Tbody > 100)
			TempFlag = 1; 
	}
	VecTemp.push_back(dummy);
	dummy.clear();


}
// ------------------------  check conservaction of energy ----------------------------------------
if (tmax <= 10.0E-6)
	dx = 0.0005/rnum;
else
	dx = 0.0005*(7.0 + log10(tmax))/rnum;

         Qkept = 0.0;

for(int j =0; j<rnum; j++){

	Qkept = Qkept + VecTemp[tnum][j] * 4.1868 * 4.0 * M_PI * dx*dx*dx * (3.0*(j+1)*(j+1) - 3.0*(j+1) + 1.0) / 3.0;
}
Qlost = Qtot - Qkept;

	if (TempFlag == 1) {
		cout<<endl;
		cout<<"////////////////////////////////////////////////////////"<<endl;
		cout<<"--------------------------------------------------------"<<endl;
		cout<<"      Warning . . . 100 degrees exceeded."<<endl;
		cout<<"--------------------------------------------------------"<<endl;
		cout<<"////////////////////////////////////////////////////////"<<endl;
	}

}

double BigTherm::BigThermTempAt(double x, double y, double z, double t){
//WARNING this function will NOT give outputs for arbitrary t, simply the closest t to the tsize given away from the min t.
if((t<tmin)||(t>tmax))
	return -1;
//fail in a misterious way
int l = (t-tmin)/tsize;
double dist, ilook, Temper, SumTemp;
SumTemp=0.0;

//!-----------------Summing effects of all granules---------------
omp_set_num_threads(numthreads);
#pragma omp parallel for reduction(+:SumTemp) private(dist, ilook, Temper)
for(int m = 0; m<melnum; m++){

	dist = sqrt((Pos[m][0]-x)*(Pos[m][0]-x) + (Pos[m][1]-y)*(Pos[m][1]-y) + (Pos[m][2]-z)*(Pos[m][2]-z));
//!-----------------Individual Temps Interpolated From VecTemp Table------
	ilook = dist/rsize[l];
	if (ilook > rnum-1) {
		Temper=0.0;
	}
	else{
		Temper = VecTemp[l][ilook] + (dist - ilook * rsize[l]) * (VecTemp[l][ilook+1] - VecTemp[l][ilook]) / rsize[l];
//							cout<<"stillwerkds"<<endl;
	}
	
	SumTemp += weight[m] * Temper;

} 
//Slight issues with noise, method only accurate to about a hundreth of a degree anyway so throwing out temperature rises on this magnitude cleans up the profiles a bit.
if(SumTemp < 1e-6)
	return 0.0;

return SumTemp;

}

//This function is a copy of the constructor if you want to run it again with the different beam perameters

void BigTherm::ReConstruct(void){
	double thermcon;
	vector<double> dummy;
	Pos.clear();
	srand(time(NULL));
	for( int j = 0; j<melnum; j++) {
		if(j== 0)
			Pos.clear();
		ran1 = (rand() % 101)/100.0;
		ran2 = (rand() % 101)/100.0; 
		ran3 = (rand() % 101)/100.0; 
		thermcon = (ran1 * BoxLength) - (BoxLength/2.0);
		dummy.push_back(thermcon);
		thermcon = (ran2 * BoxLength) - (BoxLength/2.0);
		dummy.push_back(thermcon);
		thermcon = (ran3 * BoxDepth) - (BoxDepth/2.0);
		dummy.push_back(thermcon);
		Pos.push_back(dummy);
		dummy.clear();
	}

PulseVec.clear();

PulseVec.push_back(tmin);

for(int i =1; i<PulseNum; i++){
	PulseVec.push_back(PulseVec[i-1] + granny.pulsedur + PulseSep);
}
granny.PulseVec = PulseVec;

weight.clear();
for( int j = 0; j<melnum; j++) {
	r = sqrt(Pos[j][0]*Pos[j][0] + Pos[j][1]*Pos[j][1]);
	if (improf == 0)
               weight.push_back(exp(-sigma * r*r));
	else if (improf == 3){
		if (r <=  spotsize/2.0){
			weight.push_back(1.0);
		}
		else{
			weight.push_back(0.0);
		}
	}
	else if (improf == 4){
		if( r > spotsize/2.0)
			weight.push_back(0.0);
		else { 
			if(r < dobs/2.0)
				weight.push_back(0.0);
			else
				weight.push_back(1.0);
		}
	}
}
VecTemp.clear();
for(int i = 0; i<=tnum; i++){
	t = tmin + i * tsize;
	rsize.push_back(0.0);
        if (t <= 0.00001) 
            rsize[i] = (0.0005/rnum);
         else
            rsize[i] = (0.0005*(7.0 + log10(t))/rnum);
	
	for(int j = 0; j<=rnum; j++){

		r = j * rsize[i];
       		dummy.push_back(granny.MelTemp(r,t)); 
       		if(dummy[j] + Tbody > 100)
			TempFlag = 1; 
	}
	VecTemp.push_back(dummy);
	dummy.clear();


}
// ------------------------  check conservaction of energy ----------------------------------------
if (tmax <= 10.0E-6)
	dx = 0.0005/rnum;
else
	dx = 0.0005*(7.0 + log10(tmax))/rnum;

         Qkept = 0.0;

for(int j =0; j<rnum; j++){

	Qkept = Qkept + VecTemp[tnum][j] * 4.1868 * 4.0 * M_PI * dx*dx*dx * (3.0*(j+1)*(j+1) - 3.0*(j+1) + 1.0) / 3.0;
}
Qlost = Qtot - Qkept;

if (TempFlag == 1) {
	cout<<endl;
	cout<<"////////////////////////////////////////////////////////"<<endl;
	cout<<"--------------------------------------------------------"<<endl;
	cout<<"      Warning . . . 100 degrees exceeded."<<endl;
	cout<<"--------------------------------------------------------"<<endl;
	cout<<"////////////////////////////////////////////////////////"<<endl;
}
}

void BigTherm::BigTherm2(void){
//This will run the equivalent of the Original bigtherm2, and write it to the file listed for output in the config file. it writes the entirety of the box and the entirety of the given time.
ofstream writefile(TempOutputFile.c_str(),ios::out);
writefile<<"t\t"<<"x\t"<<"y\t"<<"z\t"<<"Temp"<<"\n";
for(int i = 0; i<xnum+1; i++){
	x = xmin + i * xsize;
	for(int ii = 0; ii<ynum+1; ii++){
		y = ymin + ii*ysize;
		for(int iii = 0; iii<znum+1; iii++){
			z = zmin + iii*zsize;
			for(int iv = 0; iv < tnum+1; iv ++){
				t = tmin + iv * tsize;
				double thing = BigThermTempAt(x,y,z,t);
				writefile<<t<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<thing<<endl;
			}
		}
	}
}
}

//This runs BigTherm2 for a given subset of the box. this can be used so that one does not output or calculate an unneeded data, ie points far out enough from the beam that there is no significant temperature rise
void BigTherm::BigTherm2(double outxmin, double outxmax, int outxnum, double outymin, double outymax, int outynum, double outzmin, double outzmax, int outznum){
ofstream writefile(TempOutputFile.c_str(),ios::out);
writefile<<"t\t"<<"z\t"<<"x\t"<<"y\t"<<"Temp"<<"\n";
double outysize, outxsize, outzsize;
outysize =(outymax - outymin)/outynum;
outxsize =(outxmax - outxmin)/outxnum;
outzsize =(outzmax - outzmin)/outznum;
//using dummy variables in case one wants to use bigtherm2(void) afterwards
for(int i = 0; i<tnum; i++){
        t = tmin + i * tsize;
        for(int ii = 0; ii<outznum; ii++){
                z = outzmin + ii * outzsize;
                for(int iii = 0; iii<outxnum; iii++){
                        x = outxmin + iii*outxsize;
                        for(int vi = 0; vi<outynum; vi++) {
                                y = outymin + vi*outysize;
                                double thing = BigThermTempAt(x,y,z,t);
                                writefile<<t<<"\t"<<z<<"\t"<<x<<"\t"<<y<<"\t"<<thing<<"\n";
                        }
                }
        writefile<<"\n";
        }
}

}
void takatadamage(char* configfilename){

	int DamFlag;

	int xnum, ynum, znum, tnum;

	double Pi, BoxLength, BoxDepth, a, MelDens;

	double xmin, xmax, xsize;

	double ymin, ymax, ysize;

	double zmin, zmax, zsize;

	double tmin, tmax, tsize;

	vector< double > Temp;

	double c11, c21, c12, c22;

	double irr,pulsedur,E;

	double x, y, z, t, Tbody, rnum, absorb;

	double TempA;

	double spotsize,spotsizex,spotsizey;

	float Term1,Term2;

	float Sum1, Sum2;

	float Integrand;

	float Damage;

	int old;

	string inputfile;

	string outputfile;

	boost::property_tree::ptree configfile;

	boost::property_tree::ini_parser::read_ini(configfilename, configfile);

	string OutputFile;

//-----initializing the output file

	outputfile = configfile.get<string>("Bigtherm.DamageOutputFile");

	fstream out(outputfile.c_str());

	out.open(outputfile.c_str());

	out<<"t\t"<<"x\t"<<"y\t"<<"z\t"<<"Damage\t"<<endl;

//!-----Initializing constants:

	Tbody=37.0;

	DamFlag=0;

	c11=149.0;

	c12=50000.0;

	c21=242.0;

	c22=80000.0;

/*

!-----Reading input variables-------

*/

	BoxLength=configfile.get<double>("Bigtherm.BoxLength");

	BoxDepth=configfile.get<double>("Bigtherm.BoxDepth");

	xnum=configfile.get<int>("Bigtherm.xnum");

	xmin=configfile.get<double>("Bigtherm.xmin");

	xmax=configfile.get<double>("Bigtherm.xmax");

	ynum=configfile.get<int>("Bigtherm.ynum");

	ymin=configfile.get<double>("Bigtherm.ymin");

	ymax=configfile.get<double>("Bigtherm.ymax");

	znum=configfile.get<int>("Bigtherm.znum");

	zmin=configfile.get<double>("Bigtherm.zmin");

	zmax=configfile.get<double>("Bigtherm.zmax");

	tnum=configfile.get<int>("Bigtherm.tnum");

	tmin=configfile.get<double>("Bigtherm.tmin");

	tmax=configfile.get<double>("Bigtherm.tmax");

	MelDens=configfile.get<double>("Bigtherm.MelDens");

	rnum=configfile.get<int>("Bigtherm.rnum");

	old=configfile.get<int>("Bigtherm.old");

	absorb=configfile.get<double>("Bigtherm.absorb");

//!-----Calculating necessary parameters--------

	xsize = (xmax - xmin)/xnum;

	ysize = (ymax - ymin)/ynum;

	zsize = (zmax - zmin)/znum;

	tsize = (tmax - tmin)/tnum;
/*
!***********************************************************************

!
!     Calculate Damage(x,y,z,t)
!

!***********************************************************************
*/

// initializing input file

	inputfile = configfile.get<string>("Bigtherm.tempfile");

	ifstream damagefile(inputfile.c_str());

	double dummy;

	string line;
	
	int dam;
//cout<<"x\t"<<"y\t"<<"z\t"<<"Damage?\t"<<endl;
if( old ==0)
	getline(damagefile, line);

//!------Stepping in x direction------------

	for(int i=0;i<xnum+1;i++){

		x = xmin + i * xsize;

//!---------Stepping in y direction------------
         
		for(int ii = 0; ii<ynum+1; ii++){

			y = ymin + ii * ysize;

//!------------Stepping in z direction-------------

			for(int iii = 0; iii<znum+1; iii++){

				z = zmin + iii * zsize;

//!---------------Stepping in time---------------------

				for(int iv = 0; iv<tnum+1; iv++){

					t = tmin + iv * tsize;

					if(old ==1){
						damagefile >> dummy;
						Temp.push_back(dummy);
					}

					if(old == 0){
						damagefile >>dummy>>dummy>>dummy>>dummy>>dummy;
						Temp.push_back(dummy);
					}
				}
//!------------Simpson's Rule for Damage Integral----------------

				Sum1=0.0;

				Sum2=0.0;

				if(Tbody + Temp[0] < 50.0) 
					Term1=exp(c11 - c12/(273 + Tbody + Temp[0]));
				else
					Term1=exp(c21 - c22/(273 + Tbody + Temp[0]));


				if(Tbody + Temp[tnum-1] < 50.0)
					Term2=exp(c11 - c12/(273 + Tbody + Temp[tnum-1]));
				else
					Term2=exp(c21 - c22/(273 + Tbody + Temp[tnum-1]));


				for(int iv =1; iv<tnum; iv+=2){

					if ((Tbody + Temp[iv]) < 50.0)
						Integrand=exp(c11 - (c12/(273 + Tbody + Temp[iv])));
					else
						Integrand=exp(c21 - c22/(273 + Tbody + Temp[iv]));

					Sum1=Sum1+Integrand;
				}

				Sum1*=4.0;

				for(int iv =2; iv<tnum-1; iv +=2){

					if ((Tbody + Temp[0]) < 50.0)
						Integrand = exp(c11 - (c12/(273 + Tbody + Temp[iv])));
					else
						Integrand=exp(c21 - c22/(273 + Tbody + Temp[iv]));
					
					Sum2+=Integrand;
				}

				Sum2*=2.0;

				Damage=tsize*(Term1 + Term2 + Sum1 + Sum2)/3.0;

				if(Damage >= 1.0)
					dam = 1;
				else
					dam = 0;

					cout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<Damage<<endl;
				
				if (Damage >=1.0)
					DamFlag=1;

				Temp.clear();
		}
	}
//remember to uncomment this once proven effective, as it makes comparison easier but plotting a bitch
//	cout<<endl;
}
	out.close();

	damagefile.close();
}
void BigTherm::TakataDamage(void){

	int DamFlag;

	vector< double > Temp;

	double c11, c21, c12, c22;


	double x, y, z, t;

	double TempA;

	double spotsize,spotsizex,spotsizey;

	float Term1,Term2;

	float Sum1, Sum2;

	float Integrand;

	float Damage;

	int old;


//-----initializing the output file


	ofstream DamageOut(DamageOutputFile.c_str(),ios::out);
	ofstream TempOut(TempOutputFile.c_str(),ios::out);

	DamageOut << "x\t" << "y\t" << "z\t" <<"Damage\t" <<endl;
	if(DamageInfo == 1){
		TempOut<<"t\t"<<"x\t"<<"y\t"<<"z\t"<<"Temp\t"<<endl;
	}

//!-----Initializing constants:

	Tbody=37.0;

	DamFlag=0;

	c11=149.0;

	c12=50000.0;

	c21=242.0;

	c22=80000.0;

/*

!-----Reading input variables-------

*/

//!-----Calculating necessary parameters--------

	xsize = (xmax - xmin)/xnum;

	ysize = (ymax - ymin)/ynum;

	zsize = (zmax - zmin)/znum;

	tsize = (tmax - tmin)/tnum;
/*
!***********************************************************************

!
!     Calculate Damage(x,y,z,t)
!

!***********************************************************************
*/

// initializing input file



	double dummy;
	
	int dam;

//!------Stepping in x direction------------

	for(int i=0;i<xnum+1;i++){

		x = xmin + i * xsize;

//!---------Stepping in y direction------------
         
		for(int ii = 0; ii<ynum+1; ii++){

			y = ymin + ii * ysize;

//!------------Stepping in z direction-------------

			for(int iii = 0; iii<znum+1; iii++){

				z = zmin + iii * zsize;

//!---------------Stepping in time---------------------

				for(int iv = 0; iv<tnum+1; iv++){

					t = tmin + iv * tsize;

					dummy = BigThermTempAt(x,y,z,t);
					Temp.push_back(dummy);
					if(DamageInfo == 1)
						TempOut<<t<<"\t"<<z<<"\t"<<x<<"\t"<<y<<"\t"<<dummy<<endl;
				}
//!------------Simpson's Rule for Damage Integral----------------

				Sum1=0.0;

				Sum2=0.0;

				if(Tbody + Temp[0] < 50.0) 
					Term1=exp(c11 - c12/(273 + Tbody + Temp[0]));
				else
					Term1=exp(c21 - c22/(273 + Tbody + Temp[0]));


				if(Tbody + Temp[tnum-1] < 50.0)
					Term2=exp(c11 - c12/(273 + Tbody + Temp[tnum-1]));
				else
					Term2=exp(c21 - c22/(273 + Tbody + Temp[tnum-1]));


				for(int iv =1; iv<tnum; iv+=2){

					if ((Tbody + Temp[iv]) < 50.0)
						Integrand=exp(c11 - (c12/(273 + Tbody + Temp[iv])));
					else
						Integrand=exp(c21 - c22/(273 + Tbody + Temp[iv]));

					Sum1=Sum1+Integrand;
				}

				Sum1*=4.0;

				for(int iv =2; iv<tnum-1; iv +=2){

					if ((Tbody + Temp[0]) < 50.0)
						Integrand = exp(c11 - (c12/(273 + Tbody + Temp[iv])));
					else
						Integrand=exp(c21 - c22/(273 + Tbody + Temp[iv]));
					
					Sum2+=Integrand;
				}

				Sum2*=2.0;

				Damage=tsize*(Term1 + Term2 + Sum1 + Sum2)/3.0;

				if(Damage >= 1.0)
					dam = 1;
				else
					dam = 0;

					DamageOut<<x<<"\t"<<y<<"\t"<<z<<"\t"<<dam<<endl;
				
				if (Damage >=1.0)
					DamFlag=1;

				Temp.clear();
		}
	}
//remember to uncomment this once proven effective, as it makes comparison easier but plotting a bitch
//	cout<<endl;
}
	DamageOut.close();
	TempOut.close();

}
