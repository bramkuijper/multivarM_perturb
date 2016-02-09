
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "bramauxiliary.h"

#define NDEBUG

using namespace std;

/////////////////// parameters ///////////////////
const int Npop = 5000; // populaiton size
const int Clutch = 10; // clutch size
const int NumGen = 50000; // total number of generations
const int skip = 100; // data output each nskip-th generation

double mu_g = 0.05;            // mutation rate
double sdmu_g = 0.4;			 // standard deviation mutation size
double mu_m = 0.05;            // mutation rate
double sdmu_m = 0.4;			 // standard deviation mutation size
double mu_b = 0.05;            // mutation rate
double sdmu_b = 0.4;			 // standard deviation mutation size

bool diagonal_only = false; // only m11 and m22 allowed to evolve
double sigma_p = sqrt(0.1); // phenotypic variance

// initial values
double int1_t0 = 0;
double int2_t0 = 0;
double rate1_t0 = 0;
double rate2_t0 = 0;
double ampl1_t0 = 0;
double ampl2_t0 = 0;

// values after perturbation
double rate1ptb = 0;
double rate2ptb = 0;
double int1ptb = 0;
double int2ptb = 0;
double ampl1ptb = 0;
double ampl2ptb = 0;
double phiptb = 0;

// current values
double int1 = 0;
double int2 = 0;
double rate1 = 0;
double rate2 = 0;
double ampl1 = 0;
double ampl2 = 0;



const double pi = M_PI;
double omega1= 0;
double c = 0;           // strength of selection
double phi = 0;         // time-advance of theta_0


/////////////////// counters ///////////////////
int NSurv = 0;          // count # of survivors
bool do_stats = 0;      // aux variable whether we should do stats or not
int generation = 0;     // counter of the number of generations
unsigned seed = 0;      // keep track of the random number seed

double zopt[2] = {0,0}; // keep track of current optimum
double omega2[1000];    // strength of selection

double epsilon[2] = {0,0}; // current environment

double m[2][2] = {{0,0},{0,0}}; // initial values of matrix M
double m0 = 0;			 // initial value of the maternal effect
double meanpsurvi = 0;  // track mean survival prob


/////////////////// define the population ///////////////////
struct Individual
{
    double g[2][2]; // genetic loci (elevation)
    double m[2][2][2]; // maternal effects loci
    double b[2][2][2]; // phenotypic plasticity loci
    double phen[2]; // resulting phenotypes

    double wi; //fitness value
};

typedef Individual Population[Npop];
Population Pop;
Population Survivors;
typedef Individual NewPopulation[Npop*2];
NewPopulation NewPop;



/////////////////// random-number generators ///////////////////
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 



/////////////////// output files ///////////////////
string filename("sim_multimat");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

#ifdef DISTRIBUTION
string filename_new2(create_filename("sim_multimat_dist"));
ofstream distfile(filename_new2.c_str());
#endif //DISTRIBUTION

void initArguments(int argc, char *argv[])
{
	c = atof(argv[1]);
	mu_g= atof(argv[2]);
	sdmu_g = atof(argv[3]);
	mu_m= atof(argv[4]);
	sdmu_m = atof(argv[5]);
	mu_b= atof(argv[6]);
	sdmu_b = atof(argv[7]);
	phi = atof(argv[8]);
    m[0][0] = atof(argv[9]);
    m[0][1] = atof(argv[10]);
    m[1][0] = atof(argv[11]);
    m[1][1] = atof(argv[12]);
    sigma_p = sqrt(atof(argv[13]));
    rate1_t0 = atof(argv[14]);
    rate2_t0 = atof(argv[15]);
    rate1ptb = atof(argv[16]);
    rate2ptb = atof(argv[17]);
    phiptb = atof(argv[18]);
    int1ptb = atof(argv[19]);
    int2ptb = atof(argv[20]);
    ampl1_t0 = atof(argv[21]);
    ampl2_t0 = atof(argv[22]);
    ampl1ptb = atof(argv[23]);
    ampl2ptb = atof(argv[24]);
    diagonal_only = atoi(argv[25]);
}

void MutateG(double &G)
{
	G += gsl_rng_uniform(r)<mu_g ? gsl_ran_gaussian(r,sdmu_g) : 0;
}

void MutateM(double &M)
{
	M += gsl_rng_uniform(r)<mu_m ? gsl_ran_gaussian(r,sdmu_m) : 0;
}

void MutateB(double &B)
{
	B += gsl_rng_uniform(r)<mu_b ? gsl_ran_gaussian(r,sdmu_b) : 0;
}
void WriteParameters()
{
	DataFile << endl
		<< endl
		<< "type:;" << "multivar M" << ";" << endl
		<< "diagonalM:;" << diagonal_only << ";" << endl
		<< "c:;" <<  c << ";"<< endl
		<< "sigma_p:;" <<  sigma_p << ";"<< endl
		<< "rate1:;" <<  rate1_t0 << ";"<< endl
		<< "rate2:;" <<  rate2_t0 << ";"<< endl
		<< "ampl1:;" <<  ampl1_t0 << ";"<< endl
		<< "ampl2:;" <<  ampl2_t0 << ";"<< endl
		<< "mu_g:;" <<  mu_g << ";"<< endl
		<< "mu_b:;" <<  mu_b << ";"<< endl
		<< "mu_m:;" <<  mu_m << ";"<< endl
		<< "mu_std_m:;" <<  sdmu_m << ";"<< endl
		<< "mu_std_g:;" <<  sdmu_g << ";"<< endl
		<< "mu_std_b:;" <<  sdmu_b << ";"<< endl
        << "rate1ptb:;" << rate1ptb << ";" << endl
        << "rate2ptb:;" << rate2ptb << ";" << endl
        << "int1:;" << int1_t0 << ";" << endl
        << "int2:;" << int2_t0 << ";" << endl
        << "int1ptb:;" << int1ptb << ";" << endl
        << "int2ptb:;" << int2ptb << ";" << endl
        << "ampl1ptb:;" << ampl1ptb << ";" << endl
        << "ampl2ptb:;" << ampl2ptb << ";" << endl
        << "phi:;" << phi << ";" << endl
        << "phiptb:;" << phiptb << ";" << endl
		<< "seed:;" << seed << ";"<< endl
		<< "initm11:;" << m[0][0] << ";"<< endl
		<< "initm12:;" << m[0][1] << ";"<< endl
		<< "initm21:;" << m[1][0] << ";"<< endl
		<< "initm22:;" << m[1][1] << ";"<< endl;
}

// initialize all the phenotypes
void Init()
{
    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set the seed to the random number generator
    // stupidly enough, this can only be done by setting
    // a shell environment parameter
    stringstream s;
    s << "GSL_RNG_SEED=" << setprecision(10) << seed;
    putenv(const_cast<char *>(s.str().c_str()));

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);


	// initialize the whole populatin
	for (int i = 0; i < Npop; ++i)
	{
        Pop[i].phen[0] = 0;
        Pop[i].phen[1] = 0;

		for (int j = 0; j < 2; ++j)
		{
            // go through the number of diploid loci
            // of an individual 
            Pop[i].g[j][0] = Pop[i].phen[j];
            Pop[i].g[j][1] = Pop[i].phen[j];

            for (int k = 0; k < 2; ++k)
            {
                Pop[i].m[j][k][0] = m[j][k];
                Pop[i].m[j][k][1] = m[j][k];
            }
        }
	}

    // initialize environment
    ampl1 = ampl1_t0;
    ampl2 = ampl2_t0;
    int1 = int1_t0;
    int2 = int2_t0;
    rate1 = rate1_t0;
    rate2 = rate2_t0;

}

void Create_Kid(int mother, int father, Individual &kid)
{
    // first the first phenotype
    kid.g[0][0] = Survivors[mother].g[0][gsl_rng_uniform_int(r, 2)];
    MutateG(kid.g[0][0]);
    kid.g[0][1] = Survivors[father].g[0][gsl_rng_uniform_int(r,2)];
    MutateG(kid.g[0][1]);

    // now the other phenotype
    kid.g[1][0] = Survivors[mother].g[1][gsl_rng_uniform_int(r,2)];
    MutateG(kid.g[1][0]);
    kid.g[1][1] = Survivors[father].g[1][gsl_rng_uniform_int(r,2)];
    MutateG(kid.g[1][1]);

    // environmental variance
    kid.phen[0] = gsl_ran_gaussian(r,sigma_p);
    kid.phen[1] = gsl_ran_gaussian(r,sigma_p);

    // Kid's phenotype according to Kirkpatrick & Lande 1989, eq.4
    for (int i = 0; i < 2; ++i)
    {
        // check whether genetic loci are ok
        assert(isnan(kid.g[i][0])  == 0);
        assert(isnan(kid.g[i][1])  == 0);

        kid.phen[i] += 0.5 * (kid.g[i][0] + kid.g[i][1]);

        for (int j = 0; j < 2; ++j)
        {
            if (!diagonal_only || i == j)
            {
                assert(isnan(Survivors[mother].phen[j]) == 0);
                // inherit and mutate m
                
                kid.m[i][j][0] = Survivors[mother].m[i][j][gsl_rng_uniform_int(r,2)];
                MutateM(kid.m[i][j][0]);

                kid.m[i][j][1] = Survivors[father].m[i][j][gsl_rng_uniform_int(r,2)];
                MutateM(kid.m[i][j][1]);

                kid.phen[i] += 0.5 * (kid.m[i][j][0] + kid.m[i][j][1]) * Survivors[mother].phen[j];
            }
            else
            {
                kid.m[i][j][0] = 0;
                kid.m[i][j][1] = 0;
            }

            kid.b[i][j][0] = Survivors[mother].b[i][j][gsl_rng_uniform_int(r,2)];
            MutateB(kid.b[i][j][0]);
            
            kid.b[i][j][1] = Survivors[father].b[i][j][gsl_rng_uniform_int(r,2)];
            MutateB(kid.b[i][j][0]);

            kid.phen[i] += .5 * (kid.b[i][j][0] + kid.b[i][j][1]) * epsilon[i];
        }

        assert(isnan(kid.phen[i]) == 0);
    }
}


void Survive()
{
    double psurvi;
    NSurv = 0;

    if (generation == rint(NumGen / 2))
    {
        rate1 = rate1ptb;
        rate2 = rate2ptb;
        int1 = int1ptb;
        int2 = int2ptb;
        ampl1 = ampl1ptb;
        ampl2 = ampl2ptb;
        phi = phiptb;
    }


    // normally fluctuating equilibria
    zopt[0] = int1 + ampl1 * epsilon[0];
    zopt[1] = int2 + ampl2 * epsilon[1];
 
    meanpsurvi = 0;

    for (int i = 0; i < Npop; ++i)
    {
        psurvi = 1.0;

        for (int j = 0; j < 2; ++j)
        {
            psurvi *= exp(-pow(Pop[i].phen[j] - zopt[j],2.0)/(2.0 * pow(c,2.0)));
        }

        meanpsurvi += psurvi;

        if (gsl_rng_uniform(r) < psurvi)
        {
            Survivors[NSurv++] = Pop[i];
        }
    }

    meanpsurvi /= Npop;
}


// do the stats 
void WriteData()
{
    double meang[2] = {0,0};
    double ssg[2] = {0,0};
    double meanm[2][2] = {{0,0},{0,0}};
    double ssm[2][2] = {{0,0},{0,0}};
    double meanb[2][2] = {{0,0},{0,0}};
    double ssb[2][2] = {{0,0},{0,0}};

    double meanphen[2] = {0,0};
    double ssphen[3] = {0,0,0};

    double ssg1g2 =  0;

    for (int i =  0; i < Npop; ++i)
    {
        ssg1g2 += (Pop[i].g[0][0] + Pop[i].g[0][1])*(Pop[i].g[1][0] + Pop[i].g[1][1]);

        for (int j = 0; j < 2; ++j)
        {
            meang[j] += 0.5 * (Pop[i].g[j][0] + Pop[i].g[j][1]);
            ssg[j] += pow(0.5 * (Pop[i].g[j][0] + Pop[i].g[j][1]),2.0);

            meanphen[j] += Pop[i].phen[j];
            ssphen[j] += pow(Pop[i].phen[j],2.0);

            for (int k = 0; k < 2; ++k)
            {
                meanm[j][k] += 0.5 * (Pop[i].m[j][k][0] + Pop[i].m[j][k][1]);
                ssm[j][k] += pow(0.5 * (Pop[i].m[j][k][0] + Pop[i].m[j][k][1]),2.0);

                meanb[j][k] += 0.5 * (Pop[i].b[j][k][0] + Pop[i].b[j][k][1]);
                ssb[j][k] += pow(0.5 * (Pop[i].b[j][k][0] + Pop[i].b[j][k][1]),2.0);
            }

        }
        ssphen[2] += Pop[i].phen[0] * Pop[i].phen[1];
    }

    double meangsurv[2] = {0,0};
    double ssgsurv[2] = {0,0};
    double meanmsurv[2][2] = {{0,0},{0,0}};
    double ssmsurv[2][2] = {{0,0},{0,0}};
    double meanbsurv[2][2] = {{0,0},{0,0}};
    double ssbsurv[2][2] = {{0,0},{0,0}};

    double meanphensurv[2] = {0,0};
    double ssphensurv[2] = {0,0};

    for (int i = 0; i < NSurv; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            meangsurv[j] += 0.5 * (Survivors[i].g[j][0] + Survivors[i].g[j][1]);
            ssgsurv[j] += pow(0.5 * (Survivors[i].g[j][0] + Survivors[i].g[j][1]),2.0);

            meanphensurv[j] += Survivors[i].phen[j];
            ssphensurv[j] += pow(Survivors[i].phen[j],2.0);

            for (int k = 0; k < 2; ++k)
            {
                meanmsurv[j][k] += 0.5 * (Survivors[i].m[j][k][0] + Survivors[i].m[j][k][1]);
                ssmsurv[j][k] += pow(0.5 * (Survivors[i].m[j][k][0] + Survivors[i].m[j][k][1]),2.0);
                meanbsurv[j][k] += 0.5 * (Survivors[i].b[j][k][0] + Survivors[i].b[j][k][1]);
                ssbsurv[j][k] += pow(0.5 * (Survivors[i].b[j][k][0] + Survivors[i].b[j][k][1]),2.0);
            }
        }
    }

    DataFile << generation << ";" << NSurv << ";" << zopt[0] << ";" << zopt[1] << ";" << epsilon[0] << ";" << epsilon[1] << ";";

    for (int j = 0; j < 2; ++j)
    {
        DataFile 
            << (meanphen[j]/Npop) << ";"
        << (ssphen[j]/Npop- pow(meanphen[j]/Npop,2.0)) << ";"
        << (meang[j]/(Npop)) << ";"
        << (ssg[j]/(Npop) - pow(meang[j]/Npop,2.0)) << ";";

        for (int k = 0; k < 2; ++k)
        {
            DataFile << (meanm[j][k]/Npop) << ";";
            DataFile << (ssm[j][k]/Npop - pow(meanm[j][k]/Npop,2.0)) << ";";
            DataFile << (meanb[j][k]/Npop) << ";";
            DataFile << (ssb[j][k]/Npop - pow(meanb[j][k]/Npop,2.0)) << ";";
        }

        DataFile    << (meanphensurv[j]/NSurv) << ";"
        << (ssphensurv[j]/NSurv - pow(meanphensurv[j]/NSurv,2.0)) << ";"
        << (meangsurv[j]/(NSurv)) << ";"
        << (ssgsurv[j]/(NSurv) - pow(meangsurv[j]/NSurv,2.0)) << ";";

        for (int k = 0; k < 2; ++k)
        {
            DataFile << (meanmsurv[j][k]/NSurv) << ";";
            DataFile << (ssmsurv[j][k]/NSurv - pow(meanmsurv[j][k]/NSurv,2.0)) << ";";
            DataFile << (meanbsurv[j][k]/NSurv) << ";";
            DataFile << (ssbsurv[j][k]/NSurv - pow(meanbsurv[j][k]/NSurv,2.0)) << ";";
        }
    }

    DataFile << ssphen[2]/Npop - (meanphen[0]/Npop)*(meanphen[1]/Npop) << ";" << ssg1g2/Npop - (meang[0]/Npop)*(meang[1]/Npop) << ";" << meanpsurvi << ";" << endl;
}

void Reproduce()
{
    // first change the environment
    epsilon[0] = sin(rate1 * (generation + phi));
    epsilon[1] = sin(rate2 * generation);

    int Nkids = 0;

    if (NSurv < 2)
    {
        WriteData();
        WriteParameters();
        exit(1);
    }

    // individuals that have a  large share in the cumulative
    // fitness distribution produce a larger amount of 
    // offspring
    for (int i = 0; i < NSurv; ++i)
    {
        int mother = i;
        int father = gsl_rng_uniform_int(r,NSurv);

        Individual Kid;
        Create_Kid(mother, father, Kid);
        NewPop[Nkids++] = Kid;

        Individual Kid2;
        Create_Kid(mother, father, Kid2);
        NewPop[Nkids++] = Kid2;
    }

    for (int i = 0; i < Npop; ++i)
    {
        Pop[i] = NewPop[gsl_rng_uniform_int(r,Nkids)];
    }
}

void WriteDataHeaders()
{
    DataFile << "generation;NSurv;zopt1;zopt2;epsilon1;epsilon2"; 
    for (int j = 0; j < 2; ++j)
    {
        DataFile << ";meanphen" << (j+1) << ";ssphen" << (j+1) << ";meang" << (j+1) << ";varg" << (j+1);

        for (int k = 0; k < 2; ++k)
        {
            DataFile << ";meanm" << (j+1) << (k+1) << ";varm" << (j+1) << (k+1);
            DataFile << ";meanb" << (j+1) << (k+1) << ";varb" << (j+1) << (k+1);
        }
        
        DataFile << ";meanphensurv" << (j+1) << ";ssphensurv" << (j+1) << ";meangsurv" << (j+1) << ";vargsurv" << (j+1);

        for (int k = 0; k < 2; ++k)
        {
            DataFile << ";meanmsurv" << (j+1) << (k+1) << ";varmsurv" << (j+1) << (k+1);
            DataFile << ";meanbsurv" << (j+1) << (k+1) << ";varbsurv" << (j+1) << (k+1);
        }
    }

    DataFile << ";phencov;covg1g2;meanfitness;" << endl;
}


int main(int argc, char ** argv)
{
	initArguments(argc, argv);
	WriteDataHeaders();
	Init();

	for (generation = 0; generation <= NumGen; ++generation)
	{
        do_stats = generation % skip == 0;

		Survive();

        Reproduce();
        
        if (do_stats)
		{
			WriteData();
		}
	}

	WriteParameters();
}
