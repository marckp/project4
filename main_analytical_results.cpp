#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>



using namespace std;
using namespace arma;

ofstream ofile;

//function for periodic bondary conditions
inline int PeriodicBoundary(int i, int limit, int add){
    return (i+limit+add) % (limit);
}

//initialize energy and magnetization
void InitializeLattice(int, mat &, double&, double&);

//metroplis with monte carlo cycle loops
void MetropolisSampling(int, int, double, vec &);

//prints results to file
void WriteResultsToFile(int, int, double, vec, int counter);

//main program begins here
int main(int argc, char* argv[])
{
    string filename;
    int NSpins, MCcycles;
    double InitialTemp, FinalTemp, TempStep;
    if (argc <= 5){
        cout << "Bad Usage: " << argv[0] << " read outpute file, number of spins, MC cycles, final temperature step." << endl;
        exit(1);
    }
    else
    {
        filename            = argv[1];
        NSpins              = atoi(argv[2]);
        MCcycles            = atoi(argv[3]);
        InitialTemp         = atof(argv[4]);
        FinalTemp           = atof(argv[5]);
        TempStep            = atof(argv[6]);
    }
    //Declare new file name and add lattice size to file name
    string fileout      = filename;
    string Nstring      = to_string(NSpins);
    string MCstring     = to_string(MCcycles);
    string Tempstring   = to_string(FinalTemp);

    fileout.append("_spins");
    fileout.append(Nstring);
    fileout.append("_MC");
    fileout.append(MCstring);
    fileout.append("_temp");
    fileout.append(Tempstring);
    fileout.append(".txt");
    ofile.open(fileout);
    cout << filename << " " << fileout << endl;
    ofile << fileout << " " << Nstring << " " << MCstring << " " << Tempstring << endl;

    //Start Monte Carlo sampling by looping over temperatures
    for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature += TempStep){
        vec ExpectationValues = zeros<mat>(5);
        //Start MonteCarlo computation
        int counter = 1;
        MetropolisSampling(NSpins, MCcycles, Temperature, ExpectationValues);
        WriteResultsToFile(NSpins, MCcycles, Temperature, ExpectationValues, counter);
    }
    ofile.close();
    return 0;
}

//Make Monte Carlo and Metropolis function
void MetropolisSampling(int NSpins, int MCcycles, double Temperature, vec &ExpectationValues)
{
    //initialize
    std::random_device rd;
    std::mt19937_64 gen(rd());
    //set up uniform distirbution for x in [0,1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);
    //initialize lattice spin values
    mat SpinMatrix = zeros<mat>(NSpins, NSpins);
    //initialize energy and magnetization
    double Energy = 0.; double MagneticMoment = 0.;
    //initialize array for expectation values
    InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);
    //setup array for possible energy changes
    vec EnergyDifference = zeros<mat>(17);
    for(int de =- 8; de <= 8; de+= 4) EnergyDifference(de+8) = exp(-de/Temperature);
    //start Monte Carlo cycles
    int counter = 0;
    for(int cycles = 1; cycles <= MCcycles; cycles++)
    {
        //sweep over lattice, looping over spin sites
        for(int x = 0; x < NSpins; x++)
        {
            for(int y = 0; y < NSpins; y++)
            {
                // select a random spin element in the lattice
                int ix = (int)(RandomNumberGenerator(gen)*(double)NSpins);
                int iy = (int)(RandomNumberGenerator(gen)*(double)NSpins);
                // compute the change in energy as a result of fliping the randomly selected spin
                int deltaE = 2*SpinMatrix(ix, iy)*
                        (SpinMatrix(ix, PeriodicBoundary(iy, NSpins, -1)) +
                         SpinMatrix(PeriodicBoundary(ix, NSpins, -1), iy) +
                         SpinMatrix(ix, PeriodicBoundary(iy, NSpins, 1)) +
                         SpinMatrix(PeriodicBoundary(ix, NSpins, 1), iy));
                // determine if the new state is to be accepted or rejected according to the metropolis algorithm
                if (RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) )
                {
                    SpinMatrix(ix, iy) *= -1.0; //flip one spin and accept new spin config
                    // update magnetization and energy values for the new configuration and accepted state counter
                    MagneticMoment += (double) 2*SpinMatrix(ix, iy);
                    Energy += (double) deltaE;
                    counter += 1;
                    cout << counter << endl;
                }
            }
        }
        //update expecation values for local node
        ExpectationValues(0) += Energy;
        ExpectationValues(1) += Energy*Energy;
        ExpectationValues(2) += MagneticMoment;
        ExpectationValues(3) += MagneticMoment*MagneticMoment;
        ExpectationValues(4) += fabs(MagneticMoment);
        WriteResultsToFile(NSpins, cycles, Temperature, ExpectationValues, counter);
    }
}//end of Metropolis sampling over spins

//functon to initialize energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix, double& Energy, double& MagneticMoment){
    //setup initial energy (change the code here to initialize with an ordered or random configuration
    // setup spin matrix and initial magnetization
    srand(time(NULL));
    for(int x = 0; x < NSpins; x++){
        for(int y = 0; y < NSpins; y++){
            //SpinMatrix(x,y) = 1.0;
            double invers_period = 1.0/RAND_MAX;
            if (double(rand())*invers_period > 0.5) {
                SpinMatrix(x,y) = 1.0;
            }
            else
            {
                SpinMatrix(x,y) = -1.0;
            }

            MagneticMoment += (double)SpinMatrix(x,y);
        }
    }
    for(int x = 0; x < NSpins; x++)
    {
        for(int y = 0; y < NSpins; y++)
        {
            //SpinMatrix(x,y) = 1;
            Energy -= (double) SpinMatrix(x,y)*
                    (SpinMatrix(PeriodicBoundary(x, NSpins, -1),y) +
                     SpinMatrix(x, PeriodicBoundary(y, NSpins, -1)));
        }
    }
}//end function
// Write the results of each MC cycle to the output file
// Can be change to write different values as needed
void WriteResultsToFile(int NSpins, int MCcycles, double temperature, vec ExpectationValues , int counter)
{
    double norm                     = 1.0/((double) (MCcycles));  // divided by  number of cycles
    double E_ExpectationValues      = ExpectationValues(0)*norm;
    double E2_ExpectationValues     = ExpectationValues(1)*norm;
    double M_ExpectationValues      = ExpectationValues(2)*norm;
    double M2_ExpectationValues     = ExpectationValues(3)*norm;
    double Mabs_ExpectationValues   = ExpectationValues(4)*norm;
    // all expectation values are per spin, divide by 1/NSpins/NSpins

    double Evariance                = E2_ExpectationValues - E_ExpectationValues*E_ExpectationValues;
    double Mvariance                = M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << E_ExpectationValues;
    ofile << setw(15) << setprecision(8) << E2_ExpectationValues;
    ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues;
    ofile << setw(15) << setprecision(8) << M2_ExpectationValues;
    ofile << setw(15) << setprecision(8) << Evariance;
    ofile << setw(15) << setprecision(8) << Mvariance;
    ofile << setw(15) << setprecision(8) << counter << endl;
//    ofile << setw(15) << setprecision(8) << temperature;
//    ofile << setw(15) << setprecision(8) << E_ExpectationValues/NSpins/NSpins;
//    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
//    ofile << setw(15) << setprecision(8) << M_ExpectationValues/NSpins/NSpins;
//    ofile << setw(15) << setprecision(8) << Mvariance/temperature;
//    ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins;
//    ofile << setw(15) << setprecision(8) << counter << endl;
}
