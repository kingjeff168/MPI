/*
Last Date Modified: 1/9/2024
*/
#include <iostream>
#include <cmath>
#include <random>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <chrono>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <mpi.h>


using namespace std;


bool isPositiveInteger(const string& str) 
{
    // This for-loop is to check if the input is an integer.
    for (char c : str) 
    {
        if (!isdigit(c)) 
        {
            return false;
        }
    }
    return true;
}

// This function is for calculating the second definite integral.
double exponential(double x) 
{
    return exp(-(x * x));
}

// This function is for calculating the first definite integral.
double square(double x) 
{
    return x * x;
}

// This function estimates the integral of f(x) over the interval [0, 1].
double estimateIntegral(unsigned long numSamples, int seed, int type) 
{
    double *samples = new double[numSamples];
    double sum;
    double estimateIntegral;

    // Use seed to generate random numbers.
    srand(seed);
    for (int i = 0; i < numSamples; i++)
    {
        samples[i] = (double)rand() / RAND_MAX;
    }

    // If the user wants to calculate the first definite integral, then the type is equal to 1;
    // If the user wants to calculate the second definite integral, then the else-statement will be conducted.
    sum = 0.0;
    if (type == 1)
    {
        for (int i = 0; i < numSamples; i++)
        {
            sum += square(samples[i]);
        }
    }
    else
    {
        for (int i = 0; i < numSamples; i++)
        {
            sum += exponential(samples[i]);
        }
    }
    
    // Calculate the final estimated integral from all samples.
    estimateIntegral = sum / numSamples;

    delete[] samples;
    return estimateIntegral;
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int myRank, numberProcess, seed, type, tempType;
    double localEstimateIntegral, localAverageEstimateIntegral, averageEstimateIntegral;
    unsigned long numberSamples, tempNumberSamples;

    
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberProcess);

    // Process with rank 0 determines which definite integral need to be estimated and the number of samples.
    if (myRank == 0)
    {
        for (int i = 1; i < argc; i++)
        {
            if (argc > 1)
            {
                string arg = argv[i];
                if ((arg == "-P") && ((i + 1) < argc))
                {
                    int i2 = i + 1;

                    if (isPositiveInteger(argv[i2]))
                    {
                        tempType = stoi(argv[i2]);
                    }
                }
                else if ((arg == "-N") && ((i + 1) < argc))
                {
                    int i3 = i + 1;
                    if (isPositiveInteger(argv[i3]))
                    {
                        tempNumberSamples = stol(argv[i3]);
                    }
                }
            }
        }
    }
    
    // Process with rank 0 initializes the number of samples and the type the user assign, and it broadcasts them to all other processes.
    if (myRank == 0)
    {
        type = tempType;
        numberSamples = tempNumberSamples;
    }
    MPI_Bcast(&type, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numberSamples, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    // Each process generates its local random samples and estimates the integral.
    seed = myRank * 1000;
    localEstimateIntegral = estimateIntegral(numberSamples, seed, type);

    // Calculate the sum across all processes using MPI_Reduce.
    localAverageEstimateIntegral = 0.0;
    MPI_Reduce(&localEstimateIntegral, &localAverageEstimateIntegral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // The process with rank 0 computes the global average by dividing the sum by the number of processes.
    averageEstimateIntegral = 0.0;
    if (myRank == 0)
    {
        averageEstimateIntegral = localAverageEstimateIntegral / numberProcess;
    }

    // The process with rank 0 prints the final result.
    if (myRank == 0)
    {
        cout << "The estimate for integral " << type << " is " << averageEstimateIntegral << endl;
        cout << "Bye!" << endl;
    }

    MPI_Finalize();


    return 0;
}
