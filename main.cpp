#include <fstream>
#include <iostream>
#include <filesystem>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <Eigen/Dense>
#include <deque>
#include <random>
#include <chrono>

// type definitions
typedef double REAL_SCALAR;
typedef std::complex<double> COMPLEX;
typedef Eigen::Matrix<double,3,1> VectorDimD;

void createOutputDirectory(const std::string& outPath) {
    std::filesystem::path dirPath = outPath;
    // Create a directory.
    try {
        if (std::filesystem::create_directory(dirPath)) {
            std::cout << "output data directory created successfully." << std::endl;
        } else {
            std::cout << "output data directory already exists or cannot be created." << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << e.what() << std::endl;
    }
}

void readFromFile(const std::string &filename, REAL_SCALAR *data) {
    std::ifstream file(filename);
    int number;
    int i = 0;

    while (file >> number) {
        data[i] = number;
        i++;
    }
}

Eigen::Matrix<int,2,1> Read_dimensions(const char *fname)
{
    int NX, NY, NZ;
    char line[200];
    FILE *InFile=fopen(fname,"r");

    if (InFile == NULL)
    {
        fprintf(stderr, "Can't open noise file %s\n",fname);
        exit(1);
    }
    // return the 5th line of the vtk file
    for(int i=0;i<4;i++)
    {
        fgets(line, 200, InFile);
    }
    // scan the returned line
    fscanf(InFile, "%s %d %d %d\n", line, &(NX), &(NY), &(NZ));
    return (Eigen::Matrix<int,2,1>()<<NX,NY).finished();
}

void correlationReader(const std::string& fileName_vtk, REAL_SCALAR *Rr_xy)
{
    //std::cout << "Reading stacking fault correlation" << std::endl;

    std::ifstream vtkFile(fileName_vtk); //access vtk file
    // error check
    if (!vtkFile.is_open()) {
        throw std::runtime_error("Error opening stacking fault VTK correlation file!");
    }

    std::deque<VectorDimD> rawPoints; // second pair is CELLID
    std::deque<VectorDimD> basis1;
    std::deque<VectorDimD> basis2;

    // begin parsing structured_grid vtk file for lines
    std::string line;
    const size_t numOfPointsPerLine = 3;
    while (std::getline(vtkFile, line)) 
    {
        //if the "POINTS" string is read, read the following data
        if(line.find("POINTS")!=std::string::npos) 
        {
            // get the number of points in the file
            const size_t firstSpace(line.find(' '));
            const size_t secondSpace(line.find(' ', firstSpace+1));
            const size_t numOfPoints = std::atoi(line.substr(firstSpace+1, secondSpace-firstSpace-1).c_str());
            // data structure of ID points
            // read the point coordinates
            double x1, y1, z1, x2, y2, z2, x3, y3, z3;
            for(size_t n=0; n<numOfPoints/numOfPointsPerLine; ++n)
            {
                std::getline(vtkFile, line);
                std::stringstream ss(line);
                ss >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
                rawPoints.push_back((VectorDimD()<<x1,y1,z1).finished()); // second pair is CELLID
                rawPoints.push_back((VectorDimD()<<x2,y2,z2).finished());
                rawPoints.push_back((VectorDimD()<<x3,y3,z3).finished());
            }
        }
        //if the "POINT_DATA" string is read, read the following data
        if(line.find("POINT_DATA")!=std::string::npos) 
        {
            const size_t numOfHeaders = 2;
            // get the number of points in the file
            const size_t firstSpace(line.find(' '));
            const size_t numOfPoints = std::atoi(line.substr(firstSpace+1).c_str());
            // read the point coordinates
            for(size_t n=0; n<numOfPoints+numOfHeaders; ++n)
            {
                std::getline(vtkFile, line);
                // ignore the headers that come right after point_data
                if(n<numOfHeaders)
                    continue;
                const int ind = n-numOfHeaders;
                //correlationCoeffs.push_back(std::atoi(line.c_str()));
                Rr_xy[ind] = std::atoi(line.c_str());
            }
        }

        if(line.find("VECTORS")!=std::string::npos) 
        {
            std::getline(vtkFile, line);
            std::stringstream ss(line);
            double x, y, z;
            if(basis1.empty())
            {
                ss >> x >> y >> z;
                basis1.push_back((VectorDimD()<<x,y,z).finished());
            }
            else
            {
                ss >> x >> y >> z;
                basis2.push_back((VectorDimD()<<x,y,z).finished());
            }
        }
    }
}

void sampleCorrelationWithNoise(const std::string& fileName_vtk, const std::string& outputFolder)
{
    Eigen::Matrix<int,2,1> dimensions;
    dimensions = Read_dimensions(fileName_vtk.c_str());
    int NX = dimensions(0);
    int NY = dimensions(1);
    //int NX = dimensions(1);
    //int NY = dimensions(0);
    int NZ = 1;
    int NR = NX*NY*NZ;
    int NK = NX*NY*(NZ/2+1);

    std::cout << "NX = " << NX << " NY = " << NY << " NZ = " << NZ << " NR = " << NR << " NK = " << NK << std::endl;

    // fft plans
    fftw_plan plan_R_isf_r2c, plan_R_isf_c2r, plan_Rcor_c2r;

    // allocate
    REAL_SCALAR *Rr_isf = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR); //correlation in real space
    COMPLEX *Rk_isf= (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK); //correlation in fourier space
    COMPLEX *frHat = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK); //correlation in fourier space
    REAL_SCALAR *Rr_noise_xy = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR); //correlation in real space with noise

    COMPLEX *frHatCorrelation = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK); //correlation in fourier space
    REAL_SCALAR *frCorrelation = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR); //correlation in real space with noise

    // read stacking fault correlation
    correlationReader(fileName_vtk, Rr_isf);

    //std::default_random_engine generator(seed);
    std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);

    plan_R_isf_r2c = fftw_plan_dft_r2c_2d(NY, NX, Rr_isf, reinterpret_cast<fftw_complex*>(Rk_isf), FFTW_ESTIMATE);
    plan_R_isf_c2r = fftw_plan_dft_c2r_2d(NY, NX, reinterpret_cast<fftw_complex*>(frHat), Rr_noise_xy, FFTW_ESTIMATE);
    plan_Rcor_c2r = fftw_plan_dft_c2r_2d(NY, NX, reinterpret_cast<fftw_complex*>(frHatCorrelation), frCorrelation, FFTW_ESTIMATE);

    // FFT to fourier space
    fftw_execute(plan_R_isf_r2c);

    // ensemble average
    int sampleSize = 1000;
    for (int j=0; j<sampleSize; ++j)
    {
        unsigned int randomSeed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(randomSeed);

        // introduce noise to correlation in fourier space
        for (int i=0; i<NK; ++i) 
        {
            // random numbers
            REAL_SCALAR Nk_xy = distribution(generator);
            REAL_SCALAR Mk_xy = distribution(generator);
            frHat[i] = std::sqrt(Rk_isf[i])*((Nk_xy+Mk_xy*COMPLEX(0.0,1.0))/std::sqrt(2.0));
        }

        // calculate spatial correlation of the newly created random function
        for (int i=0; i<NK; ++i) 
        {
            frHatCorrelation[i] += frHat[i]*std::conj(frHat[i]);
        }
    }

    // ensemble average
    for (int k=0; k<NK; ++k) 
    {
        frHatCorrelation[k] = frHatCorrelation[k]/static_cast<double>(sampleSize);
    }

    // take inverse fourier transform on the sampled random function (frHat -> ft)
    fftw_execute(plan_R_isf_c2r);

    // take inverse fourier transform on the sampled correlation ( frHat*frHat_conj -> C)
    fftw_execute(plan_Rcor_c2r);

    // normalize the correlation
    for (int i=0; i<NR; ++i) 
    {
        frCorrelation[i] = frCorrelation[i]/static_cast<double>(sampleSize);
    }

    //debugging
    //std::ofstream debugOut(outputFolder+"original_Correlation.txt", std::ios::app); //append to the file
    //std::ofstream debugOut2(outputFolder+"sampled_function.txt", std::ios::app);
    //std::ofstream debugOut3(outputFolder+"sampled_correlation.txt", std::ios::app);
    std::ofstream debugOut(outputFolder+"original_Correlation.txt", std::ios::out); //overwrite
    std::ofstream debugOut2(outputFolder+"sampled_random_function.txt", std::ios::out);
    std::ofstream debugOut3(outputFolder+"sampled_correlation.txt", std::ios::out);

    //debugOut << "Rr_isf = " << std::endl;
    for (int i=0; i<NR; ++i) {
        debugOut << Rr_isf[i] << std::endl;
    }

    //debugOut2 << "Rr_noise_isf = " << std::endl;
    for (int i=0; i<NR; ++i) {
        debugOut2 << Rr_noise_xy[i] << std::endl;
    }

    //debugOut3 << "frCorrelation= " << std::endl;
    for (int i=0; i<NR; ++i) {
        debugOut3 << frCorrelation[i]  << std::endl;
    }
}

int main() {
    // read the correlation file
    const std::string fileName_vtk = "../VTKfile/AlMg5_CxFFT_R100.vtk";

    // create sampled data output directory
    const std::string outputFolder = "../sampledData/";

    // create output directory
    createOutputDirectory(outputFolder);

    // sample the correlation with noise
    sampleCorrelationWithNoise(fileName_vtk, outputFolder);

    // exclusive to Linux, send a signal to the system that the program has ended successfully
    return 0;
}
