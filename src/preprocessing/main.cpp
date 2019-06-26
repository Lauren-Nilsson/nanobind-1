#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include "vector3d.h"
#include <math.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

int main(int argc, char **argv) {

    /*Replacable variables*/
    string permitivityText = "USERINPUT_E";
    string totalTimeText = "USERINPUT_TIME";
    string initFileNameText = "USERINPUT_INIT_FILENAME";
    string mpiText = "NODESIZE";

    string initFile = "infiles/P22CD40L.in";

    // Set up the ligand xyz coordinates in a vector array

    vector<VECTOR3D> ligandCoordinates;
    double x;
    double y;
    double z;

    ifstream crds;                                      //open coordinates file
    crds.open("infiles/ligandCoordinates.coords");
    if (!crds) {                                        //check to make sure file is there
        cerr << "ERR: FILE ligandCoordinates.coords NOT OPENED. Check directory and/or filename.";
        exit(1);
    }
    for (int i = 0; i < 60; i++) {
        crds >> x >> y >> z;                // Add coordinates the the vector array
        ligandCoordinates.push_back(VECTOR3D(x, y, z));
    }


    // Get important parameters from the user

    double epsilon;
    double wallSpacing;
    double concentration;
    int numberComplexes;
    double timesteps;
    int ratio;
    bool verbose;

    options_description desc("Usage:\nrandom_mesh <options>");
    desc.add_options()
            ("help,h", "print usage message")
            ("lennard jones well depth,E", value<double>(&epsilon)->default_value(1),
             "Binding strength of the ligand-recptor interaction (KbT)")
            ("[ligand - Virus] complex concentration,C", value<double>(&concentration)->default_value(9),
             "[ligand - Virus] complex concentration (nM)") // box size adjusteded for nanomolar conc.
            ("receptor spacing,w", value<double>(&wallSpacing)->default_value(100), "receptor spacing (nm)")
            ("number of vlp-ligand complexes,S", value<int>(&numberComplexes)->default_value(108),
             "number of vlp-ligand complexes")
            ("vlp-ligand ratio,R", value<int>(&ratio)->default_value(60), "vlp-ligand ratio")
            ("simulation time,T", value<double>(&timesteps)->default_value(275), "simulation time (milliseconds)")
            ("verbose,V", value<bool>(&verbose)->default_value(true), "verbose true: provides detailed output");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    //add options for ligand - virus ratio, properties of ligand / vlp, free ligand in solution

    // GENERATE P22CD40L.in FILE!!!

    // generate atoms

    ofstream datafile(initFile);

    datafile << "#Script for P22-CD40L complex. To be read into LAMMPS template. " << endl;

    datafile << "#LAMMPS  data    file " << endl;

    double vlpDiameter = 1; //reduced units
    double ligandDiameter = 0.1; //reduced units
    double wallDiameter = 0.1; //reduced units

    double boxLength;
    double wallNumberX; //number of wall particles in x direction
    double atomNumber;
    double bondNumber;
    boxLength = pow((numberComplexes * 0.001 / (6.022E23)) / (concentration * 1E-9), (1.0 / 3.0)); //in meters

    wallSpacing = wallSpacing * 1E-9; // in meters
    wallNumberX = boxLength / wallSpacing; //calculate number of mesh points
    wallNumberX = trunc(wallNumberX); // truncate this value
    cout << "adjusted spacing is " << (boxLength / wallNumberX) * 1E9 << " nm" << endl; //this is the adjusted spacing
    atomNumber = (numberComplexes * (ratio + 1) ) + (wallNumberX * wallNumberX);
    bondNumber = numberComplexes * 180;

    double sigma = 60e-9; //sigma value currently used, diameter of P22
    double sigmaHC = 0.12; // sigma hc in reduced units
    double massSI = 3.819E-20; // mass of P22 in kg
    double epsilonSI = (1.38E-23) * (298.15); // epsilon in J
    boxLength = boxLength / sigma; // now in reduced units
    wallSpacing = boxLength / wallNumberX; // now in reduced units

    //convert simulation time to reduced units
    double tau = sigma * sqrt(massSI / epsilonSI); //1 MD timestep in seconds
    timesteps = ceil(timesteps * 1e-3 / tau);
    cout << "lammps will run for " << timesteps << " MD timesteps" << endl;

    if (boxLength < (10 * vlpDiameter)) {
        cout << "Uh oh! The z-direction is too small!" << endl;
    }


    datafile << atomNumber << " atoms" << endl;

    datafile << "3 atom types" << endl;

    datafile << "0" << " " << boxLength << " xlo xhi" << endl; // in reduced units...

    datafile << "0" << " " << boxLength << " ylo yhi" << endl;

    datafile << "0" << " " << boxLength << " zlo zhi" << endl;

    datafile << endl;

    //atoms section

    datafile << "Atoms" << endl << endl;

    double atomType;
    double vlpX;
    double vlpY;
    double vlpZ;
    double index = 0;
    vector<VECTOR3D> vlpXYZ;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);                    //setting up random seed
    //unsigned long int Seed = 23410982;
    gsl_rng_set(r, time(NULL));        //seed with time
    VECTOR3D dist;

    ///////////////////////////////...To give vlp random positions use these...////////////////////////////////
    bool add;
    int testLigand;
    vector<int> ligandsAdded;
    bool ligandAdd = true;

    while (vlpXYZ.size() < numberComplexes) {
        add = true;
        vlpX = (gsl_rng_uniform(r)) * boxLength; // guess some coordinates
        vlpY = (gsl_rng_uniform(r)) * boxLength;
        vlpZ = (gsl_rng_uniform(r)) * boxLength;

        if (vlpZ < 2.5) add = false; // if too close to the cell wall, reject the point
        if (vlpZ > (boxLength - (vlpDiameter + ligandDiameter)))
            add = false; // if overlaps with the upper z boundary reject


        if (vlpXYZ.size() == 0) {        //if no vlps have been added yet, add it
            vlpXYZ.push_back(VECTOR3D(vlpX, vlpY, vlpX));
            index += 1;
            datafile << index << " " << vlpXYZ.size() << " 1" << " " << vlpX << " " << vlpY << " " << vlpZ << " " << endl;
            while (ligandsAdded.size() < ratio) {
                ligandAdd = true;
                // Guess a ligand to add
                testLigand = floor(gsl_rng_uniform(r) * 60);
                //Make sure it hasn't been added already
                for (int j = 0; j < ligandsAdded.size(); j++) {
                   if (ligandsAdded[j] == testLigand) {
                      ligandAdd = false;
                      break;
                   }
               }
                //If it hasn't been added, push it to the vector and add to file
                if (ligandAdd == true) {
                   index += 1;
                   ligandsAdded.push_back(testLigand);
                   x = vlpX + ligandCoordinates[testLigand].x;
                   y = vlpY + ligandCoordinates[testLigand].y;
                   z = vlpZ + ligandCoordinates[testLigand].z;
                   datafile << index <<  " " << vlpXYZ.size() << " 2" << " " << x << " " << y << " " << z << " " << endl;
               }
                
            }
            //At the end of the loop, clear the ligandsAdded vector
            ligandsAdded.clear();
            add = false;
        }

        for (int j = 0; j < vlpXYZ.size(); j++) { //make sure they don't overlap with other particles
            dist.x = vlpXYZ[j].x - vlpX;
            dist.y = vlpXYZ[j].y - vlpY;
            dist.z = vlpXYZ[j].z - vlpZ;
            if (dist.x > boxLength / 2) dist.x -= boxLength; //account for periodic boundaries in x & y
            if (dist.x < -boxLength / 2) dist.x += boxLength;
            if (dist.y > boxLength / 2) dist.y -= boxLength;
            if (dist.y < -boxLength / 2) dist.y += boxLength;
            //       if (dist.z>boxLength/2) dist.z -= boxLength; 	// z is not periodic
            //       if (dist.z<-boxLength/2) dist.z += boxLength;
            //   cout << "distance is " << dist.GetMagnitude() << " btw " << j << " ( " << vlpXYZ[j].x << " , " << vlpXYZ[j].y << " , " << vlpXYZ[j].z << " ) " << " and " << index/61 << " ( " << vlpX << " , " << vlpY << " , " << vlpZ << " ) " << endl;

            if (dist.GetMagnitude() < 2.5) {  //flag it if it intersects with anything
                //cout << "too close!" << endl;
                add = false;
                break;
            }
        }


        if (add == true) {        //if it doesn't intersect, add it and its ligands to the list
            vlpXYZ.push_back(VECTOR3D(vlpX, vlpY, vlpZ));
            index += 1;
            datafile << index <<  " " << vlpXYZ.size() << " 1" << "  " << vlpX << "  " << vlpY << "  " << vlpZ << endl;
            while (ligandsAdded.size() < ratio) {
               ligandAdd = true;
               // Guess a ligand to add
               testLigand = floor(gsl_rng_uniform(r) * 60);
               //Make sure it hasn't been added already
               for (int j = 0; j < ligandsAdded.size(); j++) {
                  if (ligandsAdded[j] == testLigand) {
                     ligandAdd = false;
                     break;
                  }
               }
               //If it hasn't been added, push it to the vector and add to file
               if (ligandAdd == true) {
                  index += 1;
                  ligandsAdded.push_back(testLigand);
                  x = vlpX + ligandCoordinates[testLigand].x;
                  y = vlpY + ligandCoordinates[testLigand].y;
                  z = vlpZ + ligandCoordinates[testLigand].z;
                  datafile << index <<  " " << vlpXYZ.size() << " 2" << " " << x << " " << y << " " << z << " " << endl;
               }
               
            }
            //At the end of the loop, clear the ligandsAdded vector
            ligandsAdded.clear();
        } //else cout << "rejected!" << endl;

    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////   

//Making wall particles
    for (unsigned int i = 0; i < wallNumberX; i++) {
        for (unsigned int j = 0; j < wallNumberX; j++) {
            index += 1;
            datafile << index << " " << (vlpXYZ.size() + 1) << " 3" << "  " << i * wallSpacing << "  " << j * wallSpacing << "  " << "0.858"
                     << endl;
        }
    }

    cout << "Done making atoms!" << endl;


/*************************Generating the LAMMPS file*************************/

    ofstream inputScript("in.lammps", ios::trunc);
    if (inputScript.is_open()) {

        /*Open the template file*/
        string line;
        ifstream inputTemplate("infiles/in.lammps.template", ios::in);
        if (inputTemplate.is_open()) {
            while (getline(inputTemplate, line)) {

                std::size_t found = line.find(permitivityText);
                if (found != std::string::npos)
                    line.replace(found, permitivityText.length(),  std::to_string(epsilon));

                found = line.find(totalTimeText);
                if (found != std::string::npos)
                    line.replace(found, totalTimeText.length(), std::to_string(int(timesteps)));

                found = line.find(initFileNameText);
                if (found != std::string::npos)
                    line.replace(found, initFileNameText.length(), initFile);

                inputScript << line << endl;
            }
            inputTemplate.close();
        } else cout << "Unable to open the template input script" << endl;
        inputScript.close();
    } else cout << "Unable create a input Script" << endl;

    gsl_rng_free(r);

    return 0;

}







