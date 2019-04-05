#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include "vector3d.h"
#include <math.h>
#include <vector>
#include<gsl/gsl_rng.h>



using namespace std;




int main() {
  
  // Set up the ligand xyz coordinates in a vector array
  
  vector<VECTOR3D> ligand_xyz;
  double x,y,z;
  
  ifstream crds;                                      //open coordinates file
    crds.open("ligand_xyz.coords");
    if (!crds) {                                        //check to make sure file is there
        cerr << "ERR: FILE ligand_xyz.coords NOT OPENED. Check directory and/or filename.";
        exit(1);
    }
    for (int i = 0; i < 60; i++){
      crds >> x >> y >> z;				// Add coordinates the the vector array
      ligand_xyz.push_back(VECTOR3D(x,y,z));
    }
    
  
  // Get important parameters from the user
  
  double epsilon, wall_dist, nano_conc;
  int num_particles;
  
  cout << "What is the [ligand - Virus] complex concentration (nM)?" << endl;
  cin >> nano_conc;
  
  cout << "What is the receptor spacing (nm)?" << endl;
  cin >> wall_dist;
  
  cout << "What is the binding strength of the ligand-receptor interaction (KbT)?" << endl;
  cin >> epsilon;
  
  cout << "How many particles?" << endl;
  cin >> num_particles;
  
						//add options for ligand - virus ratio, properties of ligand / vlp, free ligand in solution
  
  // GENERATE P22CD40L.in FILE!!!
  
  // generate atoms
  
  ofstream datafile("P22CD40L.in");
  
  datafile << "#Script for P22-CD40L complex. To be read into LAMMPS template. " << endl;
  
  datafile << "#LAMMPS  data    file " << endl;
  
  double vlp_diameter = 1; //reduced units
  double ligand_diameter = 0.1; //reduced units
  double wall_diameter = 0.1; //reduced units
  
  double box_x, num_wall_x, num_atoms, num_bonds;
  box_x = pow ( (num_particles * 0.001 / (6.022E23))/(nano_conc * 1E-9) , (1.0/3.0) ); //in meters
  
   wall_dist = wall_dist * 1E-9; // in meters
  num_wall_x = box_x / wall_dist; //calculate number of mesh points
  num_wall_x = trunc(num_wall_x); // truncate this value
  cout << "adjusted spacing is " << (box_x/num_wall_x)*1E9 << " nm" << endl; //this is the adjusted spacing
  num_atoms = (num_particles * 61) + (num_wall_x * num_wall_x);
  num_bonds = num_particles * 180;
  
  double sigma = 60e-9; //sigma value currently used, diameter of P22
  double sigma_hc = 0.12; // sigma hc in reduced units
  double mass_SI = 3.819E-14; // mass of P22 in kg
  double epsilon_SI = (1.38E-23) * (298.15); // epsilon in J
  box_x = box_x / sigma; // now in reduced units
  wall_dist = box_x / num_wall_x; // now in reduced units
  
   if (box_x < (10 * vlp_diameter)){
    cout << "Uh oh! The z-direction is too small!" << endl;
  }
  
  
  datafile << num_atoms << " atoms" << endl;
  
  datafile << num_bonds << " bonds" << endl;
  
  datafile << "3 atom types" << endl;
  
  datafile << "3 bond types" << endl;
  
  datafile << "0" << " " << box_x << " xlo xhi" << endl; // in reduced units...
  
  datafile << "0" << " " << box_x << " ylo yhi" << endl;
  
  datafile << "0" << " " << box_x << " zlo zhi" << endl;
  
  datafile << endl;
  
  //atoms section
  
  datafile << "Atoms" << endl << endl;
  
  double atom_type, x_vlp, y_vlp, z_vlp;
  double index = 0;
  
  
    
   ////////////////////////////...To start on a lattice, use these lines.../////////////////////////////////////////
  
  unsigned int num_fill = int(ceil(pow((double(num_particles)), 1.0 / 3.0)));
  
//     for (unsigned int i = 0; i < num_fill; i++) {
//       for (unsigned int j = 0; j < num_fill; j++) {
// 	for (unsigned int k = 0; k < num_fill; k++) {
// 	  index += 1;
// 	  if ( (i+j+k) < num_particles ) {
// 	    x_vlp = (((double)i*box_x*(1/(double)num_fill)));
// 	    y_vlp = (((double)j*box_x*(1/(double)num_fill)));
// 	    z_vlp = (((double)k*(box_x - (vlp_diameter*2.5/2))*(1/(double)num_fill))) + (vlp_diameter*2.5);
// 	    datafile << index << "  " << "1" << "  " << x_vlp << "  " << y_vlp << "  " << z_vlp << "  " << endl;
// 	    for (unsigned int m = 0; m < ligand_xyz.size(); m++){
// 	      index += 1;
// 	      x = x_vlp + ligand_xyz[m].x;
// 	      y = y_vlp + ligand_xyz[m].y;
// 	      z = z_vlp + ligand_xyz[m].z;
// 	      datafile << index << "  " << "2" << "  " << x << "  " << y << "  " << z << "  " << endl;
// 	    }
// 	  }
// 	}
//       }
//     }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  
  ///////////////////////////////...To give vlp random positions use these.../////////////////////////////
  vector<VECTOR3D> vlp_xyz;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);                    //setting up random seed for brownian
  unsigned long int Seed = 23410981;
  gsl_rng_set(r, Seed);
  VECTOR3D dist;
  bool add;
  
  while (vlp_xyz.size() < num_particles) {
    add = true;
    x_vlp = (gsl_rng_uniform(r)) * box_x; // guess some coordinates
    y_vlp = (gsl_rng_uniform(r)) * box_x;
    z_vlp = (gsl_rng_uniform(r)) * box_x;
    
    if (z_vlp < 2.5) add = false; // if too close to the cell wall, reject the point
    if (z_vlp > (box_x - (vlp_diameter/2)) ) add = false; // if overlaps with the upper z boundary reject
    
    
    if (vlp_xyz.size() == 0) {		//if no vlps have been added yet, add it
      vlp_xyz.push_back(VECTOR3D(x_vlp, y_vlp, x_vlp));
      index += 1;
      datafile << index << " 1 " << "1" << " " << x_vlp << " " << y_vlp << " " << z_vlp << " " << endl;
      for (unsigned int m = 0; m < ligand_xyz.size(); m++){
	index += 1;
	x = x_vlp + ligand_xyz[m].x;
	y = y_vlp + ligand_xyz[m].y;
	z = z_vlp + ligand_xyz[m].z;
	datafile << index << " 1 " << "2" << " " << x << " " << y << " " << z << " " << endl;
      }
    }
    
    
    
    for (int j = 0; j < vlp_xyz.size(); j++){ //make sure they don't overlap with other particles
      dist.x = vlp_xyz[j].x - x_vlp;
      dist.y = vlp_xyz[j].y - y_vlp;
      dist.z = vlp_xyz[j].z - z_vlp;
      if (dist.x>box_x/2) dist.x -= box_x; //account for periodic bounaries in x & y  
      if (dist.x<-box_x/2) dist.x += box_x;
      if (dist.y>box_x/2) dist.y -= box_x;
      if (dist.y<-box_x/2) dist.y += box_x;
//       if (dist.z>box_x/2) dist.z -= box_x; 	// z is not periodic
//       if (dist.z<-box_x/2) dist.z += box_x;
   //   cout << "distance is " << dist.GetMagnitude() << " btw " << j << " ( " << vlp_xyz[j].x << " , " << vlp_xyz[j].y << " , " << vlp_xyz[j].z << " ) " << " and " << index/61 << " ( " << x_vlp << " , " << y_vlp << " , " << z_vlp << " ) " << endl;
      
      if (dist.GetMagnitude() < 2.5) {  //flag it if it intersects with anything
	//cout << "too close!" << endl;
	add = false;
	break;
    }  
  }
  
  if (add == true){		//if it doesn't intersect, add it and its ligands to the list
    vlp_xyz.push_back(VECTOR3D(x_vlp, y_vlp, z_vlp));
    index += 1;
    datafile << index << "  1 " << "1" << "  " << x_vlp << "  " << y_vlp << "  " << z_vlp << "  " << endl;
    for (unsigned int m = 0; m < ligand_xyz.size(); m++){
      index += 1;
      x = x_vlp + ligand_xyz[m].x;
      y = y_vlp + ligand_xyz[m].y;
      z = z_vlp + ligand_xyz[m].z;
      datafile << index << "  1 " << "2" << "  " << x << "  " << y << "  " << z << "  " << endl;
    }
  } //else cout << "rejected!" << endl;
  
}
    
    
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    
    //Making wall particles
    for (unsigned int i = 0; i < num_wall_x; i++){
      for (unsigned int j = 0; j < num_wall_x; j++){
	index += 1;
	datafile << index << " 1 " << "3" << "  " << i * wall_dist << "  " << j * wall_dist << "  " << "0" << endl;
      }
    }
    
    cout << "Done making atoms!" << endl;
    
    // Generate bonds
      
    vector<VECTOR3D> ligand_bond;
    int type, b1, b2;
    
      ifstream bonds;                                      //open coordinates file
    bonds.open("ligand_bond.coords");
    if (!bonds) {                                        //check to make sure file is there
        cerr << "ERR: FILE ligand_bond.coords NOT OPENED. Check directory and/or filename.";
        exit(1);
    }
    for (int i = 0; i < 180; i++){
      bonds >> index >> type >> b1 >> b2;				// Add coordinates the the vector array
      ligand_bond.push_back(VECTOR3D(type,b1,b2));
    }

    datafile << endl << "Bonds" << endl << endl;
    
     index = 0;
    
    for (int i = 0; i < num_particles; i++){
      for (int j = 0; j < ligand_bond.size(); j++){
	index += 1;	
	datafile << index <<  "     " << ligand_bond[j].x <<  "     " << ligand_bond[j].y + (i * 61) << "     " << ligand_bond[j].z + (i * 61) << endl;	
      }
      
    }
    
    //bond type 1 length
    dist = VECTOR3D(0,0,0) - ligand_xyz[0];
    double type1_length = dist.GetMagnitude();
    //bond type 2 length
    dist = ligand_xyz[0] - ligand_xyz[1];
    double type2_length = dist.GetMagnitude();
    //bond type 3 length
    dist = ligand_xyz[0] - ligand_xyz[20];
    double type3_length = dist.GetMagnitude();
    
    
    cout << "Done making bonds!" << endl;

  
  
  
  
  

  
  
  
  
  
  
  
  //GENERATE LAMMPS INPUT SCRIPT
	
  ofstream outfile("lammps_script.in");
  
  outfile << "#Lammps input script for nanoparticle binding with multivalency #" << endl << endl;
  
  outfile << "units lj" << endl << "boundary p p f" << endl << "atom_style bond " << endl << "neighbor 0.3 bin" << endl << "neigh_modify every 1 delay 0 check yes" << endl << endl;
  
  //Creating simulation box
  
  outfile << "## Create Simulation Box, Atoms ##" << endl;
  
  outfile << "read_data	   P22CD40L.in" << endl << endl;
  
  outfile << "group vlp type 1" << endl << "group ligand type 2" << endl << "group wall type 3" << endl << endl;
  
  //Defining particle properties
  
  outfile << "## Defining Particle/Medium Properties ##" << endl;
  
  outfile << "mass	1	590	# reduced mass of vlp P22" << endl;
  
  outfile << "mass	2	1	# reduced mass of ligand " << endl;
  
  outfile << "mass	3	1	# mass of wall, irrelevant " << endl << endl;
  
  outfile << "## Ascribing Initial Velocities ##" << endl;
  
  outfile << "velocity	all	create	1.0	4928459	rot	yes	dist	gaussian	units	box" << endl << endl;
  
  outfile << "## Fixing Wall particles ##" << endl;
  
  outfile << "fix zwalls all wall/reflect zhi EDGE zlo EDGE" << endl << endl;
  
  //Interparticle potentials
  
  outfile << "## Ascribing interparticle potentials: ##" << endl << endl;	
  
  outfile << "pair_style      hybrid lj/expand        2.5     lj/cut 2.5" << endl;
  
  outfile << "pair_coeff      1       1       lj/cut           1      " << vlp_diameter << "	" << vlp_diameter*1.12246   << "       # epsilon, sigma, cut-off  V-V" << endl;
  
  outfile << "pair_coeff      1       2       lj/expand        1      " << sigma_hc << "     " << (vlp_diameter + ligand_diameter)/2 - sigma_hc << "	" << (vlp_diameter + ligand_diameter)/2 * 1.12246 << "  # epsilon, sigma_hc, delta_V-L, cut-off" << endl;
  
  outfile << "pair_coeff      2       2       lj/cut           1      " << ligand_diameter << "     " << ligand_diameter*1.12246 << "       # L-L" << endl;
  
  outfile << "pair_coeff      1       3       lj/expand        1      " << sigma_hc << "     " << (vlp_diameter + wall_diameter)/2 - sigma_hc << "	" << (vlp_diameter + wall_diameter)/2 * 1.12246 << "  # V-W" << endl;
  
  outfile << "pair_coeff      2       3       lj/cut       " << epsilon << "	   " << ligand_diameter << "	" << "2.5" << "  # L-W" << endl;
  
  outfile << "pair_coeff      3       3       lj/cut           0      0               # W-W" << endl << endl;
  
  outfile << "#pair_modify  shift   yes     # the additive e_LJ for repulsion-only" << endl << endl;
  
  //Defining bond properties
  
  outfile << "## Making bonds ##" << endl;
  
  outfile << "bond_style harmonic" << endl;
  
  outfile << "bond_coeff 1    1000.0  " << type1_length << "  # type, K, distance" << endl;
  
  outfile << "bond_coeff 2    1000.0  " << type2_length << endl;
  
  outfile << "bond_coeff 3    1000.0  " << type3_length << endl << endl;
  
  //Ensemble fixes
  
  outfile << "## Ensemble Fixes (+ for output) ##" << endl;
  
  outfile << "variable        myTStep equal   0.001 # timestep where 1 MD step is " << sigma * sqrt( (mass_SI / 590) / (epsilon_SI) )  << " seconds" << endl;
  
  outfile << "timestep        ${myTStep}" << endl;
  
  outfile << "variable        myDStep equal   1000" << endl << endl;
  
  outfile << "fix     ens     vlp     nvt     temp    1.      1.      0.1" << endl;
  
  outfile << "fix     ens2    ligand  nvt     temp    1.      1.      0.1  # T_start, T_stop, T_damp=100*timestep" << endl;
  
  outfile << "fix_modify ens energy yes  # adds thermostat energy to potential" << endl;
  
  outfile << "fix_modify ens2 energy yes" << endl << endl;
  
  //Output dump
  
  outfile << "## Initial Image Dump ##" << endl;
  
  outfile << "dump mymovie all custom 100 ovitomovie.melt id type x y z" << endl << endl;
  
  outfile << "print \"OUTPUT COLUMNS: SIMULATION STEP NUMBER | TEMPERATURE | NP POTENTIAL ENERGY | SYSTEM VOLUME\" " << endl << endl;
  
  outfile << "thermo_style   custom  step    temp   etotal  ke      pe     #print info to file" << endl;
  
  outfile << "thermo  50000" << endl << endl;
  
  outfile << "run     1500000" << endl << endl;
  
  outfile << "## Defining Output Information ##" << endl;
  
  outfile << "dump    posD    all     custom  ${myDStep}      outfiles/dump.melt      id      type    x       y       z       # c_atomPot     c_atomKin" << endl << endl;
  
  outfile << " unfix   ens" << endl;
  
  outfile << "unfix   ens2" << endl;
  
  outfile << "undump  posD" << endl;
  
  outfile << "undump mymovie" << endl;
  
  outfile << "#shell          echo \"Lammps Simulation Ended\" " << endl;
  
  
  gsl_rng_free (r);
  
	return 0 ;
			
}







