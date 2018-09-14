/* Monte Carlo code to estimate the Lennard-Jones Potential and radial Distribution Function */
//=======================================================================
#include<iostream>
#include<math.h>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<iomanip>
#include<fstream>
#include <stdio.h>
#include <iomanip>
#include<sstream>
#include<string>

using namespace std;
#define _USE_MATH_DEFINES


double periodic_boundary_condition(double pair_distance, double box_length){   // Periodic boundary condition implemented
 if (pair_distance < - box_length/2){
  pair_distance = pair_distance + box_length;
 }
 else if (pair_distance > box_length/2){
  pair_distance = pair_distance - box_length;
 }
 return pair_distance;
}

//==========================================================================

void read_file(int n, double *x, double *y, double *z){   // function for reading a data file in xyz format
  double A[n+3][4];
  string line;
  char str[67];
  ifstream infile;
  infile.open("a.xyz");
  int counter = 0;
  while((counter < 2) && (infile.good())){   // reading the first two lines of xyz file
      getline (infile, line);
      counter++;
  }
  counter = 2;
  while((counter < (n+3)) && (infile.good())){   // reading the rest of lines from xyz file
      for (int j = 0; j < 4; j++){
        if (j == 0){
          infile >> str;
          //cout << str << "   \t   ";
        }
        else {
          infile >> A[counter][j];
        }
      }

      int j = 3;
      // incorporating the values inside x, y, z
      x[counter - 2] = A[counter][j-2];
      y[counter - 2] = A[counter][j-1];
      z[counter - 2] = A[counter][j];
      counter++;
  }
  infile.close();
}

//==========================================================================
double system_total_energy(int n, double box_length, double epsilon, double sigma, double rc, double *x, double *y, double *z){   // function for calculating the total system energy
  double energy_total = 0.0, E_stepwise;
  for (int i = 0; i < (n-1); i++){
    double W_stepwise = 0.0;
    for (int j = i+1; j < n; j++){
        double x_diff = periodic_boundary_condition((x[i] - x[j]), box_length);
        double y_diff = periodic_boundary_condition((y[i] - y[j]), box_length);
        double z_diff = periodic_boundary_condition((z[i] - z[j]), box_length);
        double R2 = (x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);   // calculating the distance between a pair

        if (R2 <= rc*rc)   // condition for checking the pair distance within the cut-off distance
        {
          E_stepwise = 4*epsilon*(pow(sigma, 12)/pow(R2, 6) - pow(sigma, 6)/pow(R2, 3));  // calculation of Lennard-Jones Potential
        }
        else {
          E_stepwise = 0;
        }
        W_stepwise = W_stepwise + E_stepwise;
    }
    energy_total = energy_total + W_stepwise;
  }
  return energy_total;
}
//==========================================================================

void rdf(int n, int binNo, double binSize, double box_length, double *count, double *x, double *y, double *z){   // bin-wise rdf calculation
  double dR[n][n];
  // Calculating the distances among all particles
  for (int i = 0; i < (n-1); i++){
    for (int j = i+1; j < n; j++){
      double dx = periodic_boundary_condition((x[i] - x[j]), box_length);
      double dy = periodic_boundary_condition((y[i] - y[j]), box_length);
      double dz = periodic_boundary_condition((z[i] - z[j]), box_length);
      dR[i][j] = pow((dx*dx + dy*dy + dz*dz), 0.5);
    }
  }

  // count the no. of particles in each bin
  for (int i = 0; i < (n-1); i++){
    for (int j = i+1; j < n; j++){
      double initial = 0.0;
      for (int k = 0; k < binNo; k++){
        if (dR[i][j] <= box_length/2 && dR[i][j] > initial && dR[i][j] <= (initial + binSize)){
          count[k] = count[k] + 1.0;
        }
        initial = initial + binSize;
      }
    }
  }
}

//==========================================================================

int main(){
  int n = 999, binNo = 500, mcMAX = 1000000;
  double rho = 0.85, K, box_length, binSize, kb = 1.0, T = 0.723, acc = 0.0, delta = 0.2, epsilon = 1.0, sigma = 1.0, rc = 2.5;

  box_length =  cbrt(n/rho);   // determination of Box length, L
  binSize = (box_length*0.5)/double(binNo);   // binsize calculation

  // declearing x, y, z as pointer array
  double *x = new double[n];
  double *y = new double[n];
  double *z = new double[n];
  double *count = new double[binNo];
  
  srand(time(0));
  for (int i = 0; i < binNo; i++){
     count[i] = 0;
  }

  read_file(n , x, y, z);   // reading the data set initially
  double energy_before_disp = system_total_energy(n, box_length, epsilon, sigma, rc, x, y, z);   // calculating the total system energy before monte carlo move

  ofstream outfile;
  outfile.open("energy_vs_mc_time_step_0.2_temp_0.723.dat");   //opening a file
  //outfile << "0" << " \t   " << energy_before_disp << endl;   // print the initial total energy

  for (int mc = 0; mc < mcMAX; mc++){
    int k = rand() % 1000;   // choosing a random particle
    
    double ranf = double(rand())/double(RAND_MAX);

    double dx = delta*(ranf - 0.5);
    double dy = delta*(ranf - 0.5);
    double dz = delta*(ranf - 0.5); 

    // update the x, y, z co-ordinates after the move
    x[k] = x[k] + dx;
    y[k] = y[k] + dy;
    z[k] = z[k] + dz;

    double energy_after_disp = system_total_energy(n, box_length, epsilon, sigma, rc, x, y, z);   // calculating the total system energy after monte carlo move

    if (energy_after_disp <= energy_before_disp){   // accept the move
      outfile << mc << "   \t    " << energy_after_disp << "   \t   " << k << "   \t   " << "1" << endl;
      acc = acc + 1.0;
      energy_before_disp = energy_after_disp;
    }
    else {
      double random = double(rand())/double(RAND_MAX);   //generating a random no between 0 and 1
      double energy_diff = energy_after_disp - energy_before_disp;
      double P = exp(-(energy_diff)/(kb*T));
      if (random < P){   // accept the move
        outfile << mc << "   \t    " << energy_after_disp << "   \t   " << k << "   \t   " << "1" << endl;
        acc = acc + random;
        energy_before_disp = energy_after_disp;
      }

      else {   // reject the move
        outfile << mc << "   \t    " << energy_before_disp << "   \t   " << k << "   \t   " << "0" << endl;
        x[k] = x[k] - dx;
        y[k] = y[k] - dy;
        z[k] = z[k] - dz;
      }
    }
    if(mc % 1000 == 0){
      rdf(n, binNo, binSize, box_length, count, x, y, z);   //here I have calculated the rdf in each 1000 step
    }
  }

  outfile.close();   // closing a file
  
  outfile.open("acceptance_prob_0.2_temp_0.723.dat");   // opening a file
  outfile << acc/double(mcMAX) << endl;
  outfile.close();

  double r[binNo], g[binNo];

  r[0] = 0.0;

  outfile.open("rdf_0.2_temp_0.4.dat");
  for (int i = 0; i < binNo; i++){
    g[i] = count[i];
    g[i] = g[i]/(((4.0*M_PI)/3.0)*(pow((r[i]+binSize),3) - pow(r[i],3))*rho);  // Normalization of g(r)
    outfile << r[i] << "   \t   " << (g[i]*2)/(n*1000) << endl;
    r[i+1] = r[i] + binSize;
  }

  outfile.close();   //closing a file


  return 0;
}
