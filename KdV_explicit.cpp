//  KdV - Program to solve the Korteweg-de Vries equation
//  by using the explicit scheme
#include <iostream>
#include <fstream>
#include <assert.h>
#include <cmath>
#include "Matrix.h"
using namespace std;

int main(){

  //* Initialize parameters (grid spacing, time step, etc.)
  int ini_offset = 4;
  cout << "Enter number of grid points: "; int N = 80; cout << N << endl;
  double L = 20;          // System extends from -L/2 to L/2
  double h = L/(N-1);     // Grid size
  cout << "Enter time step: "; double tau = 1e-4; cout << tau << endl;
  Matrix x(N);
  for (int i = 1; i <= N; i++){
    x(i) = h*(i-1) - L/2;  // Coordinates of grid points
  }

  //* Initialize the wavefunction
  Matrix rho(N);  rho.set(0.0);
  for (int i = 2; i <= (N-1); i++){
    rho(i) = 2 * pow(cosh(x(i) + ini_offset), -2);
  }

  //* Initializa loop and plot variables
  int nStep = 1 / tau;
  int nplots = 20;                      // Number of plots to record
  double plotStep = nStep / nplots;     // Iterations between plots

  Matrix rho_plot(N, nplots+2), rho_theo(N, nplots+2);
  for (int i = 1; i <= N; i++){   // Record initial condition
    rho_plot(i, 1) = rho(i);
    rho_theo(i, 1) = 2 * pow(cosh(x(i) + ini_offset), -2);
  }
  int iplot = 1;

  //* Loop over desired number of steps
  Matrix New_rho(N);
  int iStep;
  for (iStep = 1; iStep <= nStep; iStep++){

    //* Compute new rho
    New_rho(1) = rho(1) - 6 * tau * (rho(2) - 0) / 2 / h * rho(1) - tau / 2 / h / h / h * (rho(1+2) - 2 * rho(1+1) + 2 * 0 - 0);
    New_rho(2) = rho(2) - 6 * tau * (rho(3) - rho(1)) / 2 / h * rho(2) - tau / 2 / h / h / h * (rho(2+2) - 2 * rho(2+1) + 2 * rho(2-1) - 0);
    for (int i = 3; i <= (N-2); i++){
      New_rho(i) = rho(i) - 6 * tau * (rho(i+1) - rho(i-1)) / 2 / h * rho(i) - tau / 2 / h / h / h * (rho(i+2) - 2 * rho(i+1) + 2 * rho(i-1) - rho(i-2));
    }
    New_rho(N-1) = rho(N-1) - 6 * tau * (rho(N) - rho(N-2)) / 2 / h * rho(N-1) - tau / 2 / h / h / h * (0 - 2 * rho(N-1+1) + 2 * rho(N-1-1) - rho(N-1-2));
    New_rho(N) = rho(N) - 6 * tau * (0 - rho(N-1)) / 2 / h * rho(N) - tau / 2 / h / h / h * (0 - 2 * 0 + 2 * rho(N-1) - rho(N-2));
    rho = New_rho;
  
    //* Periodically record values for plotting
    if( fmod(iStep, plotStep) < 1 ){
      iplot++;
      for (int i = 1; i <= N; i++){
        rho_plot(i, iplot) = rho(i);
	rho_theo(i, iplot) = 2 * pow(cosh(x(i) - 4 * tau * iStep + ini_offset), -2);
      }
      cout << "Finish " << iStep << " of " << nStep << " steps" << endl;
    }
  }

  //* Record final probability density
  iplot++;
  for (int i = 1; i <= N; i++){
    rho_plot(i, iplot) = rho(i);
    rho_theo(i, iplot) = 2 * pow(cosh(x(i) - 4 * tau * iStep + ini_offset), -2);
  }
  nplots = iplot;

  //* Print out the plotting variables: x, rho_plot
  ofstream xOut("x_explicit.txt"), rho_plotOut("rho_explicit_plot.txt"), rho_theoOut("rho_theo_explicit_plot.txt");
  for (int i = 1; i <= N; i++){
    xOut << x(i) << endl;
    for (int j = 1; j < nplots; j++){
      rho_plotOut << rho_plot(i, j) << ", ";
    }
    rho_plotOut << rho_plot(i, nplots) << endl;
    for ( int j = 1; j < nplots; j++){
      rho_theoOut << rho_theo(i, j) << ", ";
    }
    rho_theoOut << rho_theo(i, nplots) << endl;
  } 
}

