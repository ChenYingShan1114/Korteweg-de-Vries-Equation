//  KdV - Program to solve the Korteweg-de Vries equation
//  by using the implicit scheme (Crank-Nicolson)
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

  //* Set up the M and W operator matrix
  Matrix W_tmp(N, N), W(N, N);  W_tmp.set(0.0);  W.set(0.0);
  double coeff = tau / 2 / h / h / h;
  for (int i = 3; i <= (N-2); i++){
    W_tmp(i, i-2) = coeff;
    W_tmp(i, i-1) = -2 * coeff;
    W_tmp(i, i)   = 1;
    W_tmp(i, i+1) = 2 * coeff;
    W_tmp(i, i+2) = -coeff;
  }
  W_tmp(1, 1) = 1; W_tmp(1, 2) = 2 * coeff; W_tmp(1, 3) = -coeff;
  W_tmp(2, 1) = -2 * coeff; W_tmp(2, 2) = 1; W_tmp(2, 3) = 2 * coeff; W_tmp(2, 4) = -coeff;
  W_tmp(N-1, N-3) = coeff; W_tmp(N-1, N-2) = -2 * coeff; W_tmp(N-1, N-1) = 1; W_tmp(N-1, N) = 2 * coeff; 
  W_tmp(N, N-2) = coeff; W_tmp(N, N-1) = -2 * coeff; W_tmp(N, N) = 1;
  W = W_tmp;
 
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
    //* update W matrix by rho
    W(1, 1) = W_tmp(1, 1) - 6 * tau / 2 / h * (rho(2) - rho(1));
    for (int i = 2; i <= (N-1); i++){
      W(i, i) = W_tmp(i, i) - 6 * tau / 2 / h * (rho(i+1) - rho(i-1));
    }
    W(N, N) = W_tmp(N, N) - 6 * tau / 2 / h * (rho(N) - rho(N-1));

    //* Compute new rho using New_rho = W * rho
    for (int i = 1; i <= N; i++){
      double tmp = 0;
      for (int j = 1; j <= N; j++){
        tmp = tmp + W(i, j) * rho(j);
      }
      New_rho(i) = tmp;
    }
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
  ofstream xOut("x_explicit_matrix.txt"), rho_plotOut("rho_explicit_matrix_plot.txt"), rho_theoOut("rho_theo_explicit_matrix_plot.txt");
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

