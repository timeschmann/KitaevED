////////////////////////////////////////////////////////////////////////////////
// Exact Diagonalization: Kitaev Model, Honeycomb Lattice, 2x2 unit cells     //
// cylindrical boundary conditions
// Copyright (C) 2016 by Tim Eschmann
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
//#define ARMA_64BIT_WORD
#include <armadillo>

using namespace arma;
using namespace std;

// Flip k1-th and k2-th bit of integer number
int bitflip(int *v, int k1, int k2)
{
    int erg;
    erg = *v;
    // Flip bits:
    erg ^= 1 << k1;
    erg ^= 1 << k2; 
    return erg;    
}

// Calculate diagonal element (periodic boundary conditions)
int diagel(int *w)
{
    int b, b0, b1, b2, b3, b4, b5, b6, b7;
    signed int s1, s2, s3, s4;
    
    b = *w;
    // Check bits of states 
    b0 = (b >> 0) & 1; // 0th bit and so on ... 
    b1 = (b >> 1) & 1;
    b2 = (b >> 2) & 1;
    b3 = (b >> 3) & 1;
    b4 = (b >> 4) & 1;
    b5 = (b >> 5) & 1;
    b6 = (b >> 6) & 1;
    b7 = (b >> 7) & 1;
    
    if (b0 == b4)
        s1 = 1;
    else
        s1 = -1;

    if (b5 == b1)
        s2 = 1;
    else
        s2 = -1;

    if (b2 == b6)
        s3 = 1;
    else
        s3 = -1;

    if (b7 == b3)
        s4 = 1;
    else
        s4 = -1;

    return s1 + s2 + s3 + s4;
}

// Sign check for yy interactions
int signcheck(int *x, int var1, int var2)
{
    int c, c1, c2;
    signed int s;
    c = *x;
    c1 = (c >> var1) & 1;
    c2 = (c >> var2) & 1;
    
    if (c1 == c2)
        s = 1;
    else
        s = -1;
    
    return s;    
}


// Main function 
int main(void)
{
    int i,j,k,l,m,p; // running indices 
    int steps = 100; // How many measurement steps? 
    int Jx, Jy, Jz; // coupling constants 
    int N = 8; //  number of sites 
    int size = 256; // matrix size
    double Z, E; // partition sum, energy value 
    
    int xflips[3], yflips[3]; // Arrays for spin flips, size depending on # of bonds 
    signed int si0, si1, si2;
    
    vec ev;
    double temp[steps];
    double energy[steps];

    mat ham = mat(size,size, fill::zeros); // Initialize matrix 
    
    Jx = 1;
    Jy = 1;
    Jz = 1;
    
    cout << "Calculating matrix ..." << endl;
    
    // Fill Hamiltonian matrix with elements
    for (j = 0; j < size; j++)
    {
        // Diagonal element
        ham(j,j) = -Jz * diagel(&j);
        // Calculate interaction states via spin flips
        // xx:
        xflips[0] = bitflip(&j, 1, 2);
        xflips[1] = bitflip(&j, 4, 5);
        xflips[2] = bitflip(&j, 6, 7);

        //yy:
        yflips[0] = bitflip(&j, 0, 1);
        yflips[1] = bitflip(&j, 2, 3);
        yflips[2] = bitflip(&j, 5, 6);
        
        // signs for yy-interaction matrix elements
        si0 = signcheck(&j, 0, 1);
        si1 = signcheck(&j, 2, 3);
        si2 = signcheck(&j, 5, 6);
        
        ham(j,xflips[0]) = Jx;
        ham(j,xflips[1]) = Jx;
        ham(j,xflips[2]) = Jx;

        ham(xflips[0], j) = -Jx;
        ham(xflips[1], j) = -Jx;
        ham(xflips[2], j) = -Jx;

        ham(j,yflips[0]) = Jy*si0;
        ham(j,yflips[1]) = Jy*si1;
        ham(j,yflips[2]) = Jy*si2;

        ham(yflips[0], j) = -Jy*si0;
        ham(yflips[1], j) = -Jy*si1;
        ham(yflips[2], j) = -Jy*si2;

    }
    
    cout << "Diagonalizing ... " << endl;
    
    // Compute eigenvalues of ham:
    ev = eig_sym(ham);
    ev.save("eigenvalues.mat", csv_ascii);
    
    cout << "Calculating energy per temperature ... " << endl;
    
    // Calculating energy per temperature (with Boltzmann weights):
    for (l = 0; l < steps; l++)
    {
        Z = 0;
        E = 0;
        //temp[l] = 0.000001 + l*0.01; 
        temp[l] = pow(10,-2+(3*l/float(steps))); //logspace
        
        for (m = 0; m < size; m++)
        {
            Z = Z + exp(-ev[m]/temp[l]);
            E = E + ev[m]*exp(-ev[m]/temp[l]);        
        }
        energy[l] = E/(N*Z);
    } 
     
    // Write measured data to text file:    
    std::fstream f("dataexd.txt", std::ios::out);
    for (int i = 0; i < steps; ++i)
    {
        if (i == 0 || i == steps - 1)
            f << temp[i] << " " << energy[i] << " " << 0 << "\n";
        else 
            f << temp[i] << " " << energy[i] << " " << (energy[i+1] - energy[i-1])/(2*(temp[i+1] - temp[i-1])) << "\n";
    }
    f.close();
    
    return 0;
 
}