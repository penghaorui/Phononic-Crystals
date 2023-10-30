/*************************************************************************************************************/
/*  
	Author information: Haorui Peng, written during his Bachelors study in China University of Geosciences(Wuhan), in 2013.
	This program computes the banggap structrues of 2D phononic crystals of periodic cylinder imbeddings of Nickel in an Aluminum base.

	reference: F. Wu et al., Band structure of elastic waves in the two dimensional periodic composites, ACTA ACUSTICA, 26(4), 2001.
			   M. Sigalas and E.N. Economou, Band structure of elastic waves in two dimensional systems, Solid State Communications, 86(3), 1993.
*/
/*************************************************************************************************************/


#include <stdio.h>
#include <math.h>

typedef struct G {
  	long double b1; // the x component 
	long double b2; // the y component 
} G;

/*total number of normalized frequencies*/
#define Nw 121
/*total number of reduced wave vectors, must be the divisible by 3*/
#define Nk 60

#define PI 3.141592653

/*************************************************************************************************************/
/* Material properties*/
/*************************************************************************************************************/

/* the cell constant of the crystal lattice*/
long double const a = 0.06;  

/* the radius of the cylinder*/
long double const r = 0.35 * a;

/* the filling ratio */
long double F = PI* r * r / a / a; 

/* the density of imbedding and base materials*/
long double const density_A = 8800;
long double const density_B = 2700;

/* the shear moduli of imbedding and base materials*/
long double const shear_modulus_A = 7.71 * pow(10, 10);
long double const shear_modulus_B = 2.561 * pow(10, 10);

#include "functions.h"


/*************************************************************************************************************/
// main function
void main()
{   
	char* filename = "wk--121NiAl.mat";
	/*shear velocity of the base materials*/
	long double shear_velocity_B = sqrt(shear_modulus_B/density_B);
	
	/*the reduced wave vector*/
	G K[Nk];
  
	/*counting indices used in for loop*/
	int i = 0;
	int j = 0;
	int k = 0;
	
	/*the matrix of the final results: normalized frequency for each Block wave vector
	of three areas, each has the sampling of 20*/
	static long double w[Nk][Nw] = {0};

    /*the matrix of density*/
    static long double matrix_density[Nw][Nw] = {0};

    /*the matrix of inverse density*/
	static long double inverse_density[Nw][Nw] = {0}; 
    
	/*the matrix of shear modulus*/
	static long double matrix_shear_modulus[Nw][Nw] = {0};
	
	/*the Up Heisenberg matrix used in the computation later*/
	static  long double Up_Heisenberg[Nw][Nw] = {0};
    
    /* x component of K*/
	long double k1 = 0;
    
	/* y component of K*/
	long double k2 = 0;

	/* the cell constant of the reciprocal lattice*/
	long double b = 2 * PI / a;    

	/*compute the density matrix*/
    calculate_density_matrix(b, matrix_density);		

	/*copy matrix_density to inverse_density */ 
	copy_matrix(Nw, matrix_density, inverse_density);

	/*compute the inverse matrix of the density matrix*/
	inverse_matrix(Nw, inverse_density); 

	/*************************************************************************************************************/
	/*compute three different wave vector directions, each sampled by Nk//3 */    
	for(i = (int)(Nk/3); i > 0; i--)
	{   
		k1 = i * 0.05;
		k2 = i * 0.05;
		K[k].b1 = k1 * b * 0.5;
		K[k].b2 = k2 * b * 0.5;

		k++;
	}

	for(i = 0; i < (int)(Nk/3); i++)
	{   
		k1 = i * 0.05;
		k2 = 0;		

		K[k].b1 = k1 * b * 0.5;
		K[k].b2 = k2 * b * 0.5;

		k++;
	}
    
	for(i = 0; i < (int)(Nk/3); i++)
	{   
		k1 = 1;
		k2 = i * 0.05;		
		K[k].b1 = k1 * b * 0.5;
		K[k].b2 = k2 * b * 0.5;

		k++;
	}


	/*************************************************************************************************************/
	/*loop over each wave vector*/
	for(i = 0; i < Nk; i++)
	{   
		printf("computing the wave vectors %d/%d\n",i,Nk);

		/*compute the shear modulus matrix*/
	    calculate_shear_modulus_matrix(b, K[i], matrix_shear_modulus);

		/*multiply the inverse of density matrix with the shear modulus matrix, and save it to Up_Heisenberg*/
	    times_matrix(Nw,inverse_density,matrix_shear_modulus,Up_Heisenberg);

        /*convert a matrix to Up Heisenberg matrix*/
	    calculate_Up_Heisenberg( Up_Heisenberg );

		/*apply QR decomposition to the Up Heisenberg matrix to compute the complex eigen values*/
	    QR_algorithm( Up_Heisenberg );
		
		/*extract the egen values and scale it properly*/
		for(j = 0; j < Nw; j++)
		{
			w[i][j] = sqrt( Up_Heisenberg[j][j] ) * 0.04 / PI / shear_velocity_B / 2;
		}

	}
	
	/*************************************************************************************************************/
	/*save to file*/		
	write_to_file(filename,w);

	    
}

	

