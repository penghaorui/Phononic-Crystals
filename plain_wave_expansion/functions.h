long double func_modulus_density( long double fA, long double fB, struct G *p1, struct G *p2);
void calculate_density_matrix( long double b, long double (* t)[Nw] );
void calculate_shear_modulus_matrix( long double b, struct G K, long double (* t)[Nw] );
void print_matrix(int m, long double (* t)[Nw]);
void copy_matrix(int m, long double (* t1)[Nw], long double (* t2)[Nw]);
void inverse_matrix( int m, long double (* t)[Nw] );
void transpose_matrix( int m, long double (* t)[Nw] );
void times_matrix(int m, long double (* c)[Nw], long double (* z)[Nw], long double (* d)[Nw]);
long double fact( long double n );
long double fun_J1( long double x );
int sgn( long double x);
void calculate_Householder( int n, long double x[Nw], long double (* h)[Nw]);
void calculate_Up_Heisenberg( long double (* U)[Nw]);
void calculate_Givens(int k, long double cosx[Nw], long double sinx[Nw],  long double h[Nw][Nw]);
void QR_decompose(int n, long double U[Nw][Nw]);
void QR_algorithm(long double U[Nw][Nw]);
void write_to_file(char *filename, long double w[Nk][Nw]);
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

/*save the result to a file, which will be read later for plotting with Matlab */
void write_to_file(char *filename, long double w[Nk][Nw])
{	
	int i = 0;
	int j = 0;

	/*pointer for writing result to a file*/
	FILE *fp = fopen(filename,"w");

	if(NULL == fp)
	{
		printf("cannot open file\n");
		return;
	}

	for(i = 0; i < Nk; i++)
	{   
		for(j = 0; j < Nw; j++)
		{
			fprintf(fp,"%Lf ", w[i][j]);

		}
		fprintf(fp,"\n");
	}
	
	fclose(fp);	
}

/*compute the factorial of an integer n*/
long double fact( long double n )
{   
	register long double product = 1;

	while( n > 0 )
	{
       product = product * n;
	   n = n - 1;
	}
	return product;
}

/*compute the first type bessel function*/
long double fun_J1( long double x )
{   
	// count
	register long double k = 0;
	
	register long double sum = 0;
	register long n = 0;
	
	if(x >= 35 && x < 40)
	{
		n = 60;
	}
	else if(x >= 30 && x < 35)
	{
		n = 55;
	}
	else if(x >= 25 && x < 30)
	{
		n = 45;
	}
	else if(x >= 20 && x < 25)
	{
		n = 40;
	}
	else if(x >= 15 && x < 20)
	{
		n = 35;
	}
	else if(x >= 10 && x < 15)
	{
		n = 25;
	}
	else if(x >= 5 && x < 10)
	{
		n = 20;
	}
	else if(x >= 0 && x < 5)
	{
		n = 10;
	}

	while(k <= n)
	{
		sum = sum + pow(-1, k) * pow(x, 2 * k + 1) / pow(2, 2 * k + 1) / fact( k ) / fact(k + 1);
		k = k + 1;
	}
	
	return sum;
}

/*function used to construct reciprocal lattice vector for density and shear modulus */
long double func_modulus_density(long double fA, long double fB, struct G *p1, struct G *p2)
{   
	/* the shear modulus or density when G1 equals G2 */
    long double f_G1_G2; 

	/* the length of the reciprocal lattice vector, 
	   also the variable of the stucture function */
	long double G;

	/* the value of the structure function*/
	long double PG;

	
    // when G1 is not equal to G2
    if( ((*p1).b1 == (*p2).b1 ) && ( (*p1).b2 == (*p2).b2) )
	{
        f_G1_G2 = fA * F + fB * (1 - F); 
		return f_G1_G2;
	}
    // when G1 equals G2
   
	else
	{    
		G = sqrt( (p1->b1 - p2->b1) * (p1->b1 - p2->b1) + (p1->b2 - p2->b2) * (p1->b2 - p2->b2) );
		PG = 2 * F * fun_J1(G * r) / (G * r);
 
		return (PG * (fA - fB));		
	}
}

/*matrix multiplication*/
void times_matrix(int m, long double (* c)[Nw], long double (* z)[Nw], long double (* d)[Nw])
{   
	register int i;
	register int j;
	register int k;
	register long double f = 0;
	
	for(i = 0; i < m; i++ )
	{   
		for(j = 0 ; j < m; j++)
		{   
			f = 0;
			for(k = 0 ; k < m; k++)
			{
			    f = f + * ( * (c + i) + k) * * ( * (z + k) + j);	    
			}
			* ( * (d + i) + j) = f;
		}
	}

}
/*compute the density matrix*/
void calculate_density_matrix( long double b, long double (* t)[Nw] )
{   
	
    int i = 0;
	int j = 0;

    int m1;
	int m2;
	int n1;
	int n2;

	/* the reciprocal lattices G1 and G2 */
	G G1, G2;

    for(m2 = -5; m2 <= 5; m2++)	
    {   
		G1.b2 = m2 * b;

        for(m1 = -5; m1 <= 5; m1++)
		{   
			G1.b1 = m1 * b;

		    for(n2 = -5; n2 <= 5; n2++)
			{   
	            G2.b2 = n2 * b;

	            for(n1 = -5; n1 <= 5; n1++)
				{   
	                G2.b1 = n1 * b;
			        * ( * (t + i) + j) = func_modulus_density( density_A, density_B, &G1, &G2);
				    j++;
				}
			}

			j = 0;

			i++;
		}
	}
}

/*compute the shear modulus matrix*/
void calculate_shear_modulus_matrix( long double b, struct G K, long double (* t)[Nw] )
{   
    int i = 0;
	int j = 0;

    int m1;
	int m2;
	int n1;
	int n2;

	/* the reciprocal lattices G1 and G2 */
	G G1, G2;
	
	register long double times_vector = 0;

    for(m2 = -5; m2 <= 5; m2++)	
    {   G1.b2 = m2 * b;

        for(m1 = -5; m1 <= 5; m1++)
		{   
			G1.b1 = m1 * b;

		    for(n2 = -5; n2 <= 5; n2++)
			{   
	            G2.b2 = n2 * b;

	            for(n1 = -5; n1 <= 5; n1++)
				{   
	                G2.b1 = n1 * b;
                   
					times_vector = (K.b1 + G1.b1) * (K.b1 + G2.b1) + (K.b2 + G1.b2) * (K.b2 + G2.b2);
					
			        * ( * (t + i) + j) = func_modulus_density( shear_modulus_A, shear_modulus_B, &G1, &G2) * times_vector;
				    j++;

				}
			}

			j = 0;

			i++;
		}
	}
}

/*print out small matrix for examination*/
void print_matrix(int m, long double t[Nw][Nw])
{   
	int i;
	int j;
	int n;
	int k;

    for(i = 0; i < m; i++)
	{   k = 1;
        n = 0;
	    printf("\n\n******************\n\n");
        printf(" %d ",i);
		printf("\n\n******************\n\n");
	    for(j = 0; j < m; j++)
		{   
		    
			printf("%Lf  ",* ( * (t + i) + j));
		    n++;

		    if(n == 10)
			{   
				printf("   ****  %d  ****    ",k/10);	
		        printf("\n\n");
				n = 0;
			}
			
			k++;
		}

	} 
}

/*make a copy of a matrix*/
void copy_matrix(int m, long double (* t1)[Nw], long double (* t2)[Nw])
{
	int i;
	int j;

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < m; j++)
		{
			    
	        * ( * (t2 + i) + j) = * ( * (t1 + i) + j);
		}
	}
}

/*compute inverse matrix by the approach in the Numerical Analysis textbook, i.e., by diagonizing a matrix
that lies at the left side of a bigger matrix with double column numbers */
void inverse_matrix( int m, long double (* t)[Nw] )
{    
	register int i ;
	register int j ;
	int k ;
	long double n ;
	

	long double x = 0;
    

	long double v[Nw][Nw*2] = {0};
    

	for(i = 0; i < Nw; i++)
	{
		for(j = 0; j < Nw; j++)
		{
			    
	        v[i][j] = * ( * (t + i) + j);
		}
	}


    for( i = 0; i < m; i++)
	{
		v[i][m + i] = 1;
	}

    
    for(j = 0; j < m - 1; j++)           
	{	
		for(i = j + 1; i < m; i++)       
        {   
			x = - v[i][j] / v[j][j];

			for(k = j; k < 2 * m; k++)   
            {
			    v[i][k] = v[j][k] * x + v[i][k];
			}
		}	
	}

    
    for(j = m - 1; j > 0; j--)
	{
		for(i = j - 1; i >= 0; i--)
		{
			x = - v[i][j] / v[j][j];

			for(k = j; k < 2 * m ; k++)
			{
				v[i][k] = v[j][k] * x + v[i][k];
			}
		}
	}
	
	for(i = 0; i < m; i++)
	{	
		n = v[i][i];

		for( j = 0; j < 2 * m; j++)
		{
			v[i][j] = v[i][j] / n; 
		}
	}
    
	
	for(j = 0; j <  m; j++)
	{
		for(i = 0; i < m; i++)
		{
			v[i][j] = v[i][j + m];
		}
	}


    for(i = 0; i < Nw; i++)
	{
		for(j = 0; j < Nw; j++)
		{
			    
	        t[i][j] = v[i][j];
		}
	}
}

/*transpose a matrix*/
void transpose_matrix( int m, long double (* t)[Nw] )
{   
	int i = 0;
	int j = 0;
    long double c[Nw][Nw] = {0};
	long double (* p)[Nw];
    p = c;

    for(i = 0; i < Nw; i++)
	{
		for(j = 0; j < Nw; j++)
		{
			    
	        * ( * (p + j) + i) = * ( * (t + i) + j);
		}
	}

	copy_matrix(Nw, p, t);
}

/* the sign function*/
int sgn( long double x )
{
	if( x > 0)
	{
		return 1;
	}

	if( x == 0)
	{
		return 0;
	}

	if( x < 0)
	{
		return -1;
	}
}

/*compute Householder matrix, the math part can be found in the book Numerical Analysis.*/
void calculate_Householder( int n, long double x[Nw], long double (* h)[Nw])
{   

	register int i = 0;
	register int j = 0;

	int k = 0;
   
	long double beta =  0;
	long double sigma = 0;
	long double u[Nw] = {0};
	long double sigma1 = 0;

	long double alpha = 0;
	long double sum = 0;

	long double t[Nw][Nw] = {0};
    long double (*p)[Nw];
	p = t;
    

	for( i = 0; i < Nw; i++)
	{
		for( j = 0; j < Nw; j++)
		{
			* ( * (h + i) + j) = 0;

			if(i == j)
			{
                * ( * (h + i) + j) = 1;
			}
		}
	}


	alpha = x[0];

	for(i = 1; i < n; i++)
	{
        if(x[i] > alpha)
		{
			alpha = x[i];
		}
	}


	for(i = 0; i < n; i++)
	{
		x[i] = x[i] / alpha;
	}


	for(i = 0; i < n; i++)
	{
		sum = sum + x[i] * x[i];
	}
    
	sigma = sgn( x[0] ) * sqrt( sum );



	u[0] = x[0] + sigma;

	for(i = 1; i < n; i++)
	{
		u[i] = x[i];
	}




	beta = sigma * (sigma + x[0]);


	sigma1 = alpha * sigma;


	for( i = 0; i < n; i++)
	{
		for( j = 0; j < n; j++)
		{
			* ( * (p + i) + j) = -1 * u[i] * u[j] / beta;
		}
	}
    
	for( i = 0; i < n; i++)
	{	
		* ( * (p + i) + i) = 1 + * ( * (p + i) + i);
	}

	k = Nw - n;

	for( i = k; i < Nw; i++ )
	{	
		for(j = k; j < Nw; j++)
		{
			* ( * (h + i) + j) = * ( * (p + i - k) + j - k);
		}
	}

}
/*compute Up Heisenberg matrix, the math part can be found in the book Numerical Analysis.*/
void calculate_Up_Heisenberg( long double (* U)[Nw])
{   
	int i = 0;
    int j = 0;
	int k = 0;
	int n = 0;
    long double x[Nw] = {0};
    long double A[Nw][Nw] = {0};
    long double (* p)[Nw];
	
	long double c[Nw][Nw] = {0};
	long double d[Nw][Nw] = {0};
	p = A;


	copy_matrix(Nw, U, A);
    
    for(j = 0; j < 120; j++)
    {   
		n = 0;

	    for(i = 0; i < Nw; i++)
		{
		    x[i] = 0;
		}

		for(i = j + 1; i < Nw; i++)
		{	
			x[n] = * ( * (p + i) + j);
			n++;
		}
      
		k = 120 - j;
		calculate_Householder( k, x, U );
		times_matrix(Nw, U, p, c );
		times_matrix(Nw, c, U, d );
		copy_matrix(Nw, d, p);

		for( i = j + 2; i < Nw; i++ )
		{				
			    * ( * (p + i) + j) = 0;
			
		}
	}

	copy_matrix(Nw, A, U);
}

/*compute Givens coefficients, the math part can be found in the book Numerical Analysis.*/
void calculate_Givens(int k, long double cosx[Nw], long double sinx[Nw],  long double h[Nw][Nw])
{
    long double f = 0;

	f = sqrt( h[k][k] * h[k][k] + h[k + 1][k] * h[k + 1][k]);
	cosx[k] = h[k][k] / f;
	sinx[k]  = h[k + 1][k] / f;
}

/*compute QR decomposition, the math part can be found in the book Numerical Analysis.*/
void QR_decompose(int n, long double h[Nw][Nw])
{
	int i = 0;
	int j = 0;
	long double cosx[Nw] = {0};
	long double sinx[Nw] = {0};
    long double s = 0;
	int k = 0;
	long double t = 0;
	s = h[n][n]; 
    
	for(i = 0; i <= n; i++)
	{	
		h[i][i] = h[i][i] - s;		
	}
   
    for(k = 0; k <= n - 1; k++)
    {   
		calculate_Givens( k, cosx, sinx, h );
		
		for(j = k; j <= n; j++)
		{   
		    t = h[k][j];
		    h[k][j] = cosx[k] * h[k][j] +  sinx[k] * h[k + 1][j];
            h[k + 1][j] = -1 * sinx[k] * t + cosx[k] * h[k + 1][j];
		}
	}
    
	for(k = 0; k <= n - 1; k++)
	{
		for(i = 0; i <= k + 1; i++)
		{   
		    t = h[i][k];
		    h[i][k] = cosx[k] * h[i][k] + sinx[k] * h[i][k + 1];
            h[i][k + 1] = -1 * sinx[k] * t + cosx[k] * h[i][k + 1];
		}  
    }
	
	for(i = 0; i <= n; i++)
	{	
		h[i][i] = h[i][i] + s;		
	}
}


/*compute QR algorithm, the math part can be found in the book Numerical Analysis.*/
void QR_algorithm(long double U[Nw][Nw])
{
	int n = 120;
   
	while( n >= 1)
	{  
		while( fabs(U[n][n-1]) > 0.000000001)
		{   
            QR_decompose( n, U );
		}
		n--;
	}
}