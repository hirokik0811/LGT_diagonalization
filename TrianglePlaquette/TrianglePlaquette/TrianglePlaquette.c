// TrianglePlaquette.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <stdbool.h>
//#include "triangle_plaquette_hamiltonian.h"
#include "conversion_tensor.h"
//#include "square_tessellation.h"

//#define N_LAYERS 2
//#define G 0.1
//#define ALPHA 1.0
//#define EMIN 0.0 // the minimum of the search range of the feast algorithm
//#define EMAX 1000.0 // the maxmum 


// the arguments must be: N_LAYERS, G, ALPHA, EMIN, EMAX, -o outputfilename or/and -i inputfilename

/*
int main(int argc, char* argv[])
{
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
		  if(function != SPARSE_STATUS_SUCCESS)             \
		  {                                                 \
		  printf(error_message); fflush(0);                 \
		  status = function;                                       \
		  goto memory_free;                                 \
		  }                                                 \
} while(0)

	int i = 0, ii = 0, j = 0;

	if (argc != 6 && argc != 8 && argc != 10) {
		printf("invalid number of arguments. Aborting...");
		return 1;
	}

	const int nLayers = atoi(argv[1]);
	const double g = atof(argv[2]), alpha = atof(argv[3]), Emin= atof(argv[4]), Emax = atof(argv[5]);

	// if input file or output file storing the Hamiltonian matrix is specified, use them. 
	FILE *omatfile = NULL;
	FILE *imatfile = NULL;
	bool save_mat = false, load_mat = false;
	if (argc >= 8) {
		for (i = 6; i < argc-1; i+=2) {
			if (strcmp(argv[i], "-o") == 0) {
				omatfile = fopen(argv[i + 1], "w");
				if (omatfile == NULL) 
					printf("Cannot open the matrix output file. Aborting...");
				else
				    save_mat = true;
			}
			else if (strcmp(argv[i], "-i") == 0) {
				imatfile = fopen(argv[i + 1], "r");
				if (imatfile == NULL)
					printf("Cannot open the matrix input file. Aborting...");
				else
					load_mat = true;
			}
			else { printf("invalid flag %s. Aborting...", argv[i]); }
		}
	}

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations.
	sparse_matrix_t H = NULL;
	sparse_matrix_t zeroH = NULL;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	MKL_Complex16* valuesP = NULL;
	MKL_INT* pointerB_P = NULL;
	MKL_INT* pointerE_P = NULL;
	MKL_INT* columns_P = NULL;
	MKL_INT n_rowsP, n_colsP;

	// Variables for feast eigensolver
	char UPLO = 'U'; // Type of matrix: (F=full matrix, L/U - lower/upper triangular part of matrix);
	MKL_INT* rows_indx = NULL;
	MKL_INT fpm[128];      // Array to pass parameters to Intel(R) MKL Extended Eigensolvers 

	double       epsout;        // Relative error of the trace 
	MKL_INT      loop;          // Number of refinement loop 
	MKL_INT      M0;            // Initial guess for subspace dimension to be used 
	MKL_INT      M;             // Total number of eigenvalues found in the interval 

	double* E = NULL;      // Eigenvalues 
	MKL_Complex16* X = NULL;        // Eigenvectors 
	double* res = NULL;       // Residual 

	MKL_INT      info = 0;          // Errors 
	//


	
	
	if (load_mat) { // Load the Hamiltonian from the file
		char line[128];
		int line_cnt = 0; // count the lines
		MKL_INT row = 0, col = 0;
		double real = 0.0, imag = 0.0;
		MKL_INT prev_row = 0;
		while ((fgets(line, sizeof line, imatfile)) != NULL) {
			printf("Retrieved line %s\n", line);
			if (line_cnt == 0) {
				valuesP = (MKL_Complex16 *)mkl_malloc(atoi(line) * sizeof(MKL_Complex16), 64);
				columns_P = (MKL_INT*)mkl_malloc(atoi(line) * sizeof(MKL_INT), 64);
			} else if (line_cnt == 1) {
				pointerB_P = (MKL_Complex16*)mkl_malloc(atoi(line) * sizeof(MKL_Complex16), 64);
				pointerE_P = (MKL_Complex16*)mkl_malloc(atoi(line) * sizeof(MKL_Complex16), 64);
				n_rowsP = atoi(line);
				n_colsP = atoi(line);
			}
			else {
				sscanf(line, " { %d , %d } -> %lf + %lf I", &row, &col, &real, &imag);
				MKL_Complex16 val = { real, imag };
				valuesP[line_cnt - 2] = val;
				if (prev_row != row) {
					prev_row = row;
					pointerB_P[row-1] = line_cnt-1;
					if (row > 1) pointerE_P[row - 2] = line_cnt - 1;
				}
				columns_P[line_cnt - 2] = col;
			}
			line_cnt++;
		}
		pointerE_P[n_rowsP - 1] = line_cnt - 1;
		fclose(imatfile);
	}
	else {// Compute the Hamiltonian matrix of the triangle model and store it to H
		//CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(&H, 20, 20, listOfPauliList, coefList), "Error during computing a triangle plaquette matrix");

		CALL_AND_CHECK_STATUS(triangle_plaquette_hamiltonian_matrix(&H, nLayers, g, alpha), "Error during computing a triangle plaquette matrix");
		//CALL_AND_CHECK_STATUS(square_tessellation_hamiltonian_matrix(&H, nLayers, g, alpha), "Error during computing 8 triangles model Hamiltonian matrix");
		// Compute the block Hamiltonian with zero flux
		CALL_AND_CHECK_STATUS(zero_gauge_block(&zeroH, H, nLayers),
			"Error during computing a block Hamiltonian corresponding to zero flux\n");

		// Check the data
		CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(zeroH, &indexing, &n_rowsP, &n_colsP, &pointerB_P, &pointerE_P, &columns_P, &valuesP),
			"Error after MKL_SPARSE_Z_EXPORT_CSR  H\n");

		// Convert zero-based indexing to one-based indexing
		for (i = 0; i < n_rowsP; ++i) {
			pointerB_P[i] ++;
		}
		pointerE_P[n_rowsP - 1]++;
		for (i = 0; i < pointerE_P[n_rowsP - 1] - 1; ++i) {
			columns_P[i]++;
		}
	}
	printf("%d\n", (pointerE_P[n_rowsP - 1] - 1)); // number of non-zero elements
	if (save_mat) fprintf(omatfile, "%d\n", (pointerE_P[n_rowsP - 1] - 1));
	printf("%d\n", n_rowsP); // number of rows
	if (save_mat) fprintf(omatfile,"%d\n", n_rowsP);
	for (i = 0; i < n_rowsP; i++)
	{
		for (j = pointerB_P[i]; j < pointerE_P[i]; j++)
		{
			printf(" { %d , %d } -> %5.3f + %5.3f I,\n", i + 1, columns_P[ii], valuesP[ii].real, valuesP[ii].imag); fflush(0);
			if (save_mat) fprintf(omatfile, " { %d , %d } -> %5.3f + %5.3f I,\n", i + 1, columns_P[ii], valuesP[ii].real, valuesP[ii].imag);
			ii++;
		}
	}
	if (save_mat) fclose(omatfile);

	// Solve eigenvalue problem
	printf("Computing eigenvalues...\n");
	rows_indx = (MKL_INT*)mkl_malloc((n_rowsP+1) * sizeof(MKL_INT), 64);
	for (i = 0; i < n_rowsP; ++i)
		rows_indx[i] = pointerB_P[i];
	rows_indx[n_rowsP] = pointerE_P[n_rowsP - 1];

	M0 = n_rowsP;
	E = (double*)malloc(n_rowsP * sizeof(double));
	X = (MKL_Complex16*)malloc((n_rowsP * n_rowsP ) * sizeof(MKL_Complex16));
	res = (double*)malloc(n_rowsP * sizeof(double));

	// Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters 
	feastinit(
		fpm // OUT: Array is used to pass parameters to Intel(R) MKL Extended Eigensolvers 
	);
	
	fpm[0] = 1; // Extended Eigensolver routines print runtime status to the screen. 

	zfeast_hcsrev(
		&UPLO,   // IN: UPLO = 'U', stores the upper triangle part of the matrix 
		&n_rowsP,      // IN: Size of the problem 
		valuesP,     // IN: CSR matrix A, values of non-zero elements 
		rows_indx,    // IN: CSR matrix A, index of the first non-zero element in row 
		columns_P,    // IN: CSR matrix A, columns indices for each non-zero element 
		fpm,     // IN/OUT: Array is used to pass parameters to Intel(R) MKL Extended Eigensolvers 
		&epsout, // OUT: Relative error of on the trace 
		&loop,   // OUT: Contains the number of refinement loop executed 
		&Emin,   // IN: Lower bound of search interval 
		&Emax,   // IN: Upper bound of search interval 
		&M0,     // IN: The initial guess for subspace dimension to be used. 
		E,       // OUT: The first M entries of Eigenvalues 
		X,       // IN/OUT: The first M entries of Eigenvectors 
		&M,      // OUT: The total number of eigenvalues found in the interval 
		res,     // OUT: The first M components contain the relative residual vector 
		&info    // OUT: Error code 
	);

	printf("FEAST OUTPUT INFO %d \n", info);
	if (info != 0)
	{
		printf("Routine zfeast_hcsrev returns code of ERROR: %i", (int)info);
		return 1;
	}
	
	double        Eig[64] = { 0.97904, 4.96887, 4.96887, 4.96887, 4.96887, 4.96887, 4.96887, 5., \
                               5., 5., 8.82186, 8.93845, 8.93845, 8.93845, 8.93845, 8.93845, \
                               8.93845, 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., \
                               9., 9., 9., 9., 9., 13., 13., 13., 13., 13., 13., 13., 13., 13., 13., \
                               13., 13., 13., 13., 13.0311, 13.0311, 13.0311, 13.0311, 13.0311, \
                               13.0311, 13.0616, 13.0616, 13.0616, 13.0616, 13.0616, 13.0616, \
                                13.1991 };       // Eig - array for storing exact eigenvalues 
	double        R[64];         // R = |E-Eig| 

	
	// Step 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
	// are the expected eigenvalues and E(i) are eigenvalues computed by CFEAST_HCSREV(). 
	//
	
	printf("Number of eigenvalues found %d \n", M);
	printf("Computed\n");
	printf("Eigenvalues\n");
	double eigabs = 0.0;
	for (i = 0; i < M; i++)
	{
		printf("%.7e\n", E[i]);
	}

memory_free:
	free(E);
	free(X);
	free(res);
	if (!load_mat) {
		status = mkl_sparse_destroy(H);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error after MKL_SPARSE_DESTROY(H) \n"); fflush(0);
		}
		status = mkl_sparse_destroy(zeroH);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error after MKL_SPARSE_DESTROY(zeroH) \n"); fflush(0);
		}
	}
	mkl_free_buffers();

	return 0;
	
}

*/

int main(int argc, char* argv[])
{
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
		  if(function != SPARSE_STATUS_SUCCESS)             \
		  {                                                 \
		  printf(error_message); fflush(0);                 \
		  status = function;                                       \
		  goto memory_free;                                 \
		  }                                                 \
} while(0)
	// Test trace
	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations.  

	sparse_matrix_t P = NULL;
	
	CALL_AND_CHECK_STATUS(conversion_tensor(&P, 1), "Error constructing a sigma tensor \n");
	MKL_Complex16* valuesP = NULL;
	MKL_INT* pointerB_P = NULL;
	MKL_INT* pointerE_P = NULL;
	MKL_INT* columns_P = NULL;
	MKL_INT n_rowsP, n_colsP;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
	
	// Check the data
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(P, &indexing, &n_rowsP, &n_colsP, &pointerB_P, &pointerE_P, &columns_P, &valuesP),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  P\n");


	int i, j, ii=0;
	// Convert zero-based indexing to one-based indexing
	
	for (i = 0; i < n_rowsP; ++i) {
		pointerB_P[i] ++;
	}
	pointerE_P[n_rowsP - 1]++;
	for (i = 0; i < pointerE_P[n_rowsP - 1] - 1; ++i) {
		columns_P[i]++;
	}
	for (i = 0; i < n_rowsP; i++)
	{
		for (j = pointerB_P[i]; j < pointerE_P[i]; j++)
		{
			printf(" { %d , %d } -> %5.3f + %5.3f I,\n", i + 1, columns_P[ii], valuesP[ii].real, valuesP[ii].imag); fflush(0);
			ii++;
		}
	}
memory_free:
	//Release matrix handle and deallocate arrays for which we allocate memory ourselves.
	status = mkl_sparse_destroy(P);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(P) \n"); fflush(0);
	}
}