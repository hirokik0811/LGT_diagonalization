// TrianglePlaquette.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <stdio.h>
#include <math.h>
#include "triangle_plaquette_hamiltonian.h"
#define N_LAYERS 2
#define G 1.0
#define ALPHA 1.0

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

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations.
	sparse_matrix_t P = NULL;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	MKL_Complex16* valuesP = NULL;
	MKL_INT* pointerB_P = NULL;
	MKL_INT* pointerE_P = NULL;
	MKL_INT* columns_P = NULL;
	MKL_INT n_rowsP, n_colsP;

	// Variables for feast eigensolver
	char UPLO = 'F'; // Type of matrix: (F=full matrix, L/U - lower/upper triangular part of matrix);
	MKL_INT* rows_indx = NULL;
	MKL_INT fpm[128];      /* Array to pass parameters to Intel(R) MKL Extended Eigensolvers */
	double        Emin = 0.9, Emax = 14.0;    /* Lower/upper bound of search interval [Emin,Emax] */

	double       epsout;        /* Relative error of the trace */
	MKL_INT      loop;          /* Number of refinement loop */
	MKL_INT      M0 = pow(2, 3*N_LAYERS);            /* Initial guess for subspace dimension to be used */
	MKL_INT      M;             /* Total number of eigenvalues found in the interval */

	double* E = NULL;      /* Eigenvalues */
	MKL_Complex16* X = NULL;        /* Eigenvectors */
	double* res = NULL;       /* Residual */

	MKL_INT      info = 0;          /* Errors */
	//


	int i = 0, ii = 0;


	// Compute the Hamiltonian matrix of the triangle model and store it to P
	CALL_AND_CHECK_STATUS(triangle_plaquette_hamiltonian_matrix(&P, N_LAYERS, G, ALPHA), "Error during computing a triangle plaquette matrix");

	// Check the data
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(P, &indexing, &n_rowsP, &n_colsP, &pointerB_P, &pointerE_P, &columns_P, &valuesP),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  P\n");

	// Convert zero-based indexing to one-based indexing
	for (i = 0; i < n_rowsP; ++i) {
		pointerB_P[i] ++;
	}
	pointerE_P[n_rowsP-1] ++;
	for (i = 0; i < pointerE_P[n_rowsP - 1] - 1; ++i) {
		columns_P[i]++;
	}

	printf("Hamiltonian: \n");
	printf("Shape: %d x %d \n", n_rowsP, n_colsP);
	for (i = 0; i < n_rowsP; i++)
	{
		for (int j = pointerB_P[i]; j < pointerE_P[i]; j++)
		{
			printf(" {%d, %d} -> %5.3f + %5.3f I,", i+1, columns_P[ii], valuesP[ii].real, valuesP[ii].imag); fflush(0);
			ii++;
		}
		printf("\n");
	}


	// Solve eigenvalue problem
	printf("Computing eigenvalues...\n");
	rows_indx = (MKL_INT*)mkl_malloc((n_rowsP+1) * sizeof(MKL_INT), 64);
	for (i = 0; i < n_rowsP; ++i)
		rows_indx[i] = pointerB_P[i];
	rows_indx[n_rowsP] = pointerE_P[n_rowsP - 1];

	E = (double*)malloc(n_rowsP * sizeof(double));
	X = (MKL_Complex16*)malloc((n_rowsP * n_rowsP ) * sizeof(MKL_Complex16));
	res = (double*)malloc(n_rowsP * sizeof(double));

	/* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
	feastinit(
		fpm /* OUT: Array is used to pass parameters to Intel(R) MKL Extended Eigensolvers */
	);
	
	fpm[0] = 1; /* Extended Eigensolver routines print runtime status to the screen. */

	zfeast_hcsrev(
		&UPLO,   /* IN: UPLO = 'F', stores the full matrix */
		&n_rowsP,      /* IN: Size of the problem */
		valuesP,     /* IN: CSR matrix A, values of non-zero elements */
		pointerB_P,    /* IN: CSR matrix A, index of the first non-zero element in row */
		columns_P,    /* IN: CSR matrix A, columns indices for each non-zero element */
		fpm,     /* IN/OUT: Array is used to pass parameters to Intel(R) MKL Extended Eigensolvers */
		&epsout, /* OUT: Relative error of on the trace */
		&loop,   /* OUT: Contains the number of refinement loop executed */
		&Emin,   /* IN: Lower bound of search interval */
		&Emax,   /* IN: Upper bound of search interval */
		&M0,     /* IN: The initial guess for subspace dimension to be used. */
		E,       /* OUT: The first M entries of Eigenvalues */
		X,       /* IN/OUT: The first M entries of Eigenvectors */
		&M,      /* OUT: The total number of eigenvalues found in the interval */
		res,     /* OUT: The first M components contain the relative residual vector */
		&info    /* OUT: Error code */
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
                                13.1991 };       /* Eig - array for storing exact eigenvalues */
	double        R[64];         /* R = |E-Eig| */

	/* Step 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
	* are the expected eigenvalues and E(i) are eigenvalues computed by CFEAST_HCSREV(). */
	printf("Number of eigenvalues found %d \n", M);
	printf("Computed      |    Expected    \n");
	printf("Eigenvalues   |    Eigenvalues \n");
	double eigabs = 0.0;
	for (i = 0; i < M; i++)
	{
		R[i] = fabs(E[i] - Eig[i]);
		eigabs = max(eigabs, R[i]);
		printf("%.7e  %.7e \n", E[i], Eig[i]);
	}
	printf("Max value of | computed eigenvalue(i) - expected eigenvalues(i) | %.7e \n", eigabs);


memory_free:
	free(E);
	free(X);
	free(res);
	if (mkl_sparse_destroy(P) != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(P) \n"); fflush(0); status = mkl_sparse_destroy(P);
	}

	return 0;
	
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file