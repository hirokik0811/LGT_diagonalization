// TrianglePlaquette.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <stdio.h>
#include "triangle_plaquette_hamiltonian.h"


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


	int ii = 0;

	CALL_AND_CHECK_STATUS(triangle_plaquette_hamiltonian_matrix(&P, 2, 1.0, 1.0), "Error during computing a triangle plaquette matrix");

	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(P, &indexing, &n_rowsP, &n_colsP, &pointerB_P, &pointerE_P, &columns_P, &valuesP),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  P\n");

	printf("P \n");
	printf("Shape: %d x %d \n", n_rowsP, n_colsP);
	for (int i = 0; i < n_rowsP; i++)
	{
		printf("row#%d:", i + 1); fflush(0);
		for (int j = pointerB_P[i]; j < pointerE_P[i]; j++)
		{
			printf(" %5.3f + %5.3f I at col#%6d", valuesP[ii].real, valuesP[ii].imag, columns_P[ii] + 1); fflush(0);
			ii++;
		}
		printf("\n");
	}

memory_free:
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
