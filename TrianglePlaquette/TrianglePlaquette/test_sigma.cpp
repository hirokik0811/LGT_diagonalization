// TrianglePlaquette.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
/*
#include <iostream>
#include "sigma_matrix.h"


int main()
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
	sparse_matrix_t I = NULL, X = NULL, Y = NULL, Z = NULL;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	MKL_Complex8* valuesI = NULL;
	MKL_INT* pointerB_I = NULL;
	MKL_INT* pointerE_I = NULL;
	MKL_INT* columns_I = NULL;
	MKL_INT n_rowsI, n_colsI;

	MKL_Complex8* valuesX = NULL;
	MKL_INT* pointerB_X = NULL;
	MKL_INT* pointerE_X = NULL;
	MKL_INT* columns_X = NULL;
	MKL_INT n_rowsX, n_colsX;

	MKL_Complex8* valuesY = NULL;
	MKL_INT* pointerB_Y = NULL;
	MKL_INT* pointerE_Y = NULL;
	MKL_INT* columns_Y = NULL;
	MKL_INT n_rowsY, n_colsY;

	MKL_Complex8* valuesZ = NULL;
	MKL_INT* pointerB_Z = NULL;
	MKL_INT* pointerE_Z = NULL;
	MKL_INT* columns_Z = NULL;
	MKL_INT n_rowsZ, n_colsZ;

	int ii = 0;

	CALL_AND_CHECK_STATUS(sigma_matrix(&I, 0), "Error in simga_matrix I \n");
	CALL_AND_CHECK_STATUS(sigma_matrix(&X, 1), "Error in simga_matrix X \n");
	CALL_AND_CHECK_STATUS(sigma_matrix(&Y, 2), "Error in simga_matrix Y \n");
	CALL_AND_CHECK_STATUS(sigma_matrix(&Z, 3), "Error in simga_matrix Z \n");

	CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr(I, &indexing, &n_rowsI, &n_colsI, &pointerB_I, &pointerE_I, &columns_I, &valuesI),
		"Error after MKL_SPARSE_C_EXPORT_CSR  I \n");
	CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr(X, &indexing, &n_rowsX, &n_colsX, &pointerB_X, &pointerE_X, &columns_X, &valuesX),
		"Error after MKL_SPARSE_C_EXPORT_CSR  X \n");
	CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr(Y, &indexing, &n_rowsY, &n_colsY, &pointerB_Y, &pointerE_Y, &columns_Y, &valuesY),
		"Error after MKL_SPARSE_C_EXPORT_CSR  Y \n");
	CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr(Z, &indexing, &n_rowsZ, &n_colsZ, &pointerB_Z, &pointerE_Z, &columns_Z, &valuesZ),
		"Error after MKL_SPARSE_C_EXPORT_CSR  Z \n");

	printf("I \n");
	printf("Shape: %d x %d \n", n_rowsI, n_colsI);
	for (int i = 0; i < n_rowsI; i++)
	{
		printf("row#%d:", i + 1); fflush(0);
		for (int j = pointerB_I[i]; j < pointerE_I[i]; j++)
		{
			printf(" %5.3f + %5.3f I at col#%6d", valuesI[ii].real, valuesI[ii].imag, columns_I[ii] + 1); fflush(0);
			ii++;
		}
		printf("\n");
	}
	ii = 0;
	printf("X \n");
	printf("Shape: %d x %d \n", n_rowsX, n_colsX);
	for (int i = 0; i < n_rowsX; i++)
	{
		printf("row#%d:", i + 1); fflush(0);
		for (int j = pointerB_X[i]; j < pointerE_X[i]; j++)
		{
			printf(" %5.3f + %5.3f I at col#%6d", valuesX[ii].real, valuesX[ii].imag, columns_X[ii] + 1); fflush(0);
			ii++;
		}
		printf("\n");
	}
	ii = 0;
	printf("Y \n");
	printf("Shape: %d x %d \n", n_rowsY, n_colsY);
	for (int i = 0; i < n_rowsY; i++)
	{
		printf("row#%d:", i + 1); fflush(0);
		for (int j = pointerB_Y[i]; j < pointerE_Y[i]; j++)
		{
			printf(" %5.3f + %5.3f I at col#%6d", valuesY[ii].real, valuesY[ii].imag, columns_Y[ii] + 1); fflush(0);
			ii++;
		}
		printf("\n");
	}
	ii = 0;
	printf("Z \n");
	printf("Shape: %d x %d \n", n_rowsZ, n_colsZ);
	for (int i = 0; i < n_rowsZ; i++)
	{
		printf("row#%d:", i + 1); fflush(0);
		for (int j = pointerB_Z[i]; j < pointerE_Z[i]; j++)
		{
			printf(" %5.3f + %5.3f I at col#%6d", valuesZ[ii].real, valuesZ[ii].imag, columns_Z[ii] + 1); fflush(0);
			ii++;
		}
		printf("\n");
	}
	return 0;

memory_free:
	if (mkl_sparse_destroy(I) != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(I) \n"); fflush(0); status = mkl_sparse_destroy(I);
	}
	if (mkl_sparse_destroy(X) != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(X) \n"); fflush(0); status = mkl_sparse_destroy(X);
	}
	if (mkl_sparse_destroy(Y) != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(Y) \n"); fflush(0); status = mkl_sparse_destroy(Y);
	}
	if (mkl_sparse_destroy(Y) != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(Y) \n"); fflush(0); status = mkl_sparse_destroy(Z);
	}
}*/

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
