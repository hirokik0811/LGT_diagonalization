// TrianglePlaquette.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#include <stdio.h>
//#include "kronecker.h"

/*
void test_kron()
{
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                       \
          goto memory_free;                                 \
          }                                                 \
    } while(0)

	// kron(I, X)
	MKL_INT status;
	MKL_Complex8* valuesA = (MKL_Complex8*)mkl_malloc(sizeof(MKL_Complex8) * 2, 64);
	MKL_Complex8* valuesB = (MKL_Complex8*)mkl_malloc(sizeof(MKL_Complex8) * 2, 64);
	MKL_Complex8* valuesC = NULL;
	MKL_INT* columns_A = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * 2, 64);
	MKL_INT* rowIndex_A = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * 3, 64);
	MKL_INT* columns_B = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * 2, 64);
	MKL_INT* rowIndex_B = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * 3, 64);
	MKL_INT* pointerB_C = NULL;
	MKL_INT* pointerE_C = NULL;
	MKL_INT* columns_C = NULL;
	MKL_INT n_rowsC, n_colsC;
	sparse_matrix_t A = NULL, B = NULL, C = NULL;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	int ii = 0;

	valuesA[0] = { 1, 0 }; valuesA[1] = { 1, 0 };
	columns_A[0] = 0; columns_A[1] = 1;
	rowIndex_A[0] = 0; rowIndex_A[1] = 1; rowIndex_A[2] = 2;
	valuesB[0] = { 1, 0 }; valuesB[1] = { 1, 0 };
	columns_B[0] = 1; columns_B[1] = 0;
	rowIndex_B[0] = 0; rowIndex_B[1] = 1; rowIndex_B[2] = 2;

	CALL_AND_CHECK_STATUS(mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 2, 2, rowIndex_A, rowIndex_A + 1, columns_A, valuesA),
		"Error in MKL_SPARSE_D_CREATE_CSR, A \n");
	CALL_AND_CHECK_STATUS(mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, 2, 2, rowIndex_B, rowIndex_B + 1, columns_B, valuesB),
		"Error in MKL_SPARSE_D_CREATE_CSR, B \n");
	CALL_AND_CHECK_STATUS(kronecker_sparse_c_csr(&C, A, B),
		"Error in kronecker_sparse_c_csr \n");

	CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr(C, &indexing, &n_rowsC, &n_colsC, &pointerB_C, &pointerE_C, &columns_C, &valuesC),
		"Error after MKL_SPARSE_D_EXPORT_CSR  \n");
	printf("result of kron (I, X) \n");
	printf("Shape: %d x %d \n", n_rowsC, n_colsC);
	for (int i = 0; i < n_rowsC; i++)
	{
		printf("row#%d:", i + 1); fflush(0);
		for (int j = pointerB_C[i]; j < pointerE_C[i]; j++)
		{
			printf(" %5.3f + %5.3f I at col#%6d", valuesC[ii].real, valuesC[ii].imag, columns_C[ii] + 1); fflush(0);
			ii++;
		}
		printf("\n");
	}

memory_free:
	mkl_free(valuesA); mkl_free(rowIndex_A); mkl_free(columns_A);
	mkl_free(valuesB); mkl_free(rowIndex_B); mkl_free(columns_B);
	if (mkl_sparse_destroy(A) != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(A) \n"); fflush(0); status = 1;
	}
	if (mkl_sparse_destroy(B) != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(B) \n"); fflush(0); status = 1;
	}
}
*/
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
