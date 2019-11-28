#pragma once
#ifndef SIGMA_INCLUDED
#define SIGMA_INCLUDED
#include <stdio.h>
#include <mkl.h>

sparse_status_t sigma_matrix(sparse_matrix_t* const S, int const id, double const coef) {
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)
	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	MKL_Complex16* values = (MKL_Complex16*)mkl_malloc(sizeof(MKL_Complex16) * 2, 64);
	MKL_INT* columns = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * 2, 64);
	MKL_INT* rowIndex = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * 3, 64);


	if (id == 1) {
		MKL_Complex16 val0 = { coef, 0 }, val1 = { coef, 0 };
		values[0] = val0; values[1] = val1;
		columns[0] = 1; columns[1] = 0;
		rowIndex[0] = 0; rowIndex[1] = 1; rowIndex[2] = 2;
	}
	else if (id == 2) {
		MKL_Complex16 val0 = { 0, -coef }, val1 = { 0, coef };
		values[0] = val0; values[1] = val1;
		columns[0] = 1; columns[1] = 0;
		rowIndex[0] = 0; rowIndex[1] = 1; rowIndex[2] = 2;
	}
	else if (id == 3) {
		MKL_Complex16 val0 = { coef, 0 }, val1 = { -coef, 0 };
		values[0] = val0; values[1] = val1;
		columns[0] = 0; columns[1] = 1;
		rowIndex[0] = 0; rowIndex[1] = 1; rowIndex[2] = 2;
	} else {
		MKL_Complex16 val0 = { coef, 0 }, val1 = { coef, 0 };
		values[0] = val0; values[1] = val1;
		columns[0] = 0; columns[1] = 1;
		rowIndex[0] = 0; rowIndex[1] = 1; rowIndex[2] = 2;
	}

	CALL_AND_CHECK_STATUS(mkl_sparse_z_create_csr(S, SPARSE_INDEX_BASE_ZERO, 2, 2, rowIndex, rowIndex + 1, columns, values),
		("Error in MKL_SPARSE_Z_CREATE_CSR, Sigma matrix %d \n", id));

memory_free:
	//mkl_free(values); mkl_free(columns); mkl_free(rowIndex);

	// return status
	return status;
}
#endif