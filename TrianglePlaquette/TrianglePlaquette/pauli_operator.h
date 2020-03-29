#pragma once
#ifndef PAULIOPERATOR_INCLUDED
#define PAULIOPERATOR_INCLUDED
#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <stdbool.h>
#include "sigma_matrix.h"
#include "kronecker.h"

#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)

sparse_status_t pauli_operator_matrix(sparse_matrix_t* const dest, int const nQubits, int* const pauliList, double const coef) {

	sparse_matrix_t S = NULL; // a sigma matrix to be inserted at each site
	sparse_matrix_t temp = NULL; // a matrix to store the kronecker product

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	bool is_S_destroyed = true, is_temp_destroyed = true;
	CALL_AND_CHECK_STATUS(sigma_matrix(dest, pauliList[0], coef), "Error constructing a sigma matrix \n");

	int i;
	for (i = 1; i < nQubits; ++i) {

		// Compute Sigma
		CALL_AND_CHECK_STATUS(sigma_matrix(&S, pauliList[i], 1.0), "Error constructing a sigma matrix\n");
		is_S_destroyed = false;

		// Insert the Sigma and store the product to temp
		CALL_AND_CHECK_STATUS(kronecker_sparse_z_csr_h(&temp, *dest, S), "Error during kronecker_sparse_z_csr\n");
		is_temp_destroyed = false;
		// Destroy the old dest
		status = mkl_sparse_destroy(*dest);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(*dest) \n"); fflush(0);
		}
		
		// Copy temp to dest
		struct matrix_descr descr;
		descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
		descr.mode = SPARSE_FILL_MODE_UPPER;
		descr.diag = SPARSE_DIAG_NON_UNIT;
		CALL_AND_CHECK_STATUS(mkl_sparse_copy(temp, descr, dest), "Error during kronecker_sparse_z_csr\n");

		// Destroy temp
		status = mkl_sparse_destroy(temp);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(temp) \n"); fflush(0);
		}
		else {
			is_temp_destroyed = true;
		}

		// Destroy S
		status = mkl_sparse_destroy(S);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(S) \n"); fflush(0);
		}
		else {
			is_S_destroyed = true;
		}
	}

memory_free:
	if (!is_S_destroyed) {
		status = mkl_sparse_destroy(S);

		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(S) \n"); fflush(0);
		}
	}
	if (!is_temp_destroyed) {
		status = mkl_sparse_destroy(temp);

		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(temp) \n"); fflush(0);
		}
	}
	//mkl_free_buffers();
	return status;
	
}

#endif