#pragma once
#ifndef PAULIOPERATOR_INCLUDED
#define PAULIOPERATOR_INCLUDED
#include <stdio.h>
#include <math.h>
#include <mkl.h>
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

	sparse_matrix_t Prod = NULL; // a matrix to store the kronecker product
	sparse_matrix_t S = NULL; // a sigma matrix to be inserted at each site

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 

	// Create matrix descriptor
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	CALL_AND_CHECK_STATUS(sigma_matrix(&S, pauliList[0], coef), "Error constructing a sigma matrix \n");
	CALL_AND_CHECK_STATUS(mkl_sparse_copy(S, descr, &Prod), "Error copying a sigma matrix \n");

	int i;
	for (i = 1; i < nQubits; ++i) {
		CALL_AND_CHECK_STATUS(sigma_matrix(&S, pauliList[i], 1.0), "Error constructing a sigma matrix\n");
		CALL_AND_CHECK_STATUS(kronecker_sparse_z_csr(&Prod, Prod, S), ("Error during kronecker_sparse_z_csr at site %d \n", i));
	}
	CALL_AND_CHECK_STATUS(mkl_sparse_copy(Prod, descr, dest), "Error copying the kronecker result \n");

memory_free:
	status = mkl_sparse_destroy(Prod);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(Prod) \n"); fflush(0);
	}
	status = mkl_sparse_destroy(S);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(S) \n"); fflush(0);
	}

	return status;

}

#endif