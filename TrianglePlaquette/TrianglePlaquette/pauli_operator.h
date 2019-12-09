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

	sparse_matrix_t S = NULL; // a sigma matrix to be inserted at each site

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 

	CALL_AND_CHECK_STATUS(sigma_matrix(dest, pauliList[0], coef), "Error constructing a sigma matrix \n");

	int i;
	for (i = 1; i < nQubits; ++i) {
		printf("%d th tensor product\n", i);
		CALL_AND_CHECK_STATUS(sigma_matrix(&S, pauliList[i], 1.0), "Error constructing a sigma matrix\n");
		CALL_AND_CHECK_STATUS(kronecker_sparse_z_csr_h(dest, *dest, S), "Error during kronecker_sparse_z_csr\n");
	}

memory_free:
	status = mkl_sparse_destroy(S);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(S) \n"); fflush(0);
	}

	return status;

}

#endif