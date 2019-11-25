#include <stdlib.h>
#include <mkl.h>
#include "pauli_operator.h"
#include "kronecker.h"
#include "sigma_matrix.h"

#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)

PauliOperator::PauliOperator(int const nQubits, int* const pauliList, float const alpha) {
	this->coef = alpha;
	this->nQubits = nQubits;
	this->pauliList = pauliList;

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	sparse_matrix_t S = NULL; // a sigma matrix to be inserted at each site

	CALL_AND_CHECK_STATUS(sigma_matrix(& (this->P), pauliList[0], alpha), "Error constructing a sigma matrix %d \n", pauliList[0]);
	for (int i = 1; i < nQubits; ++i) {
		CALL_AND_CHECK_STATUS(sigma_matrix(&S, pauliList[i], 1.0), "Error constructing a sigma matrix %d \n", pauliList[i]);
		CALL_AND_CHECK_STATUS(kronecker_sparse_c_csr(&(this -> P), (this->P), S), "Error during kronecker_sparse_c_csr at site %d \n", i);
	}
memory_free:
	if (mkl_sparse_destroy(S) != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(S) \n"); fflush(0); status = mkl_sparse_destroy(S);
	}
}

PauliOperator::~PauliOperator() {
}
sparse_matrix_t* PauliOperator::matrixLocation() {
	return &(this->P);
}

sparse_status_t PauliOperator::copyMatrix(sparse_matrix_t* dest) {
	
	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 

	// Create matrix descriptor
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	CALL_AND_CHECK_STATUS(mkl_sparse_copy(this->P, descr, dest), "Error during copying a pauli matrix\n");
memory_free:
	return status;
}