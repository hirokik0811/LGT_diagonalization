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

PauliOperator::PauliOperator() {};

PauliOperator::PauliOperator(int const nQubits, int* const pauliList, float const alpha) {
	this->coef = alpha;
	this->nQubits = nQubits;
	this->pauliList = pauliList;

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	sparse_matrix_t S = NULL; // a sigma matrix to be inserted at each site

	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	CALL_AND_CHECK_STATUS(sigma_matrix(& (this->P), pauliList[0], alpha), "Error constructing a sigma matrix \n");
	for (int i = 1; i < nQubits; ++i) {
		CALL_AND_CHECK_STATUS(sigma_matrix(&S, pauliList[i], 1.0), "Error constructing a sigma matrix\n");
		CALL_AND_CHECK_STATUS(kronecker_sparse_c_csr(&(this->P), (this->P), S), "Error during kronecker_sparse_c_csr at site %d \n", i);
	}

memory_free:
	status = mkl_sparse_destroy(S);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(S) \n"); fflush(0);
	}
}

PauliOperator::~PauliOperator() {
}

sparse_status_t PauliOperator::copyMatrix(sparse_matrix_t* dest) {
	
	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 

	// Create matrix descriptor
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	MKL_Complex8* valuesP = NULL;
	MKL_INT* pointerB_P = NULL;
	MKL_INT* pointerE_P = NULL;
	MKL_INT* columns_P = NULL;
	MKL_INT n_rowsP, n_colsP;
	int ii = 0;

	CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr(this->P, &indexing, &n_rowsP, &n_colsP, &pointerB_P, &pointerE_P, &columns_P, &valuesP),
		"Errpr exporting");
	
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


	CALL_AND_CHECK_STATUS(mkl_sparse_copy(this->P, descr, dest), "Error during copying a pauli matrix\n ");
memory_free:
	return status;
}
