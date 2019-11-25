#include <stdio.h>
#include "pauli_hamiltonian.h"

#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)
PauliHamiltonian::PauliHamiltonian(int const nQubits, int const PauliLength, int ** const listOfPauliList, float* const coefs) {

	this->nQubits = nQubits;
	this->pauliLength = pauliLength;
	this->pauliOpList = new PauliOperator[PauliLength];
	
	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	sparse_matrix_t Next = NULL; // a matrix to store the next pauli operator
	sparse_matrix_t Sum = NULL; // a matrix to store the sum

	int* pauliList = new int[nQubits];
	for (int j = 0; j < nQubits; ++j) {
		pauliList[j] = listOfPauliList[0][j];
	}
	this->pauliOpList[0] = PauliOperator(nQubits, pauliList, coefs[0]); // define the first pauli term
	CALL_AND_CHECK_STATUS(this->pauliOpList[0].copyMatrix(&Sum),
		"Error during copying a pauli operator matrix at 0th site\n");

	for (int i = 1; i < PauliLength; ++i) {
		for (int j = 0; j < nQubits; ++j) {
			pauliList[j] = listOfPauliList[i][j];
		}
		this->pauliOpList[i] = PauliOperator(nQubits, pauliList, coefs[i]);
		CALL_AND_CHECK_STATUS(this->pauliOpList[i].copyMatrix(&Next),
			"Error during copying a pauli operator matrix at %dth site\n", i);
		CALL_AND_CHECK_STATUS(mkl_sparse_c_add(SPARSE_OPERATION_NON_TRANSPOSE, Sum, { 1.0, 0.0 }, Next, &Sum),
			"Error during addition of two pauli terms at %d th site\n", i);
	}

	// Create matrix descriptor
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	CALL_AND_CHECK_STATUS(mkl_sparse_copy(Sum, descr, &(this->P)), "Error during copying a pauli matrix\n");


memory_free:
	delete[] pauliList;
	status = mkl_sparse_destroy(Sum);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(Sum) \n"); fflush(0);
	}
	status = mkl_sparse_destroy(Next);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(Next) \n"); fflush(0);
	}
}

PauliHamiltonian::~PauliHamiltonian() {
	delete[] this->pauliOpList;
}

sparse_status_t PauliHamiltonian::copyMatrix(sparse_matrix_t* dest) {

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 

	// Create matrix descriptor
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	CALL_AND_CHECK_STATUS(mkl_sparse_copy(this->P, descr, dest), "Error during copying a pauli hamiltonian\n");
memory_free:
	return status;
}