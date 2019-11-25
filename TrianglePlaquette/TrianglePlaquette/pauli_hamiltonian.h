#pragma once
#ifndef PAULIHAMILTONIAN_INCLUDED
#define PAULIHAMILTONIAN_INCLUDED
#include <mkl.h>
#include "pauli_operator.h"

class PauliHamiltonian
{
private:
	int nQubits; // number of qubits (spin sites)
	PauliOperator* pauliOpList; // list of pauli operators
	int pauliLength; // the number of the pauli terms

    sparse_matrix_t P = NULL; // store the matrix in the sparse form

public:
	PauliHamiltonian(int const nQubits, int const pauliLength, int ** const listOfPauliList, float* const coefs);
	~PauliHamiltonian(); // destructor 
	sparse_status_t copyMatrix(sparse_matrix_t* dest); // Copy the matrix form P to the destination. 
};

#endif