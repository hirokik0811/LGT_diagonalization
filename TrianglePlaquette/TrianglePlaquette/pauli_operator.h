#ifndef PAULIOPERATOR_INCLUDED
#define PAULIOPERATOR_INCLUDED
#pragma once
#include <stdlib.h>
#include <math.h>
#include <mkl.h>


class PauliOperator
{
private:
	float coef; // coefficient of the operator
	int nQubits; // number of qubits (spin sites)
	int* pauliList; // list of pauli matrices: 0: I, 1: X, 2: Y, 3: Z

	sparse_matrix_t P = NULL; // store the matrix in the sparse form

public:
	PauliOperator(); // Constructor. Compute the matrix form P of this pauli operator.
	PauliOperator(int const nQubits, int* const pauliList, float const alpha); // Constructor. Compute the matrix form P of this pauli operator.
	~PauliOperator(); // destructor 
	sparse_status_t copyMatrix(sparse_matrix_t *dest); // Copy the matrix form P to the destination. 

};

#endif