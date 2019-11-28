# LGT_diagonalization
Codes for diagonalization of a large scale spin Hamiltonians for lattice gauge theories. 

# What can it do?
Currently this returns the eigenvalues in a search range which can be specified by a user of the Hamiltonian of the triangle plaquette model in spin basis. 
- <img src="https://latex.codecogs.com/gif.latex?O_t=\hat{H} = \frac{g^2}{2} " /> 
You may also specify the coefficients of the Hamiltonian. 

# How to use
Compile TrianglePlaquette/TrianglePlaquette.c with -mkl. Checked it works properly with Intel C Compiler, but you may need to modify some codes for other compilers. 
