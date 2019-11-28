# LGT_diagonalization
Codes for diagonalization of a large scale spin Hamiltonians for lattice gauge theories. 

# What can it do?
Currently this returns the eigenvalues in a search range which can be specified by a user of the Hamiltonian of the triangle plaquette model in spin basis. 
![equation](https://quicklatex.com/cache3/a4/ql_050cffb910751f6935fa2d395463f1a4_l3.png)
You may also specify the coefficients of the Hamiltonian. 

# How to use
Compile TrianglePlaquette/TrianglePlaquette.c with -mkl. Checked it works properly with Intel C Compiler, but you may need to modify some codes for other compilers. 
