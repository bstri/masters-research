----Name Considerations----
Matrix Interpretations and Tools for Investigating Even Functionals

----Future work----
Some more complicated matrix expressions appear more than once in a set of expressions used for counting graphs.
Only the basic expressions A# and A_# are precomputed, but some time could be saved if expressions like A2*A2 are used multiple times.

Test performance of a parallel matrix multiplication algorithm using gmp and mpi. 

----Libraries----
Python script for parsing matrix expressions and turning them into C++ code using GMP.

Solving systems of equations in C++ using GMP.
Matrix class using mpq_t, mpz_t, or float types.
    Supports reading matrices from files in a variety of formats
    static methods for creating a variety of matrices
    row echelon form
    determinant
    hadamard product
    matrix multiplication
Faddeev-Leverrier method for finding the characteristic equation (and the determinant, inverse) of a matrix

Python notebook for
    drawing graphs with graphviz
    checking isomorphism with networkx
    creating common graph types with networkx
    converting to/from numpy adjacency matrix, networkx, and graphviz
    parsing hollow symmetric strings into graphs
    solving systems with sympy and its Matrix class and rref method
    brute force counting subgraphs
    faddeev-leverrier

Matlab scripts

----Findings----
determinant is not a fill-in for a matrix expression for K4

There may be multiple ways to write a matrix expression for a given (even) graph, some of which require less computation.
These other ways may not always be "trace equivalent". E.g. the 4-edge even graph can be written as A_2@A_2 or A_3@A depending on how it is converted to a series-parallel graph. These expressions are only equivalent for matrices of 0's and 1's. 

Using increasing sizes of complete graphs to create a system of equations for finding the coefficients to count a subgraph doesn't work for subgraphs with more than 5 edges. It results in linear dependence, meaning the system doesn't have enough unique information. This is because complete graphs contain ALL even subgraphs. Some even subgraphs are subgraphs themselves of their fellow basis graphs, and therefore it is necessary to use some host graphs where one but not the other occurs. 
This gives rise to a method of creating a system by taking the graph to be counted and attaching it to each even subgraph to create each row of the system. The number of occurrences of the graph to be counted is always 1 except when it itself is an even subgraph, in which case when it is combined with itself the number of occurrences will of course be 2.

Thought - If a basis is missing for even graphs of size n, there exists a graph with n edges that the set of coefficients and matrix expressions will fail to correctly count a subgraph of size n in.

Thought - If a given set of coefficients and matrix expressions for counting a graph with n edges is incorrect, it will fail for some (possibly weighted) subgraph of Kn. 

Parallel processing with GMP matrices (using MPI)
    -- must copy matrices to all processes
    -- convert gmp elements to string, send them, on the receiving process convert back