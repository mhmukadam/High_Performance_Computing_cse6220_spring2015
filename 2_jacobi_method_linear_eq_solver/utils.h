/**
 * @file    utils.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements common utility/helper functions.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

/*********************************************************************
 *             You can add new functions to this header.             *
 *********************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

/*********************************************************************
 * DO NOT CHANGE THE FUNCTION SIGNATURE OF THE FOLLOWING 3 FUNCTIONS *
 *********************************************************************/

inline int block_decompose(const int n, const int p, const int rank)
{
    return n / p + ((rank < n % p) ? 1 : 0);
}

inline int block_decompose(const int n, MPI_Comm comm)
{
    int rank, p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);
    return block_decompose(n, p, rank);
}

inline int block_decompose_by_dim(const int n, MPI_Comm comm, int dim)
{
    // get dimensions
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    return block_decompose(n, dims[dim], coords[dim]);
}


/*********************************************************************
 *                  DECLARE YOUR OWN FUNCTIONS HERE                  *
 *********************************************************************/

// ...

/**
 * @brief   Performs the matrix vector product: y = A*x
 *
 * @param n     The size of the dimensions.
 * @param A     A n-by-n matrix represented in row-major order.
 * @param x     The input vector of length n.
 * @param y     The output vector of length n.
 */
void diagonal(const int n, const double* A, double* D);

void nonDiagonal(const int n, const double* A, double* R);

void inverseDiagonal(const int n, const double* D, double* invD);

void vectorSub(const int n, const double* a, const double* b, double* res);

void vectorMult(const int n, const double* a, const double* b, double* res);

void init(const int n, double* x);

double norm(const int n, const double* x);

#endif // UTILS_H
