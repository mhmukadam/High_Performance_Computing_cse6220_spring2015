/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>

void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
	int rank, coords[2], blockSize;
	MPI_Comm_rank(comm,&rank);
	MPI_Cart_coords(comm,rank,2,coords);
	MPI_Comm col_comm;
	MPI_Comm_split(comm,coords[1],coords[0],&col_comm);
	
	MPI_Barrier(comm);
	if (coords[1] == 0)
	{
		int s;
		MPI_Comm_size(col_comm,&s);
		int sc[s], d[s];
		for (int i=0;i<s;i++)
		{
			sc[i] = block_decompose(n,s,i);
			if (i == 0)
				d[i] = 0;
			else
				d[i] = d[i-1]+sc[i-1];
		}
		
		blockSize = block_decompose(n,col_comm);
		*local_vector = new double[blockSize];
		MPI_Scatterv(&input_vector[0],sc,d,MPI_DOUBLE,*local_vector,blockSize,MPI_DOUBLE,0,col_comm);
	}
	
	MPI_Comm_free(&col_comm);
	MPI_Barrier(comm);
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
	int rank, coords[2], blockSize;
	MPI_Comm_rank(comm,&rank);
	MPI_Cart_coords(comm,rank,2,coords);
	MPI_Comm col_comm;
	MPI_Comm_split(comm,coords[1],coords[0],&col_comm);
	
	MPI_Barrier(comm);
	if (coords[1] == 0)
	{
		int s;
		MPI_Comm_size(col_comm,&s);
		int rc[s], d[s];
		for (int i=0;i<s;i++)
		{
			rc[i] = block_decompose(n,s,i);
			if (i == 0)
				d[i] = 0;
			else
				d[i] = d[i-1]+rc[i-1];
		}
		
		blockSize = block_decompose(n,col_comm);
		MPI_Gatherv(local_vector,blockSize,MPI_DOUBLE,output_vector,rc,d,MPI_DOUBLE,0,col_comm);
	}
	
	MPI_Comm_free(&col_comm);
	MPI_Barrier(comm);
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
	int rank, coords[2], blockRow, blockCol;
	MPI_Comm_rank(comm,&rank);
	MPI_Cart_coords(comm,rank,2,coords);
	MPI_Comm row_comm, col_comm;
	MPI_Comm_split(comm,coords[0],coords[1],&row_comm);
	MPI_Comm_split(comm,coords[1],coords[0],&col_comm);
	blockRow = block_decompose(n,row_comm);
	blockCol = block_decompose(n,col_comm);
	double *temp =  new double[n*blockCol];
	
	MPI_Barrier(comm);
	// first: scatter data at 0,0 accross the first column
	if (coords[1] == 0)
	{
		int s;
		MPI_Comm_size(col_comm,&s);
		int sc[s], d[s];
		for (int i=0;i<s;i++)
		{
			sc[i] = n*block_decompose(n,s,i);
			if (i == 0)
				d[i] = 0;
			else
				d[i] = d[i-1]+sc[i-1];
		}
		MPI_Scatterv(input_matrix,sc,d,MPI_DOUBLE,&temp[0],n*blockCol,MPI_DOUBLE,0,col_comm);
	}
	
	MPI_Barrier(comm);
	// second: scatter each temp in the first column accross their respective row
	int size;
	MPI_Comm_size(row_comm,&size);
	int sendCount[size], disp[size];
	for (int i=0;i<size;i++)
	{
		sendCount[i] = block_decompose(n,size,i);
		if (i == 0)
			disp[i] = 0;
		else
			disp[i] = disp[i-1]+sendCount[i-1];
	}
	*local_matrix = new double[blockRow*blockCol];
	for (int i=0;i<blockCol;i++)
		MPI_Scatterv(&temp[i*n],sendCount,disp,MPI_DOUBLE,&(*local_matrix)[i*blockRow],blockRow,MPI_DOUBLE,0,row_comm);
	
	delete[] temp;
	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&col_comm);
	MPI_Barrier(comm);
}

void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
	int rank, coords[2], blockRow, blockCol;
	MPI_Comm_rank(comm,&rank);
	MPI_Cart_coords(comm,rank,2,coords);
	int row = coords[0];
	int col = coords[1];
	MPI_Comm row_comm, col_comm;
	MPI_Comm_split(comm,row,col,&row_comm);
	MPI_Comm_split(comm,col,row,&col_comm);
	blockRow = block_decompose(n,row_comm);
	blockCol = block_decompose(n,col_comm);

	MPI_Barrier(comm);
	// send data from column zero to their respective diagonal locations in the same row
	if(col == 0)
	{
		int des, des_coords[2];
		des_coords[0]=row;
		des_coords[1]=row;
		MPI_Cart_rank(comm,des_coords,&des);
		MPI_Send(&col_vector[0],blockCol,MPI_DOUBLE,des,1,comm);
	}
	if(row == col)
    {
		int src, src_coords[2];
		src_coords[0]=row;
		src_coords[1]=0;
		MPI_Cart_rank(comm,src_coords,&src);
		MPI_Recv(&row_vector[0],blockCol,MPI_DOUBLE,src,1,comm,MPI_STATUS_IGNORE);
	}
	
	MPI_Barrier(comm);
	// Broadcast diagonal data to each respective column
	MPI_Bcast(&row_vector[0],blockRow,MPI_DOUBLE,col,col_comm);
	
	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&col_comm);
	MPI_Barrier(comm);
}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
	int rank, coords[2], blockRow, blockCol;
	MPI_Comm_rank(comm,&rank);
	MPI_Cart_coords(comm,rank,2,coords);
	MPI_Comm row_comm, col_comm;
	MPI_Comm_split(comm,coords[0],coords[1],&row_comm);
	MPI_Comm_split(comm,coords[1],coords[0],&col_comm);
	blockRow = block_decompose(n,row_comm);
	blockCol = block_decompose(n,col_comm);
	
	MPI_Barrier(comm);
	// distribute local_x on column zero 
	double *dist_x = new double[blockRow];
	transpose_bcast_vector(n,local_x,dist_x,comm);

	MPI_Barrier(comm);
	// compute local product pre_y from local_A and distribued x
	double *pre_y= new double[blockCol];
	for (int i=0;i<blockCol;i++)
	{
		pre_y[i] = 0;
		for (int j=0;j<blockRow;j++)
			pre_y[i] += (*(local_A + blockRow*i + j)) * (*(dist_x+j));
	}
	
	MPI_Barrier(comm);
	// reduce (sum) pre_y on all processors to local_y on column zero
	MPI_Reduce(&pre_y[0],&local_y[0],blockCol,MPI_DOUBLE,MPI_SUM,0,row_comm);
	
	delete[] dist_x;
	delete[] pre_y;
	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&col_comm);
	MPI_Barrier(comm);
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
	int rank, row, col, coords[2], blockRow, blockCol;
	MPI_Comm_rank(comm,&rank);
	MPI_Cart_coords(comm,rank,2,coords);
	MPI_Comm row_comm, col_comm;
	row = coords[0];
	col = coords[1];
	MPI_Comm_split(comm,row,col,&row_comm);
	MPI_Comm_split(comm,col,row,&col_comm);
	blockRow = block_decompose(n,row_comm);
	blockCol = block_decompose(n,col_comm);
	
	MPI_Barrier(comm);
	// collect D on first column
	double* D = new double[blockCol];
	double* invD = new double[blockCol];
	
	if(row == col)
    {
		diagonal(blockCol,local_A,&D[0]);
		int des, des_coords[2];
		des_coords[0]=row;
		des_coords[1]=0;
		MPI_Cart_rank(comm,des_coords,&des);
		MPI_Send(&D[0],blockCol,MPI_DOUBLE,des,1,comm);
	}
	if(col == 0)
	{
		int src, src_coords[2];
		src_coords[0]=row;
		src_coords[1]=row;
		MPI_Cart_rank(comm,src_coords,&src);
		MPI_Recv(&D[0],blockCol,MPI_DOUBLE,src,1,comm,MPI_STATUS_IGNORE);
		inverseDiagonal(blockCol,&D[0],&invD[0]);
	}
	
	MPI_Barrier(comm);
	// find R = A - D
	double* R = new double[blockRow*blockCol]; 
	if(row == col)
	{
		nonDiagonal(blockCol,local_A,&R[0]);
	}
	else
	{
		for (int i=0;i<blockRow*blockCol;i++)
			R[i] = *(local_A+i);
	}
	
	// initialize local_x
	init(blockCol,local_x);
	
	MPI_Barrier(comm);
	double  l2, temp;
	double* w = new double[blockCol];
	double* s = new double[blockCol];
	int rank00, rank00_coords[2] = {0,0};
	MPI_Cart_rank(comm,rank00_coords,&rank00);
	
	for (int iter=0;iter<max_iter;iter++)
	{
		distributed_matrix_vector_mult(n,&R[0],&local_x[0],&w[0],comm);
		if (col == 0)
		{
			vectorSub(blockCol,&local_b[0],&w[0],&s[0]);
			vectorMult(blockCol,&invD[0],&s[0],&local_x[0]);
		}
		MPI_Barrier(comm);
		distributed_matrix_vector_mult(n,&local_A[0],&local_x[0],&w[0],comm);
		if (col == 0)
		{
			vectorSub(blockCol,&w[0],&local_b[0],&s[0]);
		 	temp = 0;
		 	for(int i=0;i<blockCol;i++)
		 		temp += s[i]*s[i];
		 	MPI_Reduce(&temp,&l2,1,MPI_DOUBLE,MPI_SUM,0,col_comm);
		}
		if (row == 0 && col == 0)
			l2 = sqrt(l2);
		MPI_Barrier(comm);
		MPI_Bcast(&l2,1,MPI_DOUBLE,rank00,comm);
		if (l2 <= l2_termination)
		{
			delete[] D;
			delete[] invD;
			delete[] R;
			delete[] w;
			delete[] s;
			MPI_Comm_free(&row_comm);
			MPI_Comm_free(&col_comm);
		 	return;
		}
	}
}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
                            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);
	
    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);
	
    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}