/**
 * @file    utils.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements common utility/helper functions.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "utils.h"
#include <iostream>
#include <cmath>

/*********************************************************************
 *                 Implement your own functions here                 *
 *********************************************************************/

// ...

void diagonal(const int n, const double* A, double* D)
{
	for(int i=0;i<n;i++)
		D[i]=A[i*n+i];
}

void nonDiagonal(const int n, const double* A, double* R)
{
	for (int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(i==j)
				R[i*n+j]= 0;
			else
				R[i*n+j]= A[i*n+j];
		}
	}
}

void inverseDiagonal(const int n, const double* D, double* invD)
{
	for(int i=0;i<n;i++)
		invD[i]=1/D[i];
}

void vectorSub(const int n, const double* a, const double* b, double* res)
{
	for(int i=0;i<n;i++)
		res[i]=a[i]-b[i];
}

void vectorMult(const int n, const double* a, const double* b, double* res)
{
	for(int i =0;i<n;i++)
		res[i]=a[i]*b[i];
}

void init(const int n, double* x)
{
	for(int i=0;i<n;i++)
		x[i]=0;
}

double norm(const int n, const double* x)
{
	double sum=0;
	for(int i=0;i<n;i++)
		sum += x[i]*x[i];
	return sqrt(sum);
}
