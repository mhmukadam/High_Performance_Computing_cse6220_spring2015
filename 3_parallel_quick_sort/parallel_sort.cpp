/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "parallel_sort.h"
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include "io.h"
#include "utils.h"
// implementation of your parallel sorting
struct local {
	int size;
	int * elm;

};

local my_parallel_sort(int * begin, int* end, MPI_Comm comm, int random) {
	int rank;
	MPI_Comm_rank(comm, &rank);
	int p;
	MPI_Comm_size(comm, &p);
	int local_m = end - begin;
	int m;
	MPI_Allreduce(&local_m, &m, 1, MPI_INT, MPI_SUM, comm);
//Recursion Termination Case when we have just one processor
	if (p == 1) {

		struct local lp;
		lp.size = local_m;
		std::sort(begin, end);
		lp.elm = begin;
		return lp;

	}
	int pivot;
	std::srand(p);

	//Random Number Generator with fixed seed
	srand(p + random);
	pivot = rand() % m;
	//Finding processor with pivot index
	int start_range = 0;
	int end_range = -1;
	int i = 0;
	int prev = 0;
	for (; i < p; i++) {
		start_range = start_range + prev;
		prev = block_decompose(m, p, i);
		end_range = end_range + prev;
		if (start_range <= pivot && end_range >= pivot) {
			break;
		}

	}

	int pivot_element;
	if (rank == i) {
		int iter = pivot - start_range;
		pivot_element = begin[0 + iter];

	}
	//Broadcast Pivot element
	MPI_Bcast(&pivot_element, 1, MPI_INT, i, comm);

	int counter = 0;
	int less = 0;

	//	Quicksort as per the pivot
	for (; counter < local_m; counter++) {
		if (begin[counter] <= pivot_element) {

			int temp = begin[counter];
			begin[counter] = begin[less];
			begin[less] = temp;
			less++;
		}
	}
	int greater = local_m - less;

	//All Reduce and all Gather to get total elements less than pivot
	int m_less;
	int* size_vector = new int[2];
	int* recv_vector = new int[2 * p];
	size_vector[0] = less;
	size_vector[1] = greater;
	MPI_Allreduce(&less, &m_less, 1, MPI_INT, MPI_SUM, comm);
	//   std::cerr <<"Less  "<< less << "at"<<rank<<std::endl;
	MPI_Allgather(&size_vector[0], 2, MPI_INT, &recv_vector[0], 2, MPI_INT,
			comm);

	//Assignment of processors based on the split of less-greater elements than pivot. Make sure no one gets 0 zero if elements are thee
	int m_greater = m - m_less;
	int less_processors = (m_less * p) / m;
	int greater_processors = (m_greater * p) / m;
	if (less_processors == 0 && m_less > 0) {
		less_processors = 1;
	} else if (greater_processors == 0 && m_greater > 0) {
		greater_processors = 1;
	} else {
		if (((m_less * p) % m == (m_greater * p) % m)
				&& (m_greater * p) % m == 0) {
			//Nothing
		} else if ((m_less * p) % m > (m_greater * p) % m) {
			less_processors++;
		} else {
			greater_processors++;
		}

	}
	if (p == less_processors || p == greater_processors) {
		MPI_Barrier(comm);
		return my_parallel_sort(begin, end, comm, random + 1);

	}
	//Exclusive Prefix Scan for number of elements less and greater before the processor
	int* prefix_vector = new int[2];
	MPI_Exscan(&size_vector[0], &prefix_vector[0], 2, MPI_INT, MPI_SUM, comm);
	MPI_Barrier(comm);
	if (rank == 0) {
		counter = 0;
		for (; counter < 2; counter++) {
			prefix_vector[counter] = 0;
		}
	}
	//Send Count Start
	int sendcnts[p];
	int sdispls[p];
	int prefix_min = prefix_vector[0];
	int left = recv_vector[2 * rank];
	start_range = 0;
	end_range = -1;
	prev = 0;

	for (i = 0; i < less_processors; i++) {

		start_range = start_range + prev;
		prev = block_decompose(m_less, less_processors, i);
		end_range = end_range + prev;
		if (start_range <= prefix_min && end_range >= prefix_min) {
			sdispls[i] = 0;
			if (left < (end_range - prefix_min + 1)) {
				sendcnts[i] = left;
				left = 0;
			} else {
				sendcnts[i] = end_range - prefix_min + 1;
				left = left - sendcnts[i];
			}
			i = i + 1;
			break;
		} else {
			sendcnts[i] = 0;
			sdispls[i] = 0;
		}
	}

	while (i < less_processors) {
		sdispls[i] = sdispls[i - 1] + sendcnts[i - 1];
		int r = left - block_decompose(m_less, less_processors, i);
		if (r <= 0) {
			sendcnts[i] = left;
			left = 0;
		} else {
			sendcnts[i] = block_decompose(m_less, less_processors, i);
			left = left - sendcnts[i];
		}
		i++;
	}

	left = recv_vector[2 * rank + 1];
	start_range = 0;
	end_range = -1;
	prev = 0;
	prefix_min = prefix_vector[1];

	for (; i < p; i++) {
		start_range = start_range + prev;
		int real_rank = i - less_processors;
		prev = block_decompose(m_greater, greater_processors, real_rank);

		end_range = end_range + prev;
		if (start_range <= prefix_min && end_range >= prefix_min) {

			sdispls[i] = sdispls[i - 1] + sendcnts[i - 1];
			if (left < (end_range - prefix_min + 1)) {
				sendcnts[i] = left;
				left = 0;

			} else {
				sendcnts[i] = end_range - prefix_min + 1;
				left = left - sendcnts[i];
			}
			i = i + 1;
			break;

		} else {
			sendcnts[i] = 0;
			sdispls[i] = sdispls[i - 1] + sendcnts[i - 1];
		}

	}

	while (i < p) {

		sdispls[i] = sdispls[i - 1] + sendcnts[i - 1];
		int real_rank = i - less_processors;
		int r = left
				- block_decompose(m_greater, greater_processors, real_rank);
		if (r <= 0) {
			sendcnts[i] = left;
			left = 0;
		} else {
			sendcnts[i] = block_decompose(m_greater, greater_processors,
					real_rank);
			left = left - sendcnts[i];
		}
		i++;
	}

	//Send Complete
	//Start Receive
	MPI_Barrier(comm);

	int* prefix_cap = new int[2];
	int* cap = new int[2];

	if (rank < less_processors) {
		cap[0] = block_decompose(m_less, less_processors, rank);
		cap[1] = 0;
		MPI_Exscan(&cap[0], &prefix_cap[0], 2, MPI_INT, MPI_SUM, comm);
	} else {
		cap[0] = 0;
		cap[1] = block_decompose(m_greater, greater_processors,
				rank - less_processors);
		MPI_Exscan(&cap[0], &prefix_cap[0], 2, MPI_INT, MPI_SUM, comm);
	}

	//	MPI_Barrier(comm);
	if (rank == 0) {
		counter = 0;
		for (; counter < 2; counter++) {
			prefix_cap[counter] = 0;
		}
	}

	int recvcnts[p];	   // = new int[p];
	int rdispls[p]; //= new int[p];
	int fcap;

	if (rank < less_processors) {
		prefix_min = prefix_cap[0];
		prev = 0;
		//    std::cerr <<"ReceiveCounts  "<< prefix_min<<"rank"<<rank<<std::endl;
		int capacity = block_decompose(m_less, less_processors, rank);
		fcap = capacity;
		//    std::cerr <<"Prefix  "<< prefix_cap[0] << " from rank"<<rank<<std::endl;
		for (i = 0; i < p; i++) {

			prev = recv_vector[2 * i];
			if (prev < prefix_min) {
				recvcnts[i] = 0;
				rdispls[i] = 0;
				prefix_min = prefix_min - recv_vector[2 * i];
			} else {
				rdispls[i] = 0;
				recv_vector[2 * i] = recv_vector[2 * i] - prefix_min;
				int remaining = capacity - recv_vector[2 * i];

				if (remaining <= 0) {
					recvcnts[i] = capacity;
					capacity = 0;
				} else {
					recvcnts[i] = recv_vector[2 * i];
					capacity = capacity - recvcnts[i];
				}
				i++;
				break;
			}
		}

		while (i < p) {
			int remaining = capacity - recv_vector[2 * i];
			rdispls[i] = rdispls[i - 1] + recvcnts[i - 1];

			if (remaining <= 0) {
				recvcnts[i] = capacity;
				capacity = 0;
			} else {
				recvcnts[i] = recv_vector[2 * i];
				capacity = capacity - recvcnts[i];
			}
			i++;
		}
	}

	if (rank >= less_processors) {
		prefix_min = prefix_cap[1];
		prev = 0;
		int capacity = block_decompose(m_greater, greater_processors,
				rank - less_processors);
		fcap = capacity;
		for (i = 0; i < p; i++) {

			prev = recv_vector[2 * i + 1];
			if (prev < prefix_min) {
				recvcnts[i] = 0;
				rdispls[i] = 0;
				prefix_min = prefix_min - recv_vector[2 * i + 1];
			} else {

				rdispls[i] = 0;
				recv_vector[2 * i + 1] = recv_vector[2 * i + 1] - prefix_min;
				int remaining = capacity - recv_vector[2 * i + 1];
				if (remaining <= 0) {
					recvcnts[i] = capacity;
					capacity = 0;
				} else {
					recvcnts[i] = recv_vector[2 * i + 1];
					capacity = capacity - recvcnts[i];
				}
				i++;
				break;
			}
		}

		while (i < p) {
			int remaining = capacity - recv_vector[2 * i + 1];
			rdispls[i] = rdispls[i - 1] + recvcnts[i - 1];
			if (remaining <= 0) {
				recvcnts[i] = capacity;
				capacity = 0;
			} else {
				recvcnts[i] = recv_vector[2 * i + 1];
				capacity = capacity - recvcnts[i];
			}
			i++;
		}
	}

	int* recv = new int[fcap];

	MPI_Barrier(comm);
	MPI_Alltoallv(&begin[0], &sendcnts[0], &sdispls[0], MPI_INT, &recv[0],
			&recvcnts[0], &rdispls[0], MPI_INT, comm);
	MPI_Barrier(comm);

	MPI_Comm low_comm;
	MPI_Comm high_comm;
	struct local lp;

	//Recursive function on lesser element and greater lements processors with diff communicators
	if (rank < less_processors) {
		MPI_Comm_split(comm, 0, rank, &low_comm);
		lp = my_parallel_sort(&recv[0], &recv[0] + fcap, low_comm, 0);
		MPI_Barrier(low_comm);

	} else {
		MPI_Comm_split(comm, 1, rank, &high_comm);
		lp = my_parallel_sort(&recv[0], &recv[0] + fcap, high_comm, 0);
		MPI_Barrier(high_comm);

	}
	// struct local lp;
	return lp;

}

void parallel_sort(int * begin, int* end, MPI_Comm comm) {
	int rank;
	MPI_Comm_rank(comm, &rank);
	double* local_start = NULL;
	double* local_end = NULL;
	//  struct local glp ;
	struct local glp = my_parallel_sort(begin, end, comm, 0);

	int local_m = end - begin;
	int m;
	MPI_Allreduce(&local_m, &m, 1, MPI_INT, MPI_SUM, comm);
	int p;
	MPI_Comm_size(comm, &p);

	int pre = 0;
	int pre_cap = 0;
	MPI_Exscan(&glp.size, &pre, 1, MPI_INT, MPI_SUM, comm);
	MPI_Exscan(&local_m, &pre_cap, 1, MPI_INT, MPI_SUM, comm);

	int sendcnts[p];
	int sdispls[p];
	int prefix_min = pre;
	int start_range = 0;
	int end_range = -1;
	int prev = 0;
	int left = glp.size;
	int i = 0;

	for (i = 0; i < p; i++) {
		start_range = start_range + prev;
		prev = block_decompose(m, p, i);
		end_range = end_range + prev;
		if (start_range <= prefix_min && end_range >= prefix_min) {
			sdispls[i] = 0;
			if (left < (end_range - prefix_min + 1)) {
				sendcnts[i] = left;
				left = 0;
			} else {
				sendcnts[i] = end_range - prefix_min + 1;
				left = left - sendcnts[i];
			}
			i = i + 1;
			break;
		} else {
			sendcnts[i] = 0;
			sdispls[i] = 0;
		}
	}
	while (i < p) {
		sdispls[i] = sdispls[i - 1] + sendcnts[i - 1];
		int r = left - block_decompose(m, p, i);
		if (r <= 0) {
			sendcnts[i] = left;
			left = 0;
		} else {
			sendcnts[i] = block_decompose(m, p, i);
			left = left - sendcnts[i];
		}
		i++;
	}

	//Receive
	int* recv_vector = new int[p];
	MPI_Barrier(comm);
	MPI_Allgather(&glp.size, 1, MPI_INT, &recv_vector[0], 1, MPI_INT, comm);

	int recvcnts[p];
	int rdispls[p];
	int fcap;
	prefix_min = pre_cap;
	prev = 0;

	int capacity = block_decompose(m, p, rank);
	fcap = capacity;

	for (i = 0; i < p; i++) {
		prev = recv_vector[i];
		if (prev < prefix_min) {
			recvcnts[i] = 0;
			rdispls[i] = 0;
			prefix_min = prefix_min - recv_vector[i];
		} else {
			rdispls[i] = 0;
			recv_vector[i] = recv_vector[i] - prefix_min;
			int remaining = capacity - recv_vector[i];
			if (remaining <= 0) {
				recvcnts[i] = capacity;
				capacity = 0;
			} else {
				recvcnts[i] = recv_vector[i];
				capacity = capacity - recvcnts[i];
			}
			i++;
			break;
		}
	}

	while (i < p) {
		int remaining = capacity - recv_vector[i];
		rdispls[i] = rdispls[i - 1] + recvcnts[i - 1];

		if (remaining <= 0) {
			recvcnts[i] = capacity;
			capacity = 0;
		} else {
			recvcnts[i] = recv_vector[i];
			capacity = capacity - recvcnts[i];
		}
		i++;
	}

	MPI_Alltoallv(&glp.elm[0], &sendcnts[0], &sdispls[0], MPI_INT, &begin[0],
			&recvcnts[0], &rdispls[0], MPI_INT, comm);
}
