// Implement your solutions in this file
#include "findmotifs.h"
#include "hamming.h"
#include <mpi.h>

unsigned int mpi_n_bit_position(bits_t flipper, unsigned int position)
{
    return (flipper >> (position))%2;
}

bits_t getNextFlipper(bits_t flipper, unsigned int position, int isMaster)
{
	bits_t test = 1;
	if ((flipper == (test<<(position+1))-1) && isMaster)
		flipper = -1;
	else if (mpi_n_bit_position(flipper,position) == 0)
		flipper = (test<<position) | flipper;
	else
	{
		int result = mpi_n_bit_position(flipper,position);
		while (result!=0)
		{
			bits_t newflipper = (test<<position) - 1;
			position = position-1;
			flipper = flipper & newflipper;
			result = mpi_n_bit_position(flipper,position);
		}
		flipper = test<<position | flipper;
	}
    return flipper;
}

std::vector<bits_t> findmotifs_worker(const unsigned int n,
									  const unsigned int l,
									  const unsigned int d,
									  const bits_t* input,
									  const unsigned int startbitpos,
									  bits_t start_value)
{
    std::vector<bits_t> results;
	bits_t next_value = start_value;
	bits_t test = 1;
	bits_t end = (test<<d)-1;
	bits_t mask = (test<<startbitpos)-1; 

	while ((next_value != end) && (((start_value^next_value)&mask) == 0))
	{
		if (hamming(next_value,0) <= d)
		{
			bits_t candidate = next_value^input[0];	
			int flag = 1;
			for (int i=1; i<n; i++)
			{
				if (hamming(candidate,input[i]) > d)
				{
					flag = 0;
					break;
				}
			}
			if (flag == 1) 
			{
				
				results.push_back(candidate);
			}
		}
	
		next_value = getNextFlipper(next_value,l-1,0);
	}
    return results;
}

void worker_main()
{
	// get communicator size and my rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
	
	// receive input from master (including n, l, d, input, master-depth)
	unsigned int start_data[4];
	MPI_Bcast(&start_data,4,MPI_UNSIGNED,0,comm);
	unsigned int n = start_data[0];
	unsigned int l = start_data[1];
	unsigned int d = start_data[2];
	unsigned int master_depth = start_data[3];
	bits_t input[n];
	MPI_Bcast(&input,n,MPI_UINT64_T,0,comm);
	
	// declare variables for communication with master
	int terminate = 0;
	
	bits_t flipper;
	MPI_Request reqR, reqS;
	MPI_Status stat;
	int flag = 0;
	
	// work until master says to terminate
	while (!terminate)
	{
		MPI_Iprobe(0,MPI_ANY_TAG,comm,&flag,&stat);
		if (flag) // if new work has arrived
		{
			if (stat.MPI_TAG) // terminate if tag = 1 received from master
				terminate = 1; 
			else // receive new work, do the work and send results back
			{
				std::vector<bits_t> results;
				MPI_Irecv(&flipper,1,MPI_UINT64_T,0,MPI_ANY_TAG,comm,&reqR);
				results = findmotifs_worker(n,l,d,&input[0],master_depth,flipper);
				
				int i = 0;
				int count = 400;
				while (i<=results.size())
				{
					if (i+count > results.size())
					{
						count  = results.size() - i;
						MPI_Isend(&results[i],count,MPI_UINT64_T,0,1,comm,&reqS);
						break;
					}
					else
					{
						MPI_Isend(&results[i],count,MPI_UINT64_T,0,0,comm,&reqS);
					}
					i = i+ count;
				}
								
			}
		}
	}
	
}

std::vector<bits_t> findmotifs_master(const unsigned int n,
                                      const unsigned int l,
                                      const unsigned int d,
                                      const bits_t* input,
                                      const unsigned int till_depth)
{
    std::vector<bits_t> results;


    return results;
}

std::vector<bits_t> master_main(unsigned int n, unsigned int l, unsigned int d,
                                const bits_t* input, unsigned int master_depth)
{
	// get communicator size and my rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
	
	// broadcast initial data to workers (including n, l, d, input, master-depth) 
	unsigned int start_data[4] = {n,l,d,master_depth};
	MPI_Bcast(&start_data,4,MPI_UNSIGNED,0,comm);
	
	bits_t input_new[n];
	for (int i=0; i<n; i++) input_new[i] = input[i];
	MPI_Bcast(&input_new,n,MPI_UINT64_T,0,comm);
	
	// solve upto master_depth and send to workers
	std::vector<bits_t> results;
	
	// declare variables for sending work to workers
	bits_t start = 0;
	MPI_Request reqS, reqR;
	MPI_Status stat;
	int flag;
	int done = 0;
	
	// initially send two rounds of work to each worker
	for (int worker=1; worker<2*p; worker++)
	{
		if (worker%p == 0) continue;
		MPI_Isend(&start,1,MPI_UINT64_T,worker%p,0,comm,&reqS);
		start = getNextFlipper(start,master_depth-1,1);
		done++;
		if (start == -1) break;
	}
	
	// send next work when solution to old work is received, until all work is done
	while (start != -1 || done > 0)
	{
		for (int worker=1; worker<2*p; worker++)
		{
			if (worker%p == 0) continue;
			MPI_Iprobe(worker%p,MPI_ANY_TAG,comm,&flag,&stat);
			if (flag)
			{
				int count = 0;
				MPI_Get_count(&stat,MPI_UINT64_T,&count);
				std::vector<bits_t> result(count);				
				MPI_Irecv(&result[0],count,MPI_UINT64_T,worker%p,MPI_ANY_TAG,comm,&reqR);
				for (int i=0; i<count; i++)
				{
					results.push_back(result[i]);
				}
				
				if (start != -1  && stat.MPI_TAG == 1)
				{
					MPI_Isend(&start,1,MPI_UINT64_T,worker%p,0,comm,&reqS);
					start = getNextFlipper(start,master_depth-1,1);
				}
				else
				{
					if(stat.MPI_TAG == 1)
						done--;
				}
			}
		}
	}
	
	// all work done - send terminate signal to all workers
	for (int worker=1; worker<p; worker++)
	{
		MPI_Isend(&start,1,MPI_UINT64_T,worker,1,comm,&reqS);
	}

    return results;
}

