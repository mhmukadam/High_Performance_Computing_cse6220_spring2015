// Implement your solutions in this file
#include "findmotifs.h"
#include "hamming.h"
#include <iostream>

unsigned int mySetBits(bits_t value)
{
    int dist = 0;
    while (value != 0)
	{
        dist++;
        value &= value - 1;
    }
    return dist;
}

unsigned int n_bit_position(bits_t flipper, unsigned int position)
{
    return (flipper >> (position))%2;
}


bits_t find_next(bits_t flipper, unsigned int position)
{
	int result = n_bit_position(flipper,position);
	bits_t test = 1;
    while (result!=0)
	{
        bits_t newflipper = (test<<position) - 1;
        position=position-1;
        flipper = flipper & newflipper;
        result=n_bit_position(flipper,position);
    }
    flipper = test<<position | flipper;
    return flipper;
}

// implements the sequential findmotifs function
std::vector<bits_t> findmotifs(unsigned int n, unsigned int l,
                               unsigned int d, const bits_t* input)
{
    // If you are not familiar with C++ (using std::vector):
    // For the output (return value) `result`:
    //                  The function asks you to return all values which are
    //                  of a hamming distance `d` from all input values. You
    //                  should return all these values in the return value
    //                  `result`, which is a std::vector.
    //                  For each valid value that you find (i.e., each output
    //                  value) you add it to the output by doing:
    //                      result.push_back(value);
    //                  Note: No other functionality of std::vector is needed.
    // You can get the size of a vector (flipper of elements) using:
    //                      result.size()

    // create an empty vector
    std::vector<bits_t> result;
    
    l--;
    bits_t flipper = 0;
	bits_t test = 1;
    bits_t end = (test<<d)-1;
    
    while(flipper != end)
	{
		if (n_bit_position(flipper,l) == 0)
            flipper = (test<<l) | flipper;
		else 
            flipper = find_next(flipper,l);
		if (mySetBits(flipper) <= d)
		{
			bits_t candidate = flipper^input[0];	
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
				result.push_back(candidate);
		}
	}
	
    return result;
}