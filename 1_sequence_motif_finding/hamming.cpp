// Implement your solutions in this file
#include "hamming.h"

unsigned int setBits(uint64_t value)
{
    unsigned int dist = 0;
    while (value != 0)
	{
        dist++;
        value &= value - 1;
    }
    return dist;
}

unsigned int hamming(uint64_t x, uint64_t y)
{
    uint64_t val = x^y;
    return setBits(val);
}