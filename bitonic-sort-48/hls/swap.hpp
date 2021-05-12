#ifndef __SWAP_HPP__
#define __SWAP_HPP__

#include "../../firmware/data.h"

void swap1(PFOutputObj &data1, PFOutputObj &data2);
void swap2(PFOutputObj &data1, PFOutputObj &data2);

// this is instead of std::swap, since that's unsupported
void swap(PFOutputObj &data1, PFOutputObj &data2);


#endif
