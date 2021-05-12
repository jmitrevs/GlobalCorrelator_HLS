#include "swap.hpp"

// replaces std::swap since that's not supported in synthesis
void swap(PFOutputObj &data1, PFOutputObj &data2) {
  auto temp = data2;
  data2 = data1;
  data1 = temp;
}

void swap1(PFOutputObj &data1, PFOutputObj &data2) {
	//if (data1.range(BIT_LOW,BIT_HIGH) < data2.range(BIT_LOW,BIT_HIGH)) {
	if (data1.hwPt < data2.hwPt) {
		swap(data1, data2);
	}
}

void swap2(PFOutputObj &data1, PFOutputObj &data2) {
	//if (data1.range(BIT_LOW,BIT_HIGH) > data2.range(BIT_LOW,BIT_HIGH)) {
	if (data1.hwPt > data2.hwPt) {
		swap(data1, data2);
	}
}
