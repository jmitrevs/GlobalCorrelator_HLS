#include "../firmware/data.h"
#include "../firmware/simple_fullpfalgo.h"
#include "firmware/simple_puppi.h"
#include <cmath>
#include <algorithm>
#include <iostream>

void simple_puppi_ref(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0) {

    z0_t DZMAX = 256; // placeholder, big window (z0 is 10-bits)

    // input to the table is a 12-bit number, 2^12 = 4096
    ap_uint<8> puppiweight_table[PUPPI_TABLE_SIZE];
    lut_puppiweight_init< ap_uint<8>, PUPPI_TABLE_SIZE >( puppiweight_table );

    // compute alpha
    const int DR2MAX = PFPUPPI_DR2MAX; // 0.4 cone
    for (int in = 0; in < NNEUTRALS; ++in) {
      std::cout << "Running over " << in << ": pt = " << pfallne[in].hwPt << std::endl;
      //if (pfallne[in].hwPt == 0) continue;
        int sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if ((Z0 - pfch[it].hwZ0 > DZMAX) || (Z0 - pfch[it].hwZ0 < -DZMAX)) continue; // if track is PV
            int dr2 = dr2_int(pfch[it].hwEta, pfch[it].hwPhi, pfallne[in].hwEta, pfallne[in].hwPhi); // if dr is inside puppi cone
            if (dr2 <= DR2MAX) {
	      std::cout << "(ref) Looking at ch " << it << " with pt = " << pfch[it].hwPt 
			<< ", dr2 = " << dr2 << ", DR2MAX = "  << DR2MAX << std::endl;
                ap_uint<9> dr2short = dr2 >> 5; // why?
                int pt2 = pfch[it].hwPt*pfch[it].hwPt;

		std::cout << "pt2 = " << pt2 << ", pt2 >> 5 = " << (pt2 >> 5) << ", dr2short = " << dr2short;
                //int term = (std::min(pt2 >> 5, 131071) << 15)/(dr2short > 0 ? dr2short : ap_uint<9>(1));
                int term = 32768/(dr2short > 0 ? dr2short : ap_uint<9>(1)).to_int();
                term = std::min(pt2 >> 5, 131071)*term;
                sum += term;
		std::cout << ", sum = " << sum << std::endl;
            }
        }
      
        int weight = 0;
        int eToAlpha = sum >> 10;
	std::cout << "eToAlpha = " << eToAlpha << std::endl;
        if (eToAlpha > 0) { // get the weight if e^alpha is not 0
            if (eToAlpha <= 0) weight = 0; // < e^10 where that is the median
            else if (eToAlpha >= 1174) weight = 256; // e^14 where that is the median
            else{
                // bitshift down the sum by 10 bits to get it in the right range
                int sum_short = eToAlpha;
                // -- compute weight --
                // sum is the current e^alpha
                // med of alpha = 10, rms of alpha = 2
                // pass in sum, alphamed, alpharms into a LUT
                int index = sum_short >= PUPPI_TABLE_SIZE ? PUPPI_TABLE_SIZE : sum_short; // 2^16
                // ap_uint<8> weight_short = puppiweight_table[index];
                weight = (int) puppiweight_table[index];
            }
        }

	std::cout << "weight = " << weight << std::endl;

        int pT_new = pfallne[in].hwPt * weight;
        pT_new =  pT_new >> 8;
        pfallne[in].hwPtPuppi = (pt_t) pT_new;
        //std::cout<<in<<" "<<weight<<" "<<pT_new<<" "<<pfallne[in].hwPtPuppi<<std::endl;
    }

}

