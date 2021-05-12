#include <cstdio>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include "firmware/simple_fullpfalgo.h"
#include "vertexing/firmware/simple_vtx.h"
#include "puppi/firmware/simple_puppi.h"
#include "utils/random_inputs.h"
#include "utils/DiscretePFInputs_IO.h"
#include "utils/pattern_serializer.h"
//#include "utils/test_utils.h"
//#include "bitonic-sort-48/hls/sorting_network_corr.hpp"

bool compare_hwPt(PFOutputObj i1, PFOutputObj i2) 
{ 
    return (i1.hwPt > i2.hwPt); 
} 

#define NTEST 10
#define NLINKS_APX_GEN0 48
#define NFRAMES_APX_GEN0 3

template <class T> 
void print3vec(std::string name, size_t n, T col) {
  std::cout << name << ":" << std::endl;
  for (size_t i = 0; i < n; i++) {
    std::cout << "  " << i << ": pT = " << col[i].hwPt << ", eta = " 
	      << col[i].hwEta << ", phi = " << col[i].hwPhi << std::endl;
  }
}

// template <class T> 
// void print3vecId(std::string name, size_t n, T col) {
//   std::cout << name << ":" << std::endl;
//   for (size_t i = 0; i < n; i++) {
//     std::cout << "  " << i << ": pT = " << col[i].hwPt << ", eta = " 
// 	      << col[i].hwEta << ", phi = " << col[i].hwPhi << ", id = " << col[i].hwId << std::endl;
//   }
// }

// template <> 
void print3vecId(std::string name, size_t n, PFNeutralObj col[]) {
  std::cout << name << ":" << std::endl;
  for (size_t i = 0; i < n; i++) {
    std::cout << "  " << i << ": pT = " << col[i].hwPt << ", eta = " 
	      << col[i].hwEta << ", phi = " << col[i].hwPhi << ", id = " << col[i].hwId 
	      << ", ptPuppy = " << col[i].hwPtPuppi
	      << std::endl;
  }
}
//template <> 
void print3vecId(std::string name, size_t n, PFChargedObj col[]) {
  std::cout << name << ":" << std::endl;
  for (size_t i = 0; i < n; i++) {
    std::cout << "  " << i << ": pT = " << col[i].hwPt << ", eta = " 
	      << col[i].hwEta << ", phi = " << col[i].hwPhi << ", id = " << col[i].hwId 
	      << ", hwZ0 = " << col[i].hwZ0
	      << std::endl;
  }
}


bool cmp(const PFOutputObj& a, const PFOutputObj& b) {
  return a.hwPt > b.hwPt;
}

void sort_ref_alg(APxDataWord input[NALL], APxDataWord output[NOUT_SORT]) {
  
    PFOutputObj pf_comb[NALL];
    apxwrapped_unpack_in_comb<NALL>(input, pf_comb);

    // actually the sort needs to be anti-stable
    std::stable_sort(std::begin(pf_comb), std::end(pf_comb), cmp);

    // make anti-stable
    auto it = std::begin(pf_comb);
    while (it < std::end(pf_comb)) {
      auto it2 = it + 1;
      while (it2 < std::end(pf_comb) && it2->hwPt == it->hwPt) {
	++it2;
      }
      std::reverse(it, it2);
      it = it2;
    }
    
    apxwrapped_pack_out_comb<NOUT_SORT>(pf_comb, output);
}

int main() {

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    DiscretePFInputs inputs("regions_TTbar_PU140.dump");
    
    // input TP objects
    HadCaloObj calo[NCALO]; EmCaloObj emcalo[NEMCALO]; TkObj track[NTRACK]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO], calo_subem_ref[NCALO]; 
    MuObj mu[NMU];

    std::cout << "Here" << std::endl;

    // output PF objects
    PFChargedObj outch[NTRACK];
    PFNeutralObj outpho[NPHOTON];
    PFNeutralObj outne[NSELCALO];
    PFChargedObj outmupf[NMU];
    PFOutputObj outpf[NALL];
#if defined(TESTMP7)
    //MP7PatternSerializer serInPatterns( "mp7_input_patterns.txt", HLS_pipeline_II,HLS_pipeline_II-1); // mux each event into HLS_pipeline_II frames
    //MP7PatternSerializer serOutPatterns("mp7_output_patterns.txt",HLS_pipeline_II,HLS_pipeline_II-1); // assume only one PF core running per chip,
    //MP7PatternSerializer serInPatterns2( "mp7_input_patterns_magic.txt", HLS_pipeline_II,-HLS_pipeline_II+1); // mux each event into HLS_pipeline_II frames
    //MP7PatternSerializer serOutPatterns2("mp7_output_patterns_magic.txt",HLS_pipeline_II,-HLS_pipeline_II+1); // assume only one PF core running per chip,
    MP7PatternSerializer serInPatterns3( "mp7_input_patterns_nomux.txt");  // 
    MP7PatternSerializer serOutPatterns3("mp7_output_patterns_nomux.txt"); // ,
    std::cout << "Here 2" << std::endl;

#endif
#if defined(TESTCTP7)
    CTP7PatternSerializer serInPatterns4( "ctp7_input_patterns_nomux.txt",CTP7_NCHANN_IN, true);  // 
    CTP7PatternSerializer serOutPatterns4("ctp7_output_patterns_nomux.txt",CTP7_NCHANN_OUT, false); // fill the rest of the lines with empty events for now
#endif
    HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
    //HumanReadablePatternSerializer serHR("-");
    HumanReadablePatternSerializer debugHR("-"); // this will print on stdout, we'll use it for errors

    std::cout << "Here 3" << std::endl;

    // -----------------------------------------
    // run multiple tests
    for (int test = 1; test <= NTEST; ++test) {

        // initialize TP objects
        for (int i = 0; i < NTRACK; ++i) {
            track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwZ0 = 0; 
        }
        for (int i = 0; i < NCALO; ++i) {
            calo[i].hwPt = 0; calo[i].hwEmPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0; calo[i].hwIsEM = 0; 
        }
        for (int i = 0; i < NEMCALO; ++i) {
            emcalo[i].hwPt = 0; emcalo[i].hwPtErr = 0;  emcalo[i].hwEta = 0; emcalo[i].hwPhi = 0;
        }
        for (int i = 0; i < NMU; ++i) {
            mu[i].hwPt = 0; mu[i].hwPtErr = 0; mu[i].hwEta = 0; mu[i].hwPhi = 0;
        }

	std::cout << "Here 4" << std::endl;
        // get the inputs from the input object
        if (!inputs.nextRegion(calo, emcalo, track, mu, hwZPV)) break;

	std::cout << "Here 5" << std::endl;
        VtxObj curvtx;    
        simple_vtx_ref(track,&curvtx);
        printf("Vertex Z   %i\n",(int)(curvtx.hwZ0));

	std::cout << "New region" << std::endl;
 
	// print3vec("tracks", NTRACK, track);
	// print3vec("calo", NCALO, calo);
	// print3vec("emcalo", NEMCALO, emcalo);
	// print3vec("mu", NMU, mu);

        // VtxObj curvtx;    
        // simple_vtx_ref(track,&curvtx);
        // printf("Vertex Z   %i\n",(int)(curvtx.hwZ0));

        MP7DataWord data_in[MP7_NCHANN], data_out[MP7_NCHANN];
        // initialize
        for (int i = 0; i < MP7_NCHANN; ++i) {
            data_in[i] = 0;
            data_out[i] = 0;
        }
        mp7wrapped_pack_in(emcalo, calo, track, mu, data_in);

        PFChargedObj pfch_out_internal[NTRACK]; PFNeutralObj pfne_all_out_internal[NNEUTRALS]; PFChargedObj pfmu_out_internal[NMU];

        mp7wrapped_pfalgo3_only(data_in, pfch_out_internal, pfne_all_out_internal, pfmu_out_internal);

	APxDataWord output_inter[NALL];
	simple_puppi_hw_apxoutput(pfch_out_internal, pfne_all_out_internal, pfmu_out_internal, z0_t(curvtx.hwZ0), output_inter);

	std::cout << "output_inter: " << std::endl;
	for (int i = 0; i < NALL; i++) {
	  std::cout << "Entry " << i << ",: " << std::hex << output_inter[i] << std::dec << std::endl;
	}

	APxDataWord sorted_out[NOUT_SORT];
	APxDataWord sorted_ref[NOUT_SORT];
	sort_output_apxpack(output_inter, sorted_out);
	sort_ref_alg(output_inter, sorted_ref);

        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0; int ntot = 0;

        // check res
        for (int i = 0; i < NOUT_SORT; ++i) {
	  std::cout << "Result " << i << ",: " << std::hex << sorted_ref[i] << ", " <<  sorted_out[i] << std::dec << std::endl;
	  if (sorted_ref[i] != sorted_out[i] && (sorted_ref[i] & 0xFF) !=  (sorted_out[i] & 0xFF)) {
	    std::cout << "Error in entry " << i << ",: " << std::hex << sorted_ref[i] << ", " <<  sorted_out[i] << std::dec << std::endl;
	    errors++;
	  }
        }

        if (errors != 0) {
            printf("Error in pf test %d (%d)\n", test, errors);
            printf("Inputs: \n"); debugHR.dump_inputs(emcalo, calo, track, mu);
            return 1;
        } else {
            printf("Passed sort test %d (%d)\n", test, NALL);
        }

        std::cout << "end of test ---- " << test << std::endl;

    }
    return 0;
}
