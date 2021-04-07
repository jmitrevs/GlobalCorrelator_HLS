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
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFNeutralObj outpho[NPHOTON], outpho_ref[NPHOTON];
    PFNeutralObj outne[NSELCALO], outne_ref[NSELCALO];
    PFChargedObj outmupf[NMU], outmupf_ref[NMU];
    PFOutputObj outpf[NALL], outpf_ref[NALL];
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
 
	print3vec("tracks", NTRACK, track);
	print3vec("calo", NCALO, calo);
	print3vec("emcalo", NEMCALO, emcalo);
	print3vec("mu", NMU, mu);

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

	// since puppi changes the PF candidates, copy for reference
        PFChargedObj pfch_out_ref[NTRACK]; PFNeutralObj pfne_all_out_ref[NNEUTRALS]; PFChargedObj pfmu_out_ref[NMU];
	std::copy(std::begin(pfch_out_internal), std::end(pfch_out_internal), std::begin(pfch_out_ref));
	std::copy(std::begin(pfne_all_out_internal), std::end(pfne_all_out_internal), std::begin(pfne_all_out_ref));
	std::copy(std::begin(pfmu_out_internal), std::end(pfmu_out_internal), std::begin(pfmu_out_ref));

	APxDataWord output_inter[NALL];
	simple_puppi_hw_apxoutput(pfch_out_internal, pfne_all_out_internal, pfmu_out_internal, z0_t(curvtx.hwZ0), output_inter);

	simple_puppi_ref(pfch_out_ref, pfne_all_out_ref, z0_t(curvtx.hwZ0));

	// postprocess the reference
	PFOutputObj pf_comb_ref[NALL];
	for (unsigned int id = 0; id < NTRACK; id++) {
	  pf_comb_ref[id].hwPt = pfch_out_ref[id].hwPt;
	  pf_comb_ref[id].hwEta = pfch_out_ref[id].hwEta;
	  pf_comb_ref[id].hwPhi = pfch_out_ref[id].hwPhi;
	  pf_comb_ref[id].hwId = pfch_out_ref[id].hwId;
	  pf_comb_ref[id].hwZ0Pup = pfch_out_ref[id].hwZ0;
	}
	for (unsigned int id = 0; id < NPHOTON; id++) {
	  pf_comb_ref[id+NTRACK].hwPt = pfne_all_out_ref[id].hwPtPuppi;
	  pf_comb_ref[id+NTRACK].hwEta = pfne_all_out_ref[id].hwEta;
	  pf_comb_ref[id+NTRACK].hwPhi = pfne_all_out_ref[id].hwPhi;
	  pf_comb_ref[id+NTRACK].hwId = pfne_all_out_ref[id].hwId;
	  pf_comb_ref[id+NTRACK].hwZ0Pup = pfne_all_out_ref[id].hwPt;
	}
	for (unsigned int id = 0; id < NSELCALO; id++) {
	  pf_comb_ref[id+NTRACK+NPHOTON].hwPt = pfne_all_out_ref[id+NPHOTON].hwPtPuppi;
	  pf_comb_ref[id+NTRACK+NPHOTON].hwEta = pfne_all_out_ref[id+NPHOTON].hwEta;
	  pf_comb_ref[id+NTRACK+NPHOTON].hwPhi = pfne_all_out_ref[id+NPHOTON].hwPhi;
	  pf_comb_ref[id+NTRACK+NPHOTON].hwId = pfne_all_out_ref[id+NPHOTON].hwId;
	  pf_comb_ref[id+NTRACK+NPHOTON].hwZ0Pup = pfne_all_out_ref[id+NPHOTON].hwPt;
	}
	for (unsigned int id = 0; id < NMU; id++) {
	  pf_comb_ref[id+NTRACK+NPHOTON+NSELCALO].hwPt = pfmu_out_ref[id].hwPt;
	  pf_comb_ref[id+NTRACK+NPHOTON+NSELCALO].hwEta = pfmu_out_ref[id].hwEta;
	  pf_comb_ref[id+NTRACK+NPHOTON+NSELCALO].hwPhi = pfmu_out_ref[id].hwPhi;
	  pf_comb_ref[id+NTRACK+NPHOTON+NSELCALO].hwId = pfmu_out_ref[id].hwId;
	  pf_comb_ref[id+NTRACK+NPHOTON+NSELCALO].hwZ0Pup = pfmu_out_ref[id].hwZ0;
	}

	APxDataWord pf_comb_apx_ref[NALL];
	apxwrapped_pack_out_comb<NALL>(pf_comb_ref, pf_comb_apx_ref);

	print3vecId("pfch_out", NTRACK, pfch_out_internal);
	print3vecId("pfne_all_out", NNEUTRALS, pfne_all_out_internal);
	print3vecId("pfmu_out", NMU, pfmu_out_internal);

	print3vecId("pfch_ref", NTRACK, pfch_out_ref);
	print3vecId("pfne_all_ref", NNEUTRALS, pfne_all_out_ref);
	print3vecId("pfmu_ref", NMU, pfmu_out_ref);


        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0; int ntot = 0;

        // check res
        for (int i = 0; i < NALL; ++i) {
	  std::cout << "Result " << i << ",: " << std::hex << pf_comb_apx_ref[i] << ", " <<  output_inter[i] << std::dec << std::endl;
	  if (pf_comb_apx_ref[i] != output_inter[i]) {
	    std::cout << "Error in entry " << i << ",: " << std::hex << pf_comb_apx_ref[i] << ", " <<  output_inter[i] << std::dec << std::endl;
	    errors++;
	  }
        }

        if (errors != 0) {
            printf("Error in pf test %d (%d)\n", test, errors);
            printf("Inputs: \n"); debugHR.dump_inputs(emcalo, calo, track, mu);
            return 1;
        } else {
            printf("Passed pf test %d (%d)\n", test, NALL);
        }

        std::cout << "end of test ---- " << test << std::endl;

    }
    return 0;
}
