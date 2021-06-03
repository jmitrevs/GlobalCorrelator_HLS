#include <cstdio>
#include <iomanip>
#include <algorithm>
#include "firmware/simple_fullpfalgo.h"
//#include "vertexing/firmware/simple_vtx.h"
//#include "puppi/firmware/simple_puppi.h"
#include "utils/random_inputs.h"
#include "utils/DiscretePFInputs_IO.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"
//#include "bitonic-sort-48/hls/sorting_network_corr.hpp"

template<class T>
bool compare_hwPt(T i1, T i2) 
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

template <class T> 
void print3vecId(std::string name, size_t n, T col) {
  std::cout << name << ":" << std::endl;
  for (size_t i = 0; i < n; i++) {
    std::cout << "  " << i << ": pT = " << col[i].hwPt << ", eta = " 
	      << col[i].hwEta << ", phi = " << col[i].hwPhi << ", id = " << col[i].hwId << std::endl;
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
#endif
#if defined(TESTCTP7)
    CTP7PatternSerializer serInPatterns4( "ctp7_input_patterns_nomux.txt",CTP7_NCHANN_IN, true);  // 
    CTP7PatternSerializer serOutPatterns4("ctp7_output_patterns_nomux.txt",CTP7_NCHANN_OUT, false); // fill the rest of the lines with empty events for now
#endif
    HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
    //HumanReadablePatternSerializer serHR("-");
    HumanReadablePatternSerializer debugHR("-"); // this will print on stdout, we'll use it for errors

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

        // get the inputs from the input object
        if (!inputs.nextRegion(calo, emcalo, track, mu, hwZPV)) break;

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

        //PFChargedObj pfch_out_internal[NTRACK]; PFNeutralObj pfne_all_out_internal[NNEUTRALS]; PFChargedObj pfmu_out_internal[NMU];
        PFChargedObj pfch_out_internal[NTRACK]; PFNeutralObj pfpho_out_internal[NPHOTON]; PFNeutralObj pfne_out_internal[NCALO]; PFChargedObj pfmu_out_internal[NMU];

        //mp7wrapped_pfalgo3_only(data_in, pfch_out_internal, pfne_all_out_internal, pfmu_out_internal);
        mp7wrapped_pfalgo3_only(data_in, pfch_out_internal, pfpho_out_internal, pfne_out_internal, pfmu_out_internal);

	std::sort(pfne_out_internal, pfne_out_internal+NCALO, compare_hwPt<PFNeutralObj>);

        pfalgo3_full_ref(emcalo, calo, track, mu, outch_ref, outpho_ref, outne_ref, outmupf_ref);


	// // the reference puts neutrals separately; let's combine
	// PFNeutralObj outne_all_ref[NNEUTRALS];
	// for (int i = 0; i < NNEUTRALS; i++) {
	//   if (i < NPHOTON) {
	//     outne_all_ref[i] = outpho_ref[i];
	//   } else {
	//     outne_all_ref[i] = outne_ref[i - NPHOTON];
	//   }
	// }


	print3vecId("pfch_out", NTRACK, pfch_out_internal);
	print3vecId("pfpho_out", NPHOTON, pfne_out_internal);
	print3vecId("pfne_out", NCALO, pfne_out_internal);
	print3vecId("pfmu_out", NMU, pfmu_out_internal);

	print3vecId("outch_ref", NTRACK, outch_ref);
	print3vecId("outpho_ref", NPHOTON, outpho_ref);
	print3vecId("outne_ref", NSELCALO, outne_ref);
	print3vecId("outmupf_ref", NMU, outmupf_ref);

        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0; int ntot = 0;

        // check pf
        for (int i = 0; i < NTRACK; ++i) {
            if (!pf_equals(outch_ref[i], pfch_out_internal[i], "PFch", i)) errors++;
            if (outch_ref[i].hwPt > 0) { ntot++; }
        }
        // check pf
        for (int i = 0; i < NPHOTON; ++i) {
	    if (!pf_equals(outpho_ref[i], pfpho_out_internal[i], "PFpho", i)) errors++;
            if (outpho_ref[i].hwPt > 0) { ntot++; }
        }
        // check pf
        for (int i = 0; i < NSELCALO; ++i) {
	    if (!pf_equals(outne_ref[i], pfne_out_internal[i], "PFne", i)) errors++;
            if (outne_ref[i].hwPt > 0) { ntot++; }
        }
        // check pf
        for (int i = 0; i < NMU; ++i) {
            if (!pf_equals(outmupf_ref[i], pfmu_out_internal[i], "PFmu", i)) errors++;
            if (outmupf_ref[i].hwPt > 0) { ntot++; }
        }

        if (errors != 0) {
            printf("Error in pf test %d (%d)\n", test, errors);
            printf("Inputs: \n"); debugHR.dump_inputs(emcalo, calo, track, mu);
            return 1;
        } else {
            printf("Passed pf test %d (%d)\n", test, ntot);
        }

        std::cout << "end of test ---- " << test << std::endl;

    }
    return 0;
}
