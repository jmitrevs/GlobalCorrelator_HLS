#ifndef SIMPLE_PFALGO3_H
#define SIMPLE_PFALGO3_H

#include "data.h"
#include <hls_stream.h>

bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) ;
etaphi_t dr_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
template<int NB> ap_uint<NB>  dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) ;

void pfalgo3_full_ref(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
void pfalgo3_full_ref_set_debug(bool debug);
//void pfalgo3_full(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU], tk2calo_dr_t drvals_tk2em[NTRACK][NPHOTON], tk2calo_dr_t drvals_tk2calo[NTRACK][NSELCALO]) ;
//void pfalgo3_part1(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFChargedObj outmu[NMU], tk2em_dr_t drvals_tk2em[NTRACK][NPHOTON], HadCaloObj hadcalo_sub[NCALO], bool isMu[NTRACK], bool isEle[NTRACK]) ;
//void pfalgo3_part2(TkObj track[NTRACK], HadCaloObj hadcalo_sub[NCALO], bool isMu[NTRACK], bool isEle[NTRACK], PFChargedObj outch[NTRACK], PFNeutralObj outne[NSELCALO], tk2calo_dr_t drvals_tk2calo[NTRACK][NSELCALO]) ;
void pfalgo3_full(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
void pfalgo3_part1(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFChargedObj outmu[NMU], HadCaloObj hadcalo_sub[NCALO], bool isMu[NTRACK], bool isEle[NTRACK]) ;
void pfalgo3_part2(TkObj track[NTRACK], HadCaloObj hadcalo_sub[NCALO], bool isMu[NTRACK], bool isEle[NTRACK], PFChargedObj outch[NTRACK], PFNeutralObj outne[NSELCALO]) ;
void mp7wrapped_pack_in(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], MP7DataWord data[MP7_NCHANN]) ;
void mp7wrapped_unpack_in(MP7DataWord data[MP7_NCHANN], EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], TkObj track[NTRACK], MuObj mu[NMU]) ;
void mp7wrapped_pack_out(PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU], MP7DataWord data[MP7_NCHANN]) ;
void mp7wrapped_pack_out_necomb(PFChargedObj outch[NTRACK], PFNeutralObj outne_all[NNEUTRALS], PFChargedObj outmu[NMU], MP7DataWord data[MP7_NCHANN]) ;
void mp7wrapped_unpack_out(MP7DataWord data[MP7_NCHANN], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
void mp7wrapped_unpack_out_comb(MP7DataWord data[MP7_NCHANN], PFOutputObj outch[NALL]) ;
void mp7wrapped_unpack_out_necomb(MP7DataWord data[MP7_NCHANN], PFChargedObj outch[NTRACK], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
//void mp7wrapped_pfalgo3_full(MP7DataWord input[MP7_NCHANN], MP7DataWord output[2*NOUT_SORT], z0_t Z0) ;
void mp7wrapped_pfalgo3_full(MP7DataWord input[MP7_NCHANN], MP7DataWord output[2*NALL], z0_t Z0) ;
void mp7wrapped_pfalgo3_only(MP7DataWord input[MP7_NCHANN], PFChargedObj pfch_out[NTRACK], PFNeutralObj pfne_all_out[NNEUTRALS], PFChargedObj pfmu_out[NMU]) ;
void sort_output(PFOutputObj pf_comb[NALL], MP7DataWord output[2*NALL]) ;
void sort_output_onlycands(PFOutputObj pf_comb[NALL], PFOutputObj pf_sort[NOUT_SORT]) ;
//void sort_output_apxpack(PFOutputObj pf_comb[NALL], APxDataWord output[NOUT_SORT]) ;
void sort_output_apxpack(APxDataWord input[NALL], APxDataWord output[NOUT_SORT]) ;
template<int NOUT>
void mp7wrapped_pack_out_comb( PFOutputObj pfout[NOUT], MP7DataWord data[NALL*2]) ;
template<int NOUT>
void apxwrapped_pack_out_comb( PFOutputObj pfout[NOUT], APxDataWord data[NOUT]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfout complete
    // pack outputs
    for (unsigned int i = 0; i < NOUT; ++i) {
        data[i] = ( pfout[i].hwId, pfout[i].hwPhi, pfout[i].hwEta, pfout[i].hwZ0Pup, pfout[i].hwPt );
    }

}

template<int NOUT>
void apxwrapped_unpack_in_comb( APxDataWord data[NOUT], PFOutputObj pfout[NOUT]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfout complete
    // pack outputs
    for (unsigned int i = 0; i < NOUT; ++i) {
        pfout[i].hwPt     = data[i](15, 0);
        pfout[i].hwZ0Pup  = data[i](31, 16);
        pfout[i].hwEta    = data[i](41, 32);
        pfout[i].hwPhi    = data[i](51, 42);
        pfout[i].hwId     = data[i](54, 52);
    }

}

#endif

#ifndef DRVALSET
#define DRVALSET
//#define PFALGO3_DR2MAX_TK_CALO 756
#define PFALGO3_DR2MAX_TK_CALO 377
//#define PFALGO3_DR2MAX_EM_CALO 525
#define PFALGO3_DR2MAX_EM_CALO 262
//#define PFALGO3_DR2MAX_TK_MU   2101
#define PFALGO3_DR2MAX_TK_MU   1049
//#define PFALGO3_DR2MAX_TK_EM   84
#define PFALGO3_DR2MAX_TK_EM   42
//for demonstrator, altered encoding to handle large region. will need to modify in any case to handle actual inputs
#define PFALGO3_TK_MAXINVPT_LOOSE    40
#define PFALGO3_TK_MAXINVPT_TIGHT    80
//#define PFPUPPI_DR2MAX 8405
#define PFPUPPI_DR2MAX 4195
#endif
