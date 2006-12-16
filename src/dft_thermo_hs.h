/* This file was automatically generated.  Do not edit! */
double dmu_drho_hs_PY(double *rho);
double dp_drho_hs_PY(double *rho);
void FMT3_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define FMT3       2
void FMT2_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
#define FMT2       1
void FMT1_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
#define FMT1       0
extern int Type_func;
extern int Nrho_bar_s;
extern int Ndim;
void dphi_drb_bulk(double *rhobar,double *dphi_drb);
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
extern double Rho_coex[2];
extern int Iliq_vap;
extern int Nwall;
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
extern int Lsteady_state;
#define NMER_MAX     40
extern double Rho_seg_b[NMER_MAX];
void rhobar_icomp(double rho,int icomp,double *rhobar);
extern int Unk2Comp[NMER_MAX];
extern int Nseg_tot;
extern double Dphi_Drhobar_RTF[10];
extern double Dphi_Drhobar_LBB[10];
extern double Rhobar_b_RTF[10];
extern double Rhobar_b_LBB[10];
extern int Nrho_bar;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void chempot_PY_hs(double *rho);
double pressure_PY_hs(double *rho);
extern double Betamu_hs_ex[NCOMP_MAX];
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
void chempot_FMT_hs(double *rho);
extern double Dphi_Drhobar_b[10];
double phispt(double *rho_bar);
extern double Rhobar_b[10];
#define NDIM_MAX  3
double pressure_FMT_hs(double *rho);
double uLJ12_6_cut(double r,double sigma,double eps,double rcut);
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define VERBOSE      3 
extern int Iwrite;
double integrand_BH(double r,int icomp);
#define LJ_BH_MF  1
extern int Type_attr;
void calc_HS_diams();
extern double Inv_4pirsq[NCOMP_MAX];
extern double Inv_4pir[NCOMP_MAX];
extern double HS_diam[NCOMP_MAX];
extern double Inv_rad[NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
#define PI    M_PI
extern double Inv_4pi;
void compute_bulk_FMT_properties(char *output_file1);
void calc_InvR_params();
extern double Fac_overlap_hs[NCOMP_MAX];
extern int Ncomp;
#define WTC          3
extern int Type_poly;
void HS_thermo_precalc(char *output_file1);