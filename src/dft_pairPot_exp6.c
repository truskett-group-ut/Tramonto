/*
 *  FILE: dft_pairPot_exp6.c
 *
 *  This file contains routines specific to a strict mean field implementation of
 *  an exponential-6 attraction.
 *
 */

#include "dft_pairPot_exp6.h"

/******************************************************************************/
/* uexp6_CS: The cut and shifted exponential-6 potential                         */

//! Functional form of the expoenential-6 potential
/*!
 * \param r Separation distance.
 * \param sigma Position of potential minimum.
 * \param epsilon Energy at potential minimum.
 * \param alpha Dimensionless shape parameter.
 *
 * \returns The potential evaluated at \a r.
 */
extern double uexp6(double r, double sigma, double epsilon, double alpha)
    {
    return epsilon/(alpha-6.)*(6.*exp(alpha*(1.-r/sigma))-alpha*POW_DOUBLE_INT(sigma/r, 6));
    }

double uexp6_CS(double r,double sigma, double eps, double rcut,double yukawaK)
{
  double u,alpha;
  alpha=yukawaK;

  if (r <= rcut) {
     u = uexp6(r, sigma, eps, alpha) - uexp6(rcut, sigma, eps, alpha);
  }
  else u = 0.0;
  return (u);
}
/*******************************************************************************/
/* uexp6_CS_setparams: The parameters for the cut and shifted exponential-6 potential */
void uexp6_CS_setparams(int context, int i, int j, double *param1,double *param2,double *param3,double *param4)
{
  switch (context){
     case FLUID_FLUID:
        *param1 = Sigma_ff[i][j];
        *param2 = Eps_ff[i][j];
        *param3 = Cut_ff[i][j];
        *param4 = YukawaK_ff[i][j];
        break;
     case WALL_FLUID:
        *param1 = Sigma_wf[i][WallType[j]];
        *param2 = Eps_wf[i][WallType[j]];
        *param3 = Cut_wf[i][WallType[j]];
        *param4 = YukawaK_wf[i][WallType[j]];
        break;
     case WALL_WALL:
        *param1 = Sigma_ww[WallType[i]][WallType[j]];
        *param2 = Eps_ww[WallType[i]][WallType[j]];
        *param3 = Cut_ww[WallType[i]][WallType[j]];
        *param4 = YukawaK_ww[WallType[i]][WallType[j]];
        break;
     default:
        if (Iwrite_screen != NONE) printf("problem with potential context uexp6_CS_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* uexp6_DERIV1D: The derivative of the exponential-6 potential in the x (or y or z) direction */

double uexp6_DERIV1D(double r,double x,double sigma, double eps, double rcut,double yukawaK)
{
  double uderiv,alpha;
  alpha=yukawaK;
  
  if (r <= rcut) {
     uderiv = -(eps/sigma)*(6.*alpha/(alpha-6.))*(exp(alpha*(1.-r/sigma))-POW_DOUBLE_INT(sigma/r,7))*(x/r);
  }
  else uderiv = 0.0;
  return (uderiv);
}
/******************************************************************************/
/* uEXP_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void uexp6_InnerCore(int i, int j,double *rCore_left, double *rCore_right, double *epsCore)
{
   // rmin == sigma for this potential, so not all options are available
   switch(Type_CoreATT_R){
      case ATTCORE_SIGMA:
          *rCore_right=Sigma_ff[i][j];
          *rCore_left=0.0; break;
      case ATTCORE_UMIN:
          *rCore_right=Rmin_ff[i][j];
          *rCore_left=0.0; break;
      case ATTCORE_UCSZERO:
          *rCore_right=Rzero_ff[i][j];
          *rCore_left=0.0; break;
/* This case would cause weird equality issues, make it an error.
      case ATTCORE_SIGTOUMIN:
          *rCore_right=Rmin_ff[i][j];
          *rCore_left=Sigma_ff[i][j]; break;*/
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_R - set to %d\n",Type_CoreATT_R);
        exit(-1);
   }

   switch(Type_CoreATT_CONST){
      case CORECONST_UCONST:   *epsCore=uexp6_ATT_noCS(*rCore_right,i,j); break;
      case CORECONST_ZERO:     *epsCore=0.0; break;
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_CONST - set to %d\n",Type_CoreATT_CONST);
        exit(-1);
   }
   return;
}
/******************************************************************************/
/* uexp6_ATT_CS: the attractive part of the potential for a cut and shifted 
                   exponential-6 system */
double uexp6_ATT_CS(double r,int i, int j)
{
  double uatt,r_min,rcut,sigma,alpha,eps;
  sigma=Sigma_ff[i][j];
  rcut=Cut_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];

  switch(Type_CoreATT_R){
     case ATTCORE_SIGMA:      r_min=Sigma_ff[i][j]; break;
     case ATTCORE_UMIN:       r_min=Rmin_ff[i][j]; break;
     case ATTCORE_UCSZERO:    r_min=Rzero_ff[i][j]; break;
     default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_R - set to %d\n",Type_CoreATT_R);
        exit(-1);
  }

  if ((r < r_min && Type_CoreATT_CONST == CORECONST_ZERO)) uatt = 0.0;
  else
    {
    if (r < r_min) r = r_min;
    uatt = uexp6_CS(r, sigma, eps, rcut, alpha);
    }

  return uatt;
}
/******************************************************************************/
/* uexp6_ATT_noCS: the attractive part of the potential for an exponential-6 potential */
double uexp6_ATT_noCS(double r,int i, int j)
{
  double uatt,r_min,sigma,alpha,eps;
  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];

  switch(Type_CoreATT_R){
     case ATTCORE_SIGMA:      r_min=Sigma_ff[i][j]; break;
     case ATTCORE_UMIN:       r_min=Rmin_ff[i][j]; break;
     case ATTCORE_UCSZERO:    r_min=Rzero_ff[i][j]; break;
     default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_R - set to %d\n",Type_CoreATT_R);
        exit(-1);
  }

  if ((r < r_min && Type_CoreATT_CONST == CORECONST_ZERO)) uatt = 0.0;
  else
    {
    if (r < r_min) r = r_min;
    uatt = uexp6(r, sigma, eps, alpha);
    }

  return uatt;
}
/****************************************************************************/
/* uexp6_IntStencil:  the integral of the exponential-6 potential that is used
                        to define the magnitude of the DFTMFT UATTRACT
                        integration stencil. */

double uexp6_Integral(double r,int i, int j)
{
  double uatt_int, sigma,sigma2,sigma3,alpha,a2,a3,eps;

  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];
  sigma2 = sigma*sigma;
  sigma3 = sigma2*sigma;
  a2 = alpha*alpha;
  a3 = a2*alpha;

  uatt_int = -4.*PI*(eps/(alpha-6.))*((6.*sigma/a3)*exp(alpha*(1.-r/sigma))*(a2*r*r+2.*alpha*sigma*r+2*sigma2)-(alpha*sigma3/3.)*POW_DOUBLE_INT(sigma/r,3));

  return uatt_int;
}
/****************************************************************************/
