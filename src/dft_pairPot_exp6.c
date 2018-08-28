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

//! Functional form of the exponential-6 potential
/*!
 * \param r Separation distance.
 * \param sigma Position of potential minimum.
 * \param epsilon Energy at potential minimum.
 * \param alpha Dimensionless shape parameter.
 * \param C van der Waals dispersion energy (for 1/r^6).
 *
 * \returns The potential evaluated at \a r.
 */
extern double uexp6(double r, double sigma, double eps, double alpha, double C)
    {
    return -eps*exp((sigma-r)/alpha) - C*POW_DOUBLE_INT(1./r,6);
    }

double uexp6_CS(double r,double sigma, double eps, double rcut, double K_Y, double eps_Y)
{
  double u,alpha,C;
  alpha=K_Y;
  C=eps_Y;

  if (r <= rcut) {
     u = uexp6(r, sigma, eps, alpha, C) - uexp6(rcut, sigma, eps, alpha, C);
  }
  else u = 0.0;
  return (u);
}
/*******************************************************************************/
/* uexp6_CS_setparams: The parameters for the cut and shifted exponential-6 potential */
void uexp6_CS_setparams(int context, int i, int j, double *param1,double *param2,double *param3,double *param4, double *param5)
{
  switch (context){
     case FLUID_FLUID:
        *param1 = Sigma_ff[i][j];
        *param2 = Eps_ff[i][j];
        *param3 = Cut_ff[i][j];
        *param4 = YukawaK_ff[i][j];
        *param5 = EpsYukawa_ff[i][j];
        break;
     case WALL_FLUID:
        *param1 = Sigma_wf[i][WallType[j]];
        *param2 = Eps_wf[i][WallType[j]];
        *param3 = Cut_wf[i][WallType[j]];
        *param4 = YukawaK_wf[i][WallType[j]];
        *param5 = EpsYukawa_wf[i][WallType[j]];
        break;
     case WALL_WALL:
        *param1 = Sigma_ww[WallType[i]][WallType[j]];
        *param2 = Eps_ww[WallType[i]][WallType[j]];
        *param3 = Cut_ww[WallType[i]][WallType[j]];
        *param4 = YukawaK_ww[WallType[i]][WallType[j]];
        *param5 = EpsYukawa_ww[WallType[i]][WallType[j]];
        break;
     default:
        if (Iwrite_screen != NONE) printf("problem with potential context uexp6_CS_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* uexp6_DERIV1D: The derivative of the exponential-6 potential in the x (or y or z) direction */

double uexp6_DERIV1D(double r,double x,double sigma, double eps, double rcut,double K_Y, double eps_Y)
{
  double uderiv,alpha,C;
  alpha=K_Y;
  C=eps_Y;

  if (r <= rcut) {
     uderiv = ((eps/alpha)*exp((sigma-r)/alpha) + C/POW_DOUBLE_INT(1./r,7))*(x/r);
  }
  else uderiv = 0.0;
  return (uderiv);
}
/******************************************************************************/
/* uEXP_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void uexp6_InnerCore(int i, int j,double *rCore_left, double *rCore_right, double *epsCore)
{
   *rCore_left=0.0;
   *rCore_right=Sigma_ff[i][j];   /* monotonic potential - UMIN & UZERO options not available */

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
  double uatt,r_min,rcut,sigma,alpha,eps,C;
  sigma=Sigma_ff[i][j];
  rcut=Cut_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];
  C=EpsYukawa_ff[i][j];

  if(r<sigma){
     if (Type_CoreATT_CONST==CORECONST_ZERO) uatt = 0.0;
     else                                    uatt = uexp6_CS(sigma,sigma,eps,rcut,alpha,C);
  }
  else if (r<=rcut){
     uatt = uexp6_CS(r,sigma,eps,rcut,alpha,C);
  }
  else uatt=0.0;
  return uatt;
}
/******************************************************************************/
/* uexp6_ATT_noCS: the attractive part of the potential for an exponential-6 potential */
double uexp6_ATT_noCS(double r,int i, int j)
{
  double uatt,r_min,sigma,alpha,eps,C;
  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];
  C=EpsYukawa_ff[i][j];

  if (r<sigma){
       if (Type_CoreATT_CONST==CORECONST_ZERO) uatt = 0.0;
       else                                    uatt = uexp6(sigma,sigma,eps,alpha,C);
  }
  else   uatt = uexp6(sigma,sigma,eps,alpha,C);

  return uatt;
}
/****************************************************************************/
/* uexp6_IntStencil:  the integral of the exponential-6 potential that is used
                        to define the magnitude of the DFTMFT UATTRACT
                        integration stencil. */

double uexp6_Integral(double r,int i, int j)
{
  double uatt_int,sigma,alpha,eps,C;
  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j];
  C=EpsYukawa_ff[i][j];

  uatt_int = 4.*PI*(eps*alpha*exp((sigma-r)/alpha)*(2.*alpha*alpha+2.*alpha*r+r*r) + C/(3.*r*r*r));
  return uatt_int;
}
/****************************************************************************/
