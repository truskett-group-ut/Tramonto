/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_pairPot_Yukawa.c
 *
 *  This file contains routines specific to a strict mean field implementation of
 *  a Yukawa fluid .
 *
 */

#include "dft_pairPot_LJ12_6YukawaSum.h"

/******************************************************************************/
/* uLJandYUKAWA_CS: The cut and shifted LJ+YUKAWA potential                         */

double uLJandYUKAWA_CS(double r,double sigma, double eps, double rcut,double yukawaK,double Ayukawa)
{
  double u,alpha;
  alpha=yukawaK*sigma;

  /* note that writing the potential this way implies that the Yukawa parameters in the
     input file will be entered as K*sigma_ij */
  
  if (r <= rcut) {
     u = 4.0*fabs(eps)*(( POW_DOUBLE_INT(sigma/r,12) - POW_DOUBLE_INT(sigma/rcut,12) ) -
            ( POW_DOUBLE_INT(sigma/r,6) - POW_DOUBLE_INT(sigma/rcut,6) )) +
         Ayukawa*(exp(-alpha*(r/sigma-1.0))/(r/sigma)- exp(-alpha*(rcut/sigma-1.0))/(rcut/sigma));
  }
  else u = 0.0;
  return (u);
}
/*******************************************************************************/
/* uLJandYUKAWA_CS_setparams: The parameters for the cut and shifted LJ12_6 + Yukawa potential */
void uLJandYUKAWA_CS_setparams(int context, int i, int j, double *param1,double *param2,double *param3,double *param4,double *param5)
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
        if (Iwrite_screen != NONE) printf("problem with potential context uLJandYUKAWA_CS_setparams\n");
        exit(-1);
   }
   return;
}
/*******************************************************************************/
/* uLJandYUKAWA_DERIV1D: The derivative of the summed LJ and Yukawa potential in the x (or y or z) direction */

double uLJandYUKAWA_DERIV1D(double r,double x,double sigma, double eps, double rcut,double yukawaK,double Ayukawa)
{
  double uderiv,alpha;
  alpha=yukawaK*sigma;
  
  if (r <= rcut) {
     uderiv = (4.0*fabs(eps)/(sigma*sigma)) * (
            -12.*x*POW_DOUBLE_INT(sigma/r,14) + 6.*x*POW_DOUBLE_INT(sigma/r,8) )
            -Ayukawa*sigma*x*exp(-alpha*(r/sigma-1.0))*((1./r)+(alpha/(sigma*sigma)))/(r*r);
  }
  else uderiv = 0.0;
  return (uderiv);
}
/******************************************************************************/
/* uLJandYUKAWA_InnerCore : define the properties of the inner core of the potential based on
                  input parameters */
void uLJandYUKAWA_InnerCore(int i, int j,
       double *rCore_left, double *rCore_right, double *epsCore)
{
   /* note if the summed potential is purely repulsive then Rmin and Rzero will be set to Sigma_ff */
   switch(Type_CoreATT_R){
      case ATTCORE_SIGMA:      *(rCore_right)=Sigma_ff[i][j]; *(rCore_left)=0.0; break;
      case ATTCORE_UMIN:       *rCore_right=Rmin_ff[i][j]; *rCore_left=0.0; break;
      case ATTCORE_UCSZERO:    *rCore_right=Rzero_ff[i][j]; *rCore_left=0.0; break;
      case ATTCORE_SIGTOUMIN:  *rCore_right=Rmin_ff[i][j]; *rCore_left=Sigma_ff[i][j]; break;
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_R - set to %d\n",Type_CoreATT_R);
        exit(-1); break;
   }
   switch(Type_CoreATT_CONST){
      case CORECONST_UCONST:   *epsCore=uLJandYUKAWA_ATT_noCS(*rCore_right,i,j); break;
      case CORECONST_ZERO:     *epsCore=0.0; break;
      default:
        if (Iwrite_screen != NONE) printf("Problem with Type_CoreATT_CONST - set to %d\n",Type_CoreATT_CONST);
        exit(-1); break;
   }
   return;
}
/******************************************************************************/
/* uLJandYUKAWA_ATT_CS: the attractive part of the potential for a cut and shifted 
                        potential that sums LJ and yukawa potentials */
double uLJandYUKAWA_ATT_CS(double r,int i, int j)
{
  double uatt,r_min,rcut,sigma,alpha,eps,Ayukawa;
  double sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;
  double rc_inv,rc2_inv,rc6_inv,rc12_inv;

  sigma=Sigma_ff[i][j];
  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

  rc_inv   = 1.0/Cut_ff[i][j];
  rc2_inv  = rc_inv*rc_inv;
  rc6_inv  = rc2_inv*rc2_inv*rc2_inv;
  rc12_inv = rc6_inv*rc6_inv;

  rcut=Cut_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;
  Ayukawa=EpsYukawa_ff[i][j];

   /* note if the summed potential is purely repulsive then Rmin and Rzero will be set to Sigma_ff */
  switch(Type_CoreATT_R){
     case ATTCORE_SIGMA:      r_min=Sigma_ff[i][j]; break;
     case ATTCORE_SIGTOUMIN:       
     case ATTCORE_UMIN:       r_min=Rmin_ff[i][j]; break; 
     case ATTCORE_UCSZERO:    r_min=Rzero_ff[i][j]; break;
  }

  if ((r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) || 
      (r<sigma && Type_CoreATT_R==ATTCORE_SIGTOUMIN)){
         uatt=0.0;
  }
  else{
     if (r<=rcut){
        if (r<r_min) r=r_min;

        r_inv = 1.0/r;
        r2_inv  = r_inv*r_inv;
        r6_inv  = r2_inv*r2_inv*r2_inv;
        r12_inv = r6_inv*r6_inv;

        uatt=4.0 * fabs(eps)* sigma6 * (
                  sigma6*(r12_inv - rc12_inv)
                    - (r6_inv  - rc6_inv ) )+
             Ayukawa*exp(-alpha*(r/sigma-1.0))/(r/sigma)
           - Ayukawa*exp(-alpha*(rcut/sigma-1.0))/(rcut/sigma);
     }
     else uatt=0.0;
  }
  return uatt;
}
/******************************************************************************/
/* uLJandYUKAWA_ATT_noCS: the attractive part of the potential for a 
       combined LJ and yukawa potential */
double uLJandYUKAWA_ATT_noCS(double r,int i, int j)
{
  double uatt,sigma,alpha,r_min,eps,Ayukawa;
  double sigma2,sigma6;
  double r_inv,r2_inv,r6_inv,r12_inv;

  sigma=Sigma_ff[i][j];
  eps=Eps_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;
  Ayukawa=EpsYukawa_ff[i][j];
  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;

   /* note if the summed potential is purely repulsive then Rmin and Rzero will be set to Sigma_ff */
  switch(Type_CoreATT_R){
     case ATTCORE_SIGMA:      r_min=Sigma_ff[i][j]; break;
     case ATTCORE_SIGTOUMIN:       
     case ATTCORE_UMIN:       r_min=Rmin_ff[i][j]; break; 
     case ATTCORE_UCSZERO:    r_min=Rzero_ff[i][j]; break;
  }

  if ((r<r_min && Type_CoreATT_CONST==CORECONST_ZERO) || 
      (r<sigma && Type_CoreATT_R==ATTCORE_SIGTOUMIN)){
     uatt=0.0;
  }
  else{
     if (r<r_min) r=r_min;

     r_inv = 1.0/r;

     r2_inv  = r_inv*r_inv;
     r6_inv  = r2_inv*r2_inv*r2_inv;
     r12_inv = r6_inv*r6_inv;

     uatt= 4.0 * fabs(eps)* sigma6 * ( sigma6*r12_inv  - r6_inv)
           + Ayukawa*exp(-alpha*(r/sigma-1.0))/(r/sigma);
  }
  return uatt;
}
/****************************************************************************/
/* uLJandYUKAWA_IntStencil:  the integral of the summed LJ and Yukawa potential 
                        that is used to define the magnitude of the DFTMFT UATTRACT
                        integration stencil. */

double uLJandYUKAWA_Integral(double r,int i, int j)
{
  double uatt_int, sigma,c,alpha,eps,Ayukawa;
  double sigma2,sigma6;
  double r_inv,r3_inv,r9_inv;

  sigma=Sigma_ff[i][j];
  sigma2 = Sigma_ff[i][j]*Sigma_ff[i][j];
  sigma6 = sigma2*sigma2*sigma2;
  eps=Eps_ff[i][j];
  Ayukawa=EpsYukawa_ff[i][j];
  alpha=YukawaK_ff[i][j]*sigma;
  c=alpha/sigma;

  r_inv = 1.0/r;

  r3_inv  = r_inv*r_inv*r_inv;
  r9_inv  = r3_inv*r3_inv*r3_inv;


  uatt_int = 16*PI*fabs(eps)*sigma6 * ( - sigma6*r9_inv/9.0  + r3_inv/3.0 )+
             4*PI*Ayukawa*sigma*exp(alpha)*((exp(-c*r)/(c*c))*(-c*r-1.0));

  return uatt_int;
}
/****************************************************************************/

