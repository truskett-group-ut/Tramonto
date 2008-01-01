/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

#include "dft_energy_wjdc.h"

/****************************************************************************/
double integrand_WJDC_freen(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     int iseg,unk_field;

     iseg = iunk-Phys2Unk_first[DENSITY];
     rho_i = x[iunk][inode_box];
     unk_field=Phys2Unk_first[WJDC_FIELD]+iseg;

     if (rho_i > 1.e-9){
        integrand = -log(x[unk_field][inode_box]) + 0.5*Nbonds_SegAll[iseg]-1.0;
        integrand *= rho_i;
     }
     return(integrand);
}
/****************************************************************************/
double integrand_WJDC_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     int iseg,icomp;

     iseg = iunk-Phys2Unk_first[DENSITY];
     icomp=Unk2Comp[iseg];
     rho_i = Rho_seg_b[icomp];

     if (rho_i > 1.e-9){
        integrand = -log(Field_WJDC_b[iseg]) + 0.5*Nbonds_SegAll[iseg]-1.0;
        integrand *= rho_i;
     }

     return(integrand);
}
/****************************************************************************/
