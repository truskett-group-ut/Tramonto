
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

#ifndef DFT_EQLOADERNEWPHYS_H
#define DFT_EQLOADERNEWPHYS_H

#include <dft_eqLoader.h>

//! dft_eqLoader: Parent class for all equation (physics) loading classes

class dft_eqLoader_NewPhys: public class dft_eqLoader  {
      
 public:

  //@{ \name Constructors.
  dft_eqLoader_NewPhys(int iBlock, int numBlocks, int block2Phys[NEQ_TYPE], int block2Nunk[NEQ_TYPE]);

  //@{ \name Destructor.
    //! Destructor
   ~dft_eqLoader_NewPhys();
  //@}

  //@}
  //@{ \name Public Method
    //! Top-level call to load the entire residual/matrix for this phsics
    //! This is implemented in parent class!
  //!void loadAll(double **x, int iter, int resid_only_flag);
  //@}
  
 private:
   void loadDensity();
   void loadRhoBarRosen();
   void loadPoisson();
   void loadDiffusion();
   void loadCavityWTC();
   void loadBondWTC();
   void loadCMSField();
   void loadCMSG();
   void loadNewPhys();

   void printBlockInfo(string phys);
};

#endif
