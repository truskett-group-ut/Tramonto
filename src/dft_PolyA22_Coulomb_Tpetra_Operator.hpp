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

#ifndef DFT_POLYA22_COULOMB_TPETRA_OPERATOR_H
#define DFT_POLYA22_COULOMB_TPETRA_OPERATOR_H

#include "Tpetra_Headers.hpp"
#include "dft_PolyA22_Tpetra_Operator.hpp"
#include "dft_direct_solver_const.h"

//! dft_PolyA22_Tpetra_Operator: An implementation of the Operator class for Tramonto Schur complements with Coulomb effects.
/*! Special 2*numBeads by 2*numBeads (plus Coulomb) for Tramonto polymer problems.
*/

template<class Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
	 class Node = Kokkos::DefaultNode::DefaultNodeType>
class dft_PolyA22_Coulomb_Tpetra_Operator:
  public virtual dft_PolyA22_Tpetra_Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{

 public:
TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node);
  typedef dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> P22TO;

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads (plus Coulomb) system

  dft_PolyA22_Coulomb_Tpetra_Operator
  (const RCP<const MAP> & cmsMap,const RCP<const MAP> & densityMap,
   const RCP<const MAP> & poissonMap, const RCP<const MAP> & cmsDensMap,
   const RCP<const MAP> & block2Map, RCP<ParameterList> parameterList);

  //@}
  //@{ \name Assembly methods.
  void
  initializeProblemValues
  ();

  void
  insertMatrixValue
  (GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value, GlobalOrdinal blockColFlag);

  void
  finalizeProblemValues
  ();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_PolyA22_Coulomb_Tpetra_Operator();
  //@}

  //@{ \name Atribute get methods.

  //! Returns an Operator pointer that is actually the \e this object, since this class implements Operator.
  RCP<OP>
  getA22Inv
  ()
  {
    return(rcp(this));
  }
  //@}

  //@{ \name Mathematical functions.

    //! Returns the result of a dft_PolyA22_Coulomb_Tpetra_Operator applied to a MultiVector X in Y.
    /*!
    \param In
	   X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A MultiVector of dimension NumVectors containing result.
  */
  void
  apply
  (const MV& X, MV& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS, Scalar alpha = 1.0 , Scalar beta = 0.0) const;

  //! Returns the result of an inverse dft_PolyA22_Coulomb_Tpetra_Operator applied to a MultiVector X in Y.
  /*!
    \param In
    X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y - A MultiVector of dimension NumVectors containing result.
  */
  void
  applyInverse
  (const MV& X,MV& Y) const;

  //! Check for inconsistencies in operators.
  /* \param verbose (In) Print the residual of inv(A22)*A22*x_random.

  //! Throws an exception if the residual error is "large".
  */
  void
  Check
  (bool verbose) const;

  //@}

  //@{ \name Atribute access functions

  //! Returns a pointer to the Comm communicator associated with this operator.
  const RCP<const COMM> &
  Comm
  () const
  {
    return(block2Map_->getComm());
  };

  //! Returns the Map object associated with the domain of this operator.
  const RCP<const MAP > &
  getDomainMap
  () const
  {
    return(block2Map_);
  };

  //! Returns the Map object associated with the range of this operator.
  const RCP<const MAP > &
  getRangeMap
  () const
  {
    return(block2Map_);
  };
  //@}

protected:

  void insertRow();
  const RCP<const MAP> poissonMap_;
  const RCP<const MAP> cmsDensMap_;
  const RCP<const MAP> block2Map_;
  RCP<MAT> poissonOnPoissonMatrix_;
  RCP<MAT> cmsOnPoissonMatrix_;
  RCP<MAT> poissonOnDensityMatrix_;
  GlobalOrdinal curPoissonRow_;
  std::map<GlobalOrdinal, Scalar> curPoissonRowValues_;
  GlobalOrdinal curCPRow_;
  std::map<GlobalOrdinal, Scalar> curCPRowValues_;
  GlobalOrdinal curPDRow_;
  std::map<GlobalOrdinal, Scalar> curPDRowValues_;
  std::map<GlobalOrdinal, Scalar> curRowValuesCmsOnPoisson_, curRowValuesPoissonOnPoisson_, curRowValuesPoissonOnDensity_;
  Array<GlobalOrdinal> indicesCmsOnPoisson_, indicesPoissonOnPoisson_, indicesPoissonOnDensity_;
  Array<Scalar> valuesCmsOnPoisson_, valuesPoissonOnPoisson_, valuesPoissonOnDensity_;
  RCP<Hierarchy> H_;
  RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > > mueluPP_;
  RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps> > mueluPP;
  FactoryManager M_;
  RCP<Teuchos::ParameterList> MueluList_;

  using P22TO::isGraphStructureSet_;
  using P22TO::Label_;
  using P22TO::isLinearProblemSet_;
  using P22TO::firstTime_;
  using P22TO::curRow_;
  using P22TO::curRowValues_;
  using P22TO::curRowValuesCmsOnDensity_;
  using P22TO::curRowValuesCmsOnCms_;
  using P22TO::indices_;
  using P22TO::indicesCmsOnCms_;
  using P22TO::indicesCmsOnDensity_;
  using P22TO::values_;
  using P22TO::valuesCmsOnCms_;
  using P22TO::valuesCmsOnDensity_;
  using P22TO::cmsOnDensityMatrix_;
  using P22TO::densityMap_;
  using P22TO::densityOnDensityMatrix_;
  using P22TO::cmsOnCmsMatrix2_;
  using P22TO::densityOnCmsMatrix_;
  using P22TO::F_location_;
  using P22TO::cmsMap_;
  using P22TO::isFLinear_;
  using P22TO::parameterList_;
};
#endif /* DFT_POLYA22_COULOMB_TPETRA_OPERATOR_H */
