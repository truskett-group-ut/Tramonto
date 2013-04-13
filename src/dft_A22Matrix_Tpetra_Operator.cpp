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

#include "dft_A22Matrix_Tpetra_Operator.hpp"

//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_A22Matrix_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_A22Matrix_Tpetra_Operator
(const RCP<const MAP> & block2Map, RCP<ParameterList> parameterList)
  : block2Map_(block2Map),
    parameterList_(parameterList),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1) {

  A22Matrix_ = rcp(new MAT_P(block2Map, 0));
  Label_ = "dft_A22Matrix_Tpetra_Operator";
  A22Matrix_->setObjectLabel("dft_A22Matrix_Tpetra_Operator::A22Matrix");
}
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_A22Matrix_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_A22Matrix_Tpetra_Operator
()
{
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_A22Matrix_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{
  TEUCHOS_TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n");
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (firstTime_) {
    A22Matrix_->resumeFill();
    A22Matrix_->setAllToScalar(0.0);
  } else {
    A22MatrixStatic_->resumeFill();
    A22MatrixStatic_->setAllToScalar(0.0);
  }

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_A22Matrix_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value)
{
  Array<GlobalOrdinal> cols(1);
  Array<precScalar> vals(1);
  cols[0] = colGID;
  vals[0] = PREC_CAST(value);

  if (firstTime_) {
    if (rowGID!=curRow_) {
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
    }
    curRowValues_[colGID] += PREC_CAST(value);
  }
  else
    A22MatrixStatic_->sumIntoGlobalValues(rowGID, cols, vals);

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_A22Matrix_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
()
{
  if (curRowValues_.empty()) {
    return;
  }
  size_t numEntries = curRowValues_.size();
  if (numEntries>indices_.size()) {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  }
  GlobalOrdinal i=0;
  ITER pos;
  for (pos = curRowValues_.begin(); pos != curRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  A22Matrix_->insertGlobalValues(curRow_, indices_, values_);

  indices_.clear();
  values_.clear();
  curRowValues_.clear();

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_A22Matrix_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_) {
    return;
  }

  if (firstTime_) {
    insertRow();
    RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
    pl->set( "Preserve Local Graph", true );
    A22Matrix_->fillComplete(pl);

    ArrayRCP<size_t> numEntriesPerRow(block2Map_->getNodeNumElements());
    for (LocalOrdinal i = 0; i < block2Map_->getNodeNumElements(); ++i) {
      numEntriesPerRow[i] = A22Matrix_->getNumEntriesInLocalRow( i );
    }
    A22Graph_ = rcp(new GRAPH(block2Map_, A22Matrix_->getColMap(), numEntriesPerRow, Tpetra::StaticProfile));
    for (LocalOrdinal i = 0; i < block2Map_->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const precScalar> values;
      A22Matrix_->getLocalRowView( i, indices, values );
      A22Graph_->insertLocalIndices( i, indices );
    }
    A22Graph_->fillComplete();
    A22MatrixStatic_ = rcp(new MAT_P(A22Graph_));
    A22MatrixStatic_->setAllToScalar(0.0);

    for (LocalOrdinal i = 0; i < block2Map_->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const precScalar> values;
      A22Matrix_->getLocalRowView( i, indices, values );
      A22MatrixStatic_->sumIntoLocalValues( i, indices(), values() );
    }
    A22MatrixStatic_->fillComplete();
    A22MatrixOperator_ = rcp(new MMOP_P(A22MatrixStatic_));
  }

  if (!A22MatrixStatic_->isFillComplete()) {
    RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
    A22MatrixStatic_->fillComplete(pl);
  }

  LocalOrdinal overlapLevel = 0;
  A22Inverse_ = rcp(new PRECOND_AS(A22MatrixStatic_,overlapLevel));
  A22InverseOp_ = rcp(new PRECOND_AS_OP(A22Inverse_));

  TEUCHOS_TEST_FOR_EXCEPT(A22Inverse_==Teuchos::null);

  ParameterList ifpack2List = parameterList_->sublist("ifpack2ListA22");
  A22Inverse_->setParameters(ifpack2List);

  A22Inverse_->initialize();

  A22Inverse_->compute();

  isLinearProblemSet_ = true;
  firstTime_ = false;

}
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_A22Matrix_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const
{
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());

  A22InverseOp_->apply(X, Y);

}
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_A22Matrix_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());

  A22MatrixOperator_->apply(X, Y);

}
#if LINSOLVE_PREC == 0
// Use float
template class dft_A22Matrix_Tpetra_Operator<float, int, int>;
#elif LINSOLVE_PREC == 1
// Use double
template class dft_A22Matrix_Tpetra_Operator<double, int, int>;
#elif LINSOLVE_PREC == 2
// Use double double
template class dft_A22Matrix_Tpetra_Operator<dd_real, int, int>;
#elif LINSOLVE_PREC == 3
// Use quad double
template class dft_A22Matrix_Tpetra_Operator<qd_real, int, int>;
#endif
