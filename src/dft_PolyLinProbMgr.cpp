/*@HEADER
// ***********************************************************************
// 
//                Tramonto: Molecular Theories Modeling Code
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// Questions? Contact Laura J.D. Frink (ljfrink@sandia.gov)
// 
// ***********************************************************************
//@HEADER
*/

#include "dft_PolyLinProbMgr.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "AztecOO.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "dft_PolyA11_Epetra_Operator.hpp"
#include "dft_PolyA22_Epetra_Operator.hpp"
#include "dft_PolyA22Full_Epetra_Operator.hpp"
#include "dft_PolyA22Bsor_Epetra_Operator.hpp"
#include "dft_Schur_Epetra_Operator.hpp"


//=============================================================================
dft_PolyLinProbMgr::dft_PolyLinProbMgr(int numUnknownsPerNode, int * solverOptions, double * solverParams, MPI_Comm comm, bool debug) 
  : dft_BasicLinProbMgr(numUnknownsPerNode, solverOptions, solverParams, comm),
    debug_(debug) {
  

  return;
}
//=============================================================================
dft_PolyLinProbMgr::~dft_PolyLinProbMgr() {
   return;
}
//=============================================================================
int dft_PolyLinProbMgr::finalizeBlockStructure() {

  if (isBlockStructureSet_) return(1); // Already been here, return warning

  if (numGlobalNodes_==0 ||
      numGlobalBoxNodes_==0 ||
      gInvEquations_.Length()==0 ||
      gInvEquations_.Length()==0 ||
      cmsEquations_.Length()==0 ||
      densityEquations_.Length()==0) return(-1); // Error: One or more set methods not called
      
  
  // Fill physics ordering vector with the concatenated contents of the IDs for all physics types
  // Load Schur block mappings
  physicsOrdering_.Size(numUnknownsPerNode_);
  physicsIdToSchurBlockId_.Size(numUnknownsPerNode_);
  int * ptr = physicsOrdering_.Values();
  for (int i=0; i<gEquations_.Length(); i++) {
    *ptr++ = gEquations_[i];
    physicsIdToSchurBlockId_[gEquations_[i]] = 1;
  }
  for (int i=0; i<gInvEquations_.Length(); i++) {
    *ptr++ = gInvEquations_[i];
    physicsIdToSchurBlockId_[gInvEquations_[i]] = 1;
  }
  for (int i=0; i<cmsEquations_.Length(); i++) {
    *ptr++ = cmsEquations_[i];
    physicsIdToSchurBlockId_[cmsEquations_[i]] = 2;
  }
  for (int i=0; i<densityEquations_.Length(); i++) {
    *ptr++ = densityEquations_[i];
    physicsIdToSchurBlockId_[densityEquations_[i]] = 2;
  }

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.Size(numUnknownsPerNode_);
  for (int i=0; i<physicsOrdering_.Length(); i++) solverOrdering_[physicsOrdering_[i]]=i;

  const int numUnks = numOwnedNodes_*numUnknownsPerNode_;
  const int numUnks1 = numOwnedNodes_*(gEquations_.Length()+gInvEquations_.Length());
  const int numUnks2 = numOwnedNodes_*(cmsEquations_.Length()+densityEquations_.Length());
  assert(numUnks==(numUnks1+numUnks2));  // Sanity test
  const int numCms = numOwnedNodes_*(cmsEquations_.Length());
  const int numDensity = numOwnedNodes_*(densityEquations_.Length());
  Epetra_IntSerialDenseVector globalGIDList(numUnks);

  int * GIDs = ownedMap_->MyGlobalElements();
  int k=0;
  for (int i=0; i<numUnknownsPerNode_; i++) {
    int ii=physicsOrdering_[i];
    for (int j=0; j<numOwnedNodes_; j++) 
      globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
  }

  ptr = globalGIDList.Values();
  globalRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks, ptr, 0, comm_));
  block1RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks1, ptr, 0, comm_));
  block2RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks2, ptr+numUnks1, 0, comm_));
  cmsRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numCms, ptr+numUnks1, 0, comm_));
  densityRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numDensity, ptr+numUnks1+numCms, 0, comm_));
  /*
    std::cout << " Global Row Map" << *globalRowMap_.get() << std::endl
    << " Block 1 Row Map " << *block1RowMap_.get() << std::endl
    << " Block 2 Row Map " << *block2RowMap_.get() << std::endl
    << " CMS     Row Map " << *cmsRowMap_.get() << std::endl
    << " Density Row Map " << *densityRowMap_.get() << std::endl;
  */
  A11_ = Teuchos::rcp(new dft_PolyA11_Epetra_Operator(*(ownedMap_.get()), *(block1RowMap_.get())));
  A12_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *(block1RowMap_.get()), 0));
  A21_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *(block2RowMap_.get()), 0));
  //A22_ = Teuchos::rcp(new dft_PolyA22_Epetra_Operator(*(cmsRowMap_.get()), *(densityRowMap_.get()), *(block2RowMap_.get())));
  //A22_ = Teuchos::rcp(new dft_PolyA22Full_Epetra_Operator(*(cmsRowMap_.get()), *(densityRowMap_.get()), *(block2RowMap_.get())));
  A22_ = Teuchos::rcp(new dft_PolyA22Bsor_Epetra_Operator(*(cmsRowMap_.get()), *(densityRowMap_.get()), *(block2RowMap_.get())));
  if (debug_) 
    globalMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *globalRowMap_, 0));
  else
    globalMatrix_ = Teuchos::null; // not used by this solver
  
  globalRhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalLhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));

  rhs1_ = Teuchos::rcp(new Epetra_Vector(View, *(block1RowMap_.get()), globalRhs_->Values()));
  rhs2_ = Teuchos::rcp(new Epetra_Vector(View, *(block2RowMap_.get()), globalRhs_->Values()+numUnks1));
  rhsSchur_ = Teuchos::rcp(new Epetra_Vector(*(rhs2_.get())));
  lhs1_ = Teuchos::rcp(new Epetra_Vector(View, *(block1RowMap_.get()), globalLhs_->Values()));
  lhs2_ = Teuchos::rcp(new Epetra_Vector(View, *(block2RowMap_.get()), globalLhs_->Values()+numUnks1));

  schurOperator_ = Teuchos::rcp(new dft_Schur_Epetra_Operator(A11_.get(), A12_.get(), A21_.get(), A22_.get()));
  implicitProblem_ = Teuchos::rcp(new Epetra_LinearProblem(schurOperator_.get(), lhs2_.get(), rhsSchur_.get()));

    
  ownedToBoxImporter_ = Teuchos::rcp(new Epetra_Import(*(boxMap_.get()), *(ownedMap_.get())));

  isBlockStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    A12_->PutScalar(0.0);
    A21_->PutScalar(0.0);
    globalRhs_->PutScalar(0.0);
    globalLhs_->PutScalar(0.0);
    if (debug_) globalMatrix_->PutScalar(0.0);
  }

  A11_->initializeProblemValues();
  A22_->initializeProblemValues();
  
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value) {

  int schurBlockRow = physicsIdToSchurBlockId_[ownedPhysicsID];
  int schurBlockCol = physicsIdToSchurBlockId_[boxPhysicsID];
  int rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  int colGID = boxToSolverGID(boxPhysicsID, boxNode);
  if (schurBlockRow==1 && schurBlockCol==1) { // A11 block
    A11_->insertMatrixValue(solverOrdering_[ownedPhysicsID], ownedMap_->GID(ownedNode), rowGID, colGID, value); 
  }
  else if (schurBlockRow==2 && schurBlockCol==2) { // A22 block
    A22_->insertMatrixValue(rowGID, colGID, value); 
  }
  else if (schurBlockRow==2 && schurBlockCol==1) { // A21 block
    if (firstTime_)
      A21_->InsertGlobalValues(rowGID, 1, &value, &colGID);
    else
      A21_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
  }
  else { // A12 block
    if (firstTime_)
      A12_->InsertGlobalValues(rowGID, 1, &value, &colGID);
    else
      A12_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
  }

  if (debug_) dft_BasicLinProbMgr::insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, value);
  
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) {
    A12_->FillComplete(*(block2RowMap_.get()),*(block1RowMap_.get()));
    A12_->OptimizeStorage();
    A21_->FillComplete(*(block1RowMap_.get()),*(block2RowMap_.get()));
    A21_->OptimizeStorage();

    if (debug_) globalMatrix_->FillComplete();
  }
  //std::cout << *A12_.get() << endl 
  //          << *A21_.get() << endl;

  A11_->finalizeProblemValues();
  A22_->finalizeProblemValues();

  //Check(true);
  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::setupSolver() {

  if (!isLinearProblemSet_) return(-1);

  schurOperator_->ComputeRHS(*rhs1_.get(), *rhs2_.get(), *rhsSchur_.get());
  
  solver_ = Teuchos::rcp(new AztecOO(*(implicitProblem_.get())));

  if (solverOptions_!=0) solver_->SetAllAztecOptions(solverOptions_);
  if (solverParams_!=0) solver_->SetAllAztecParams(solverParams_);

  const int * options = solver_->GetAllAztecOptions();
  const double * params = solver_->GetAllAztecParams();

  solver_->SetAztecOption(AZ_scaling, AZ_none); 
  int maxiter = 500;
  solver_->SetAztecOption(AZ_max_iter, maxiter);
  solver_->SetAztecOption(AZ_kspace, maxiter); 
  //solver_->SetAztecParam(AZ_tol, 1.0e-12); 
  //solver_->SetAztecOption(AZ_conv, AZ_noscaled); 
  solver_->SetPrecOperator(A22_.get());
  //solver_->SetAztecParam(AZ_ill_cond_thresh, 0.0);
  //solver_->SetAztecOption(AZ_precond, AZ_none);
  //solver_->SetAztecOption(AZ_solver, AZ_bicgstab);
  //solver_->SetAztecOption(AZ_precond, AZ_dom_decomp);
  //solver_->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  //solver_->SetAztecParam(AZ_ilut_fill, 4.0);
  //solver_->SetAztecParam(AZ_drop, 0.0);

  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::solve() {
  
  //writeMatrix("2D.mm", "Small Polymer Matrix", "Global Matrix from Small Polymer Problem");
  //abort();

  const int * options = solver_->GetAllAztecOptions();
  const double * params = solver_->GetAllAztecParams();

  solver_->Iterate(options[AZ_max_iter], params[AZ_tol]); // Try to solve
  schurOperator_->ComputeX1(*rhs1_.get(), *lhs2_.get(), *lhs1_.get()); // Compute rest of solution

  if (debug_) {
    Epetra_Vector tmpRhs(*globalRowMap_);
    Epetra_Vector tmprhs1(View, *(block1RowMap_.get()), tmpRhs.Values());
    Epetra_Vector tmprhs2(View, *(block2RowMap_.get()), tmpRhs.Values()+block1RowMap_->NumMyElements());
    
    schurOperator_->ApplyGlobal(*lhs1_.get(), *lhs2_.get(), tmprhs1, tmprhs2);
    
    tmpRhs.Update(-1.0, *globalRhs_.get(), 1.0);
    double resid=0.0;
    tmpRhs.Norm2(&resid);
    std::cout << "Global Residual for solution = " << resid << std::endl;
    bool writeMatrixNow = false;
    if (writeMatrixNow) {
      writeMatrix("A.dat", "GlobalMatrix", "GlobalMatrix");
      writeLhs("x.dat");
      writeRhs("b.dat");
      writePermutation("p.dat");
      abort();
    }
  }
  //std::cout << "Global RHS = " << *globalRhs_.get() << std::endl
  //          << "Global LHS = " << *globalLhs_.get() << std::endl;

  //solver_->AdaptiveIterate(solverOptions_[AZ_max_iter], 5, solverParams_[AZ_tol]); // Try to solve
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::applyMatrix(const double** x, double** b) const {
  
  setLhs(x);
  schurOperator_->ApplyGlobal(*lhs1_.get(), *lhs2_.get(), *rhs1_.get(), *rhs2_.get());
  getRhs(b);
  
  return(0);
}
//=============================================================================
  int dft_PolyLinProbMgr::Check(bool verbose) const {

  int ierr1 = A11_->Check(verbose);
  int ierr2 = A22_->Check(verbose);
  if (ierr1!=0 || ierr2!=0) return(-1);
  return(0);
    
  }
//=============================================================================
int dft_PolyLinProbMgr::writeMatrix(const char * filename, const char * matrixName, const char * matrixDescription) const  {
  if (debug_)
    return(EpetraExt::RowMatrixToMatrixMarketFile(filename, *globalMatrix_.get(), matrixName, matrixDescription));
  else
    return(-1); // Not available if not in debug mode
}
