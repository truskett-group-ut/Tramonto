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

/*
 *  FILE: dft_output.c
 *
 *  This file contains routines that post-process the density profiles.
 *
 */

#include "dft_out_profiles.h"

/*******************************************************************************
collect_x_old: This gathers all of the densities into xold on proc 0.        */ 

void collect_x_old(double **x,double *xold)
{
  int i,iunk,loc_inode, loc_i,idim,nunk_per_proc;
  int *index=NULL;
  double *unk_global, *unk_loc;

  /* allocate temporary arrays */

  nunk_per_proc = Nnodes_per_proc*Nunk_per_node;
  unk_loc = (double *) array_alloc (1, nunk_per_proc, sizeof(double));

  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ )
     for (iunk=0; iunk<Nunk_per_node; iunk++){
     unk_loc[iunk+Nunk_per_node*loc_inode] = x[iunk][L2B_node[loc_inode]];  /* always use nodal ordering here */
  }

  if (Proc == 0) {
    unk_global = (double *) array_alloc (1, Nunknowns, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }
  else {
    unk_global=NULL;
    index=NULL;
  }

  /* collect the node numbers from all the processors */

  MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,nunk_per_proc,MPI_DOUBLE,
              unk_global,Comm_unk_proc,Comm_offset_unk,
              MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  if (Proc == 0){
     for (i=0; i<Nnodes; i++){
        for (iunk=0; iunk<Nunk_per_node; iunk++){
           xold[index[i]*Nunk_per_node+iunk] = unk_global[i*Nunk_per_node+iunk];
        }
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }
  safe_free((void *) &unk_loc);

  return;
}
/*******************************************************************************
collect_vext_old: This gathers all of the external field into Vext_old on proc 0   */ 

void collect_vext_old()
{
  int i,iunk,loc_inode,nunk_per_proc;
  int *index=NULL;
  int *comm_comp_proc_tmp,*comm_offset_comp_tmp;
  double *unk_global, *unk_loc;

  /* allocate temporary arrays */

  nunk_per_proc = Nnodes_per_proc*Ncomp;
  unk_loc = (double *) array_alloc (1, nunk_per_proc, sizeof(double));

  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ )
     for (iunk=0; iunk<Ncomp; iunk++){
     unk_loc[loc_inode*Ncomp+iunk] = Vext[loc_inode][iunk];
  }


  if (Proc == 0) {
    unk_global = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }
  else {
    unk_global=NULL;
    index=NULL;
  }

  /* collect the node numbers from all the processors */
 

  MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

  /* collect the unknowns from all the processors */

  
  comm_comp_proc_tmp = (int *) array_alloc (1, Num_Proc, sizeof(int));
  comm_offset_comp_tmp = (int *) array_alloc (1, Num_Proc, sizeof(int));
  for (i=0; i<Num_Proc; i++){ comm_comp_proc_tmp[i] = 0; comm_offset_comp_tmp[i] = 0;}
  if (Proc==0){
     for (i=0; i<Num_Proc; i++){
        comm_comp_proc_tmp[i] = Comm_unk_proc[i]*Ncomp/Nunk_per_node;
        comm_offset_comp_tmp[i] = Comm_offset_unk[i]*Ncomp/Nunk_per_node;
     }
  }


  MPI_Gatherv(unk_loc,nunk_per_proc,MPI_DOUBLE,
              unk_global,comm_comp_proc_tmp,comm_offset_comp_tmp,
              MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);
  safe_free((void *) &comm_comp_proc_tmp);
  safe_free((void *) &comm_offset_comp_tmp);

  if (Proc == 0){
     for (i=0; i<Nnodes; i++){
        for (iunk=0; iunk<Ncomp; iunk++)
           Vext_old[index[i]*Ncomp+iunk] = unk_global[i*Ncomp+iunk];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }
  safe_free((void *) &unk_loc);

  return;
}
/*******************************************************************************
print_profile_box: This routine prints out the density profile. It gathers
the solution vector on Proc 0 and then prints it out using print_profile      */ 

void print_profile_box(double **x, char *outfile)
{

  if (Proc == 0){
    X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
    Vext_old = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
  }

  collect_x_old(x,X_old);
  collect_vext_old();

  if (Proc==0) {
     print_profile(outfile,X_old);
     safe_free((void *) &X_old);
     safe_free((void *) &Vext_old);
  }
}
/*******************************************************************************
print_profile: This routine prints out the density profile.        
this routine is only ever called by Proc 0                                    */ 

void print_profile(char *output_file4,double *xold)
{
  int icomp,iunk,i,inode,ijk[3],idim,ipol,iseg,itype_mer,ibond,unk_GQ,unk_B;
  int unk_field, node_start,jcomp,pol_number;
  double kappa_sq,kappa,r,rsq,bondproduct,site_dens=0.,sumsegdens[NCOMP_MAX],flag_type_mer[NMER_MAX],scale_term,scalefac;
  char *unk_char;
    
  char gfile[20],gfile2[20];
  char compfile[20],compfile2[20];
  FILE *ifp=NULL,*fp6=NULL,*fp7=NULL;
  /* 
   *  print out the densities (and electrostatic potential)
   *  to the file dft_dens.dat or dft_dens.?.?   
   */

           /* open primary output file .... densities, electrostatic potential, and CMS fields */
     ifp = fopen(output_file4,"w");

           /* open file for G_CHAIN variables ... */
     if (Iwrite==VERBOSE &&(Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)){
       sprintf(gfile,"%sg",output_file4);
       fp6 = fopen(gfile,"w");
     } 

           /* open file for segment type densities per chain ... */
     if (Iwrite==VERBOSE &&( Type_poly == WTC || Type_poly==WJDC || Type_poly==WJDC2)){
       sprintf(compfile,"%s_comp",output_file4);
       fp7 = fopen(compfile,"w");
     } 

           /* open file for segment densities */
     if (Iwrite==VERBOSE &&(Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC3)){
       sprintf(gfile2,"%s_seg",output_file4);
       fp7 = fopen(gfile2,"w");
     }

           /* print order of unknowns at the top of the file */
     for (i=0; i<NEQ_TYPE; i++){
         switch(i){
            case DENSITY: 
               if (Lseg_densities) unk_char = "DENSSEG"; 
               else                unk_char = "DENSITY"; 
               if (Phys2Nunk[i] > 0){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); break;
               }
            case MF_EQ: 
               unk_char = "MFEQ"; 
               if (Phys2Nunk[i] > 0 && (Iwrite==VERBOSE)){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); 
               } break;
            case POISSON: 
               unk_char = "POISSON"; 
               if (Phys2Nunk[i] > 0){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); break;
               }
            case DIFFUSION: 
               unk_char = "CHEMPOT";
               if (Phys2Nunk[i] > 0){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); break;
               }
            case CMS_FIELD: 
               unk_char = "CMSFIELD";
               if (Phys2Nunk[i] > 0 && (Iwrite==VERBOSE)){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); 
               }break;
            case WJDC_FIELD: 
               unk_char = "WJDCFIELD";
               if (Phys2Nunk[i] > 0 && (Iwrite==VERBOSE|| Type_poly==WJDC3)){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); 
               } break;
            case HSRHOBAR: 
               unk_char="HSRHOBAR";
               if (Phys2Nunk[i] > 0 && (Iwrite==VERBOSE)){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); 
               }break;
            case CAVWTC:
               unk_char="CAVWTC";
               if (Phys2Nunk[i] > 0 && (Iwrite==VERBOSE)){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); 
               } break;
            case BONDWTC:
               unk_char="BONDWTC";
               if (Phys2Nunk[i] > 0 && (Iwrite==VERBOSE)){
                 fputs (unk_char,ifp); 
                 fprintf(ifp,"\n"); 
               } break;
	    case SCF_CONSTR:
		unk_char="SCF_CONSTR";
		if (Phys2Nunk[i] > 0 && (Iwrite==VERBOSE)){
			fputs (unk_char,ifp); 
			fprintf(ifp,"\n"); 
		} break;				
	    case SCF_FIELD:
		unk_char="SCFFIELD";
		if (Phys2Nunk[i] > 0 && (Iwrite==VERBOSE)){
			fputs (unk_char,ifp); 
			fprintf(ifp,"\n"); 
		} break;
	 case YW_DENS:
	   break;
         }
     }

           /* compute DeBroglie wavelength for charged systems */
     if (Npoisson >0 ){
        kappa_sq = 0.0;
        if (Ipot_ff_c == COULOMB){
        for(icomp = 0; icomp<Ncomp; icomp++)
           kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*
                        Charge_f[icomp]*Charge_f[icomp];
        }
        kappa = sqrt(kappa_sq);
     }

     if (Type_poly==WTC || Type_poly==WJDC || Type_poly==WJDC2 || ((Type_poly == CMS  || Type_poly==CMS_SCFT || Type_poly==WJDC3))) {
        if (Lseg_densities) unk_char = "DENSITY"; 
        else                unk_char = "DENSSEG"; 
        if (Iwrite==VERBOSE) fputs (unk_char,fp7); 
        if (Iwrite==VERBOSE) fprintf(fp7,"\n"); 
     }
     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        node_start = Nunk_per_node*inode;

                         /* print ijk coordinates of this node in the files */ 
        for (idim=0; idim<Ndim; idim++) {
                                    fprintf(ifp,"%9.6f\t ", ijk[idim]*Esize_x[idim]);
            if (Iwrite==VERBOSE &&(Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3))  {
                                   fprintf(fp6,"%9.6f\t ",ijk[idim]*Esize_x[idim]);
            }
            if (Iwrite==VERBOSE &&(Type_poly==WTC || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly == CMS  || Type_poly==CMS_SCFT || Type_poly==WJDC3)) 
                 fprintf(fp7,"%9.6f\t ", ijk[idim]*Esize_x[idim]);
        }

        for (iunk=0; iunk<Nunk_per_node; iunk++){
            
            switch(Unk2Phys[iunk]){
               case DENSITY:   icomp = iunk-Phys2Unk_first[DENSITY]; break;
               case MF_EQ:     icomp = iunk-Phys2Unk_first[MF_EQ]; break;
               case DIFFUSION: icomp = iunk-Phys2Unk_first[DIFFUSION]; break;
               case CMS_FIELD: icomp = iunk-Phys2Unk_first[CMS_FIELD]; break;
	       case SCF_FIELD: icomp = iunk-Phys2Unk_first[SCF_FIELD]; break;
	       case SCF_CONSTR: icomp = iunk-Phys2Unk_first[SCF_CONSTR]; break;
               case WJDC_FIELD:
                    if (Type_poly==WJDC) icomp=Unk2Comp[iunk-Phys2Unk_first[WJDC_FIELD]]; 
                    else                 icomp=iunk-Phys2Unk_first[WJDC_FIELD];
                    break;
            }
            switch(Unk2Phys[iunk]){
                case DENSITY:
/*                fprintf(ifp,"%g\t", xold[iunk+node_start]/Rho_b[icomp]);*/
                  fprintf(ifp,"%g\t", xold[iunk+node_start]);
                  break;

                case DIFFUSION:
                  if (LDeBroglie){
                       fprintf(ifp,"%g\t", xold[iunk+node_start]
                            + 3.0*log(Sigma_ff[icomp][icomp]) + 1.5*log(Mass[icomp]*Temp)  );
                  }
                  else fprintf(ifp,"%g\t", xold[iunk+node_start]);

                case POISSON:
                  fprintf(ifp,"%g\t", xold[iunk+node_start]);
                  break;

                case MF_EQ:
                case HSRHOBAR:
                case CAVWTC:
                case BONDWTC:
                  if (Iwrite==VERBOSE) fprintf(ifp,"%g\t", xold[iunk+node_start]);
                  break;

                case CMS_FIELD:
                case WJDC_FIELD:
	        case SCF_FIELD:
                   if(Iwrite==VERBOSE || Type_poly==WJDC3){
                      /*if (fabs(xold[iunk+node_start]) > 1.e-12 && -log(xold[iunk+node_start]) < VEXT_MAX){
                          fprintf(ifp,"%g\t", -log(xold[iunk+node_start]));
                      }
                      else fprintf(ifp,"%g\t", VEXT_MAX);*/

                     if (Unk2Phys[iunk]!=WJDC_FIELD){ 
                         fprintf(ifp,"%g\t", xold[iunk+node_start]);
                     }
                     else{
                       for (pol_number=0;pol_number<Npol_comp;pol_number++) {
                          if (Nseg_type_pol[pol_number][icomp] !=0) scalefac=Scale_fac_WJDC[pol_number][icomp];
                        }
                        fprintf(ifp,"%g\t", xold[iunk+node_start]/exp(scalefac));
                     }
                   }
                   break;
					
	       case SCF_CONSTR:
		    if(Iwrite==VERBOSE){
		         fprintf(ifp,"%g\t", xold[iunk+node_start]);
		    }
                   break;

                case G_CHAIN:
                   if (Iwrite==VERBOSE) fprintf(fp6,"%g\t", xold[iunk+node_start]);
                   break;
            }

        }    /* end loop over unknowns in the run */

                /* print the Poisson-Boltzmann solution based on the computed electrostatic field */
        if (Ipot_ff_c == 1 && Type_poly==NONE){
        for (icomp=0; icomp<Ncomp; icomp++)
          fprintf(ifp,"%g\t",
                  Rho_b[icomp]*exp(-Charge_f[icomp]*xold[Phys2Unk_first[POISSON]+node_start]
                                                              -Vext_old[inode*Ncomp+icomp]));
        }
 
               /* print segment densities for a CMS polymer run ... print component densities in WTC run*/
        if ((Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC3)){
              for (itype_mer=0;itype_mer<Ncomp;itype_mer++) sumsegdens[itype_mer]=0.0;
              for (ipol=0; ipol<Npol_comp; ipol++){
                 if (Type_poly==WJDC3){
                 scale_term=0.0;
                 for (jcomp=0;jcomp<Ncomp;jcomp++) {
                        scale_term-=Scale_fac_WJDC[ipol][jcomp]*Nseg_type_pol[ipol][jcomp];
                 }
                 }
                 for(iseg=0;iseg<Nmer[ipol];iseg++){
                    itype_mer=Type_mer[ipol][iseg];
                    bondproduct=1.0;
                    for(ibond=0;ibond<Nbond[ipol][iseg];ibond++){
			 unk_GQ  = Geqn_start[ipol] + Poly_to_Unk[ipol][iseg][ibond];
                         bondproduct *= xold[unk_GQ+node_start];
                    }  

					/* check this code!!! needs fixing */
		   if(Type_poly == CMS)            unk_B=Phys2Unk_first[CMS_FIELD]+itype_mer;
		   else if(Type_poly == CMS_SCFT)  unk_B=Phys2Unk_first[SCF_FIELD]+itype_mer;
                   else if (Type_poly==WJDC3)      unk_B=Phys2Unk_first[WJDC_FIELD]+itype_mer;

                   if (fabs(xold[unk_B+node_start])>1.e-12){
                      site_dens=bondproduct*POW_DOUBLE_INT(xold[unk_B+node_start],-(Nbond[ipol][iseg]-1));
                      if (Type_poly==CMS || Type_poly==CMS_SCFT) site_dens*=Rho_b[itype_mer]/Nmer_t[ipol][itype_mer];
                      else if (Type_poly==WJDC3)                 site_dens*=exp(Betamu_chain[ipol]+scale_term);
                   }
                   else site_dens=0.0;

                   sumsegdens[itype_mer]+=site_dens;
                   if (Iwrite==VERBOSE)fprintf(fp7,"%g\t", site_dens);
                 }
              }
        }
        else if (Type_poly==WTC || Type_poly==WJDC || Type_poly==WJDC2){
              for (ipol=0; ipol<Npol_comp; ipol++){
                for (itype_mer=0;itype_mer<Ncomp;itype_mer++) {
                     sumsegdens[itype_mer]=0.0;
                     flag_type_mer[itype_mer]=FALSE;
                }
                for (iseg=0; iseg<Nmer[ipol]; iseg++){
                   iunk=Phys2Unk_first[DENSITY]+SegChain2SegAll[ipol][iseg];
                   sumsegdens[Type_mer[ipol][iseg]]+=xold[iunk+node_start];
                   flag_type_mer[Type_mer[ipol][iseg]]=TRUE;
                }
                for (itype_mer=0;itype_mer<Ncomp;itype_mer++){
                   if (Iwrite==VERBOSE && flag_type_mer[itype_mer]==TRUE)  fprintf(fp7,"%g\t",sumsegdens[itype_mer]);
                }
              }
        }
 
                /* add a carriage return to the file to start a new line */
        fprintf(ifp,"\n");
        if (Iwrite==VERBOSE && (Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)) fprintf(fp6,"\n");
        if (Iwrite==VERBOSE && (Type_poly==WTC || Type_poly==WJDC || Type_poly==WJDC2 ||((Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC3)))) fprintf(fp7,"\n");

                /* add some blank lines for improved graphics in 2D and 3D gnuplot */
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");

     }    /* loop over all nodes  */

          /* close files */
     fclose(ifp);
     if (Iwrite==VERBOSE &&(Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC || Type_poly==WJDC2 | Type_poly==WJDC3)) fclose(fp6);
     if (Iwrite==VERBOSE &&(Type_poly==WTC || Type_poly==WJDC || Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC2 || Type_poly==WJDC3)) fclose(fp7);

  return;
}
/*******************************************************************************
print_gofr: This routine prints out the density profile.        
this routine is only ever called by Proc 0                                    */ 

void print_gofr(char *output_file6,double *xold)
{
  int icomp,i,inode,ijk[3],idim,nunk_print,npol=0,itype_mer,iwall,iunk,end_loop;
  double kappa_sq,kappa,r,rsq;
  FILE *ifp=NULL;
  /* 
   *  print out the densities (and electrostatic potential)
   *  to the file dft_dens.dat or dft_dens.?.?   
   */

     ifp = fopen(output_file6,"w");

     if (L_HSperturbation) nunk_print = Nunk_per_node;
     else nunk_print = 2*Ncomp;

     if(Nwall > 0) end_loop=Nwall;
     if(Nwall==0 && Nlocal_charge !=0) end_loop=Nlocal_charge;

     for (iwall=0; iwall<end_loop; iwall++){  /*compute g(r) for different atoms
                                             in one linked wall --- e.g. 
                                             could represent H2O as 3 atoms
                                             then need to compute gHH,gOO,gHO */
     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
 
        rsq=0.0; r=0.0;
        for (idim=0; idim<Ndim; idim++) {
            if (Nwall==0 && Nlocal_charge ==1){
                rsq = rsq+(ijk[idim]*Esize_x[idim]-(Charge_x[idim][iwall]+0.5*Size_x[idim]))*
                          (ijk[idim]*Esize_x[idim]-(Charge_x[idim][iwall]+0.5*Size_x[idim]));
            }
            else{
                rsq = rsq+(ijk[idim]*Esize_x[idim]-(WallPos[idim][iwall]+0.5*Size_x[idim]))*
                          (ijk[idim]*Esize_x[idim]-(WallPos[idim][iwall]+0.5*Size_x[idim]));
            }
        }
        if (rsq > 0.0) r=sqrt(rsq); 
        fprintf(ifp,"%9.6f\t ",r);

        for (iunk=0; iunk<Nunk_per_node; iunk++){
            if (Unk2Phys[iunk]==DENSITY){
                icomp = iunk-Phys2Unk_first[DENSITY];
                if (Lprint_gofr==2) {
                    if (xold[iunk+Nunk_per_node*inode]>1.e-8) fprintf(ifp,"%22.17f\t", -log(xold[iunk+Nunk_per_node*inode]/Rho_b[icomp]));
                    else fprintf(ifp,"%22.17f\t",VEXT_MAX);
                }
                else fprintf(ifp,"%22.17f\t", xold[iunk+Nunk_per_node*inode]/Rho_b[icomp]);
            }
        }

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");

     } /* loop over nodes */
     }    /* loop over all walls  */
     fclose(ifp);
  return;
}
/************************************************************************
print_zeroTF: This routine collects the zero_TF array and prints it out  */
void print_zeroTF(int **zero_TF, char *output_file)
{
  int icomp,loc_inode,inode,ijk[3],*index,idim,inode_box;
  int *unk_loc,*unk_global,**unk;
  FILE *ifp=NULL;

  if (Proc == 0) {
       ifp = fopen(output_file,"w");
       unk = (int **) array_alloc (2, Nnodes, Ncomp+1, sizeof(int));
  }

  unk_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));

  for (icomp=0; icomp<Ncomp+1; icomp++){

     /*  define the local array for each component separately */

     for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
         inode_box = L2B_node[loc_inode];
         unk_loc[loc_inode] = zero_TF[inode_box][icomp];
     }

     if (Proc ==0){
       index = (int *) array_alloc (1, Nnodes, sizeof(int));
       unk_global = (int *) array_alloc (1, Nnodes, sizeof(int));
     }

     /* collect the global indices from all processors */
     MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

     /* collect the unknowns from all the processors */

     MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_INT,
              unk_global,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

     if (Proc == 0){
        for (inode=0; inode<Nnodes; inode++){
            unk[index[inode]][icomp] = unk_global[inode];
        }
        safe_free((void *) &unk_global);
        safe_free((void *) &index);
     }
  }
  safe_free((void *) &unk_loc);
     
  /* 
   *  now print out the Zero_density_TF array to the file dft_zeroTF.dat
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        for (idim=0; idim<Ndim; idim++){
            fprintf(ifp,"%9.6f\t ", ijk[idim]*Esize_x[idim]);
        }
        for (icomp=0; icomp<Ncomp+1; icomp++){
            fprintf(ifp,"%d\t", unk[inode][icomp]);
        }

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */
     fclose(ifp);
     safe_free((void *) &unk);
  }       /* end of Proc ==0 test */


  return;
}
/************************************************************************
print_Nodes_to_zone: This routine collects and prints nodes_to_zone  */
void print_Nodes_to_zone(int *node_to_zone, char *output_file)
{
  int loc_inode,inode,ijk[3],*index,idim,inode_box;
  int *unk,*unk_loc, *unk_global;
  FILE *ifp=NULL;

  if (Proc == 0){
     ifp = fopen(output_file,"w");
     unk = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  unk_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));

  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
      inode_box = L2B_node[loc_inode];
      unk_loc[loc_inode] = node_to_zone[inode_box];
  }

  if (Proc ==0){
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
    unk_global = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  /* collect the global indices from all processors */
  MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
           index,Comm_node_proc,Comm_offset_node,
           MPI_INT,0,MPI_COMM_WORLD);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_INT,
           unk_global,Comm_node_proc,Comm_offset_node,
           MPI_INT,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  if (Proc == 0){
     for (inode=0; inode < Nnodes; inode++ ){
        unk[index[inode]] = unk_global[inode];
     }

     if (Proc==0) safe_free((void *) &index);
     if (Proc==0) safe_free((void *) &unk_global);
  }

  /* 
   * now print the array.
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t ",
            (double)ijk[idim]*Esize_x[idim]);

            fprintf(ifp,"%d \t", unk[inode]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */
     fclose(ifp);
     safe_free((void *) &unk);
  }       /* end of Proc ==0 test */

  return;
}
/************************************************************************
print_charge_surf: This routine collects and prints Charge_w_sum_els  */
void print_charge_surf(double **charge_w_sum, char *output_file)
{
  int i,icount;
  int loc_inode,inode,ijk[3],*index_loc,*index,idim;
  int reflect_flag[3];
  int *comm_icount_proc, *comm_offset_icount;
  double **unk,*unk_loc, *unk_global;
  FILE *ifp=NULL;


  if (Proc == 0){
     unk = (double **) array_alloc (2, Nnodes, Ndim, sizeof(double));
     ifp = fopen(output_file,"w");
  }
  reflect_flag[0] = reflect_flag[1] = reflect_flag[2] = FALSE;

  for (idim=0; idim<Ndim; idim++) {

  index_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));
  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  icount=0;
  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
      index_loc[icount] = loc_inode;
      unk_loc[icount++] = charge_w_sum[loc_inode][idim];
  }

  if (Proc ==0){
    unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  comm_icount_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
  comm_offset_icount = (int *) array_alloc (1, Num_Proc, sizeof(int));

  MPI_Gather(&icount,1,MPI_INT,
             comm_icount_proc,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc == 0){
     comm_offset_icount[0] = 0; 
     for (i=1; i<Num_Proc; i++){
        comm_offset_icount[i] = comm_offset_icount[i-1] + comm_icount_proc[i-1];
     }
  }

  /* collect the global indices from all processors */
  MPI_Gatherv(index_loc,Nnodes_per_proc,MPI_INT,
           index,comm_icount_proc,comm_offset_icount,
           MPI_INT,0,MPI_COMM_WORLD);
  safe_free((void *) &index_loc);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_DOUBLE,
           unk_global,comm_icount_proc,comm_offset_icount,
           MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  
  safe_free((void *) &comm_icount_proc);
  safe_free((void *) &comm_offset_icount);

  if (Proc == 0) {
     for (inode=0; inode<Nnodes; inode++){
         unk[index[inode]][idim] = unk_global[inode];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }

  }


  /* 
   *  now print out the volumetric charges
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t ", ((double)ijk[idim])*Esize_x[idim]);

        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t", unk[inode][idim]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]) fprintf(ifp,"\n");
     }    /* loop over all nodes  */

     safe_free((void *) &unk);
     fclose(ifp);
  }       /* end of Proc ==0 test */
  return;
}
/************************************************************************
print_free_energy_profile: This routine collects and prints freen_profile_1D  */
void print_freen_profile_1D(double *freen_profile_1D, char *output_file)
{
  int i,iel,iel_box,icount,logical;
  int loc_inode,inode,ijk[3],*index_loc,*index,idim,inode_box;
  int reflect_flag[3];
  int *comm_icount_proc, *comm_offset_icount;
  double *unk,*unk_loc, *unk_global,charge_total;
  FILE *ifp=NULL;

  reflect_flag[0] = reflect_flag[1] = reflect_flag[2] = FALSE;

  if (Proc == 0){
     unk = (double *) array_alloc (1, Nnodes, sizeof(double));
     ifp = fopen(output_file,"w");
  }

  index_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));
  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  icount=0;
  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
      inode_box = L2B_node[loc_inode];
      inode = L2G_node[loc_inode];
      index_loc[icount] = inode;
      unk_loc[icount++] = freen_profile_1D[loc_inode];
  }

  if (Proc ==0){
    unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  comm_icount_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
  comm_offset_icount = (int *) array_alloc (1, Num_Proc, sizeof(int));

  MPI_Gather(&icount,1,MPI_INT, comm_icount_proc,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc == 0){
     comm_offset_icount[0] = 0; 
     for (i=1; i<Num_Proc; i++){
        comm_offset_icount[i] = comm_offset_icount[i-1] + comm_icount_proc[i-1];
     }
  }

  /* collect the global indices from all processors */
  MPI_Gatherv(index_loc,icount,MPI_INT,
           index,comm_icount_proc,comm_offset_icount,
           MPI_INT,0,MPI_COMM_WORLD);
  safe_free((void *) &index_loc);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,icount,MPI_DOUBLE,
           unk_global,comm_icount_proc,comm_offset_icount,
           MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  
  safe_free((void *) &comm_icount_proc);
  safe_free((void *) &comm_offset_icount);

  if (Proc == 0) {
     for (inode=0; inode<Nnodes; inode++){
         unk[index[inode]] = unk_global[inode];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }

  /* 
   *  now print out the free energy profile
   */
  if (Proc ==0){
     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
/*        for (idim=0; idim<Ndim; idim++)*/  /* only do free energy profile in 1D so far */
            idim=0;
            fprintf(ifp,"%9.6f\t ",
            ((double)ijk[idim])*Esize_x[idim]);

            fprintf(ifp,"%9.6f\t", unk[inode]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */

     safe_free((void *) &unk);
     fclose(ifp);
  }       /* end of Proc ==0 test */
  return;
}
/************************************************************************
print_charge_vol: This routine collects and prints Charge_vol_els  */
void print_charge_vol(double *charge_els, char *output_file)
{
  int i,iel,iel_box,icount,logical;
  int loc_inode,inode,ijk[3],*index_loc,*index,idim,inode_box;
  int reflect_flag[3];
  int *comm_icount_proc, *comm_offset_icount;
  double *unk,*unk_loc, *unk_global,charge_total;
  FILE *ifp=NULL;

  reflect_flag[0] = reflect_flag[1] = reflect_flag[2] = FALSE;

  if (Proc == 0){
     unk = (double *) array_alloc (1, Nelements, sizeof(double));
     ifp = fopen(output_file,"w");
  }

  index_loc = (int *) array_alloc (1, Nnodes_per_proc, sizeof(int));
  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  icount=0;
  for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
      inode_box = L2B_node[loc_inode];
      inode = L2G_node[loc_inode];
      node_to_ijk(inode,ijk);

      logical = FALSE;
      for (idim=0; idim<Ndim; idim++) 
          if (ijk[idim] == Nodes_x[idim] -1 
              && Type_bc[idim][1] != PERIODIC) logical = TRUE;

      if (!logical) {
         iel   = node_to_elem(inode,0,reflect_flag);
         iel_box = node_box_to_elem_box_reflect(inode_box,0,reflect_flag);
         index_loc[icount] = iel;
         unk_loc[icount++] = charge_els[iel_box];
      }
  }

  if (Proc ==0){
    unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
    index = (int *) array_alloc (1, Nnodes, sizeof(int));
  }

  comm_icount_proc = (int *) array_alloc (1, Num_Proc, sizeof(int));
  comm_offset_icount = (int *) array_alloc (1, Num_Proc, sizeof(int));

  MPI_Gather(&icount,1,MPI_INT, comm_icount_proc,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc == 0){
     comm_offset_icount[0] = 0; 
     for (i=1; i<Num_Proc; i++){
        comm_offset_icount[i] = comm_offset_icount[i-1] + comm_icount_proc[i-1];
     }
  }

  /* collect the global indices from all processors */
  MPI_Gatherv(index_loc,icount,MPI_INT,
           index,comm_icount_proc,comm_offset_icount,
           MPI_INT,0,MPI_COMM_WORLD);
  safe_free((void *) &index_loc);

  /* collect the unknowns from all the processors */

  MPI_Gatherv(unk_loc,icount,MPI_DOUBLE,
           unk_global,comm_icount_proc,comm_offset_icount,
           MPI_DOUBLE,0,MPI_COMM_WORLD);
  safe_free((void *) &unk_loc);

  
  safe_free((void *) &comm_icount_proc);
  safe_free((void *) &comm_offset_icount);

  if (Proc == 0) {
     charge_total=0.0;
     for (iel=0; iel<Nelements; iel++){
         unk[index[iel]] = unk_global[iel];
         charge_total+=unk_global[iel];
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
     printf("THE NET VOLUME CHARGE IN THE DOMAIN IS : %9.6f\n",charge_total);
  }

  /* 
   *  now print out the volumetric charges
   */
  if (Proc ==0){

     for (iel=0; iel<Nelements; iel++){
        inode = element_to_node(iel);
        node_to_ijk(inode,ijk);
        for (idim=0; idim<Ndim; idim++)
            fprintf(ifp,"%9.6f\t ",
            ((double)ijk[idim]+0.5)*Esize_x[idim]);

            fprintf(ifp,"%9.6f\t", unk[iel]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */

     safe_free((void *) &unk);
     fclose(ifp);
  }       /* end of Proc ==0 test */
  return;
}
/************************************************************************
print_vext: This routine collects the vext array and prints it out  */
void print_vext(double **vext, char *output_file)
{
  int icomp,loc_inode,inode,ijk[3],*index,idim;
  double *unk_loc,*unk_global,**unk,rsq,r;
  FILE *ifp=NULL;

  if (Proc == 0) {
       ifp = fopen(output_file,"w");
       unk = (double **) array_alloc (2, Nnodes, Ncomp, sizeof(double));
  }

  unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));

  for (icomp=0; icomp<Ncomp; icomp++){

     /* 
      *  first collect all the icomp unknowns on processor 0 
      */

     for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ){
         unk_loc[loc_inode] = vext[loc_inode][icomp];
     }

     if (Proc ==0){
       index = (int *) array_alloc (1, Nnodes, sizeof(int));
       unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
     }

     /* collect the global indices from all processors */
     MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

     /* collect the unknowns from all the processors */

     MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_DOUBLE,
              unk_global,Comm_node_proc,Comm_offset_node,
              MPI_DOUBLE,0,MPI_COMM_WORLD);

     if (Proc == 0){
        for (inode=0; inode<Nnodes; inode++){
            unk[index[inode]][icomp] = unk_global[inode];
        }
        safe_free((void *) &unk_global);
        safe_free((void *) &index);
     }
  }
  safe_free((void *) &unk_loc);

     
  /* 
   *  now print out the densities (and electrostatic potential)
   *  to the file dft_dens.dat. 
   */
  if (Proc ==0){

     for (inode=0; inode<Nnodes; inode++){
        node_to_ijk(inode,ijk);
        r=0.0; rsq=0.0;
        for (idim=0; idim<Ndim; idim++){
            fprintf(ifp,"%9.6f\t ", ijk[idim]*Esize_x[idim]);
            rsq +=  ijk[idim]*Esize_x[idim]* ijk[idim]*Esize_x[idim];
        }
        r=sqrt(rsq);
        if (Lprint_gofr) fprintf(ifp,"%9.6f\t ",r);

        for (icomp=0; icomp<Ncomp; icomp++)
            fprintf(ifp,"%22.17f\t", unk[inode][icomp]);

        fprintf(ifp,"\n");
        if (ijk[0] == Nodes_x[0]-1) fprintf(ifp,"\n");
     }    /* loop over all nodes  */
     fclose(ifp);
     safe_free((void *) &unk);
  }       /* end of Proc ==0 test */
  return;
}
/************************************************************************/
