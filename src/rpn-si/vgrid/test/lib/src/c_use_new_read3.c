/* libdescrip - Vertical grid descriptor library for FORTRAN programming
 * Copyright (C) 2016  Direction du developpement des previsions nationales
 *                     Centre meteorologique canadien
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vgrid.h"
#include "rmn.h"
#include "rmn/rpnmacros.h"
#include "c_ut_report.h"

int main() {

  int ier, iun = 10, iun2 = 11;
  int quiet = 0, *i_val = NULL, in_log = 0, dpidpis = 0;
  int nl_t, nt, ni, nj, nk, ni2, nj2, nk2, k, key, ij, ijk, status;
  char filename[]="data/dm_5005_from_model_run";
  char mode[]="RND";
  char format[] = "FST";
  char name[5];
  float *f_val = NULL, *p0 = NULL, *px = NULL;
  double *a_8_t = NULL, *b_8_t = NULL, *table = NULL, *levels_8 = NULL, *p0_8 = NULL;
  vgrid_descriptor *vgd = NULL, *vgd2 = NULL;
  int   dateo, deet, npas, nbits, datyp, ip1, ip2, ip3, ig1, ig2, ig3, ig4;
  int   swa, lng, dltf, ubc, extra1, extra2, extra3;
  char  typvar[3], nomvar[5], etiket[13], grtyp[2];
  
  status = VGD_OK;

  // Initialisation of string must be done for fstprm to work!
  strcpy(nomvar, "    ");
  strcpy(typvar, "  ");
  strcpy(etiket, "            ");
  strcpy(grtyp, " ");
  
  ier = c_fnom(&iun,filename,mode,0);
  if( ier < 0 ) {
    printf("ERROR with c_fnom on iun, file %s\n", filename);
    return(1);
  }
  ier = c_fstouv(iun,"RND");  
  if( ier < 0 ) {
    printf("ERROR with c_fstouv on iun, file %s\n", filename);
    return(1);
  }
  key = c_fstinf(iun,&ni,&nj,&nk,-1," ",-1,  -1,  -1," ","!!");
  if(key < 0){
    printf("ERROR cannot find !!\n");
    return(1);
  }
  if( c_fstprm(key,
	       &dateo,  &deet,   &npas, 
	       &ni,     &nj,     &nk,
	       &nbits,  &datyp,  
	       &ip1,    &ip2,    &ip3,
	       typvar,  nomvar,  etiket,
	       grtyp,   &ig1,    &ig2,    &ig3, &ig4,
	       &swa,    &lng,    &dltf,   &ubc,
	       &extra1, &extra2, &extra3) < 0 ) {
    printf("(Cvgd) ERROR: cannot fstprm for fstkey %d\n",key);
    return(1);
  }
  if( Cvgd_new_read(&vgd, iun, ip1, ip2, 5, 5) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun\n");
    return(1);
  }
  if( Cvgd_new_read2(&vgd2, iun, ip1, ip2, 5, 5, 0) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read2 on iun\n");
    return(1);
  }
  if( Cvgd_vgdcmp(vgd, vgd2) != 0 ){
    printf("ERROR, vgd and vgd2 shouldne the same with Cvgd_new_read2\n");
    return(1);
  }

  // Test that we can read with all selection parameter and get equality with Cvgd_new_read
  if( Cvgd_new_read3(&vgd2, iun, 0, "STG_CP_GEMV4", ip1, ip2, 0, 5, 5, 0) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun\n");
    return(1);
  }
  if( Cvgd_vgdcmp(vgd, vgd2) != 0 ){
    printf("ERROR, vgd and vgd2 should be the same with Cvgd_new_read3\n");
    return(1);
  }

  Cvgd_free(&vgd2);

  // Test wrong datev
  printf("The following error message is normal since testing wrong datev\n");
  if( Cvgd_new_read3(&vgd2, iun, 354514400, etiket, ip1, ip2, 0, 5, 5, 0) == VGD_OK ) {
    printf("ERROR with Cvgd_new_read3 should not find !! and did\n");
    return(1);
  }

  // Test wrong etiket
  printf("The following error message is normal since testing wrong etiket\n");
  if( Cvgd_new_read3(&vgd2, iun, 0,"WRONG ETIKET", ip1, ip2, 0, 5, 5, 0) == VGD_OK ) {
    printf("ERROR with Cvgd_new_read3 should not find !! and did\n");
    return(1);
  }

  // Test wrong ip1
  printf("The following error message is normal since testing wrong ip1\n");
  if( Cvgd_new_read3(&vgd2, iun, 0,etiket, ip1+1, ip2, 0, 5, 5, 0) == VGD_OK ) {
    printf("ERROR with Cvgd_new_read3 should not find !! and did\n");
    return(1);
  }

  // Test wrong ip2
  printf("The following error message is normal since testing wrong ip2\n");
  if( Cvgd_new_read3(&vgd2, iun, 0,etiket, ip1, ip2+1, 0, 5, 5, 0) == VGD_OK ) {
    printf("ERROR with Cvgd_new_read3 should not find !! and did\n");
    return(1);
  }

  // Test wrong ip3
  printf("The following error message is normal since testing wrong ip3\n");
  if( Cvgd_new_read3(&vgd2, iun, 0,etiket, ip1, ip2, 1, 5, 5, 0) == VGD_OK ) {
    printf("ERROR with Cvgd_new_read3 should not find !! and did\n");
    return(1);
  }

  Cvgd_free(&vgd);
  Cvgd_free(&vgd2);

  return(c_ut_report(status,"testing Cvgd_new_read*"));
  
}
