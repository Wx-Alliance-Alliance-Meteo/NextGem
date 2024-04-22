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
#include <string.h>
#include <math.h>
#include "vgrid.h"
#include "rmn.h"
#include "c_ut_report.h"

char *filenames[] = {
    "data/dm_1001_from_model_run",
    "data/dm_1002_from_model_run",
    "data/2001_from_model_run",
    "data/dm_4001_from_model_run",
    "data/dm_5001_from_model_run",
    "data/dm_5005_from_model_run",
    "data/dm_5100_from_model_run",
    "data/dm_5999_from_model_run",
    "data/dm_21001_from_model_run_SLEVE",
    "data/dm_21001_from_model_run_NON_SLEVE",
    "data/dm_21002_from_model_run_SLEVE",
    "data/dm_21002_from_model_run_NON_SLEVE"
 };


#define n_file (sizeof (filenames) / sizeof (const char *))

int compare_values(int iun, int vcode, float *levels, double *levels_8, float *p0, int *i_val, int nl, char *ref_nomvar){
  int k, ij, ijk, ni2, nj2, nk2, key, ier, ip1;  
  char nomvar[5];
  float fact;
  // Load PX/GZ to see if pressure/height computation is OK
  strcpy(nomvar,"PX  ");  
  fact=0.01f;
  if(! strcmp(ref_nomvar,"ME  ")){
    strcpy(nomvar,"GZ  ");
    fact=.1f;
  }
  if(vcode == 4001 ){
    strcpy(nomvar,"GZ  ");
    fact=.1f;
  }

  printf("nl=%d\n",nl);
  for( k = 0, ijk = 0; k < nl; k++){
    // Il n'est pas normal que GZ 0 m ne soit pas dans le fichier pour Vcode 21001
    ip1=i_val[k];
    if(i_val[k] == 82837504 && ! strcmp(ref_nomvar,"ME  ")){
      ip1=93423364;
    }
    key = c_fstinf( iun, &ni2, &nj2, &nk2, -1, " ", ip1, -1, -1, " ", nomvar);
    if(key < 0){
      printf("Problem getting info for %s for ip1 = %d\n",nomvar,ip1);
      return(VGD_ERROR);
    }
    // To simplify, PX are assumed to be on the same grid as P0. but this should be check in an operational program!
    ier = c_fstluk((uint32_t*)p0, key, &ni2, &nj2, &nk2 );
    if( ier < 0 ){
      printf("Error with c_fstluk on level = %d\n",i_val[k]);
      return(VGD_ERROR);
    }
    
    for( ij = 0; ij < ni2*nj2; ij++, ijk++){
      if(abs(p0[ij]) < 1.e-6){
	if(fabs(p0[ij] - levels[ijk]*fact) > 1.e-6 ){
	  printf("Difference is too large (float), expected %f, got %f\n", p0[ij], levels[ijk]*fact);
	  return(VGD_ERROR);
	}
	if(fabs(p0[ij] - levels_8[ijk]*fact) > 1.e-6 ){
	  printf("Difference is too large (double), expected %f, got %f\n", p0[ij], levels[ijk]*fact);
	  return(VGD_ERROR);
	}
      } else {
	if( fabs(p0[ij] - levels[ijk]*fact)/abs(p0[ij]) >  1.e-6 ) {
	  printf("Difference is too large (float), expected %f, got %f\n", p0[ij], levels[ijk]*fact);
	  return(VGD_ERROR);
	}
	if( fabs(p0[ij] - levels_8[ijk]*fact)/abs(p0[ij]) > 1.e-6 ) {
	  printf("Difference is too large (double), expected %f, got %f\n", p0[ij], levels_8[ijk]*fact);
	  return(VGD_ERROR);
	}
      }
    }
  }
  return(VGD_OK);
}

int test_it(char *filename, char *ip1_name, int ind) {

  int ier, iun, vcode;
  int quiet=0, *i_val = NULL, in_log = 0, dpidpis = 0;
  int nl, ni, nj, nk, ni2, nj2, nk2, key, ij, ijk, kind, ref1, ref2;
  char mode[]="RND", key_name[]="1234";
  char nomvar1[5], nomvar2[5];
  float *p0 = NULL, *p0ls = NULL, *levels = NULL, fact;
  double *p0_8 = NULL, *p0ls_8 = NULL, *levels_8 = NULL;
  vgrid_descriptor *vgd = NULL;
      
  iun = 10 + ind;
  
  ier = c_fnom(&iun,filename,mode,0);
  if( ier < 0 ) {
    printf("ERROR with c_fnom on iun, file %s\n", filename);
    return(VGD_ERROR);
  }
  ier = c_fstouv(iun,"RND");  
  if( ier < 0 ) {
    printf("ERROR with c_fstouv on iun, file %s\n", filename);
    return(VGD_ERROR);
  }

  if (strcmp(ip1_name, "VIPM") == 0) {
    strcpy(key_name,"NL_M");
  } else if (strcmp(ip1_name, "VIPT") == 0) {
    strcpy(key_name,"NL_T");
  }else if (strcmp(ip1_name, "VIPW") == 0) {
    strcpy(key_name,"NL_W");
  } else {
    printf("Wrong name passed to ip1_name %s\n",key_name);
    return(VGD_ERROR);
  }
      
  if( Cvgd_new_read(&vgd, iun, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun\n");
    return(VGD_ERROR);
  }
  //ier = Cvgd_print_desc(vgd, -1, -1);

  if( Cvgd_get_int_1d(vgd, ip1_name, &i_val, NULL, quiet) ==  VGD_ERROR ) {
    printf("ERROR with Cvgd_get_int for %s\n",ip1_name);
    return(VGD_ERROR);
  }
  
  ier = Cvgd_get_int(vgd, key_name, &nl, quiet);
  if(ier == VGD_ERROR){
    printf("ERROR cannot Cvgd_get_int on %s\n",key_name);
    return(VGD_ERROR);
  }
  //for(ijk=0; ijk<nl; ijk++){
  //  printf("i_val[ijk]=%d\n",i_val[ijk]);
  //}

  // Compute 3D pressure levels or heights
  ier = Cvgd_get_char(vgd, "RFLD", nomvar1, 1);
  ref1 = strcmp(nomvar1,VGD_NO_REF_NOMVAR) ? 1 : 0;
  ier = Cvgd_get_char(vgd, "RFLS", nomvar2, 1);
  ref2 = strcmp(nomvar2,VGD_NO_REF_NOMVAR) ? 1 : 0;

  if(! ref1){
    // if no ref1 test if it kind == 2 (pressure level)
    // If is is then put nomvar1 to TT to get problem size. but the actual TT value will not be used
    ier = Cvgd_get_int(vgd, "KIND", &kind, 0);
    if(ier == VGD_ERROR){
      printf("ERROR cannot Cvgd_get_int on KIND\n");
      return(VGD_ERROR);
    }
    if( kind != 2 && kind != 4 ){
      printf("ERROR in test: Expecting kind = 2 (pressure) or 4 (HAGL) since no RFLD but got kind = %d\n",kind);
      return(VGD_ERROR);
    }
    if( kind == 2 ){
      strcpy(nomvar1,"TT  ");
    }
    if( kind == 4 ){
      strcpy(nomvar1,"GZ  ");
    }
  }

  key = c_fstinf( iun, &ni2, &nj2, &nk2, -1, " ", -1, -1, -1, " ", nomvar1);
  if(key < 0){
    printf("Problem getting info for %s", nomvar1);
    return(VGD_ERROR);
  }
  p0 = malloc(ni2*nj2 * sizeof(float));
  if(! p0){
    printf("Problem allocating surface reference field of size %d\n",ni2*nj2);
    return(VGD_ERROR);
  }
  p0_8 = malloc(ni2*nj2 * sizeof(double));
  if(! p0_8){
    printf("Problem allocating double surface reference field of size %d\n",ni2*nj2);
    return(VGD_ERROR);
  }
  ier = c_fstluk((uint32_t*)p0, key, &ni2, &nj2, &nk2 );
  if(ier < 0){
    printf("Problem with fstluk for reference field\n");
    return(VGD_ERROR);
  }
  fact=100.f;
  if(! strcmp(nomvar1,"ME  ")){
    fact=1.f;
  }
  for( ij = 0; ij < ni2*nj2; ij++, ijk++){    
    p0[ij] = p0[ij]*fact;
    p0_8[ij] = p0[ij];
  }

  if(ref2){
    key = c_fstinf( iun, &ni, &nj, &nk, -1, " ", -1, -1, -1, " ", nomvar2);
    if(key < 0){
      printf("Problem getting info for %s\n", nomvar2);
      return(VGD_ERROR);
    }
    if(ni != ni2 || nj != nj2){
      printf("Incompatible size of %s\n",nomvar2);
      return(VGD_ERROR);
    }
    p0ls = malloc(ni2*nj2 * sizeof(float));
    if(! p0ls){
      printf("Problem allocating p0ls of size %d\n",ni2*nj2);
      return(VGD_ERROR);
    }
    p0ls_8 = malloc(ni2*nj2 * sizeof(double));
    if(! p0ls_8){
      printf("Problem allocating p0ls_8 of size %d\n",ni2*nj2);
      return(VGD_ERROR);
    }
    ier = c_fstluk((uint32_t*)p0ls, key, &ni2, &nj2, &nk2 );
    if(ier < 0){
      printf("Problem with fstluk for p0\n");
      return(VGD_ERROR);
    }
    for( ij = 0; ij < ni2*nj2; ij++, ijk++){
      p0ls[ij] = p0ls[ij]*fact;
      p0ls_8[ij] = p0ls[ij];
    }    
  }

  levels = malloc(ni2*nj2*nl * sizeof(float));
  if(! levels){
    printf("Problem allocating levels of size %d\n",ni2*nj2*nl);
    return(VGD_ERROR);
  }
  levels_8 = malloc(ni2*nj2*nl * sizeof(double));
  if(! levels_8){
    printf("Problem allocating levels_8 of size %d\n",ni2*nj2*nl);
    return(VGD_ERROR);
  }

  //===================================================
  printf("   testing Cvgd_levels* float interface for ip1 list %s\n",ip1_name);
  if(ref2){
    ier = Cvgd_levels_2ref(vgd, ni2, nj2, nl, i_val, levels, p0, p0ls, in_log);
  } else {
    ier = Cvgd_levels(vgd, ni2, nj2, nl, i_val, levels, p0, in_log);
  }
  if(ier == VGD_ERROR){
    printf("Error with Cvgd_levels*\n");
    return(VGD_ERROR);
  }
  printf("   testing Cvgd_levels* double interface for ip1 list %s\n",ip1_name);
  if(ref2){
    ier = Cvgd_levels_2ref_8(vgd, ni2, nj2, nl, i_val, levels_8, p0_8, p0ls_8, in_log);
  } else {
    ier = Cvgd_levels_8(vgd, ni2, nj2, nl, i_val, levels_8, p0_8, in_log);
  }
  if(ier == VGD_ERROR){
    printf("Error with Cvgd_levels*_8\n");
    return(VGD_ERROR);
  }

  //===================================================
  printf("   testing Cvgd_diag_withref float interface\n");
  if(ref2){
    ier = Cvgd_diag_withref_2ref(vgd, ni2, nj2, nl, i_val, levels, p0, p0ls, in_log, dpidpis);
  } else {
    ier = Cvgd_diag_withref(vgd, ni2, nj2, nl, i_val, levels, p0, in_log, dpidpis);
  }
  if(ier == VGD_ERROR){
    printf("Error with Cvgd_diag_withref*\n");
    return(VGD_ERROR);
  }
  printf("   testing Cvgd_levels* double interface for ip1 list %s\n",ip1_name);
  if(ref2){
    ier = Cvgd_diag_withref_2ref_8(vgd, ni2, nj2, nl, i_val, levels_8, p0_8, p0ls_8, in_log, dpidpis);
  } else {
    ier = Cvgd_diag_withref_8(vgd, ni2, nj2, nl, i_val, levels_8, p0_8, in_log, dpidpis);
  }
  if(ier == VGD_ERROR){
    printf("Error with Cvgd_diag_withref*_8\n");
    return(VGD_ERROR);
  }
  strcpy(key_name,"VCOD");
  ier = Cvgd_get_int(vgd, key_name, &vcode, quiet);
  if( compare_values(iun, vcode, levels, levels_8, p0, i_val, nl,nomvar1) == VGD_ERROR){
    printf("ERROR comparison failed\n");
    return(VGD_ERROR);
  }

  ier = c_fstfrm(iun);
  ier = c_fclos(iun);
      
  Cvgd_free(&vgd);
  free(p0);
  free(levels);
  free(p0_8);
  free(levels_8);
  free(i_val);

  return(VGD_OK);
  
}

//========================================================================
//========================================================================

int main() {
  
  int i, ier, status = VGD_OK;

  ier = Cvgd_putopt_int("ALLOW_SIGMA",1);
  
  for (i = 0; i < (int) n_file; i++) {
    printf ("Testing %s\n", filenames[i]);
    if(test_it(filenames[i],"VIPM",i) == VGD_ERROR){
      printf("ERROR with %s\n",filenames[i]);
      status = VGD_ERROR;
      return(1);
    }
    if(test_it(filenames[i],"VIPT",i) == VGD_ERROR){
      printf("ERROR with %s\n",filenames[i]);
      status = VGD_ERROR;
      return(1);
    }
    if(test_it(filenames[i],"VIPW",i) == VGD_ERROR){
      printf("ERROR with %s\n",filenames[i]);
      status = VGD_ERROR;
      return(1);
    }
  }  
  printf("status=%d, VGD_OK=%d, VGD_ERROR=%d\n",status, VGD_OK, VGD_ERROR);
  return(c_ut_report(status,"testing Cvgd_levels"));  
}
