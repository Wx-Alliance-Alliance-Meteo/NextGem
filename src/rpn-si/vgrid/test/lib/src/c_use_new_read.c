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
#include "c_ut_report.h"

int main() {

  int ier, iun = 10, iun2 = 11;
  int quiet = 0, *i_val = NULL, in_log = 0, dpidpis = 0;
  int nl_t, nt, ni, nj, nk, ni2, nj2, nk2, k, key, ij, ijk, status;
  char filename[]="data/dm_5005_from_model_run";
  char mode[]="STD+RND";
  char format[] = "FST";
  char name[5];
  float *f_val = NULL, *p0 = NULL, *px = NULL;
  double *a_8_t = NULL, *b_8_t = NULL, *table = NULL, *levels_8 = NULL, *p0_8 = NULL;
  vgrid_descriptor *vgd = NULL, *vgd2 = NULL;

  status = VGD_OK;

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
  
  if( Cvgd_new_read(&vgd, iun, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun\n");
    return(1);
  }
  if( Cvgd_get_int_1d(vgd, "VIPT", &i_val, NULL, quiet) ==  VGD_ERROR ) {
    printf("ERROR with Cvgd_get_int for VIPT\n");
    return(1);
  }
  if( Cvgd_get_float_1d(vgd, "VCDT", &f_val, NULL , quiet) ==  VGD_ERROR ) {
    printf("ERROR with Cvgd_get_float_1d for VCDT\n");
    return(1);
  }
  if( Cvgd_get_double_1d(vgd, "CA_T", &a_8_t, NULL, quiet) ==  VGD_ERROR ) {
    printf("ERROR with Cvgd_get_double_1d for CA_T\n");
    return(1);
  }
  if( Cvgd_get_double_1d(vgd, "CB_T", &b_8_t, &nt, quiet) ==  VGD_ERROR ) {
    printf("ERROR with Cvgd_get_double_1d for CB_T\n");
    return(1);
  }

  // Size of thermo may also be obtained by this:
  ier = Cvgd_get_int(vgd, "NL_T", &nl_t, quiet);
  if(nl_t != nt ) {
    printf("ERROR: nt and nl_t should be equal, got %d, %d\n",nt, nl_t);
    return(1);
  }
  printf("nl_t = %d\n", nl_t);
  
  for( k = 0; k < nl_t; k++) {
    printf("k = %d, ip1 = %d, val = %f, A = %f, B = %f\n",
  	   k, i_val[k], f_val[k], a_8_t[k], b_8_t[k]);
  }

  // Load table (this is the actual data in fst record !! which may also be
  // obtained with fstlir, but why do it if vgd already contains it!)
  if ( Cvgd_get_double_3d(vgd, "VTBL", &table, &ni, &nj, &nk, quiet) ==  VGD_ERROR ) {
    printf("ERROR with Cvgd_get_double_3d for VTBL\n");
    return(1);
  }
  
  // Constructing new vgd with this table
  if ( Cvgd_new_from_table(&vgd2, table, ni, nj, nk) ==  VGD_ERROR ) {
    printf("ERROR with Cvgd_new_from_table for VTBL\n");
    return(1);
  }

  // Comparing new table with original table, must be the same.
  if( Cvgd_vgdcmp(vgd, vgd2) != 0 ){
    printf("ERROR, vgd and vgd2 shouldne the same\n");
    return(1);
  }

  // Write descriptor in new file  
  char command[50];
  strcpy( command, "rm -f to_erase" );
  system(command);
  ier = c_fnom(&iun2,"to_erase",mode,0);
  if( ier < 0 ) {
    printf("ERROR with c_fnom on iun2\n");
    return(1);
  }
  ier = c_fstouv(iun2,"RND");  
  if( ier < 0 ) {
    printf("ERROR with c_fstouv on iun2\n");
    return(1);
  }
  if( Cvgd_write_desc(vgd, iun2) == VGD_ERROR ){
    printf("ERROR with Cvgd_write_desc on iun2\n");
    return(1);
  }

  // Compute 3D pressure levels
  // First get P0 (the reference field
  key = c_fstinf( iun, &ni2, &nj2, &nk2, -1, " ", -1, -1, -1, " ", "P0");
  if(key < 0){
    printf("Problem getting info for P0\n");
    return(1);
  }
  p0 = malloc(ni2*nj2 * sizeof(float));
  if(! p0){
    printf("Problem allocating P0 of size %d\n",ni2*nj2);
    return(1);
  }
  p0_8 = malloc(ni2*nj2 * sizeof(double));
  if(! p0_8){
    printf("Problem allocating P0_8 of size %d\n",ni2*nj2);
    return(1);
  }
  levels_8 = malloc(ni2*nj2*nl_t * sizeof(double));
  if(! levels_8){
    printf("Problem allocating levels_8 of size %d\n",ni2*nj2);
    return(1);
  }
  ier = c_fstluk((uint32_t*)p0, key, &ni2, &nj2, &nk2 );
  if( ier < 0 ){
    printf("Problem with fstluk on p0\n");
  }
  for( k = 0; k < ni2*nj2; k++){
    p0_8[k] = p0[k]*100.;
  }
  ier = Cvgd_diag_withref_8(vgd, ni2, nj2, nl_t, i_val, levels_8, p0_8, in_log, dpidpis);
  if(ier == VGD_ERROR){
    status=VGD_ERROR;
  }

  // Load PX to see if pressure computation is OK
  for( k = 0, ijk = 0; k < nl_t; k++){
    key = c_fstinf( iun, &ni2, &nj2, &nk2, -1, " ", i_val[k], -1, -1, " ", "PX");
    if(key < 0){
      printf("Problem getting info for PX for ip1 = %d\n",i_val[k]);
      return(1);
    }
    // To simplify, PX are assumed to be on the same grid as P0. bur rhis should be check in an operational program!
    ier = c_fstluk((uint32_t*)p0, key, &ni2, &nj2, &nk2 );    
    if( ier < 0 ){
      printf("Problem with fstluk on px ip1 = %d\n", i_val[k]);
      return(1);
    }
    for( ij = 0; ij < ni2*nj2; ij++, ijk++){
      if( fabs(p0[ij] - levels_8[ijk]/100.)/p0[ij] > 1.e-6 ) {
	printf("Difference is too large, expected %f, got %f\n", p0[ij], levels_8[ijk]/100.);
	return(1);
      }
    }
  }

  ier = c_fstfrm(iun2);
  ier = c_fclos(iun2);

  // Re open, read and compare
  ier = c_fnom(&iun2,"to_erase",mode,0);
  if( ier < 0 ) {
    printf("ERROR with c_fnom on iun2\n");
    return(1);
  }
  ier = c_fstouv(iun2,"RND");  
  if( ier < 0 ) {
    printf("ERROR with c_fstouv on iun2\n");
    return(1);
  }
  if( Cvgd_new_read(&vgd2, iun, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read vgd2\n");
    return(1);
  }
  if( Cvgd_vgdcmp(vgd, vgd2) != 0 ){
    printf("ERROR, vgd and vgd2 should be the same after write in file, read from file\n");
    return(1);
  }

  Cvgd_free(&vgd);
  Cvgd_free(&vgd2);
  free(table);
  free(i_val);
  free(f_val);
  free(a_8_t);
  free(b_8_t);

  return(c_ut_report(status,"testing Cvgd_levels"));
}
