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

// Do not let write_control set to 1 since this will invalidate the test.
// Just use this to produce the validation data. To do so, run the program
//  with write_control set to 1, copy the files
// data/*_stda76_temp.txt to data dir
#define write_control 0


int main() {

  FILE *fp;
  int ier, nk, key, ij, nl, nl_c, in_log, k;
  int iun = 10, quiet=0, *i_val = NULL;
  float me, *pres, *hgts, ff;
  char mode[]="RND";
  char *filename = "data/dm_21001_from_model_run_NON_SLEVE";
  char *filename_c ="data/c_stda76_pres_from_ghts.txt";
  char buff[255];
  vgrid_descriptor *vgd = NULL;

  // Get any heights vertical descriptor
  if(c_fnom(&iun,filename,mode,0) < 0 ) {
    printf("ERROR with c_fnom on iun, file %s\n", filename);
    return(1);
  }
  if(c_fstouv(iun,"RND") < 0 ) {
    printf("ERROR with c_fstouv on iun, file %s\n", filename);
    return(1);
  }
  if( Cvgd_new_read(&vgd, iun, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read\n");
    return(1);
  }
  if( Cvgd_get_int(vgd, "NL_M", &nl, quiet) == VGD_ERROR){
    printf("ERROR cannot Cvgd_get_int on NL_M\n");
    return(1);
  }
  printf("Found %d heights\n", nl);
  
  i_val = malloc(nl * sizeof(int));
  if(! i_val){
    printf("Problem allocating i_val of size %d\n",nl);
    return(1);
  }
  if(Cvgd_get_int_1d(vgd, "VIPM", &i_val, NULL, quiet) ==  VGD_ERROR ) {
    printf("ERROR with Cvgd_get_int for VIPM\n");
    return(1);
  }
  //for(k=0; k<nl; k++){
  //  printf("i_val[k] = %d\n",i_val[k]);
  //}


  me = -100.;
  hgts = malloc(nl * sizeof(float));
  if(! hgts){
    printf("Problem allocating hgts of size %d\n",nl);
    return(1);
  }
  if(Cvgd_levels(vgd, 1, 1, nl, i_val, hgts, &me, 0) == VGD_ERROR){
    printf("Problem Computing heights");
    return(1);
  }
  //for(k=0; k<nl; k++){
  //  printf("ghts[k]. = %f\n",hgts[k]);
  //}
  pres = malloc(nl * sizeof(float));
  if(! pres){
    printf("Problem allocating pres of size %d\n",nl);
    return(1);
  }
  if(Cvgd_stda76_pres_from_hgts_list(pres, hgts, nl)
     == VGD_ERROR){
    printf("Problem Computing pres from heights value");
    return(1);
  }
  //for(k=0; k<nl; k++){
  //  printf("k = %d, pres[k] = %f, hgts[k] = %f\n", k, pres[k], hgts[k]);
  //}
  if(write_control){
    fp = fopen(filename_c, "w");
    fprintf(fp, "size %d\n",nl);
    for( k=0; k<nl; k++){
      fprintf(fp, "%f\n", pres[k]);
    }
    fclose(fp);
  }
  fp = fopen(filename_c, "r");
  if( fp == NULL ){
    printf("ERROR in test, validation file not present %s\n",filename_c);
  }
  fscanf(fp, "%s %d",buff, &nl_c);
  if( nl_c != nl){
    printf("ERROR in tests, size problem with validation file\n");
    return(1);
  }
  for( k=0; k<nl; k++){
    fscanf(fp,"%f", &ff);
    if( fabs(ff - pres[k]) > .01f ){
      printf("ERROR differences found, expecting: %f, got %f\n", ff, pres[k]);
      printf("ERROR TEST on file %s failled\n", filename_c);
      return(1);
    }
  }
  printf("The following error is expected\n");
  hgts[0]=84852. + 1.;
  if(Cvgd_stda76_pres_from_hgts_list(pres,hgts,nl)
     == VGD_OK){
    printf("Problem: call sould have produce a bound error by did not\n");
    return(1);
  }
  
  printf("All tests OK");
  return(c_ut_report(VGD_OK,"testing Cvgd_stda76_pres_from_hgts_list"));
}
