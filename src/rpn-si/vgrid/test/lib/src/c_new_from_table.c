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
#include "vgrid.h"
#include "rmn.h"
#include "c_ut_report.h"

char *filenames[] = {
    "data/dm_1001_from_model_run",
    "data/dm_1002_from_model_run",
    "data/2001_from_model_run",
    "data/dm_5001_from_model_run",
    "data/dm_5002_from_model_run",
    "data/dm_5005_from_model_run",
    "data/dm_5100_from_model_run",
    "data/dm_5999_from_model_run",
    "data/dm_21001_from_model_run_SLEVE",
    "data/dm_21001_from_model_run_NON_SLEVE",
    "data/dm_21002_from_model_run_SLEVE",
    "data/dm_21002_from_model_run_NON_SLEVE"
};

#define n_file (sizeof (filenames) / sizeof (const char *))

int test_it(char *filename, int ind) {

  int ni, nj, nk, iun, ier;
  double *table = NULL;
  vgrid_descriptor *vgd = NULL,  *vgd2 = NULL;

  iun = 10 + ind;

  if( c_fnom(&iun,filename,"RND+R/O",0) < 0 ) {
    printf("ERROR with c_fnom on iun, file %s\n", filename);
    return(1);
  }
  if( c_fstouv(iun,"RND") < 0 ) {
    printf("ERROR with c_fstouv on iun, file %s\n", filename);
    return(1);
  }  
  if( Cvgd_new_read(&vgd, iun, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun\n");
    return(1);
  }
  // Get table  
  if( Cvgd_get_double_3d(vgd, "VTBL", &table, &ni, &nj, &nk, 0) == VGD_ERROR ){
    printf("ERROR with Cvgd_get_double_3d on VTBL\n");
    return(1);
  }
  if( Cvgd_new_from_table(&vgd2, table, ni, nj, nk) == VGD_ERROR ){
    printf("ERROR with Cvgd_new_from_table\n");
    return(1);
  }
  // Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);    
    return (VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  Cvgd_free(&vgd);
  Cvgd_free(&vgd2);
  free(table);
  return (VGD_OK);
}

int main() {
  
  int i, ier, status = VGD_OK;

  ier = Cvgd_putopt_int("ALLOW_SIGMA",1);
  
  for (i = 0; i < (int) n_file; i++) {
    printf ("Testing %s\n", filenames[i]);
    if(test_it(filenames[i],i) == VGD_ERROR){
      printf("ERROR with %s\n",filenames[i]);
      status = VGD_ERROR;
      return(1);
    }
  }  
  printf("status=%d, VGD_OK=%d, VGD_ERROR=%d\n",status, VGD_OK, VGD_ERROR);
  return(c_ut_report(status,"testing Cvgd_levels"));  
  
}
