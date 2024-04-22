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

int main() {

  int ier, iun = 10;
  int quiet = 0, *i_val = NULL, in_log = 0;
  int ind, k, nl_t, status;
  char filename[]="data/dm_5005_from_model_run";
  char mode[]="RND";
  float *levels = NULL, p0;
  vgrid_descriptor *vgd = NULL;

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

  ier = Cvgd_get_int(vgd, "NL_T", &nl_t, quiet);
  if(ier == VGD_ERROR){
    status = VGD_ERROR;
  }

  levels = malloc(nl_t * sizeof(float));
  if(! levels){
    printf("Problem allocating levels of size %d\n",nl_t);
    return(1);
  }
  p0=1000*100.;
  ier = Cvgd_levels(vgd, 1, 1, nl_t, i_val, levels, &p0, in_log);
  if(ier == VGD_ERROR){
    printf("Error with Cvgd_diag_withref\n");
    status = VGD_ERROR;
  }
  // Find the first hyb level at or above 70000. hpa
  ind=-1;
  for(k=0;k<nl_t;k++){
    printf("k = %d, press = %f hPa\n",k,levels[k]);
    if(ind == -1 && levels[k] >= 70000.){
      ind = k;
    }
  }
  printf("First level at or above 70000 is at ind = %d, hyb = %d, and at pressure = %f hPa\n",ind,i_val[ind],levels[ind]);

  ier = c_fstfrm(iun);
  ier = c_fclos(iun);

  Cvgd_free(&vgd);
  free(levels);
  free(i_val);
  
  return(c_ut_report(status,"testing Cvgd_levels"));
}
