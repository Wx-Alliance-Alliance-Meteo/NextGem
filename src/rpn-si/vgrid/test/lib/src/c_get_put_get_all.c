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
// Do not let write_control set to 1 since this will invalidate the test. Just use this to produce the 
// validation data. To do so, run the program with write_control set to 1, copy the files
// data/*_stda76_temp.txt to data dir
#define write_control 0

char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the null-terminator
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

int test_it(char *filename, int ind) {
  int i, iun, ier;
  int quiet=0;
  int value=0, value2=0;
  char mode[]="RND";
  vgrid_descriptor *vgd = NULL;
  char key[11][5] = {"DATE","IG_1","IG_2","IG_3","IG_4","IP_1","IP_2","IP_3","DIPM", "DIPT", "DIPW"};

  iun = 10 + ind;
  
  if( c_fnom(&iun,filename,mode,0) < 0 ) {
    printf("ERROR in test with c_fnom on iun, file %s\n", filename);
    return(VGD_ERROR);
  }
  if( c_fstouv(iun,"RND") < 0 ) {
    printf("ERROR in test with c_fstouv on iun, file %s\n", filename);
    return(VGD_ERROR);
  }

  if( Cvgd_new_read(&vgd, iun, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR in test with Cvgd_new_read on iun\n");
    return(VGD_ERROR);
  }
  //ier = Cvgd_print_desc(vgd, -1, -1);
  for( i = 0; i < 11; i++ ){
    printf("   %s\n", key[i]);
    if (strcmp(key[i], "DIPM") == 0){
      if(Cvgd_is_valid(vgd,"dhm_valid")){
	if( Cvgd_get_int(vgd, key[i], &value, quiet) == VGD_ERROR ) {
	  printf("ERROR in test with Cvgd_get_int on key %s\n", key[i]);
	  return(VGD_ERROR);
	}
      } else {
	continue;
      }
    }
    if (strcmp(key[i], "DIPT") == 0){
      if(Cvgd_is_valid(vgd,"dht_valid")){
	if( Cvgd_get_int(vgd, key[i], &value, quiet) == VGD_ERROR ) {
	  printf("ERROR in test with Cvgd_get_int on key %s\n", key[i]);
	  return(VGD_ERROR);
	}
      } else {
	continue;
      }
    }
    if (strcmp(key[i], "DIPW") == 0){
      if(Cvgd_is_valid(vgd,"dhw_valid")){
	if( Cvgd_get_int(vgd, key[i], &value, quiet) == VGD_ERROR ) {
	  printf("ERROR in test with Cvgd_get_int on key %s\n", key[i]);
	  return(VGD_ERROR);
	}
      } else {
	continue;
      }
    }    
    value2 = value;
    value = value + 1;
    if( Cvgd_put_int(&vgd, key[i] , value) == VGD_ERROR ){
      return(VGD_ERROR);
    }
    if( Cvgd_get_int(vgd, key[i] , &value2, quiet) == VGD_ERROR ){
      return(VGD_ERROR);
    }
    if ( value2 != value ){
      printf("ERROR in test with key %s\n",key[i]);
      return(VGD_ERROR);
    }
  }
  
  ier = c_fstfrm(iun);
  ier = c_fclos(iun);
  printf("ier = %d\n",ier);
      
  Cvgd_free(&vgd);

  return(VGD_OK);
  
}

//========================================================================
//========================================================================

int main() {
  
  int ier, i, status = VGD_OK;

  if( Cvgd_putopt_int("ALLOW_SIGMA",1) == VGD_ERROR){
    printf("ERROR with Cvgd_putopt_int on ALLOW_SIGMA)\n");
    status = VGD_ERROR;
    return(1);
  }
  
  for (i = 0; i < (int) n_file; i++) {
    printf ("Testing %s\n", filenames[i]);
    if(test_it(filenames[i],i) == VGD_ERROR){
      printf("ERROR with %s\n",filenames[i]);
      status = VGD_ERROR;
      return(1);
    }
  }  

  return(c_ut_report(status,"testing Cvgd_stda76"));  
}
