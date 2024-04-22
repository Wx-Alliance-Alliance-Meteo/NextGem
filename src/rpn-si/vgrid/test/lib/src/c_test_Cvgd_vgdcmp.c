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
#include "vgrid.h"
#include "rmn.h"
#include "c_ut_report.h"

int main() {

  int ier, iun = 10, iun2 = 11, status;
  char filename[]="data/dm_5005_from_model_run";
  char filename2[]="data/dm_5001_from_model_run";
  char mode[]="RND";
  char format[] = "FST";
  char name[5];
  vgrid_descriptor *vgd = NULL, *vgd2 = NULL;
  vgrid_descriptor *vgd3 = NULL, *vgd4 = NULL;

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
  if( Cvgd_new_read(&vgd2, iun, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun\n");
    return(1);
  }

  ier = c_fnom(&iun2,filename2,mode,0);
  if( ier < 0 ) {
    printf("ERROR with c_fnom on iun2, file %s\n", filename2);
    return(1);
  }
  ier = c_fstouv(iun2,"RND");  
  if( ier < 0 ) {
    printf("ERROR with c_fstouv on iun2, file %s\n", filename2);
    return(1);
  }
  if( Cvgd_new_read(&vgd3, iun2, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun2\n");
    return(1);
  }
  if( Cvgd_new_read(&vgd4, iun2, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun2\n");
    return(1);
  }

  // Comparing all
  if( Cvgd_vgdcmp(vgd, vgd2) != 0 ){
    printf("ERROR, vgd and vgd2 should be the same\n");
    return(1);
  }

  vgd2->vcode=0;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to vcode\n");
    status = VGD_ERROR;
  }
  vgd2->vcode=vgd->vcode;
  
  vgd2->kind=0;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to kind\n");
    status = VGD_ERROR;
  }
  vgd2->kind=vgd->kind;

  vgd2->version=0;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to version\n");
    status = VGD_ERROR;
  }
  vgd2->version=vgd->version;
 
  strncpy(vgd2->ref_name,"NULL",4);
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to ref_name\n");
    status = VGD_ERROR;
  }
  strncpy(vgd2->ref_name,vgd->ref_name,4);

  vgd4->ptop_8=0.;
  if( Cvgd_vgdcmp(vgd3, vgd4) == 0 ){
    printf("ERROR, vgd3 and vgd4 should not be the same due to ptop_8\n");
    status = VGD_ERROR;
  }
  vgd4->ptop_8=vgd3->ptop_8;

  vgd4->pref_8=0.;
  if( Cvgd_vgdcmp(vgd3, vgd4) == 0 ){
    printf("ERROR, vgd3 and vgd4 should not be the same due to pref_8\n");
    status = VGD_ERROR;
  }
  vgd4->pref_8=vgd3->pref_8;

  vgd4->rcoef1=0.;
  if( Cvgd_vgdcmp(vgd3, vgd4) == 0 ){
    printf("ERROR, vgd3 and vgd4 should not be the same due to rcoef1\n");
    status = VGD_ERROR;
  }
  vgd4->rcoef1=vgd3->rcoef1;

  vgd4->rcoef2=0.;
  if( Cvgd_vgdcmp(vgd3, vgd4) == 0 ){
    printf("ERROR, vgd3 and vgd4 should not be the same due to rcoef2\n");
    status = VGD_ERROR;
  }
  vgd4->rcoef2=vgd3->rcoef2;

  vgd2->ip1_m[vgd2->nl_m-1]=vgd->ip1_m[vgd->nl_m-1]+1.;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to ip1_m\n");
    status = VGD_ERROR;
  }
  vgd2->ip1_m[vgd2->nl_m-1]=vgd->ip1_m[vgd2->nl_m-1];
  
  vgd2->ip1_t[vgd2->nl_t-1]=vgd->ip1_t[vgd->nl_t-1]+1.;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to ip1_t\n");
    status = VGD_ERROR;
  }
  vgd2->ip1_t[vgd2->nl_t-1]=vgd->ip1_t[vgd2->nl_t-1];
  
  vgd2->a_m_8[vgd2->nl_m-1]=vgd->a_m_8[vgd->nl_m-1]+1.;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to a_m_8\n");
    status = VGD_ERROR;
  }
  vgd2->a_m_8[vgd2->nl_m-1]=vgd->a_m_8[vgd2->nl_m-1];

  vgd2->b_m_8[vgd2->nl_m-1]=vgd->b_m_8[vgd->nl_m-1]+1.;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to b_m_8\n");
    status = VGD_ERROR;
  }
  vgd2->b_m_8[vgd2->nl_m-1]=vgd->b_m_8[vgd2->nl_m-1];

  vgd2->a_t_8[vgd2->nl_t-1]=vgd->a_t_8[vgd->nl_t-1]+1.;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to a_t_8\n");
    status = VGD_ERROR;
  }
  vgd2->a_t_8[vgd2->nl_t-1]=vgd->a_t_8[vgd2->nl_t-1];

  vgd2->b_t_8[vgd2->nl_t-1]=vgd->b_t_8[vgd->nl_t-1]+1.;
  if( Cvgd_vgdcmp(vgd, vgd2) == 0 ){
    printf("ERROR, vgd and vgd2 should not be the same due to b_t_8\n");
    status = VGD_ERROR;
  }
  vgd2->b_t_8[vgd2->nl_t-1]=vgd->b_t_8[vgd2->nl_t-1];

  Cvgd_free(&vgd);
  Cvgd_free(&vgd2);

  return(c_ut_report(status,"testing Cvgd_levels"));
}
