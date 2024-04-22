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

//========================================================================
int check_build_1001_2001_5999_4001(vgrid_descriptor *vgd, int vcode){
  int kind, version, nk, ier;
  int *ip1_m = NULL;
  double *a_m_8 = NULL, *b_m_8 = NULL;
  vgrid_descriptor *vgd2 = NULL;
  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"NL_M", &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_M", &a_m_8, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_M", &b_m_8, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  printf("  Testing generic interface\n");
  if( Cvgd_new_build_vert(&vgd2, kind, version, nk, -1, -1, NULL, NULL, NULL, NULL,
  			  a_m_8, b_m_8, NULL, NULL, ip1_m, NULL, nk, 0) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  // Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  Cvgd_free(&vgd2);

  printf("  Testing specific interface\n");
  switch(vcode) {
  case 1001:
    if( Cvgd_new_build_vert_1001(&vgd2, -1, -1,
				 a_m_8, b_m_8, ip1_m, nk) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 2001:
    if( Cvgd_new_build_vert_2001(&vgd2, -1, -1,
				 a_m_8, b_m_8, ip1_m, nk) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 4001:
    if( Cvgd_new_build_vert_4001(&vgd2, -1, -1,
				 a_m_8, b_m_8, ip1_m, nk) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 5999:
    if( Cvgd_new_build_vert_5999(&vgd2, -1, -1,
				 a_m_8, b_m_8, ip1_m, nk) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  default:
    printf("In test ERROR in check_build_1001_2001_5999_4001, unsupported Vcode %d\n",vcode);
    return(VGD_ERROR);
  }
  //Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  free(ip1_m);
  free(a_m_8);
  free(b_m_8);
  Cvgd_free(&vgd2);

  return(VGD_OK);
}
//========================================================================
int check_build_1002(vgrid_descriptor *vgd){
  int kind, version, nk, ier;
  int *ip1_m = NULL;
  double *a_m_8 = NULL, *b_m_8 = NULL, ptop_8;
  vgrid_descriptor *vgd2 = NULL;

  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"NL_M", &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_M", &a_m_8, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_M", &b_m_8, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double(vgd,"PTOP", &ptop_8, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  
  printf("  Testing generic interface\n");
  if( Cvgd_new_build_vert(&vgd2, kind, version, nk, -1, -1, &ptop_8, NULL, NULL, NULL,
  			  a_m_8, b_m_8, NULL, NULL, ip1_m, NULL, nk, 0) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  // Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal (generic itf), Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  Cvgd_free(&vgd2);

  if( Cvgd_new_build_vert_1002(&vgd2, -1, -1, ptop_8, a_m_8, b_m_8, ip1_m, nk) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal (specific itf), Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  free(ip1_m);
  free(a_m_8);
  free(b_m_8);
  Cvgd_free(&vgd2);  
  return(VGD_OK);
}

//========================================================================
int check_build_5001(vgrid_descriptor *vgd){
  int kind, version, nk, ier;
  int *ip1_m = NULL;
  float rc_1;
  double *a_m_8 = NULL, *b_m_8 = NULL, ptop_8, pref_8;
  vgrid_descriptor *vgd2 = NULL;

  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"NL_M", &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_M", &a_m_8, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_M", &b_m_8, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double(vgd,"PTOP", &ptop_8, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double(vgd,"PREF", &pref_8, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_1", &rc_1, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  
  printf("  Testing generic interface\n");
  if( Cvgd_new_build_vert(&vgd2, kind, version, nk, -1, -1, &ptop_8, &pref_8, &rc_1, NULL,
  			  a_m_8, b_m_8, NULL, NULL, ip1_m, NULL, nk, 0) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  // Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  Cvgd_free(&vgd2);

  printf("  Testing specific interface\n");
  if( Cvgd_new_build_vert_5001(&vgd2, -1, -1, ptop_8, pref_8, rc_1, 
			       a_m_8, b_m_8, ip1_m, nk) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  // Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  free(ip1_m);
  free(a_m_8);
  free(b_m_8);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}

//========================================================================
int check_build_5002(vgrid_descriptor *vgd){
  int kind, version, nl_m, nl_t, ier;
  int *ip1_m = NULL, *ip1_t = NULL;
  float rc_1, rc_2;
  double *a_m_8 = NULL, *b_m_8 = NULL, *a_t_8 = NULL, *b_t_8 = NULL, ptop_8, pref_8;
  vgrid_descriptor *vgd2 = NULL;

  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_M", &a_m_8, &nl_m, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_M", &b_m_8, &nl_m, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nl_m, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_T", &a_t_8, &nl_t, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_T", &b_t_8, &nl_t, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPT", &ip1_t, &nl_t, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double(vgd,"PTOP", &ptop_8, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double(vgd,"PREF", &pref_8, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_1", &rc_1, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_2", &rc_2, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  
  printf("  Testing generic interface\n");
  if( Cvgd_new_build_vert(&vgd2, kind, version, 0, -1, -1, &ptop_8, &pref_8, &rc_1, &rc_2,
  			  a_m_8, b_m_8, a_t_8, b_t_8, ip1_m, ip1_t, nl_m, nl_t) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  // Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }

  Cvgd_free(&vgd2);
  printf("  Testing specific interface\n");
  if( Cvgd_new_build_vert_5002(&vgd2, -1, -1, ptop_8, pref_8, rc_1, rc_2,
			       a_m_8, b_m_8, a_t_8, b_t_8, ip1_m, ip1_t, nl_m, nl_t) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  // Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  
  free(ip1_m);
  free(ip1_t);
  free(a_m_8);
  free(b_m_8);
  free(a_t_8);
  free(b_t_8);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}
//========================================================================
int check_build_5005(vgrid_descriptor *vgd){
  int kind, version, nl, ier;
  int *ip1_m = NULL, *ip1_t = NULL;
  float rc_1, rc_2;
  double *a_m_8 = NULL, *b_m_8 = NULL, *a_t_8 = NULL, *b_t_8 = NULL, pref_8;
  vgrid_descriptor *vgd2 = NULL;

  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_M", &a_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_M", &b_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_T", &a_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_T", &b_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPT", &ip1_t, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double(vgd,"PREF", &pref_8, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_1", &rc_1, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_2", &rc_2, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  
  printf("  Testing generic interface\n");
  if( Cvgd_new_build_vert(&vgd2, kind, version, -1, -1, -1, NULL, &pref_8, &rc_1, &rc_2,
  			  a_m_8, b_m_8, a_t_8, b_t_8, ip1_m, ip1_t, nl, nl) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  // Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }

  Cvgd_free(&vgd2);
  printf("  Testing specific interface\n");
  if( Cvgd_new_build_vert_5005(&vgd2, -1, -1, pref_8, rc_1, rc_2,
  			       a_m_8, b_m_8, a_t_8, b_t_8, ip1_m, ip1_t, nl) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  //Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  free(ip1_m);
  free(ip1_t);
  free(a_m_8);
  free(b_m_8);
  free(a_t_8);
  free(b_t_8);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}

//========================================================================
int check_build_5100(vgrid_descriptor *vgd){
  int kind, version, nl, ier;
  int *ip1_m = NULL, *ip1_t = NULL;
  float rc_1, rc_2, rc_3, rc_4;
  double *a_m_8 = NULL, *b_m_8 = NULL, *c_m_8 = NULL, *a_t_8 = NULL,
    *b_t_8 = NULL, *c_t_8 = NULL, pref_8;
  vgrid_descriptor *vgd2 = NULL;

  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_M", &a_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_M", &b_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CC_M", &c_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_T", &a_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_T", &b_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CC_T", &c_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPT", &ip1_t, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double(vgd,"PREF", &pref_8, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_1", &rc_1, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_2", &rc_2, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_3", &rc_3, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_4", &rc_4, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  
  // No generic interface for 5100
  Cvgd_free(&vgd2);
  printf("  Testing specific interface\n");
  if( Cvgd_new_build_vert_5100(&vgd2, -1, -1, pref_8, rc_1, rc_2, rc_3, rc_4,
  			       a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8, ip1_m, ip1_t, nl) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  //Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  
  free(ip1_m);
  free(ip1_t);
  free(a_m_8);
  free(b_m_8);
  free(c_m_8);
  free(a_t_8);
  free(b_t_8);
  free(c_t_8);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}

//========================================================================
int check_build_21001(vgrid_descriptor *vgd){
  int kind, version, nl, ier;
  int *ip1_m = NULL, *ip1_t = NULL;
  float rc_1, rc_2, rc_3, rc_4;
  double *a_m_8 = NULL, *b_m_8 = NULL, *c_m_8 = NULL, *a_t_8 = NULL,
    *b_t_8 = NULL, *c_t_8 = NULL, pref_8;
  vgrid_descriptor *vgd2 = NULL;

  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_M", &a_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_M", &b_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CC_M", &c_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_T", &a_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_T", &b_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CC_T", &c_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPT", &ip1_t, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_1", &rc_1, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_2", &rc_2, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_3", &rc_3, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_4", &rc_4, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  
  // No generic interface for 21001
  Cvgd_free(&vgd2);
  printf("  Testing specific interface\n");
  if( Cvgd_new_build_vert_21001(&vgd2, -1, -1, rc_1, rc_2, rc_3, rc_4,
  			       a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8, ip1_m, ip1_t, nl) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  //Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  
  free(ip1_m);
  free(ip1_t);
  free(a_m_8);
  free(b_m_8);
  free(c_m_8);
  free(a_t_8);
  free(b_t_8);
  free(c_t_8);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}

//========================================================================
int check_build_21002(vgrid_descriptor *vgd){
  int kind, version, nl, ier;
  int *ip1_m = NULL, *ip1_t = NULL, *ip1_w = NULL;
  float rc_1, rc_2, rc_3, rc_4;
  double
    *a_m_8 = NULL, *b_m_8 = NULL, *c_m_8 = NULL,
    *a_t_8 = NULL, *b_t_8 = NULL, *c_t_8 = NULL,
    *a_w_8 = NULL, *b_w_8 = NULL, *c_w_8 = NULL;
  vgrid_descriptor *vgd2 = NULL;
  
  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_M", &a_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_M", &b_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CC_M", &c_m_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_T", &a_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_T", &b_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CC_T", &c_t_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPT", &ip1_t, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CA_W", &a_w_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CB_W", &b_w_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double_1d(vgd,"CC_W", &c_w_8, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPW", &ip1_w, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_1", &rc_1, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_2", &rc_2, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_3", &rc_3, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_4", &rc_4, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  
  // No generic interface for 21001
  Cvgd_free(&vgd2);
  printf("  Testing specific interface\n");
  if( Cvgd_new_build_vert_21002(&vgd2, -1, -1, rc_1, rc_2, rc_3, rc_4,
				a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8, a_w_8, b_w_8, c_w_8,
				ip1_m, ip1_t, ip1_w, nl) == VGD_ERROR) {
    return(VGD_ERROR);
  }
  //Test equality
  ier = Cvgd_vgdcmp(vgd, vgd2);
  if( ier != 0 ){
    printf("     Descritors not equal, Cvgd_vgdcmp code is %d\n", ier);
    return(VGD_ERROR);
  } else {
    printf("     Descritors are equal.\n");
  }
  
  free(ip1_m);
  free(ip1_t);
  free(ip1_w);
  free(a_m_8);
  free(b_m_8);
  free(c_m_8);
  free(a_t_8);
  free(b_t_8);
  free(c_t_8);
  free(a_w_8);
  free(b_w_8);
  free(c_w_8);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}

//========================================================================
int test_it(char *filename, int ind) {
  int ier, iun, vcode;
  iun = 10 + ind;
  char mode[]="RND";
  vgrid_descriptor *vgd = NULL;

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
  if( Cvgd_new_read(&vgd, iun, -1, -1, -1, -1) == VGD_ERROR ) {
    printf("ERROR with Cvgd_new_read on iun\n");
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd, "VCOD", &vcode, 0) == VGD_ERROR ){
    printf("ERROR with  Cvgd_get_int on VCOD\n");
    return(VGD_ERROR);
  }
  switch(vcode) {
  case 1001:
  case 2001:
  case 4001:
  case 5999:
    if( check_build_1001_2001_5999_4001(vgd, vcode) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 1002:
    if( check_build_1002(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 5001:
    if( check_build_5001(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 5002:
    if( check_build_5002(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 5005:
    if( check_build_5005(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 5100:
    if( check_build_5100(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 21001:
    if( check_build_21001(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 21002:
    if( check_build_21002(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  default:
    printf("In test ERROR unsupported Vcode %d\n",vcode);
    return(VGD_ERROR);
  }
  Cvgd_free(&vgd);
  ier = c_fstfrm(iun);
  ier = c_fclos(iun);
  return(VGD_OK);
}

//========================================================================
int main() {
  
  int i, ier, status = VGD_OK;

  ier = Cvgd_putopt_int("ALLOW_SIGMA",1);

  for (i = 0; i < (int) n_file; i++) {
    printf ("Testing %s\n", filenames[i]);
    if(test_it(filenames[i],i) == VGD_ERROR)
      status = VGD_ERROR;
  }  

  return(c_ut_report(status,"testing new_build"));  
  
}
