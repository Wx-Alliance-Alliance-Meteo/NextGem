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
    "data/dm_21001_from_model_run_SLEVE",
    "data/dm_21002_from_model_run_SLEVE",
};

#define n_file (sizeof (filenames) / sizeof (const char *))

static double epsilon_8 = 1.e-300;

static float c_convip_IP2Level(int IP,int *kind) {

   int    mode=-1,flag=0, strglen=0;
   float  level=0.0;
   char   format;

   /*Convertir en niveau reel*/
    f77name(convip_plus)(&IP,&level,kind,&mode,&format,&flag,strglen);

   return(level);
}

static int similar_vec_r8(double *vec1, int n1, double *vec2, int n2, double precision) {
  int i;
  if(vec1) {
    if (vec2) {
      if ( n1 == n2 ) {
	for(i = 0; i < n1; i++) {
	  if( fabs(vec1[i]) < 1.e-307 ){
	    if( fabs(vec2[i]) > 1.e-307 ){
	      //printf("zero, vec1[i] = %f, vec2[i] = %f\n", vec1[i], vec2[i]);
	      return(-1);
	    }
	  } else {
	    if ( fabs(vec1[i]-vec2[i])/fabs(vec1[i]) > precision ){
	      //printf("non zero, vec1[i] = %f, vec2[i] = %f, fabs(vec1[i]-vec2[i])/vec1[i]= %f\n", vec1[i], vec2[i], fabs(vec1[i]-vec2[i])/vec1[i]);
	      return(-1);
	    }
	  }
	}
      } else {
	// Vectors are not the same size.
	return(-2);
      }
    } else {
      // vec2 not allocated
      return(-3);
    }
  }
  // Vector are the same or are not allocated.
  return(0);
}

static int c_vgdcmp_no_ip1_check(vgrid_descriptor *vgd1, vgrid_descriptor *vgd2, double precision) {
  if(! vgd1){
    return(-21);
  }
  if(! vgd2){
    return(-22);
  }
  if (vgd1->vcode != vgd2->vcode)                   return(-101);
  if (vgd1->kind != vgd2->kind)                     return(-102);
  if (vgd1->version != vgd2->version)               return(-103);
  if (strcmp(vgd1->ref_name, vgd2->ref_name) != 0 ) return(-104);
  if (strcmp(vgd1->ref_namel, vgd2->ref_namel) != 0 )return(-120);
  if (vgd1->nl_w != vgd2->nl_w) return(-23);
  if (memcmp(&(vgd1->ptop_8),&(vgd2->ptop_8), sizeof(double)/sizeof(char) ))return(-105);
  if (memcmp(&(vgd1->pref_8),&(vgd2->pref_8), sizeof(double)/sizeof(char) ))return(-106);
  if (memcmp(&(vgd1->rcoef1),&(vgd2->rcoef1), sizeof(float) /sizeof(char) ))return(-107);
  if (memcmp(&(vgd1->rcoef2),&(vgd2->rcoef2), sizeof(float) /sizeof(char) ))return(-108);
  if (memcmp(&(vgd1->rcoef3),&(vgd2->rcoef3), sizeof(float) /sizeof(char) ))return(-116);
  if (memcmp(&(vgd1->rcoef4),&(vgd2->rcoef4), sizeof(float) /sizeof(char) ))return(-117);
  if(similar_vec_r8(vgd1->a_m_8, vgd1->nl_m, vgd2->a_m_8, vgd2->nl_m, precision) != 0) return (-111);
  if(similar_vec_r8(vgd1->b_m_8, vgd1->nl_m, vgd2->b_m_8, vgd2->nl_m, precision) != 0) return (-112);
  if(similar_vec_r8(vgd1->a_t_8, vgd1->nl_t, vgd2->a_t_8, vgd2->nl_t, precision) != 0) return (-113);
  if(similar_vec_r8(vgd1->b_t_8, vgd1->nl_t, vgd2->b_t_8, vgd2->nl_t, precision) != 0) return (-114);
  if(similar_vec_r8(vgd1->a_w_8, vgd1->nl_w, vgd2->a_w_8, vgd2->nl_w, precision) != 0) return (-125);
  if(similar_vec_r8(vgd1->b_w_8, vgd1->nl_w, vgd2->b_w_8, vgd2->nl_w, precision) != 0) return (-126);

  return(0);
}


//========================================================================
int check_gen_1001_2001_4001(vgrid_descriptor *vgd, int vcode){
  int kind, kind2, version, nk, ier, k;
  int *ip1_m = NULL;
  float *hyb;
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
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  hyb = malloc( nk * sizeof(float) );
  if(! hyb ){ 
    printf("ERROR: check_gen_1001_2001_4001, problem allocating hyb\n");
    return(VGD_ERROR);
  }
  for( k=0; k < nk; k++){
    hyb[k] = c_convip_IP2Level(ip1_m[k], &kind2);
  }
  
  printf("  Testing generic interface\n");
  if( Cvgd_new_gen(&vgd2, kind, version, hyb, nk, NULL, NULL, NULL, NULL, NULL,
		   -1, -1, NULL, NULL, 0) == VGD_ERROR ) {
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
    if( Cvgd_new_gen_1001(&vgd2, hyb, nk, 0, 0) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 2001:
    if( Cvgd_new_gen_2001(&vgd2, hyb, nk, 0, 0) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 4001:
    if( Cvgd_new_gen_4001(&vgd2, hyb, nk, 0, 0) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  default:
    printf("In test ERROR in check_gen_1001_2001_4001, unsupported Vcode %d\n",vcode);
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
  free(hyb);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}
//========================================================================
int check_gen_1002(vgrid_descriptor *vgd){
  int kind, kind2, version, nk, ier, k;
  int *ip1_m = NULL;
  float *hyb;
  double ptop_8;
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
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nk, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_double(vgd,"PTOP", &ptop_8, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  hyb = malloc( nk * sizeof(float) );
  if(! hyb ){ 
    printf("ERROR: check_gen_1002, problem allocating hyb\n");
    return(VGD_ERROR);
  }
  for( k=0; k < nk; k++){
    hyb[k] = c_convip_IP2Level(ip1_m[k], &kind2);
  }
  printf("  Testing generic interface\n");
  if( Cvgd_new_gen(&vgd2, kind, version, hyb, nk, NULL, NULL, &ptop_8, NULL, NULL,
		   0, 0, NULL, NULL, 0) == VGD_ERROR) {
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
  if( Cvgd_new_gen_1002(&vgd2, hyb, nk, ptop_8, 0, 0) == VGD_ERROR) {
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
  free(hyb);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}

//========================================================================
int check_gen_5001(vgrid_descriptor *vgd){
  int kind, kind2, version, nk, ier, k;
  int *ip1_m = NULL;
  float rc_1, *hyb;
  double ptop_8, pref_8;
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
  hyb = malloc( nk * sizeof(float) );
  if(! hyb ){ 
    printf("ERROR: check_gen_1002, problem allocating hyb\n");
    return(VGD_ERROR);
  }
  for( k=0; k < nk; k++){
    hyb[k] = c_convip_IP2Level(ip1_m[k], &kind2);
  }
  
  printf("  Testing generic interface\n");
  if( Cvgd_new_gen(&vgd2, kind, version, hyb, nk, &rc_1, NULL,
	      &ptop_8, &pref_8, NULL, 0, 0, NULL, NULL, 0) == VGD_ERROR) {
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
  if( Cvgd_new_gen_5001(&vgd2, hyb, nk, ptop_8, pref_8, rc_1, 0, 0) == VGD_ERROR) {
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
    
  free(hyb);
  free(ip1_m);
  Cvgd_free(&vgd2);

  return(VGD_OK);
}

//========================================================================
int check_gen_5002(vgrid_descriptor *vgd){
  int kind, kind2, version, nl_m, ier, nk, k;
  int *ip1_m = NULL;
  float rc_1, rc_2, *hyb;
  double ptop_8, pref_8;
  vgrid_descriptor *vgd2 = NULL;

  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nl_m, 0) == VGD_ERROR ){
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
  nk = nl_m - 1;
  hyb = malloc( nk * sizeof(float) );
  if(! hyb ){ 
    printf("ERROR: check_gen_5002, problem allocating hyb\n");
    return(VGD_ERROR);
  }
  for( k=0; k < nk; k++){
    hyb[k] = c_convip_IP2Level(ip1_m[k], &kind2);
    //printf("k = %d, hyb[k] = %f\n", k, hyb[k]);
  }
 
  printf("  Testing generic interface\n");
  if( Cvgd_new_gen(&vgd2, kind, version, hyb, nk, &rc_1, &rc_2,
		   &ptop_8, &pref_8, NULL,0, 0, NULL, NULL, 0) == VGD_ERROR ) {
    printf("ERROR: check_gen_5002, problem with Cvgd_new_gen\n");
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
  if( Cvgd_new_gen_5002(&vgd2, hyb, nk, ptop_8, pref_8, rc_1, rc_2, 0, 0) == VGD_ERROR) {
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
  
  free(hyb);
  free(ip1_m);
  Cvgd_free(&vgd2);
  return(VGD_OK);
}
//========================================================================
int check_gen_5005_5100_21001_21002(vgrid_descriptor *vgd, int vcode){
  int kind, kind2, version, vcode2, nl, ier, nk, k, ok;
  int *ip1_m = NULL, *ip1_t = NULL, high_precision, test_ind=5;
  float rc_1, rc_2, rc_3, rc_4, *hyb, dhm, dht, dhw;
  double pref_8, ptop_out_8, *b_m_8 = NULL;
  vgrid_descriptor *vgd2 = NULL;

  printf("\n\n=== TESTING Vcode %d ======================================\n", vcode);

  if( Cvgd_get_int(vgd,"KIND", &kind, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd,"VERS", &version, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int_1d(vgd,"VIPM", &ip1_m, &nl, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_int(vgd, "VCOD", &vcode2, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if(vcode2 != vcode){
    printf("In test check_gen_5005_5100_21001_21002 Vcode don't matck %d, vs %d", vcode, vcode2);
  }
  if(Cvgd_is_valid(vgd,"pref_8_valid")){
    if( Cvgd_get_double(vgd,"PREF", &pref_8, 0) == VGD_ERROR ){
      return(VGD_ERROR);
    } 
  }
  if( Cvgd_get_float(vgd,"RC_1", &rc_1, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"RC_2", &rc_2, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if(Cvgd_is_valid(vgd,"rcoef3_valid")){
    if( Cvgd_get_float(vgd,"RC_3", &rc_3, 0) == VGD_ERROR ){
      return(VGD_ERROR);
    }      
    if( Cvgd_get_float(vgd,"RC_4", &rc_4, 0) == VGD_ERROR ){
      return(VGD_ERROR);
    }
  }
  if( Cvgd_get_float(vgd,"DHM ", &dhm, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }
  if( Cvgd_get_float(vgd,"DHT ", &dht, 0) == VGD_ERROR ){
    return(VGD_ERROR);
  }  
  if(Cvgd_is_valid(vgd,"dhw_valid")){
    if( Cvgd_get_float(vgd,"DHW ", &dhw, 0) == VGD_ERROR ){
      return(VGD_ERROR);
    }
  }
  
  nk = nl - 2;
  hyb = malloc( nk * sizeof(float) );
  if(! hyb ){ 
    printf("ERROR: check_gen_5005_5100_21001_21002, problem allocating hyb\n");
    return(VGD_ERROR);
  }
  for( k=0; k < nk; k++){
    hyb[k] = c_convip_IP2Level(ip1_m[k], &kind2);
    //printf("k = %d, hyb[k] = %f\n", k, hyb[k]);
  }

  //=========================
  //Testing generic interface
  //-------------------------
  printf("  Testing generic interface for Vcode %d\n", vcode);
  high_precision=1;
  switch(vcode) {
  case 5005:
    if( Cvgd_new_gen(&vgd2, kind, version, hyb, nk, &rc_1, &rc_2, 
		     NULL, &pref_8, &ptop_out_8, 0, 0, &dhm, &dht, 0) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    if( abs(ptop_out_8 - 5.611828)/5.611828 > 1.e-4){
      printf("   Error with ptop_out_8 for Vcode %d\n", vcode);
      return(VGD_ERROR);
    }
    break;
  case 5100:
    if( Cvgd_new_gen2(&vgd2, kind, version, hyb, nk, &rc_1, &rc_2, &rc_3, &rc_4,
		      NULL, &pref_8, &ptop_out_8, 0, 0, &dhm, &dht, NULL, 1) ){
      return(VGD_ERROR);
    }
    if( abs(ptop_out_8 - 5.611828)/5.611828 > 1.e-4){
      printf("   Error with ptop_out_8 for Vcode %d\n", vcode);
      return(VGD_ERROR);
    }
    break;
  case 21001:
    high_precision=0;
    // Small diffence in average hyb to compute thermo levels result is ip1_t
    // that are different for some levels. Therefore at the moment, Vcode 21001
    // cannot successfully be rebuilt. Test will not be make on ip1 and at low
    // precision.
    if( Cvgd_new_gen2(&vgd2, kind, version, hyb, nk, &rc_1, &rc_2, &rc_3, &rc_4,
    		      NULL, NULL, NULL, 0, 0, &dhm, &dht, NULL, 0) ){
      return(VGD_ERROR);
    }
    break;    
  case 21002:
    high_precision=0;
    // Small diffence in average hyb to compute thermo levels result is ip1_t
    // that are different for some levels. Therefore at the moment, Vcode 21002
    // cannot successfully be rebuilt. Test will not be make on ip1 and at low
    // precision.
    if( Cvgd_new_gen2(&vgd2, kind, version, hyb, nk, &rc_1, &rc_2, &rc_3, &rc_4,
    		      NULL, NULL, NULL, 0, 0, &dhm, &dht, &dhw, 0) ){
      return(VGD_ERROR);
    }    
    break;    
  default:
    printf("In test ERROR in check_gen_5005_5100_21001_21002, on test generic interface, unsupported Vcode %d\n",vcode);
    return(VGD_ERROR);
  }

  // Test equality
  if (high_precision){
    ier = Cvgd_vgdcmp(vgd, vgd2);
    if( ier != 0 ){
      printf("     Descritors not equal for Vcode %d, Cvgd_vgdcmp code is %d\n", vcode, ier);
      //ier = Cvgd_print_desc(vgd, -1, -1);
      //ier = Cvgd_print_desc(vgd2, -1, -1);
      return(VGD_ERROR);
    } else {
      printf("     Descritors are equal for Vcode %d.\n", vcode);
    }
    Cvgd_free(&vgd2);
  }else{
    // Note, since A, B and C are comming from hyb float value, we don't expecte a hight
    // precision on this compution even if these are stored in double arrays.
    ier = c_vgdcmp_no_ip1_check(vgd, vgd2, 1.e-4);
    if( ier != 0 ){
      printf("     Descritors not equal for Vcode %d, Cvgd_vgdcmp code is %d\n", vcode, ier);
      //ier = Cvgd_print_desc(vgd, -1, -1);
      //ier = Cvgd_print_desc(vgd2, -1, -1);
      return(VGD_ERROR);
    } else {
      printf("     Descritors are equal for Vcode %d. but without ip1 list check\n", vcode);
    }
    Cvgd_free(&vgd2);
  }
  
  //==========================
  //Testing specific interface
  //--------------------------
  printf("  Testing specific interface for Vcode %d\n", vcode);
  high_precision=1;
  switch(vcode) {
  case 5005:
    if( Cvgd_new_gen_5005(&vgd2, hyb, nk, pref_8, &ptop_out_8, rc_1, rc_2, 0, 0, dhm, dht) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 5100:
    if( Cvgd_new_gen_5100(&vgd2, hyb, nk, pref_8, &ptop_out_8, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, 1) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 21001:
    high_precision=0;
    // Small diffence in average hyb to compute thermo levels result is ip1_t
    // that are different for some levels. Therefore at the moment, Vcode 21001
    // cannot successfully be rebuilt.
    if( Cvgd_new_gen_21001(&vgd2, hyb, nk, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;    
  case 21002:
    high_precision=0;
    // Small diffence in average hyb to compute vertical velocity levels result is ip1_t
    // that are different for some levels. Therefore at the moment, Vcode 21002
    // cannot successfully be rebuilt.
    if( Cvgd_new_gen_21002(&vgd2, hyb, nk, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, dhw) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;    
  default:
    printf("In test ERROR in check_gen_5005_5100_21001_21002, unsupported Vcode for specific interface %d\n",vcode);
    return(VGD_ERROR);
  }
  //Test equality
  if(high_precision){
    ier = Cvgd_vgdcmp(vgd, vgd2);
    if( ier != 0 ){
      printf("     Descritors not equal for specific interface Vcode %d, Cvgd_vgdcmp code is %d\n", vcode, ier);
      return(VGD_ERROR);
    } else {
      printf("     Descritors are equal for specific interface Vcode %d.\n", vcode);
    }
    Cvgd_free(&vgd2);
  } else {
    // Note, since A, B and C are comming from hyb float value, we don't expecte a hight
    // precision on this compution even if these are stored in double arrays.
    ier = c_vgdcmp_no_ip1_check(vgd, vgd2, 1.e-4);
    if( ier != 0 ){
      printf("     Descritors not equal for Vcode %d, Cvgd_vgdcmp code is %d\n", vcode, ier);
      //ier = Cvgd_print_desc(vgd, -1, -1);
      //ier = Cvgd_print_desc(vgd2, -1, -1);
      return(VGD_ERROR);
    } else {
      printf("     Descritors are equal for Vcode %d. but without ip1 list check\n", vcode);
    }
    Cvgd_free(&vgd2);
  }
  
  //=====================
  // Test hyb_flat option
  //---------------------   
  // Reproduce default vgrid but with option for option hyb_flat = hyb[0]
  printf("  Testing specific interface to reproduce default but we hyb_flat = hyb[0] for Vcode %d\n", vcode);
  high_precision=1;
  switch(vcode) {
  case 5005:
    if( Cvgd_new_gen_5005_2(&vgd2, hyb, nk, pref_8, &ptop_out_8, rc_1, rc_2, 0, 0, dhm, dht,hyb[0]) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
    break;
  case 5100:
    if( Cvgd_new_gen_5100_2(&vgd2, hyb, nk, pref_8, &ptop_out_8, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, 1, hyb[0]) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 21001:
    high_precision=0;
    // Small diffence in average hyb to compute thermo levels result is ip1_t
    // that are different for some levels. Therefore at the moment, Vcode 21001
    // cannot successfully be rebuilt.
    if( Cvgd_new_gen_21001_2(&vgd2, hyb, nk, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, hyb[0]) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;    
  case 21002:
    high_precision=0;
    // Small diffence in average hyb to compute vertical velocity levels result is ip1_t
    // that are different for some levels. Therefore at the moment, Vcode 21002
    // cannot successfully be rebuilt.
    if( Cvgd_new_gen_21002_2(&vgd2, hyb, nk, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, dhw, hyb[0]) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;        
  default:
    printf("In test ERROR in check_gen_5005_5100_21001_21002, unsupported Vcode for Test default hyb_flat for option value %d\n",vcode);
    return(VGD_ERROR);
  }
    
  //Test equality
  if(high_precision){
    ier = Cvgd_vgdcmp(vgd, vgd2);
    if( ier != 0 ){
      printf("     Descriptors not equal in trying to reproduce default with option hyb_float = hyb[0] on %d, Cvgd_vgdcmp code is %d\n", vcode, ier);
      return(VGD_ERROR);
    } else {
      printf("     Descriptors are equal in reproduce default with option hyb_float = hyb[0]\n");
    }
    Cvgd_free(&vgd2);
  }else{
    // Note, since A, B and C are comming from hyb float value, we don't expecte a hight
    // precision on this compution even if these are stored in double arrays.
    ier = c_vgdcmp_no_ip1_check(vgd, vgd2, 1.e-4);
    if( ier != 0 ){
      printf("     Descritors not equal for Vcode %d, Cvgd_vgdcmp code is %d\n", vcode, ier);
      //ier = Cvgd_print_desc(vgd, -1, -1);
      //ier = Cvgd_print_desc(vgd2, -1, -1);
      return(VGD_ERROR);
    } else {
      printf("     Descritors are equal for Vcode %d. but without ip1 list check\n", vcode);
    }
    Cvgd_free(&vgd2);
  }
  
  /* //Test option hyb_flat out of range */
  printf("  In test c_new_gen_all on Vcode %d, following 4 errors lines on hyb_flat range are expected\n", vcode);
  switch(vcode) {
  case 5005:
    if(Cvgd_new_gen_5005_2(&vgd2, hyb, nk, pref_8, &ptop_out_8, rc_1, rc_2, 0, 0, dhm, dht, hyb[0]*0.99) != VGD_ERROR) {
      printf("   ERROR in test c_new_gen_all on vcode %d, hyb_flat out of low range but did not get cut\n", vcode);
      return(VGD_ERROR);
    }
    if(Cvgd_new_gen_5005_2(&vgd2, hyb, nk, pref_8, &ptop_out_8, rc_1, rc_2, 0, 0, dhm, dht, hyb[nk-1]*1.01) != VGD_ERROR) {
      printf("   ERROR in test c_new_gen_all on Vcode %d, hyb_flat out of high range but did not get cut\n", vcode);
      return(VGD_ERROR);
    }
    break;
  case 5100:
    if(Cvgd_new_gen_5100_2(&vgd2, hyb, nk, pref_8, &ptop_out_8, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, 1, hyb[0]*0.99) != VGD_ERROR) {
      printf("   ERROR in test c_new_gen_all on vcode %d, hyb_flat out of low range but did not get cut\n", vcode);
      return(VGD_ERROR);
    }
    if(Cvgd_new_gen_5100_2(&vgd2, hyb, nk, pref_8, &ptop_out_8, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, 1, hyb[nk-1]*1.01) != VGD_ERROR) {
      printf("   ERROR in test c_new_gen_all on vcode %d, hyb_flat out of high range but did not get cut\n", vcode);
      return(VGD_ERROR);
    }
    break;
  case 21001:
    if(Cvgd_new_gen_21001_2(&vgd2, hyb, nk, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, hyb[0]*1.01) != VGD_ERROR) {
      printf("   ERROR in test c_new_gen_all on vcode %d, hyb_flat out of high range but did not get cut\n", vcode);
      return(VGD_ERROR);
    }
    if(Cvgd_new_gen_21001_2(&vgd2, hyb, nk, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, hyb[nk-1]*0.99) != VGD_ERROR) {
      printf("   ERROR in test c_new_gen_all on vcode %d, hyb_flat out of low range but did not get cut\n", vcode);
      return(VGD_ERROR);
    }
    break;
  case 21002:
    if(Cvgd_new_gen_21002_2(&vgd2, hyb, nk, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, dhw, hyb[0]*1.01) != VGD_ERROR) {
      printf("   ERROR in test c_new_gen_all on vcode %d, hyb_flat out of high range but did not get cut\n", vcode);
      return(VGD_ERROR);
    }
    if(Cvgd_new_gen_21002_2(&vgd2, hyb, nk, rc_1, rc_2, rc_3, rc_4, 0, 0, dhm, dht, dhw, hyb[nk-1]*0.99) != VGD_ERROR) {
      printf("   ERROR in test c_new_gen_all on vcode %d, hyb_flat out of low range but did not get cut\n", vcode);
      return(VGD_ERROR);
    }
    break;
  default:
    printf("ERROR in check_gen_5005_5100_21001_21002, unsupported Vcode %d, for Test hyb_flat out of range \n",vcode);
    return(VGD_ERROR);
 }

  //Test option hyb_flat non default value
  // Reduce rcoefs otherwise B is zero near the top and test may not be valid.
  printf("  Testing option hyb_flast = hyb[%d-1], just checking if top %d levels have Bm = 0.\n",test_ind,test_ind);
  switch(vcode) {
  case 5005:
    if( Cvgd_new_gen_5005_2(&vgd2, hyb, nk, pref_8, &ptop_out_8, 1, 2, 0, 0, dhm, dht, hyb[test_ind-1]) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 5100:
    if( Cvgd_new_gen_5100_2(&vgd2, hyb, nk, pref_8, &ptop_out_8, 1, 1, 1, 2, 0, 0, dhm, dht, 1, hyb[test_ind-1]) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 21001:
    if( Cvgd_new_gen_21001_2(&vgd2, hyb, nk, 1, 1, 1, 2, 0, 0, dhm, dht, hyb[test_ind-1]) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  case 21002:
    if( Cvgd_new_gen_21002_2(&vgd2, hyb, nk, 1, 1, 1, 2, 0, 0, dhm, dht, dhw, hyb[test_ind-1]) == VGD_ERROR) {
      return(VGD_ERROR);
    }
    break;
  default:
    printf("In test ERROR in check_gen_5005_5100_21001_21002, unsupported Vcode %d, for Test hyb_flat = hyb[%d-1] non default value\n",vcode,test_ind);
    return(VGD_ERROR);
  }
  // Here we cannot test the equality since we have no control fst file for hyb_flat = hyb[test_ind-1].
  // But we only really care to see if there are test_ind level with B = 0 at domain top.
  // The code to generate Vcode 5100 is the same for all values of hyb_flat.
  // The rest of vgrid functions are the same for 5100 for all value of hyb_flat too.
  // Therefore further testing in Vcode 5100 and option hyb_flat are not required.
  if( Cvgd_get_double_1d(vgd2,"CB_M", &b_m_8, &nk, 0) == VGD_ERROR ){
    printf("   ERROR in test c_new_gen_all on %d, getting CB_M\n", vcode);
    return(VGD_ERROR);
  }
  Cvgd_free(&vgd2);
  ok = 1;
  for (k=0; k < test_ind; k++){
    if( b_m_8[k] > epsilon_8 ){
      ok = 0;
      printf("    ERROR in test c_new_gen_all on Vcode %d, b_m_8 should be 0.0 for k = %d, got %-# 25.15G\n",vcode,k,b_m_8[k]);
    }
  }
  if(ok){
    printf("    Bm is zero for top %d levels\n",test_ind);
  } else {
    return(VGD_ERROR);
  }

  free(hyb);
  free(ip1_m);
  free(ip1_t);
  return(VGD_OK);
}

//========================================================================
int test_it(char *filename, int ind) {
  int ier, iun, vcode;
  iun = 10 + ind;
  char mode[]="RND+R/O";
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
    if( check_gen_1001_2001_4001(vgd, vcode) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 1002:
    if( check_gen_1002(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 5001:
    if( check_gen_5001(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 5002:
    if( check_gen_5002(vgd) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  case 5005:
  case 5100:
  case 21001:
  case 21002:
    if( check_gen_5005_5100_21001_21002(vgd, vcode) == VGD_ERROR){
      return(VGD_ERROR);
    }
    break;
  default:
    printf("In test ERROR unsupported Vcode %d\n",vcode);
    return(VGD_ERROR);
  }
  Cvgd_free(&vgd);
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

  return(c_ut_report(status,"testing new_gen"));  
}
