
  if(! Cvgd_is_valid(self,"SELF")){
    Lib_Log(APP_LIBVGRID,APP_ERROR,"%s: invalid vgrid.\n",__func__);
    return(VGD_ERROR);
  }
  
  switch(self->vcode) {
  case 1:
    if( dpidpis ){
      Lib_Log(APP_LIBVGRID,APP_ERROR,"%s: dpidpis not supported for vertical coordinate 1\n",__func__);
      return(VGD_ERROR);
    }
    if(double_interface){
      if( C_compute_heights_0001_8(self, ni, nj, nk, ip1_list, levels_8) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_heights_0001(self, ni, nj, nk, ip1_list, levels) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;    
  case 1001:
    if(double_interface){
      if( C_compute_pressure_1001_1002_8(self, ni, nj, nk, ip1_list, levels_8, sfc_field_8, in_log) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_pressure_1001_1002(self, ni, nj, nk, ip1_list, levels, sfc_field, in_log) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;
  case 1002:
    if( dpidpis ){
      printf("(Cvgd) ERROR: dpidpis not implemented for vertical coordinate 1002\n");
      return(VGD_ERROR);
    }
    if(double_interface){
      if( C_compute_pressure_1001_1002_8(self, ni, nj, nk, ip1_list, levels_8, sfc_field_8, in_log) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_pressure_1001_1002(self, ni, nj, nk, ip1_list, levels, sfc_field, in_log) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;
  case 2001:
    if( dpidpis ){
      Lib_Log(APP_LIBVGRID,APP_ERROR,"%s: dpidpis not implemented for vertical coordinate 2001\n",__func__);
      return(VGD_ERROR);
    }
    if(double_interface){
      if( C_compute_pressure_2001_8(self, ni, nj, nk, ip1_list, levels_8, in_log) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_pressure_2001(self, ni, nj, nk, ip1_list, levels, in_log) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;
  case 4001:
    if( dpidpis ){
      Lib_Log(APP_LIBVGRID,APP_ERROR,"%s: dpidpis not implemented for vertical coordinate 4001\n",__func__);
      return(VGD_ERROR);
    }
    if( in_log ){
      Lib_Log(APP_LIBVGRID,APP_ERROR,"%s: option in_log not supported for vertical coordinate 4001\n",__func__);
      return(VGD_ERROR);
    }
    if(double_interface){
      if( C_compute_heights_4001_8(self, ni, nj, nk, ip1_list, levels_8) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_heights_4001(self, ni, nj, nk, ip1_list, levels) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;
  case 1003:
  case 5001:
    if(double_interface){
      if( C_compute_pressure_1003_5001_8(self, ni, nj, nk, ip1_list, levels_8, sfc_field_8, in_log, dpidpis) == VGD_ERROR )
	return(VGD_ERROR);
    } else {
      if( C_compute_pressure_1003_5001(self, ni, nj, nk, ip1_list, levels, sfc_field, in_log, dpidpis) == VGD_ERROR )
	return(VGD_ERROR);
    }
    break;
  case 5002:
  case 5003:
  case 5004:
  case 5005:
    if(double_interface){
      if( C_compute_pressure_5002_5003_5004_5005_8(self, ni, nj, nk, ip1_list, levels_8, sfc_field_8, in_log, dpidpis) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_pressure_5002_5003_5004_5005(self, ni, nj, nk, ip1_list, levels, sfc_field, in_log, dpidpis) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;
  case 5100:
    if(double_interface){
      if( C_compute_pressure_5100_8(self, ni, nj, nk, ip1_list, levels_8, sfc_field_8, sfc_field_ls_8, in_log, dpidpis) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_pressure_5100(self, ni, nj, nk, ip1_list, levels, sfc_field, sfc_field_ls, in_log, dpidpis) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;
  case 5999:
    if(double_interface){
      if( C_compute_pressure_1001_1002_8(self, ni, nj, nk, ip1_list, levels_8, sfc_field_8, in_log) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_pressure_1001_1002(self, ni, nj, nk, ip1_list, levels, sfc_field, in_log) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;
  case 21001:
  case 21002:
    if(double_interface){
      if( C_compute_heights_21001_8(self, ni, nj, nk, ip1_list, levels_8, sfc_field_8, sfc_field_ls_8) == VGD_ERROR)
	return(VGD_ERROR);
    } else {
      if( C_compute_heights_21001(self, ni, nj, nk, ip1_list, levels, sfc_field, sfc_field_ls) == VGD_ERROR)
	return(VGD_ERROR);
    }
    break;
  default:
    Lib_Log(APP_LIBVGRID,APP_ERROR,"%s: invalid kind or version: kind = %d, version = %d\n",__func__,self->kind,self->version);
    return(VGD_ERROR);
  }
  
  return(VGD_OK);

