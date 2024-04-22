#!/bin/bash

#remove the diff tests for all tests except for vgrid_sample

convert_toctoc_5002=oui
add_toctoc=oui
compute_pressure_gz=oui
print_toctoc=oui
vgrid_sample=oui

set -e
#set -x

C_RED='\033[0;31m'
C_GREEN='\033[0;32m'
C_BLUE='\033[0;34m'
C_NC='\033[0m\n' # No Color

#===========================================================================
#===========================================================================
#===========================================================================
BIN_PATH=${1}
VALIDATION_DATA=${2}
TMPDIR=$(pwd)/data_out

#===========================================================================
#===========================================================================
#===========================================================================
if [ ${convert_toctoc_5002} = oui ];then
   rm -f ${TMPDIR}/dm2007050912-00-00_001_sans_p0
   editfst -s ${VALIDATION_DATA}/data_tests/dm2007050912-00-00_001 -d ${TMPDIR}/dm2007050912-00-00_001_sans_p0 <<EOF
      exclure(-1,P0)
EOF
   if [ ${?} -ne 0 ];then
      printf "${C_GREEN}Problem with editfst${C_NC}"
      exit 1
   fi
   #===========================================================================
   if [ oui = oui ];then
      # Test to convert !! and convert all !!SF to P0 if P0 not there
      # Write result in separate file
      ${BIN_PATH}/r.convert_toctoc_5002 -s ${TMPDIR}/dm2007050912-00-00_001_sans_p0 -d ${TMPDIR}/out_convert_toctoc_5002 
      if [ ${?} -ne 0 ];then
         printf "${C_RED}TEST 1 FAILLED!${C_NC}"
         exit 1
      fi
   fi
   #===========================================================================
   if [ oui = oui ];then
      # Test to convert !! and convert all !!SF to P0 if P0 not there
      # Write result in same file and rm !!SF and old !!
      rm -f  ${TMPDIR}/dm2007050912-00-00_001_same_file
      cp ${TMPDIR}/dm2007050912-00-00_001_sans_p0 ${TMPDIR}/dm2007050912-00-00_001_same_file
      ${BIN_PATH}/r.convert_toctoc_5002 -samefile -s ${TMPDIR}/dm2007050912-00-00_001_same_file
      if [ ${?} -ne 0 ];then
         printf "${C_RED}TEST 2 FAILLED!${C_NC}"
         exit 1
      fi
   fi
   #===========================================================================
   if [ oui = oui ];then
      # Selec one time, convert and recompute PX and compaire result with model PX
      rm -f ${TMPDIR}/tempo
      editfst -s ${TMPDIR}/dm2007050912-00-00_001_sans_p0 -d ${TMPDIR}/tempo <<EOF
         desire(-1,['>>','^^','!!'])
        desire(-1,-1,-1,338995474)
EOF
      ${BIN_PATH}/r.convert_toctoc_5002 -samefile -s ${TMPDIR}/tempo
      rm -f ${TMPDIR}/px
      ${BIN_PATH}/r.compute_pressure_gz -compute pressure -s ${TMPDIR}/tempo -d ${TMPDIR}/px -levels ALL_LEVELS
      fstcomp -ne -a ${TMPDIR}/px -b ${TMPDIR}/dm2007050912-00-00_001_sans_p0 | grep PX | awk '{if($7 > 5.e-4)print "NOTOK"}' > res
#xrec -imflds ${TMPDIR}/px 
      if grep -q NOTOK res;then
         printf "${C_RED}TEST 3 FAILED${C_NC}"
         exit 1
      fi
   fi
fi

#===========================================================================
#===========================================================================
#===========================================================================

if [ ${add_toctoc} = oui ];then

   echo "==============="
   printf "${C_BLUE}Test add_toctoc${C_NC}"
   echo "---------------"

   for ITEM in east.eta glbeta glbhyb regeta reghyb dm_1001_from_model_run
   do

     echo ==============
     echo ${ITEM}
     echo ==============
      rm -f $TMPDIR/toctoc
      ${BIN_PATH}/r.add_toctoc -s ${VALIDATION_DATA}/data_tests/${ITEM} -d $TMPDIR/toctoc -allow_sigma
      if [ ${?} -ne 0 ];then
         printf "${C_RED}TEST add_toctoc with ${ITEM} FAILLED!${C_NC}"
         exit 1
      fi
   done
fi

#===========================================================================
#===========================================================================
#===========================================================================

if [ ${compute_pressure_gz} = oui ];then

    function Get_Var_List(){
        VAR_LIST=MOMENTUM
	if [ ${ITEM} = 2001_from_model_run ];then
            VAR_LIST=TT
	elif [ ${ITEM} = dm_5002_from_model_run ];then
	    VAR_LIST='MOMENTUM THERMO'
	elif [ ${ITEM} = dm_5003_from_model_run ];then
	    VAR_LIST='MOMENTUM THERMO'
	elif [ ${ITEM} = dm_5004_from_model_run ];then
	    VAR_LIST='MOMENTUM THERMO'
	elif [ ${ITEM} = dm_5005_from_model_run ];then
	    VAR_LIST='MOMENTUM THERMO'
	elif [ ${ITEM} = dm_5100_from_model_run ];then
	    VAR_LIST='MOMENTUM THERMO'
	elif [ ${ITEM} = dm_21001_from_model_run_NON_SLEVE_HU_VT ];then
	    VAR_LIST='MOMENTUM THERMO'
	elif [ ${ITEM} = dm_21001_from_model_run_SLEVE_HU_VT ];then
	    VAR_LIST='MOMENTUM THERMO'
	fi
    }
    function Test_Comp_PX_GZ(){
	set -x
	for VAR in $VAR_LIST;do	   
            rm -f $TMPDIR/px tempo
	    if [ "${2}" = "" ];then
		ln -s ${VALIDATION_DATA}/data_tests/${ITEM} tempo
		USED_TTVT=""
		MESSAGE_TTVT=""
	    else
		if [ ${2} = VT ];then
		    editfst -s ${VALIDATION_DATA}/data_tests/${ITEM} -d tempo<<EOF
   exclure(-1,['TT','HU'])
EOF
		    USED_VAR="The virtual temperature VT were used."
		    MESSAGE_TTVT=", with VT"
		elif [ ${2} = TT ];then
		    editfst -s ${VALIDATION_DATA}/data_tests/${ITEM} -d tempo<<EOF
   exclure(-1,VT)
EOF
		    USED_TTVT="The variables used to compute the virtual temperature were TT and HU."
		    MESSAGE_TTVT=", with TT and HU"
		else
		    printf "${C_RED}TEST internal error, second parameter to ${0} must be blank or TT or VT, got ${2}${C_NC}"
		    exit 1
		fi
	    fi
	    file=${VALIDATION_DATA}/data_tests/${ITEM}
	    if [ ! -f ${file} ];then
		printf "${C_RED}TEST ERROR file doesn't exist ${file}${2}${C_NC}"
		exit 1
	    fi
            printf "${C_BLUE}TESTING ${1} computations${MESSAGE_TTVT}, for ${VAR} levels on file ${file}${C_NC}"
	    printf "${C_BLUE}  Inpuit file is $PWD/tempo${C_NC}"
            ${BIN_PATH}/r.compute_pressure_gz \
		       -compute ${1}\
		       -allow_sigma\
		       -s tempo\
		       -d $TMPDIR/px\
		       -levels ${VAR}
            fstcomp -ne \
		    -a $TMPDIR/px\
		    -b ${VALIDATION_DATA}/data_tests/${ITEM} > listing
	    cat listing
	    if [ ${1} = pressure ];then
		NOMVAR=PX
	    elif [ ${1} = gz ];then
		NOMVAR=GZ
	    else
		printf "${C_RED}Internal error in ${0}, parameter one should be pressure or gz but got ${1}"
		exit 1
	    fi
	    ERROR_REL_MAX=1.e-5
	    # Note that GZ may be close to zero at some points over the sea. Therefore large relative error are
	    # very likely for the low levels with GZ. Therefore if this happens, we check if the maximum error is larget
	    # than ERROR_MAX.
	    ERROR_MAX=1.e-3
	    
	    sed -e 's/</ /' -e 's/>/ /' listing | grep ${NOMVAR} | awk '{if($6 > ERROR_REL_MAX)print $0, " NOTOK"}' ERROR_REL_MAX=${ERROR_REL_MAX} > res
            if grep -q NOTOK res;then
		OTHER_TEST_OK=NO
		if [ ${NOMVAR} = GZ ];then
		    sed -e 's/</ /' -e 's/>/ /' res | grep ${NOMVAR} | awk '{if($12 > ERROR_MAX )print $0, " NOTOK"}' ERROR_MAX=${ERROR_MAX} > res2
		    if grep -q NOTOK res2;then
			printf "${C_RED}"
			echo "TEST on ${1} computations${MESSAGE_TTVT}, for levels ${VAR} on file ${VALIDATION_DATA}/data_tests/${ITEM} FAILED"
			echo "${USED_TTVT}"
			echo "The following records were NOTOK due to E-MAX > ${ERROR_MAX}"
			grep 'NOM *ETIKET' listing
			grep NOTOK res2 | sed "s/${NOMVAR}/<${NOMVAR}>/"
			printf "${C_NC}"
			exit 1
		    else
			OTHER_TEST_OK=YES
		    fi
		fi
		if [ ${OTHER_TEST_OK} = NO ];then
		    printf "${C_RED}"
		    echo "TEST on ${1} computations${MESSAGE_TTVT}, for levels ${VAR} on file ${VALIDATION_DATA}/data_tests/${ITEM} FAILED"
		    echo "${USED_TTVT}"
		    echo "The following records were NOTOK due to E-REL-MAX > ${ERROR_REL_MAX}"
		    grep 'NOM *ETIKET' listing
		    grep NOTOK res | sed "s/${NOMVAR}/<${NOMVAR}>/"
		    printf "${C_NC}"
		    exit 1
		fi
            fi
	    if grep -q 'PAS TROUVE' listing;then
		printf "${C_RED}"
		echo "TEST on ${1} computations${MESSAGE_TTVT}, for levels ${VAR} on file ${VALIDATION_DATA}/data_tests/${ITEM} FAILED"
		echo "${USED_TTVT}"
		echo "The following records were NOT FOUND in validation data file:"
		echo ${VALIDATION_DATA}/data_tests/${ITEM}		
		grep 'PAS TROUVE' listing
		echo "Output test file is $TMPDIR/px"
		printf "${C_NC}"
		exit 1
	    fi
            printf "${C_GREEN}Test computing ${1} on ${VAR} levels for file ${ITEM} is OK${C_NC}"
	done
    }
    
    # Test vgd_levels related computations:
    #    For pressure based coordinate vgd_levels returns the pressure
    #    For height based coordinate vgd_levels returns GZ
    # Pressure based files
    for ITEM in east.eta glbeta glbhyb regeta reghyb dm_1001_from_model_run dm_1002_from_model_run dm_1003_from_pgsm_lam_east_ops dm_5001_from_model_run dm_5001_from_model_run dm_5001_from_model_run dm_5002_from_model_run dm_5003_from_model_run dm_5004_from_model_run dm_5005_from_model_run dm_5100_from_model_run dm_5999_from_model_run; do
	Get_Var_List
	for VAR in $VAR_LIST;do
	    Test_Comp_PX_GZ pressure
	done
    done
    # Height based files
    for ITEM in	dm_21001_from_model_run_NON_SLEVE_HU_VT dm_21001_from_model_run_SLEVE_HU_VT;do
	Get_Var_List
	Test_Comp_PX_GZ gz
    done
    
    # Test computations related to the hydrostatic variables:
    #    gz for pressure based coordinate 
    #    pressure height based coordinate
    # Pressure based files
    for ITEM in sigma_from_interpolation eta_from_interpolation 5001_from_interpolation 5002_from_interpolation dm_5005_from_model_run_with_gz_vt dm_5100_from_model_run_with_gz_vt; do
    	Get_Var_List
	for VAR in $VAR_LIST;do
	    Test_Comp_PX_GZ gz VT
	    Test_Comp_PX_GZ gz TT
	done
    done
    # Height based files
    for ITEM in	dm_21001_from_model_run_NON_SLEVE_HU_VT dm_21001_from_model_run_SLEVE_HU_VT;do
	Get_Var_List
	Test_Comp_PX_GZ pressure VT
	Test_Comp_PX_GZ pressure TT
    done

fi
#===========================================================================
#===========================================================================
#===========================================================================
if [ ${print_toctoc} = oui ];then
    
   echo =================
   printf "${C_BLUE}Test print_toctoc${C_NC}"
   echo -----------------

   DEBUG=non

   cat > no_ip1t.txt<<EOF
 dm_1003_from_pgsm_lam_east_ops
EOF

   for ITEM in dm_2001_5001_from_editfst dm_1001_from_model_run dm_1002_from_model_run dm_1003_from_pgsm_lam_east_ops 2001_from_model_run dm_4001_from_model_run dm_5001_from_model_run dm_5001_from_model_run dm_5001_from_model_run dm_5002_from_model_run dm_5003_from_model_run dm_5004_from_model_run dm_5005_from_model_run dm_5100_from_model_run dm_5999_from_model_run dm_21001_from_model_run_NON_SLEVE dm_21001_from_model_run_SLEVE dm_21002_from_model_run_NON_SLEVE dm_21002_from_model_run_SLEVE
   do       
     
      kinds='=-1'
      if [ ${ITEM} = dm_2001_5001_from_editfst ];then
         kinds='2 5'
      fi
      
      for kind in ${kinds}
      do

	      lable=_kind_${kind}
         if [ ${kind} = =-1 ];then
            lable=""
         fi

	 printf "${C_BLUE}Test print_toctoc with kind=${kind} on ${VALIDATION_DATA}/data_tests/${ITEM}${C_NC}"

         ${BIN_PATH}/r.print_toctoc -fst ${VALIDATION_DATA}/data_tests/${ITEM} -no_box -kind ${kind} > to_erase.txt 2>&1
         if [ $? != 0 ];then
            printf "${C_RED}ERROR 1:  r.print_toctoc on ${VALIDATION_DATA}/data/${ITEM} -no_box -kind ${kind}${C_NC}"
		      exit 1
         fi

         ${BIN_PATH}/r.print_toctoc -fst ${VALIDATION_DATA}/data_tests/${ITEM} -no_box -kind ${kind} -convip > to_erase.txt 2>&1
         if [ $? != 0 ];then
            printf "${C_RED}ERROR 1:  r.print_toctoc on ${VALIDATION_DATA}/data/${ITEM} -no_box -kind ${kind} -convip${C_NC}"
		      exit 1
         fi

         ${BIN_PATH}/r.print_toctoc -fst ${VALIDATION_DATA}/data_tests/${ITEM} -ip1m_only -kind ${kind} > to_erase.txt 2>&1
         if [ $? != 0 ];then
            printf "${C_RED}ERROR 2:  r.print_toctoc -ip1m_only on ${VALIDATION_DATA}/data_tests/${ITEM}${C_NC}"
		      exit 1
         fi

         ${BIN_PATH}/r.print_toctoc -fst ${VALIDATION_DATA}/data_tests/${ITEM} -ip1t_only -kind ${kind} > to_erase.txt 2>&1
         if [ $? != 0 ];then
	     printf "${C_RED}ERROR 3:  r.print_toctoc -ip1t_only on ${VALIDATION_DATA}/data_tests/${ITEM}${C_NC}"
	     exit 1
      	fi

         ${BIN_PATH}/r.print_toctoc -fst ${VALIDATION_DATA}/data_tests/${ITEM} -ip1m_only -out output_file.txt -kind ${kind} > to_erase.txt 2>&1
         if [ $? != 0 ];then
             printf "${C_RED}ERROR 4.2:  r.print_toctoc on ${VALIDATION_DATA}/data_tests/${ITEM}${C_NC}"
	     exit 1
         fi

      done

   done

   #rm -f to_erase.txt
fi

#===========================================================================
#===========================================================================
#===========================================================================

if [ ${vgrid_sample} = oui ];then

   echo "==================="
   printf "${C_BLUE}Test r.vgrid_sample${C_NC}"
   echo "-------------------"

   printf "${C_BLUE}Test r.vgrid_sample 1${C_NC}"
   # r.vgrid_sample test 1
   rm -rf out_dir vgrid_sample.nml
   ${BIN_PATH}/r.vgrid_sample -out_dir out_dir > listing_r.vgrid_sample    
   #mkdir -p ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}
   #rm -rf ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}/out_dir
   #cp -r out_dir ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}

   # Test if all files produced by the current program
   # are in the validation directory
   OK=1
   for file in out_dir/*;do
       if [ ! -f ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}/out_dir/${file##*/} ];then
	   printf "${C_RED}ERROR: new file \"${file##*/}\" not in validation data ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}/out_dir${C_NC}"
	   OK=0
       fi
   done
    
   _file_checked=0
    for file in $(ls ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}/out_dir/*);do
	   ((n_file_checked=n_file_checked+1))
      ${BIN_PATH}/r.print_toctoc -no_box -fst out_dir/${file##*/} > ${file##*/}_tests.txt
	   ${BIN_PATH}/r.print_toctoc -no_box -fst ${file} > ${file##*/}_ctrl.txt
	   if diff  ${file##*/}_tests.txt ${file##*/}_ctrl.txt;then
         :
      else
	      OK=0
         printf "${C_RED}ERROR with ${file##*/} on r.vgrid_sample test 1${C_NC}"
      fi
   done
   if [ ${OK} = 0 ];then
      printf "${C_RED}ERROR with r.vgrid_sample for r.vgrid_sample test 1 for some Vcode, see diff output above${C_NC}"
      exit 1
   fi    
   if [ ${n_file_checked} = 0 ];then
       printf "${C_RED}ERROR with r.vgrid_sample for r.vgrid_sample test 1 no file were checked${C_NC}"
       exit 1
   fi

   # r.vgrid_sample test 2
   printf "${C_BLUE}Test r.vgrid_sample 2${C_NC}"
   cat > vgrid_sample.nml<<EOF
   &cfg
   vc_eta%ptop_8 = 1000.
   levs_eta = 0.0 ,0.5, 1.0
/
EOF
   rm -rf out_dir_2
   ${BIN_PATH}/r.vgrid_sample -out_dir out_dir_2 > listing_r.vgrid_sample_2
   #mkdir -p ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}
   #rm -rf ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}/out_dir_2
   #cp -r out_dir_2 ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}
   # Test if all files produced by the current program
   # are in the validation directory
   OK=1
   for file in out_dir_2/*;do
       if [ ! -f ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}/out_dir_2/${file##*/} ];then
	   printf "${C_RED}ERROR: new file \"${file##*/}\" not in validation data ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}/out_dir_2${C_NC}"
	   OK=0
       fi
   done
   n_file_checked=0
   for file in $(ls ${VALIDATION_DATA}/data_tests_res/vgrid_sample/${EC_ARCH}/out_dir_2/*);do
      ((n_file_checked=n_file_checked+1))
	   ${BIN_PATH}/r.print_toctoc -no_box -fst out_dir_2/${file##*/} > ${file##*/}_tests.txt
	   ${BIN_PATH}/r.print_toctoc -no_box -fst ${file} > ${file##*/}_ctrl.txt
	   if diff  ${file##*/}_tests.txt ${file##*/}_ctrl.txt;then
         :
      else
	  OK=0
          printf "${C_RED}ERROR with ${file##*/} on r.vgrid_sample test 2${C_NC}"
      fi
   done
   if [ ${OK} = 0 ];then
      printf "${C_RED}ERROR with r.vgrid_sample for r.vgrid_sample test 2 for some Vcode, see diff output above${C_NC}"
      exit 1
   fi
   if [ ${n_file_checked} = 0 ];then
	   printf "${C_RED}ERROR with r.vgrid_sample for r.vgrid_sample test 2 no file were checked${C_NC}"
	   exit 1
   fi

fi

echo 
printf "${C_GREEN}TEST(S) OK${C_NC}"
echo "${C_BLUE}Reminder: remove WORK to cleanup${C_NC}"
printf "\n=====>  testutilvgrid.sh ends: $(date) ###########\n"
echo
