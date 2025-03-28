#!/bin/bash
## This script aims to gather the queried variables from a simulation output files,
## into monthly files ready for being treated by the make_interpol_iagos script.

#MSUB -q skylake
#MSUB -A gen2201
#MSUB -r extr_Oslo
#MSUB -e $CCCSCRATCHDIR/out_extract_model_Oslo
#MSUB -o $CCCSCRATCHDIR/out_extract_model_Oslo
#MSUB -n 1
#MSUB -c 1
#MSUB -T 86400
#MSUB -m work,scratch,store
#MSUB -Q normal
#MSUB -x
set +x

module purge
source ~igcmg/MachineEnvironment/irene/env_irene
module load r
module load nco

## First of all: please indicate here the paths for the following tools.
## 1/ NetCDF libraries.
NETCDF_INCDIR="-I$NETCDF_INCDIR -I$NETCDFFORTRAN_INCDIR"
NETCDF_LIBDIR="-L$NETCDF_LIBDIR -L$NETCDFFORTRAN_LIBDIR"
NETCDF_LIB="-lnetcdf -lnetcdff -lstdc++"
## 2/ The directory where the Shell command 'gfortran' is defined.
gfortran_dir=/usr/bin
## 3/ The directory where the Shell command 'rscript' is defined.
rscript_dir=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/bin/
### Current and parent directories
# pdir=$base_dir/Pack_code/routines/grid_IAGOS_model
pdir=`pwd`
parent_dir=$base_dir/Pack_code/routines

## Enter your arguments.

yyyymm1=200101 ; yyyymm2=201712

model=INCA
config=LMDZORINCA_144_142_39
experiment=Rescaled_Nudged
time_resol=daily
sci_program=idk

model=Oslo-CTM3
config=160_80_60
experiment=ACACIA
time_resol=daily
sci_program=ACACIA

# model=MOZART3
# config=128_64_60
# experiment=ACACIA
# time_resol=daily
# sci_program=ACACIA

# model=EMAC
# config=128_64_90
# experiment=ACACIA
# sci_program=ACACIA

# model=UKESM1.1
# config=192_144_85
# experiment=ACACIA
# sci_program=ACACIA

bin_extract_model_variables=TRUE

list_NOy_variables_ACACIA='vmrno vmrno2 vmrhno3 vmrpan'
if [ "$model" = INCA ] ; then
    list_model_variables="Ozone Carbon_monoxide Water_vapour Nitric_acid Peroxyacetyl_nitrate Stratospheric_ozone Nitrogen_monoxide Nitrogen_dioxide Production_rate_ozone Loss_rate_ozone Temperature Inert_strat_ozone"
    extract_NOy=TRUE
    extract_P=TRUE
    extract_PV=TRUE
    list_NOy_variables=${list_NOy_variables_ACACIA}' vmrn vmrno3 vmrhno2 vmrhno4 vmrn2o5 vmrmpan vmrisopno3 vmrapinpan vmrpco3pan vmronitu vmronitr vmrclono2 vmrclno2 vmrbrono2'
elif [ "$model" = "Oslo-CTM3" ] ; then
    list_model_variables="Ozone Carbon_monoxide Water_vapour Nitric_acid Peroxyacetyl_nitrate Nitrogen_monoxide Nitrogen_dioxide Temperature"
    extract_NOy=TRUE
    extract_P=TRUE
    extract_PV=TRUE
    list_NOy_variables=${list_NOy_variables_ACACIA}
elif [ "$model" = EMAC ] ; then
    list_model_variables="Ozone Carbon_monoxide Water_vapour Nitric_acid Peroxyacetyl_nitrate Nitrogen_monoxide Nitrogen_dioxide Temperature"
    extract_NOy=TRUE
    extract_P=TRUE
    extract_PV=TRUE
    list_NOy_variables=${list_NOy_variables_ACACIA}
elif [ "$model" = MOZART3 ] ; then
    list_model_variables="Ozone Carbon_monoxide Water_vapour Nitric_acid Peroxyacetyl_nitrate Nitrogen_monoxide Nitrogen_dioxide Temperature"
    extract_NOy=TRUE
    extract_P=FALSE
    extract_PV=FALSE
    list_NOy_variables=${list_NOy_variables_ACACIA}
elif [ "$model" = "UKESM1.1" ] ; then
    list_model_variables="Ozone Carbon_monoxide Water_vapour Nitric_acid Peroxyacetyl_nitrate Nitrogen_monoxide Nitrogen_dioxide Temperature"
    extract_NOy=TRUE
    extract_P=TRUE
    extract_PV=TRUE
    list_NOy_variables=${list_NOy_variables_ACACIA}
fi


echo -e '\e[00;31m Extracting variables from the simulation outputs.\e[00m'

## PV, P_S and P. To be in separated files.
list_other_variables=""
if [ "$extract_PV" = TRUE ] ; then
    list_other_variables=$list_other_variables" "Potential_vorticity
fi
if [ "$extract_P" = TRUE ] ; then
    list_other_variables=$list_other_variables" "Pressure
else
    list_other_variables=$list_other_variables" "Surface_pressure
fi
########### Finding the variables under their name in the current model. #############
base_dir=$CCCWORKDIR
pdir=./
model_var_header=$(head -1 $pdir/table_variable_names_in_models.txt)
icolmodel=1
current_column_in_header=$(echo $model_var_header | awk -v icolmodel=$icolmodel '{print $icolmodel}')
## Finds the right column in the variables tabular for the input name of the variable (the one in the initial model file).
if [ "$sci_program" = ACACIA ] ; then column_name=ACACIA ; else column_name=$model ; fi
while [ "$current_column_in_header" != "$column_name" ] ; do
    icolmodel=$((icolmodel+1))
    current_column_in_header=$(echo $model_var_header | awk -v icolmodel=$icolmodel '{print $icolmodel}')
done
list_model_variables_extracted=""
list_other_variables_extracted=""
## Converts the input list for variables, by creating the same list with their names in the current model.
for var in $list_model_variables ; do
    var_model=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w $var | awk -v icolmodel=$icolmodel '{print $icolmodel}')
    list_model_variables_extracted=$list_model_variables_extracted" "$var_model
done
for var in $list_other_variables ; do
    var_model=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w $var | awk -v icolmodel=$icolmodel '{print $icolmodel}')
    list_other_variables_extracted=$list_other_variables_extracted" "$var_model
done
list_all_variables_extracted=$list_model_variables_extracted" "$list_other_variables_extracted


case "$config" in
    LMDZORINCA_144_142_39)
	config_simul=LMDZORINCA
	horiz_resol=144x142
	vertic_resol=39 
	;;
    LMDZORINCA_144_142_79)
	config_simul=LMDZORINCA
	horiz_resol=144x142
	vertic_resol=79 
	;;
    160_80_60)
	horiz_resol=160x80
	vertic_resol=60
	;;
    128_64_90)
	horiz_resol=128x64
	vertic_resol=90
	;;
    128_64_60)
	horiz_resol=128x64
	vertic_resol=60
	;;
    192_144_85)
	horiz_resol=192x144
	vertic_resol=85
	;;
esac
DA_suffix=''
if [ "$experiment" = AMIP_Free ] ; then
    # experiment_simul=Free$horiz_resol
    experiment_simul=Bis-free
    if [ $vertic_resol -eq 79 ] ; then
	experiment_simul=Free-NewPhys-79
	DA_suffix='_mth'
    fi
elif [ "$experiment" = AMIP_Nudged ] ; then
    # experiment_simul=Nudge$horiz_resol
    experiment_simul=Nudge${horiz_resol}x$vertic_resol
    if [ $vertic_resol -eq 79 ] ; then
	experiment_simul=test-nudge
    fi
elif [ "$experiment" = Rescaled_Nudged ] ; then
    experiment_simul=Rescaled-${horiz_resol}x$vertic_resol
    if [ $vertic_resol -eq 79 ] ; then
	experiment_simul=Rescaled-$vertic_resol
    fi
    DA_suffix='_mth'
elif [ "$experiment" = Rescaled_Nudged_washout ] ; then
    if [ $vertic_resol -eq 79 ] ; then
	experiment_simul=Rescaled-${vertic_resol}-washout
    fi
    DA_suffix='_mth'
elif [ "$experiment" = Aircraft_0_pct ] ; then
    experiment_simul=Aircraft-0-pct
    DA_suffix='_mth'
elif [ "$experiment" = Rescaled_no_BB ] ; then
    experiment_simul=Zero-BB-rescaled-${horiz_resol}x$vertic_resol
    DA_suffix='_mth'
elif [ "$experiment" = AMIP_MEIC_Nudged ] ; then
    experiment_simul=MEIC-Nudge$horiz_resol
elif [ "$experiment" = AMIP_Nudged_no_LNOx ] ; then
    experiment_simul=Nudge${horiz_resol}x${vertic_resol}-no-LNOx
elif [ "$experiment" = AMIP_Nudged_no_BB ] ; then
    experiment_simul=Nudge${horiz_resol}x${vertic_resol}-no-wb
elif [ "$experiment" = AMIP_free_washout ] ; then
    experiment_simul=washout-free
    DA_suffix='_mth'
elif [ "$experiment" = Prev_MF_no_LNOx_79 ] ; then
    experiment_simul=Prev-mass-flux-no-LNOx-79
    DA_suffix='_mth'
elif [ "$experiment" = ACACIA ] ; then
    experiment_simul=ACACIA
    DA_suffix='_mth'
fi
if [ $vertic_resol -eq 39 ] ; then
    dir_prod_devt=PROD
    file_type_PV=species
elif [ $vertic_resol -eq 79 ] ; then
    dir_prod_devt=DEVT
    file_type_PV=chem
fi
# Compared to this program, simul_output_dir is the very first input directory.
# It leads to the monthly files containing the simulations output.
if [ "$model" = "INCA" ] ; then
    simul_output_dir=$CCCSTOREDIR/IGCM_OUT/LMDZORINCA/$dir_prod_devt/NMHC_AER_S/$experiment_simul/CHM/Output/DA$DA_suffix
else
    simul_output_dir=$CCCSTOREDIR/data/models_output/$model/$config/$sci_program/$time_resol/DA$DA_suffix
fi
extraction_output_dir=$CCCSCRATCHDIR/Pack_code/data/models_output/$model/$config/$experiment/$time_resol
mkdir -p $CCCSCRATCHDIR/Pack_code/data/models_output/$model
mkdir -p $CCCSCRATCHDIR/Pack_code/data/models_output/$model/$config
mkdir -p $CCCSCRATCHDIR/Pack_code/data/models_output/$model/$config/$experiment
mkdir -p $extraction_output_dir
tmpdir=$CCCSCRATCHDIR/tmp_files/$config/$experiment/${yyyymm1}_${yyyymm2}
mkdir -p $CCCSCRATCHDIR/tmp_files
mkdir -p $CCCSCRATCHDIR/tmp_files/$config
mkdir -p $CCCSCRATCHDIR/tmp_files/$config/$experiment
mkdir -p $tmpdir

yyyymm=$yyyymm1
listmonths=$yyyymm
while [ $yyyymm -lt $yyyymm2 ] ; do
    yyyymm=`date -d "${yyyymm}01 1 months" +%Y%m`
    listmonths=${listmonths}" "${yyyymm}
done
if [ ! -e $extraction_output_dir/monthly_files_status.txt ] ; then
    for yyyymm in $listmonths ; do
	echo $yyyymm raw >> $extraction_output_dir/monthly_files_status.txt
    done
fi

if [ "$bin_extract_model_variables" = TRUE ] ; then
    for yyyymm in $listmonths ; do
	## We want the downwind script (make_interpol_iagos.sh) to know which model monthly files are still to be adjusted.
	## If one of these files has been extracted again and thus has to be adjusted one more time,
	## then the current script modifies its status in the monthly_files_status.txt file. It sets it at 'raw'.
	if [ $(cat $extraction_output_dir/monthly_files_status.txt | grep $yyyymm | wc -l) -eq 1 ] ; then
	    line_status_1=${yyyymm}' adjusted'
	    line_status_2=${yyyymm}' raw'
	    line_status_fatal=${yyyymm}' fatal'
	    sed -r -i -e "s/$line_status_1/$line_status_2/g" $extraction_output_dir/monthly_files_status.txt
	    adj2raw_status=$?
	    while [ $adj2raw_status -ne 0 ] ; do
		sleep 1
		sed -r -i -e "s/$line_status_1/$line_status_2/g" $extraction_output_dir/monthly_files_status.txt
		adj2raw_status=$?
	    done
	    sed -r -i -e "s/$line_status_fatal/$line_status_2/g" $extraction_output_dir/monthly_files_status.txt
	    fatal2raw_status=$?
	    while [ $fatal2raw_status -ne 0 ] ; do
		sleep 1
		sed -r -i -e "s/$line_status_fatal/$line_status_2/g" $extraction_output_dir/monthly_files_status.txt
		fatal2raw_status=$?
	    done
	else
	    echo $yyyymm raw >> $extraction_output_dir/monthly_files_status.txt
	fi
	yyyy=${yyyymm:0:4}
	echo extract $yyyymm
	mkdir -p $extraction_output_dir/$yyyy
	file_vmr=$extraction_output_dir/$yyyy/${experiment}_${yyyymm}.nc
	file_PV=$extraction_output_dir/$yyyy/PV_${yyyymm}.nc
	file_P_S=$extraction_output_dir/$yyyy/P_S_${yyyymm}.nc
	file_P=$extraction_output_dir/$yyyy/P_${yyyymm}.nc
	if [ -f "$file_vmr" ] ; then
	    rm $file_vmr
	fi
	if [ "$model" = INCA ] ; then
	    for var in $list_model_variables_extracted ; do
		var_ref_name=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w $var | awk '{print $1}')
		case "$var_ref_name" in 
		    Production_rate_ozone | Loss_rate_ozone)
			file_simul=$simul_output_dir/${experiment_simul}_${yyyymm}01_*_chem.nc
			;;
		    Ozone | Carbon_monoxide | Water_vapour | Nitric_acid | Peroxyacetyl_nitrate | Stratospheric_ozone | Nitrogen_monoxide | Nitrogen_dioxide | Temperature)
			file_simul=$simul_output_dir/${experiment_simul}_${yyyymm}01_*_species.nc
			;;
		esac
		ncks -A -h -v $var $file_simul -o $file_vmr
	    done
	    for var in $list_other_variables_extracted ; do
		var_ref_name=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w $var | awk '{print $1}')
		case $var_ref_name in 
		    Potential_vorticity)
			file_simul=$simul_output_dir/${experiment_simul}_${yyyymm}01_*_${file_type_PV}.nc
			ncks -O -h -v $var $file_simul -o $file_PV
			;;
		    Surface_pressure)
			file_simul=$simul_output_dir/${experiment_simul}_${yyyymm}01_*_species.nc
			ncks -O -h -v $var $file_simul -o $file_P_S
			;;
		    Pressure)
			file_simul=$simul_output_dir/${experiment_simul}_${yyyymm}01_*_species.nc
			ncks -O -h -v $var $file_simul -o $file_P
			;;
		esac
	    done
	else
	    file_simul=$simul_output_dir/$yyyy/${experiment_simul}_${yyyymm}*.nc
	    for var in $list_all_variables_extracted ; do
		var_ref_name=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w $var | awk '{print $1}')
		if [ "$model" = EMAC ] ; then
		    var_short=$var
		    if [ "$var_ref_name" = Surface_pressure ] ; then
			var=${var}_mm
		    else
			var=${var}_dm
		    fi
		fi
		case "$var_ref_name" in 
		    Potential_vorticity)
			ncks -O -h -v $var $file_simul -o $file_PV
			;;
		    Surface_pressure)
			ncks -O -h -v $var $file_simul -o $file_P_S
			;;
		    Pressure)
			ncks -O -h -v $var $file_simul -o $file_P
			;;
		    *)
			ncks -A -h -v $var $file_simul -o $file_vmr
			## Temporary correction, should be deleted for the next UKESM release.
			## The mixing ratios are given in ppm, but they should be in mol/mol.
			## Except water vapour, already in mol/mol.
			if [ "$model" = "UKESM1.1" ] && [ ${var:0:3} = "vmr" ] && [ "$var" != "vmrh2o" ] ; then
			    ncap2 -A -h -s "${var}=${var}*1E-6" $file_vmr -o $file_vmr
			fi
			;;
		esac
		if [ "$model" = EMAC ] ; then # Now standardizing the name of the variables.
		    case "$var_ref_name" in 
			Potential_vorticity)
			    ncrename -v $var,$var_short $file_PV
			    ;;
			Surface_pressure)
			    ncrename -v $var,$var_short $file_P_S
			    ;;
			Pressure)
			    ncrename -v $var,$var_short $file_P
			    ;;
			*)
			    ncrename -v $var,$var_short $file_vmr
			    ;;
		    esac
		fi
	    done
	fi
	last_task=$?
	if [ "$extract_NOy" = TRUE ] ; then
	    echo -e '\e[00;31m Now preparing the NOy variable \e[00m'$yyyymm
	    if [ "$model" = INCA ] ; then
		file_simul=$simul_output_dir/${experiment_simul}_${yyyymm}01_*_species.nc
	    else
		file_simul=$simul_output_dir/$yyyy/${experiment_simul}_${yyyymm}*.nc
	    fi
	    file_vmr=$extraction_output_dir/$yyyy/${experiment}_${yyyymm}.nc
	    ## One first loop with a limited list of NOy species.
	    ivar=0
	    echo -e '\e[00;31m This one with the ACACIA definition: NOx + HNO3 + PAN. \e[00m'
	    for var in $list_NOy_variables_ACACIA ; do
		if [ "$model" = EMAC ] ; then var=${var}_dm ; fi
		echo $var
		file_var=$tmpdir/${var}.nc
		ivar=$((ivar+1))
		file_sum=$tmpdir/vmrnoy_ACACIA.nc
		if [ $ivar -eq 1 ] ; then
		    ncks -O -h -v $var $file_vmr -o $file_sum
		    ncrename -v $var,vmrnoy_ACACIA $file_sum
		else
		    ncks -A -h -v $var $file_vmr -o $file_var
		    ncrename -v $var,vmrnoy_ACACIA $file_var
		    echo 'file_var = '$file_var
		    echo 'file_sum = '$file_sum
		    ncbo -A --op_typ='addition' $file_var $file_sum $file_sum
		    rm $file_var
		fi
	    done
	    ncks -A -h -v vmrnoy_ACACIA $file_sum -o $file_vmr
	    last_task=$?
	    rm $file_sum
	    ## The next loop is in case we need two definitions for NOy.
	    ## Actually, it concerns the INCA model: treated individually, we can take all the NOy species
	    ## present in the model output. Treated in the context of the ACACIA project, we would better
	    ## focus on the species present commonly through all the participating models, i.e. NOx, HNO3 and PAN.
	    ivar=0
	    if [ "$list_NOy_variables" != "$list_NOy_variables_ACACIA" ] ; then
		echo -e '\e[00;31m This one with the model own definition: a lot of N species. \e[00m'
		for var in $list_NOy_variables ; do
		    echo $var
		    file_var=$tmpdir/${var}.nc
		    ivar=$((ivar+1))
		    file_sum=$tmpdir/vmrnoy.nc
		    if [ $ivar -eq 1 ] ; then
			ncks -O -h -v $var $file_simul -o $file_sum
			ncrename -v $var,vmrnoy $file_sum
		    else
			ncks -A -h -v $var $file_simul -o $file_var
			if [ "$var" = 'vmrn2o5' ] ; then ## 2 nitrogen atoms: counts double in the NOy sum.
			    ncbo -A --op_typ='addition' $file_var $file_var $file_var
			    echo 'Doubling N2O5 in the NOy count'
			fi
			ncrename -v $var,vmrnoy $file_var
			echo 'file_var = '$file_var
			echo 'file_sum = '$file_sum
			ncbo -A --op_typ='addition' $file_var $file_sum $file_sum
			rm $file_var
		    fi
		done
		ncks -A -h -v vmrnoy $file_sum -o $file_vmr
		last_task=$?
		rm $file_sum
	    fi
	fi
	if [ $last_task -ne 0 ] ; then
	    echo 'Incomplete extraction for the current month: '$yyyymm
	    sed -r -i -e "s/$line_status_1/$line_status_fatal/g" $extraction_output_dir/monthly_files_status.txt
	    sed -r -i -e "s/$line_status_2/$line_status_fatal/g" $extraction_output_dir/monthly_files_status.txt
	    exit
	fi
    done
fi
if [ -n "$tmpdir" ] ; then ## Just to make sure that I do not launch a rm -rf command with no argument.
    rm -rf $tmpdir
fi
echo -e '\e[00;31m Extraction: done.'
echo -e 'Destination directory: \e[00m'$extraction_output_dir
