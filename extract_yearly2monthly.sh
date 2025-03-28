#!/bin/bash
## This script aims to gather the queried variables from a simulation output files,
## into monthly files ready for being treated by the make_interpol_iagos script.

#MSUB -q skylake
#MSUB -A gen2201
#MSUB -r yr2mth
#MSUB -e $CCCSCRATCHDIR/out_yearly2monthly
#MSUB -o $CCCSCRATCHDIR/out_yearly2monthly
#MSUB -n 1
#MSUB -c 1
#MSUB -T 86400
#MSUB -m work,scratch,store
#MSUB -Q normal
#MSUB -x
set +x

## Enter your arguments.
# experiment_simul=Rescaled-144x142x39
# experiment_simul=TOAR-ANT-1995
experiment_simul=ACACIA
model=MOZART3
# devt_prod=DEVT # or PROD
yyyy1=1997
yyyy2=2007

yyyymm1=${yyyy1}01 ; yyyymm2=${yyyy2}12 ; 
yyyymm=$yyyymm1
listmonths=$yyyymm
while [ $yyyymm -lt $yyyymm2 ] ; do
    yyyymm=`date -d "${yyyymm}01 1 months" +%Y%m`
    listmonths=$listmonths" "$yyyymm
done

if [ "$model" = "INCA" ] ; then
    list_files_CHM=_species' '_chem' '_dep
    var_time=time_counter

    case $experiment_simul in
	Aircraft-0-pct | Free-NewPhys-79 | No-LNOx-79 | Prev-mass-flux-no-LNOx-79 | Rescaled-79 | Rescaled-79-washout | test-nudge | washout-free)
	    devt_prod=DEVT ;;
	*)
	    devt_prod=PROD ;;
    esac
    yearly_dir=$CCCSTOREDIR/IGCM_OUT/LMDZORINCA/$devt_prod/NMHC_AER_S/$experiment_simul/CHM/Output/DA
    monthly_dir=$CCCSTOREDIR/IGCM_OUT/LMDZORINCA/$devt_prod/NMHC_AER_S/$experiment_simul/CHM/Output/DA_mth
    mkdir -p $monthly_dir
else
    list_files_CHM='no_need'
    var_time=time
    yearly_dir=$CCCSTOREDIR/models_output/$model/$experiment_simul
    monthly_dir=$CCCSTOREDIR/models_output/$model/${experiment_simul}_mth
    mkdir -p $monthly_dir
fi


for yyyymm in $listmonths ; do
    yyyy=${yyyymm:0:4}
    if [ ${yyyymm:4:2} = 01 ] ; then
	day_1=1
	day_N=0
    fi
    yyyymm01_next_month=`date -d "${yyyymm}01 1 months" +%Y%m%d`
    yyyymmdd_end_month=`date -d "${yyyymm01_next_month} -1 days" +%Y%m%d`
    ndays_this_month=${yyyymmdd_end_month:6:2}
    day_N=$(($day_N+$ndays_this_month))
    for file_CHM in $list_files_CHM ; do
	echo $yyyymm $day_1 'to' $day_N 'cause ndays=' $ndays_this_month 'in' $file_CHM
	if [ "$model" = "INCA" ] ; then
	    yearly_file=$yearly_dir/${experiment_simul}_${yyyy}0101_${yyyy}1231_1D_inca${file_CHM}.nc
	    monthly_file=$monthly_dir/${experiment_simul}_${yyyymm}01_${yyyymm}${ndays_this_month}_1D_inca${file_CHM}.nc
	else
	    yearly_file=$yearly_dir/${model}_A1_${yyyy}.nc
	    monthly_file=$monthly_dir/${experiment_simul}_${yyyymm}.nc
	fi
	ncks -O -F -d $var_time,$day_1,$day_N $yearly_file -o $monthly_file
    done
    day_1=$(($day_N+1))
done
