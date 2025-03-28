#!/bin/sh
## The latitude and vertical coordinates I used are indexed northwards and downwards.
## It corresponds to the MOCAGE configuration, the model I first used.
## The LMDZ-OR-INCA configuration, on the contrary, has a reverse indexing for latitude and vertical levels.
## For consistence, this script aims at imposing the MOCAGE configuration to the INCA outputs.

## Reading arguments.
model=$1
config=$2
experiment=$3
sci_program=$4
yyyymm1=$5
yyyymm2=$6
model_output_dir=$7
lon_name=$8
lat_name=$9
lev_name=${10}
extract_PV=${11}
extract_P=${12}
shift_lon_axis=${13}
reverse_lat_axis=${14}
reverse_P_axis=${15}
echo '###########################################'
echo 'Longitude name in model files:' $lon_name
echo 'Latitude name in model files:' $lat_name
echo 'Vertical level name in model files:' $lev_name
echo 'PV extraction:' $extract_PV
echo 'P extraction:' $extract_P

yyyymm=$yyyymm1
listmonths=$yyyymm
while [ $yyyymm -lt $yyyymm2 ] ; do
    yyyymm=`date -d "${yyyymm}01 1 months" +%Y%m`
    listmonths=$listmonths" "$yyyymm
done

file_do_it_if_raw=$model_output_dir/monthly_files_status.txt
if [ ! -f "$file_do_it_if_raw" ] ; then
    for yyyymm in $listmonths ; do echo $yyyymm raw >> $file_do_it_if_raw ; done
fi
## Making sure we modify the axes once only.
## Careful: in the CMIP6 program, the latitude axis has already the right direction.
## However, the longitude dimension is in the American coordinate system (from 0 to 360)
## and must be brought back to reas... to the European system (from -180 to 180).
for yyyymm in $listmonths ; do
    yyyy=${yyyymm:0:4}
    file_vmr=$model_output_dir/$yyyy/${experiment}_${yyyymm}.nc
    file_PV=$model_output_dir/$yyyy/PV_${yyyymm}.nc
    file_P_S=$model_output_dir/$yyyy/P_S_${yyyymm}.nc
    file_P=$model_output_dir/$yyyy/P_${yyyymm}.nc
    if [ $(cat $file_do_it_if_raw | grep $yyyymm | wc -l) -eq 0 ] ; then
	echo $yyyymm raw >> $file_do_it_if_raw
    fi
    status_overall=0 ; n_adj_done=0
    if [ $(cat $file_do_it_if_raw | grep ${yyyymm}' raw' | wc -l) -eq 1 ] ; then
	echo 'Adjusting the model file for the month' $yyyymm
	if [ "$reverse_P_axis" = TRUE ] ; then ## The vertical axis in the EMAC and MOZART3 output is already indexed downward.
	    ### Vertical grid levels
	    ncpdq -O -a '-'$lev_name $file_vmr $file_vmr
	    status_z_vmr=$? ; status_overall=$((status_overall+status_z_vmr)) ; n_adj_done=$((n_adj_done+1))
	    if [ "$extract_PV" = "yes" ] ; then
		ncpdq -O -a '-'$lev_name $file_PV $file_PV
		status_z_PV=$? ; status_overall=$((status_overall+status_z_PV)) ; n_adj_done=$((n_adj_done+1))
	    fi
	    if [ "$extract_P" = "yes" ] ; then
		ncpdq -O -a '-'$lev_name $file_P $file_P
		status_z_P=$? ; status_overall=$((status_overall+status_z_P)) ; n_adj_done=$((n_adj_done+1))
	    fi
	    ### P_S files: no vertical grid.
	    ### Latitude or longitude, depending on the file format.
	fi
	if [ "$sci_program" = "CMIP6" ] ; then
	    ### VMR files
	    ncks -O --msa -d $lon_name,180.,360. -d $lon_name,0.,179.0 $file_vmr $file_vmr 
	    txt_command='ncap2 -O -s '"'"'where('${lon_name}' >= 180) '${lon_name}'='${lon_name}'-360'"'"' $file_vmr $file_vmr'
	    eval $txt_command
	    ### P files
	    if [ "$extract_P" = "yes" ] ; then
		ncks -O --msa -d $lon_name,180.,360. -d $lon_name,0.,179.0 $file_P $file_P
		txt_command='ncap2 -O -s '"'"'where('${lon_name}' >= 180) '${lon_name}'='${lon_name}'-360'"'"' $file_P $file_P'
		eval $txt_command
	    fi
	    ### P_S files
	    ncks -O --msa -d $lon_name,180.,360. -d $lon_name,0.,179.0 $file_P_S $file_P_S
	    txt_command='ncap2 -O -s '"'"'where('${lon_name}' >= 180) '${lon_name}'='${lon_name}'-360'"'"' $file_P_S $file_P_S'
	    eval $txt_command
	    ### No need for PV, because not in CMIP6 anyway.
	if [ "$reverse_lat_axis" = TRUE ] ; then
	    ncpdq -O -a '-'$lat_name $file_vmr $file_vmr
	    status_lat_vmr=$? ; status_overall=$((status_overall+status_lat_vmr)) ; n_adj_done=$((n_adj_done+1))
	    if [ "$extract_PV" = "yes" ] ; then
		ncpdq -O -a '-'$lat_name $file_PV $file_PV
		status_lat_PV=$? ; status_overall=$((status_overall+status_lat_PV)) ; n_adj_done=$((n_adj_done+1))
	    fi
	    if [ "$extract_P" = "yes" ] ; then
		ncpdq -O -a '-'$lat_name $file_P $file_P
		status_lat_P=$? ; status_overall=$((status_overall+status_lat_P)) ; n_adj_done=$((n_adj_done+1))
	    else
		ncpdq -O -a '-'$lat_name $file_P_S $file_P_S
		status_lat_P_S=$? ; status_overall=$((status_overall+status_lat_P_S)) ; n_adj_done=$((n_adj_done+1))
	    fi
	fi
	if [ "$shift_lon_axis" = TRUE ] ; then
	    ### VMR files
	    ncks -O --msa -d $lon_name,180.,360. -d $lon_name,0.,179.0 $file_vmr $file_vmr 
	    txt_command='ncap2 -O -s '"'"'where('${lon_name}' >= 180) '${lon_name}'='${lon_name}'-360'"'"' $file_vmr $file_vmr'
	    eval $txt_command
	    status_lon_vmr=$? ; status_overall=$((status_overall+status_lon_vmr)) ; n_adj_done=$((n_adj_done+1))
	    ### PV files
	    if [ "$extract_PV" = "yes" ] ; then
		ncks -O --msa -d $lon_name,180.,360. -d $lon_name,0.,179.0 $file_PV $file_PV
		txt_command='ncap2 -O -s '"'"'where('${lon_name}' >= 180) '${lon_name}'='${lon_name}'-360'"'"' $file_PV $file_PV'
		eval $txt_command
		status_lon_PV=$? ; status_overall=$((status_overall+status_lon_PV)) ; n_adj_done=$((n_adj_done+1))
	    fi
	    ### P files
	    if [ "$extract_P" = "yes" ] ; then
		ncks -O --msa -d $lon_name,180.,360. -d $lon_name,0.,179.0 $file_P $file_P
		txt_command='ncap2 -O -s '"'"'where('${lon_name}' >= 180) '${lon_name}'='${lon_name}'-360'"'"' $file_P $file_P'
		eval $txt_command
		status_lon_P=$? ; status_overall=$((status_overall+status_lon_P)) ; n_adj_done=$((n_adj_done+1))
	    fi
	    ### P_S files
	    ncks -O --msa -d $lon_name,180.,360. -d $lon_name,0.,179.0 $file_P_S $file_P_S
	    txt_command='ncap2 -O -s '"'"'where('${lon_name}' >= 180) '${lon_name}'='${lon_name}'-360'"'"' $file_P_S $file_P_S'
	    eval $txt_command
	    status_lon_P_S=$? ; status_overall=$((status_overall+status_lon_P_S)) ; n_adj_done=$((n_adj_done+1))
	fi
	line_status_raw=$yyyymm' raw'
	line_status_adj=$yyyymm' adjusted'
	line_status_par=$yyyymm' partial'
	## Writing the new status in the info file.
	if [ $status_overall -eq 0 ] & [ $n_adj_done -gt 0 ] ; then
	    sed -r -i -e "s/$line_status_raw/$line_status_adj/g" $file_do_it_if_raw
	elif [ $status_overall -ne 0 ] & [ $n_adj_done -gt 0 ] ; then
	    sed -r -i -e "s/$line_status_raw/$line_status_par/g" $file_do_it_if_raw
	fi
	raw2adj_status=$?
	while [ $raw2adj_status -ne 0 ] ; do
	    sleep 1
	    sed -r -i -e "s/$line_status_raw/$line_status_adj/g" $file_do_it_if_raw
	    raw2adj_status=$?
	done
    elif [ $(cat $file_do_it_if_raw | grep ${yyyymm}' partial' | wc -l) -eq 1 ] ; then
	echo $yyyymm' partially adjusted: I skip'
    elif [ $(cat $file_do_it_if_raw | grep ${yyyymm}' fatal' | wc -l) -eq 1 ] ; then
	echo $yyyymm' previously fatal: I skip'
    fi
done

echo 'End of adjust_model_output_files.sh script.'
echo '###########################################'
