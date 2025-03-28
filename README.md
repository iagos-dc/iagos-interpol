This package aims to assess model simulations with IAGOS measurements.
The main task is a shell script (make_interpol_iagos.sh)
that treats the data, both simulation outputs and observations, to build several different products:
- seasonal climatologies averaged vertically (by default: through the cruise altitudes only)
- seasonal climatologies averaged vertically and in longitude (by default: through the cruise altitudes only, and through predefined longitude intervals that correspond to South America-Atlantic, West Africa, and South Asia)
- time series above a set of predefined regions, with one record dimension for time, and another dimension
for the season. The latter is five element-sized: the four seasons and the whole year.
All the routines are gathered in the path Pack_code/routines/grid_IAGOS_model.
###############
How to use
###############
#
Before launching, choose your parameters.
#
For each of the tasks described above, the main program starts with all the options that the user can change, with their explanations.
For example, the simulation characteristics (the model name, the horizontal grid resolution,
the time resolution of the output), the variables to process or the start and end of the period to process.
The end of this 'user' zone is clearly specified. Everything that comes after is automatic.
### Please read the comments given in this 'user' zone before running its program. ### 
They explain the program behaviour and the utility/consequences of each parameter you fix.

### Advanced specifications for bold users ###
If the user wants to increase the spectrum of variables or models that can be asked to the program,
then he has to add the corresponding information into the tables at this effect (in the folder grid_IAGOS_model):
table_variable_names_in_... and table_variable_units_in_...
There are 4 of them, two for observations and two for the simulation output. They contain information
on modelled/measured variables, necessary to read and convert automatically their fields from the data files.
- Concerning the simulation, the table_variable_name file contains the name of each variable, one column per model.
It also contains an output name column for homogenisation between the models, and a standard names column made to
match with the user request. The variable names in the latter are as clear as possible.
Concerning the observations, the table_variable_name file follows the same organisation, with one column per IAGOS package
because the variable names differ from one package to another. Unless the database is updated, the user is not
concerned by this file.
- The table_variable_units file contains two columns per model/package: one for the units for each model/package,
one for the conversion factor to apply to homogenise the outputs. For variables which conversion does not consist in a product,
 the factor has to be set at 1 and the conversion has to be added manually in the code, as it is done for temperature in IAGOS files.

###########
Description
###########

1/
The script make_interpol_iagos.sh starts with a preliminary part, extracting additional model variables,
as the potential vorticity, the surface pressure and the main pressure if it exists.
This first part also calculates new variables from the model output, like the O3/CO ratio, or the NOx mixing ratio.

To activate it, the user has to specify bin_extract_last_variables=TRUE

2/
A short second step consists of changing the coordinate axes, in order to make sure that
there is no artefact linked to a wrong axis direction during the routines execution.
Here, the user has to indicate which axes have to be modified to fit the required format.
This required format is defined as follow:
- the longitude axis starts at 180 W and stops at 180 E.
- the latitude axis starts at 90 S and stops at 90 N.
- the vertical axis starts at the top of atmosphere and stops at the surface.

Three logical bins are proposed to the user then, describing the axis setup in the model output as is submitted: 
- lon_axis_0_360
- lat_axis_90N_90S
- P_axis_bottom_to_top
If at least one of them is TRUE, then the corresponding coordinate will be automatically modified.

3/
A fortran routine is called if the user wants to interpolate the model output onto a new grid,
with constant pressure levels.
In this case, one will have to specify this grid characteristics, and follow the instructions
in the make_interpol_iagos header.
To activate it, the user has to specify bin_change_grid=TRUE. Then one will have to indicate
which horizontal and vertical resolution one wants.
The required subroutines and modules (in the same repertory) are:
- global_var_interpol_mod.f90
- functions_mod.f90
- ncdf.f90 (in the following repertory: ./eumetsat_ropp/ropp_io/ncdf/)
- log_P_interp.f90
- vertical_interp.f90

4/
The main purpose of this script: the Fortran routine interpol_iagos.f90,
which projects the IAGOS data (in NetCDF files) onto the model grid,
and computes monthly means for several measured variables, either with a daily or monthly resolution.
It can read the surface pressure fields and the vertical coefficients to compute pressure levels,
but it can also read the pressure variable directly if it is available.
Then, it reads each IAGOS file and projects linearly each observation onto the grid.
It also computes monthly means on each grid cell, using spatial weighting factors based on distances.
interpol_iagos.f90 generates monthly gridded data, both from IAGOS and from the model, with one file
per layer and per month.
The required subroutines and modules (in the same repertory) are:
- global_var_interpol_mod.f90
- functions_mod.f90
- ncdf.f90 (in the following repertory: ./eumetsat_ropp/ropp_io/ncdf/)
- interpol_iagos.f90
- mean_values_interpol.f90

5/
make_interpol_iagos.sh then generates regional output files on the specified UTLS
grid levels. For each region, there are:
- one IAGOS file, consisting in monthly-mean IAGOS data extrapolated
on the model grid. (IAGOS-DM)
- one model file, consisting in the results from the simulation after applying a mask, depending on the IAGOS sampling (model-M).
The required subroutines and modules (in the same repertory) are:
- global_var_interpol_mod.f90
- functions_mod.f90
- ncdf.f90 (in the following repertory: ./eumetsat_ropp/ropp_io/ncdf/)
- regional_average_IAGOS.f90

6/
make_interpol_iagos.sh then calculates mean climatologies for each season, and over the whole year. The output format is 3D.

The last two steps are made to prepare not-voluminous files ready for the analysis.

7/ make_interpol_iagos.sh calls the Fortran routine column_average_IAGOS.f90 to average the 3D climatological output vertically, through predefined pressure boundaries. The output format is then 2D.
The required subroutines and modules (in the same repertory) are:
- global_var_interpol_mod.f90
- functions_mod.f90
- ncdf.f90 (in the following repertory: ./eumetsat_ropp/ropp_io/ncdf/)
- column_average_IAGOS.f90

8/ make_interpol_iagos.sh calls the Fortran routine zonal_cross_sections_IAGOS.f90 to average the 3D climatological output vertically through predefined pressure boundaries, but also through longitudes over predefined longitude domains. The output format is then 1D, for each longitude domain.
The required subroutines and modules (in the same repertory) are:
- global_var_interpol_mod.f90
- functions_mod.f90
- ncdf.f90 (in the following repertory: ./eumetsat_ropp/ropp_io/ncdf/)
- zonal_cross_sections_IAGOS.f90

######################
In make_interpol_iagos.sh, the header shows you the dates delimitating
the period you want to study. You also have to pick the model/simulation you
want to compare to the IAGOS data, and the location of the output
files. The header also shows several binaries. Each one of them
that is TRUE activates one part of the script, coded in the file content.

1/ Call of interpol_iagos.f90:
The subroutines you need in the same repertory are:
- INT2CHAR.f90
- REAL2CHAR.f90
- compute_pressures.f90
- mean_values_interpol.f90
- standard_deviations_interpol.f90
- tonetcdf.F90
with the modules:
- global_var_interpol_mod.f90
- functions_mod.f90
The input repertories are:
- data/reference_values/this_model/, containing notably the two coefficient tables 
zvamh_nlev.txt and zvbmh_nlev.txt, to be able to generate the 3D pressure fields.
(nlev stands for the number of vertical grid levels)
- data/IAGOS_database/this_file_format/, where the IAGOS files (one file per flight)
are available until December 2017 for now. this_file_format can be either NetCDF or ASCII.
Each repertory corresponds to a month, they are gathered in yearly repertories.
- data/models_output/this_model/this_configuration/this_experiment/this_time_resolution

2/ Preparing regional time series, ready for the analysis.
2.1/ Using an R routine to derive the regions coordinates with the grid resolution, both
in lon/lat degrees and in their corresponding indexes in the model grid.
2.2/ Using the latter as input for nco commands to extract the grid cells for each region.
2.3/ Call of average_IAGOS.f90:
The subroutines you need in the same repertory are:
- INT2CHAR.f90
- REAL2CHAR.f90
- tonetcdf_time_series.F90
with the modules:
- global_var_interpol_mod.f90
- functions_mod.f90

3/ Preparing climatologies, ready for the analysis.
Using nco commands to concatenate throughout the years and to compute seasonal means,
and yearly means also. The first and last months defined by the user are readapted
in a way that allows an equivalent period sampling between the seasons.
######################

Author: Yann Cohen
yann.cohen.09@gmail.com
