# This is a placeholder/template showing how to load ISIMIP Water Quality data.
# Do not run this without valid local NetCDF files.

using TwoTimescaleResilience
using NCDatasets
using CairoMakie

# Edit these paths and variable names to match your local ISIMIP Water Quality files.
WT_file       = raw"C:\path\to\watertemp_monthlyAvg_1980_2019.nc"
BOD_file      = raw"C:\path\to\bod_monthlyAvg_1980_2019.nc"
TDS_file      = raw"C:\path\to\tds_monthlyAvg_1980_2019.nc"
FC_file       = raw"C:\path\to\fc_monthlyAvg_1980_2019.nc"
Nutrient_file = raw"C:\path\to\nutrient_monthlyAvg_1980_2019.nc"
Chemical_file = raw"C:\path\to\chemical_monthlyAvg_1980_2019.nc"
Plastic_file  = raw"C:\path\to\plastic_monthlyAvg_1980_2019.nc"

WT_var       = "watertemp"
BOD_var      = "bod"
TDS_var      = "tds"
FC_var       = "fc"
Nutrient_var = "nutrient"
Chemical_var = "chemical"
Plastic_var  = "plastic"

fill_missing_variables_with_zero = true

println("This template is not meant to be run directly. See nc_monthly_longterm_isimip_moa_deb_inspection.jl for a full example.")
