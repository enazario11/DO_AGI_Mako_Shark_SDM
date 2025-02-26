
# DO_AGI_Mako_Shark_SDM

(add DOI after publication)

This repository contains the code for our paper:

> Nazario EC, Lezama-Ochoa N, Czapanskiy MF, Dewar HH, Preti A, Fredston AL, Pinsky ML, Pozo Buil M, Hazen EL (2025). _`Dissolved oxygen and metabolic parameters improve species distribution models for a marine predator_. insert journal here (DOI link coming soon)

### How to cite

Please cite code from this page as:

> Nazario EC, Lezama-Ochoa N, Czapanskiy MF, Dewar HH, Preti A, Fredston AL, Pinsky ML, Pozo Buil M, Hazen EL (2025). _R code for Dissolved oxygen and metabolic parameters improve species distribution models for a marine predator_. Online at (https://github.com/enazario11/DO_AGI_Mako_Shark_SDM)

## Contents

  - [:file\_folder: 1_species_data](/1_species_data): Species data processing
  - [:file\_folder: 2_environmental_data](/2_environmental_data): Environmental data processing
  - [:file\_folder: 3_extract_environmental_data](/3_extract_environmental_data): Environmental data extracted for each observed and pseudo-absence location
  - [:file\_folder: 4_calculate_AGI](/4_calculate_AGI): Aerobic Growth Index (AGI) calculations
  - [:file\_folder: 5_BRT_fit_and_validate](/5_BRT_fit_and_validate): BRT model fitting and validation, including the spatiotemporal analysis
  - [:file\_folder: 6_suitability_map_predictions](/6_suitability_map_predictions): Final BRT models used to generate habitat suitability index map predictions
  - [:file\_folder: 7_diet_data](/7_diet_data): Diet data processing and GII calculations
  - [:file\_folder: 8_ms_figures](/8_ms_figures): Figures that appear in the manuscript "Dissolved oxygen and metabolic parameters improve species distribution models for a marine predator"
  - [:file\_folder: functions](/functions): Functions called throughout the rest of the data processing and analysis folders. 

Folders do not include all code used for model and data exploration. If there is an analysis you are interested in that is not included here, please don't hesitate to reach out. 

## Data

The raw mako shark data is available for [download on Dryad](https://doi.org/10.5061/dryad.31zcrjdxz). Code for the mako shark location data processing can be found in the [1_species_data folder](https://github.com/nazarioe11/DO_AGI_Mako_Shark_SDM/tree/main/1_species_data).

Environmental data is publicly accessible and available for download using [Copernicus Marine Service dashboard](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description), [Mercator Ocean International](https://www.mercator-ocean.eu/solutions-expertises/acceder-aux-donnees-numeriques/details-produit/?offer=4217979b-2662-329a-907c-602fdc69c3a3&system=89090a76-f3e8-36ea-252d-5bb624b22b67), and [GEBCO Bathymetry](https://www.gebco.net/data_and_products/gridded_bathymetry_data/#area). 
