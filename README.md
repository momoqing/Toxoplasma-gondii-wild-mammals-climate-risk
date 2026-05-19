# T. gondii infection risk in wild mammals under climate change

This repository contains the study-level data and R analysis scripts supporting the manuscript:

**Persistent tropical hotspots and regional redistribution of *Toxoplasma gondii* infection risk in wild mammals under climate change**

This repository is prepared for research transparency, reproducibility, and post-publication data sharing.

## Repository status

This repository accompanies a manuscript under submission. A permanent archived version is available through Zenodo:
        
        
        
https://doi.org/10.5281/zenodo.20226412
        
        .

Recommended data-sharing statement:

> Extracted study-level data and analysis code will be made available in this public repository upon publication. The repository includes code for meta-analysis, Bayesian spatiotemporal modeling, machine-learning prediction, future climate-risk projection, climate novelty analysis, positive-change spatial extent analysis, response-surface analysis, SHAP-dependence diagnostics, and figure generation. A permanent archived version with a DOI will be provided upon publication.

## Contents

Suggested repository structure:

```text
tgondii-wild-mammals-climate-risk/
│
├── README.md
├── LICENSE
│
├── data/
│   └── ml.xlsx
│
├── code/
│   ├── meta-analysis R script.R
│   ├── Bayesian spatiotemporal model code.R
│   ├── machine-learning scripts.R
│   ├── Spatial block CV.R
│   ├── future projection scripts.R
│   ├── climate novelty analys.R
│   ├── Training space.R
│   └── response surface and SHAP-dependence diagnostics.R
│
└── outputs/
    ├── tables/
    ├── figures/
    └── model_outputs/
```

File names can be renamed for clarity before public release. For example:

```text
01_meta_analysis.R
02_bayesian_spatiotemporal_model.R
03_machine_learning_model_comparison.R
04_spatial_block_cross_validation.R
05_future_climate_projection.R
06_climate_novelty_analysis.R
07_positive_change_extent_analysis.R
08_response_surface_and_shap_diagnostics.R
09_figure_generation.R
```

## Data

The primary dataset is `ml.xlsx`, which contains extracted study-level information from published studies of *Toxoplasma gondii* infection in wild mammals.

The workbook contains two sheets:

### Sheet1

Study-level analytical dataset with the following variables:

| Variable | Description |
|---|---|
| `Author` | First author and year label of the included study |
| `Year` | Publication or sampling year used in analysis |
| `Event` | Number of positive samples |
| `Total` | Total number of tested samples |
| `Species` | Host species name |
| `Order` | Host mammalian order |
| `Family` | Host family |
| `Diet` | Dietary guild |
| `Habitat` | Host habitat category |
| `Sample_Type` | Biological sample type |
| `Methods` | Diagnostic method |
| `Continent` | Continent of sampling location |
| `Country` | Country of sampling location |
| `Rainfall` | Mean daily precipitation used in the machine-learning model |
| `Temperature` | Annual mean temperature |
| `Longitude` | Longitude of sampling location |
| `Latitude` | Latitude of sampling location |
| `Prevalence` | Observed prevalence, calculated as `Event / Total` |

### Sheet2

Host trait lookup table with the following variables:

| Variable | Description |
|---|---|
| `Species` | Host species name |
| `Order` | Host mammalian order |
| `Family` | Host family |
| `Diet_Type` | Dietary guild assignment |

## Analysis workflow

The analysis workflow includes the following steps:

1. **Meta-analysis**  
   Estimate pooled prevalence, subgroup prevalence, publication bias, and heterogeneity.

2. **Bayesian spatiotemporal modeling**  
   Reconstruct reported hotspot dynamics while accounting for sampling effort.

3. **Machine-learning model comparison**  
   Compare regression algorithms and retain the best-performing random forest model for interpretation and projection.

4. **Spatial validation**  
   Evaluate geographic transferability using spatial block cross-validation and leave-one-continent-out validation.

5. **Future climate-risk projection**  
   Predict contemporary and future infection risk using WorldClim 2.1 baseline and CMIP6 bioclimatic layers under four SSP scenarios and two future periods.

6. **Climate novelty analysis**  
   Assess whether future projections fall outside the observed training climate space.

7. **Positive-change spatial extent analysis**  
   Quantify area-weighted spatial extent of positive relative change across SSP scenarios, periods, and latitude bands.

8. **Temperature–precipitation response surface**  
   Visualize the joint response of predicted risk to annual mean temperature and mean daily precipitation.

9. **SHAP-dependence diagnostics**  
   Examine targeted model-based interaction diagnostics for precipitation × habitat, precipitation × dietary guild, and temperature × latitude band.

10. **Figure generation**  
   Generate main-text and supplementary figures.

## External data sources

This repository does not redistribute large external raster datasets. Users should download these datasets from their original providers and update file paths in the R scripts before running the projection workflow.

External datasets include:

- WorldClim 2.1 baseline bioclimatic variables
- WorldClim 2.1 CMIP6 downscaled bioclimatic variables
- BCC-CSM2-MR and EC-Earth3-Veg climate model outputs
- Shared Socioeconomic Pathway scenarios SSP1-2.6, SSP2-4.5, SSP3-7.0, and SSP5-8.5

In the projection workflow, WorldClim bio1 was used as annual mean temperature. WorldClim bio12, which represents annual precipitation, was converted to mean daily precipitation by dividing by 365 to match the rainfall scale used in the machine-learning dataset.

## Software requirements

Analyses were conducted in R. The following R packages are commonly required across scripts:

```r
dplyr
tidyr
ggplot2
metafor
brms
rstan
caret
randomForest
fastshap
terra
sf
RANN
readxl
RColorBrewer
scales
```

Package versions should be recorded before final release. Users are encouraged to run the scripts in the order listed above.

## Reproducibility notes

- File paths in scripts should be updated to match the user's local directory structure.
- Large raster files are not included in this repository.
- Some scripts require model objects generated by earlier scripts.
- Future projections depend on the corrected precipitation scale, with bio12 divided by 365.
- Results in the manuscript were interpreted as broad-scale, scenario-based estimates rather than precise local forecasts.

## Citation

If you use these data or code, please cite the associated manuscript after publication.

Suggested citation placeholder:

> Zhu X-K, Wang J-L, Yang T, et al. Persistent tropical hotspots and regional redistribution of *Toxoplasma gondii* infection risk in wild mammals under climate change. Submitted.

A permanent DOI will be added after publication and repository archiving.

## Contact

For questions about the data or code, please contact:

- Xin-Kun Zhu: [zxk19971222@163.com]
- Wei Cong: [congwei2016@sdu.edu.cn]
- Hany M. Elsheikha: [elsheikha@nottingham.ac.uk]

## License

Code is released under the MIT License. Extracted study-level data are provided for academic transparency and reproducibility; users should cite the associated manuscript and the original studies where appropriate.
