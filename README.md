# ecRhtoGEM
### Brief repository description:
* This repository contains the work on metabolic engineering of oleaginous, non-conventional yeast *R.toruloides* using genome-scale modelling.

* Repository consists of code, data, model files and result files that can be also used for visualising metabolic networks.

#### Model files (xml, xlsx, mat, json)
* As a starting point for this work, a metabolic genome-scale model, named **rhto-GEM**, developed by [Tiukova et al. 2019](https://onlinelibrary.wiley.com/doi/abs/10.1002/bit.27162), was used.
* The initial model files used in this work are available [here](https://github.com/SysBioChalmers/rhto-GEM).

#### Scripts (ipynb, m)
* The scripts were written in Matlab R2018b and Jupyter Notebook (Python) from [Miniconda](https://docs.conda.io/en/latest/miniconda.html) software (as an alternative to full [Anaconda](https://docs.anaconda.com/anaconda/install/)).
* The modelling was done with the [COBRApy](https://github.com/opencobra/cobrapy), RAVEN and COBRA toolboxes.
* The documentation of COBRApy is browseable online at [readthedocs](https://cobrapy.readthedocs.io/en/stable/).

#### Visualisations (json)
* Metabolic network maps were constructed using [Escher](https://escher.github.io/#/) software.
* The results were directly visualised in Jupyter Notebook using Escher Python package.
* The documentation of Escher is browseable [here](https://escher.readthedocs.io/en/latest/).


### Instructions:
* `edit_rhtoGEM.m`: run this script to prepare the initial rhto model for upgrading with enzymatic constraints.
* `reconstruct_ecRhtoGEM.m`: generate carbon source specific enzyme-constrained models with proteomics data integration. Models are reconstructed based on optimized parameters, including manually curated enzymatic kcat values as provided in `manualModifications.m`. Model reconstruction involves the following steps: 1) `geckomat/enhanceGEM.m` pipeline which creates the first version of ec-model with an enzyme pool used for screening required kcat modifications; 2) `gecomat/utilities/integrate_proteomics/generate_protModels.m` pipeline which sets constraints on high-quality measured individual enzymes from the provided proteomics dataset; 3) adding ribosomal subunits to the ec-models. In `results/generate_protModels_pipeline` is written which enzyme abundances were automatically modified in order to reach the experimental conditions. In `results/ribosome_integration` is written distribution of average ribosomal subunit abundances and also which subunit abundances were flexibilized in order to reach experimental conditions.
* `analyzeUsage.m`: runs `enzymeUsage.m` to generate enzyme usage reports that are stored in `results/model_simulation`.
* `boxplotEnzymeUsage.R`: takes enzyme usage data and plots in various ways.
* `geckopy_data_prepare.ipynb` and `geckopy_json.ipynb`: take flux balance analysis and enzyme usage data and saves for visualisation in metabolic maps.

##### Last update:
2021-01-21
