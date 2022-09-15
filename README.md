# ecRhtoGEM
### Brief repository description:
* This repository contains the work on enzyme-constrained genome-scale model of oleaginous, non-conventional yeast *Rhodotorula toruloides*.

* Repository consists of code, data, model files and result files. Result visualisation is not provided.

#### Model files (xml, yml, mat, json)
* As a starting point for this work, a metabolic genome-scale model, named **rhto-GEM**, developed by [Tiukova et al. 2019](https://onlinelibrary.wiley.com/doi/abs/10.1002/bit.27162), was used.
* The initial model files used in this work are available [here](https://github.com/SysBioChalmers/rhto-GEM).

#### Scripts (m)
* The scripts were written in Matlab.
* Enzymatic constraints were applied and models were generated using the [GECKO](https://github.com/SysBioChalmers/GECKO) toolbox.
* The modelling was done with the [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox.

#### Visualisations (json)
* Metabolic network maps are available upon request. They were constructed for internal use using [Escher](https://escher.github.io/#/) software.


### Instructions:
* `edit_rhtoGEM.m`: run this script to introduce reactions for alternative xylose assimilation pathway.
* `reconstruct_ecRhtoGEM.m`: generate condition-specific enzyme-constrained models with integrated absolute proteomics data. Models are reconstructed based on optimized parameters, including manually curated enzymatic kcat values as provided in `manualModifications.m`. Model reconstruction involves the following steps: 1) `geckomat/enhanceGEM.m` pipeline which creates the first version of ec-model with an enzyme pool used for screening required kcat modifications; 2) `gecomat/utilities/integrate_proteomics/generate_protModels.m` pipeline which sets constraints on high-quality measured individual enzymes from the provided proteomics dataset; 3) adding ribosomal subunits to the ec-models. In `results/generate_protModels_pipeline` is written which enzyme abundances were automatically modified in order to reach the experimental conditions. In `results/ribosome_integration` is written distribution of average ribosomal subunit abundances and also which subunit abundances were flexibilized in order to reach experimental conditions.
* `analyze_ecRhtoGEM.m`: runs random sampling, calculates flux mean value and standard deviation; calculate enzyme usage; extracts ATP, NADPH, and NADH production balances from flux results; converts fluxes to original, non-ecModel version.

##### Last update:
09-15-2022
