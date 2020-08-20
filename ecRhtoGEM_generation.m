%
% ecRhtoGEM_generation
%
%   Generation of enzyme-constrained R.toruloides genome-scale metabolic
%   model (ecRhtoGEM). Firstly, enzyme-constrained model for growth on
%   xylose at 0.064 [h-1] with total protein pool is generated. Then
%   secondly, absolute proteomics data are incorporated, allowing the
%   modelled but not individually measured enzymes be drawn from the
%   enzymatic protein pool. Thus, generation of ecRhtoGEM employs two GECKO
%   pipelines: enhanceGEM and generate_protModels
%
%   Alina Rekena.   Last modified: 2020-08-20
%

%% Load model:
model    = load('rhto_edit_v1_P1.mat');
model    = model.model;
modelVer = model.description(strfind(model.description,'_v')+1:end);

%% Checklist in the GECKO folder for enhanceGEM pipeline:
%
%   /databases: uniprot.tab, relative_proteomics.txt, chemostatData.tsv;
%   /geckomat/getModelParameters.m: sigma, Ptot, gR_exp, c_source,
%   exch_name{2};
%   /geckomat/enhanceGEM.m: disp(['Sigma factor (not fitted)...
%   /geckomat/limit_proteins/getConstrainedModel.m: OptSigma=sigma;
%   ecModel_batch=...currentEnzymeUB;
%   /geckomat/kcat_sensitivity_analysis/changeMedia_batch.m: %block xylose
%   and oxygen production;
%   /geckomat/change_model/manualModifications.m: %growth limitng Kcats
%   section;

%% Run GECKO/enhanceGEM pipeline:
cd GECKO
GECKOver = git('describe --tags');
cd geckomat/get_enzyme_data
updateDatabases;
cd ..
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA','ecYeastGEM',modelVer);
cd ../..

%% ecModel_batch enhanceGEM pipeline:
% Print all constraints of the model:
printConstraints(ecModel_batch,-1000,1000);
printConstraints(ecModel,-1000,1000);
tempModel = setParam(ecModel_batch,'ub','r_1718_REV',1.74)  %xylose uptake
tempModel = setParam(tempModel,'eq','r_2104',0.228)         %xylitol exchange
tempModel = setParam(tempModel,'eq','r_4340',0.372)         %D-arabinitol
printConstraints(tempModel,-1000,1000);

%% Explore model solution:
% Print exchange reactions using RAVEN:
solveLP(tempModel,1)
printFluxVector(tempModel, ans.x, 'true', 'true');
topUsedEnzymes(ans.x,tempModel,{''},{''},false)

%    If the solution shows that the ecRhtoGEM model can reach the experimental
% parameters and the majority of topUsedEnzymes are less than 10% of
% protein pool, proceed with the generate_protModels pipeline for the
% proteomics integration.

%% Checklist in the GECKO folder for generate_protModels pipeline:
%
%   /databases: abs_proteomics.txt, fermentationData.txt,
%   chemostatData.tsv and relative_proteomics.txt (for the protein pool);
%   /geckomat/getModelParameters.m: sigma=1, Ptot, gR_exp, c_source,
%   exch_names{2};
%   /geckomat/kcat_sensitivity_analysis/changeMedia_batch.m: %block glucose
%   and oxygen production;
%   

%% Incorporate proteomics
grouping   = [3]; %Our dataset contains three replicates per condition
flexFactor = 1.05;  %Allowable flexibilization factor for fixing carbon uptake rate
%Use GECKo utilities for proteomics integration
cd geckomat/utilities/integrate_proteomics
generate_protModels(ecModel,grouping,'ecYeastGEM',ecModel_batch);

ecModelP    = load('ecYeastGEM_P1.mat');
ecModelP    = ecModelP.ecModelP;
checkObjective(ecModelP);
ecModelP=changeObjective(ecModelP,'r_2111');
printConstraints(ecModelP,-1000,1000);

tempModel = setParam(ecModelP,'ub','r_1718_REV',1.74)  %xylose uptake
tempModel = setParam(tempModel,'lb','r_1718_REV',0)
tempModel = setParam(tempModel,'lb','r_2111',0)        %growth
tempModel = setParam(tempModel,'lb','r_4041',0)        %biomass pseudoreaction
tempModel = setParam(tempModel,'eq','r_2104',0.228)    %xylitol exchange
tempModel = setParam(tempModel,'eq','r_4340',0.372)    %D-arabinitol

printConstraints(tempModel,-1000,1000);

% Explore model solutions:
solveLP(tempModel,1)
printFluxVector(tempModel, ans.x, 'true', 'true');
topUsedEnzymes(ans.x,tempModel,{''},{''},false)
[capUsage,absUsage] = enzymeUsage(tempModel,ans.x)

% Save
exportToExcelFormat(tempModel,'ecModelP_XP1_v2_C_tempModel.xlsx')
% and save matlab version of tempModel