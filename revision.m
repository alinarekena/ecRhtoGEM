% revision_ecRhtoGEM.m

% This script contains workflow for:
%
% I. phosphoketolase pathway;
% I.a acetate kinase (ACK) - S5 Dataset;
% I.b blocked phosphoketolase (XPK) - S6 Dataset;
% II. new xylose uptake pathway;
% II.a NAD-dependent D-arabinitol 2-dehydrogenase - S7 Dataset.

% 2023-02-08
% Alina Rekena

% prepare software:
initCobraToolbox

% set path:
cd ../../Github/ecRhtoGEM/
root = pwd; % Get the root directory of the folder

%% I. Phosphoketolase pathway:

% Comment:
% the authors need to provide evidence from their own data or from the 
% literature on phosphoketolase pathway activity. At least, the authors 
% should comment and follow up on: [..] (ii) FBA simulation without 
% phosphoketolase being active to confirm no alternative pathways exist 
% to functionally replace it [..]

% Response:
% [..] that we have evaluated the possibility of acetate kinase (ACK) 
% being present in the pathway, as reported in R. toruloides models for 
% strain IFO 0880 (Dinh et al. 2019, Kim et al. 2021), while 
% phosphotransacetylase (PTA) reaction was included for strain NP11 in 
% rhtoGEM (Tiukova et al. 2019). 

% "As the results were highly similar, further flux analysis was carried 
% out based on a metabolic route where PTA is active.”

% “The pyruvate decarboxylase and ACL, which exist as alternative pathways
% for producing cytosolic acetyl-CoA during lipid accumulation, were 
% activated only when we blocked the XPK pathway (S6 Dataset).”

%% I.a acetate kinase (ACK) - S5 Dataset:

% Load model:
model    = importModel(fullfile(root,'models','rhto_edit.xml'));
modelVer = model.description(strfind(model.description,'_v')+1:end);

% Remove reactions:
model=removeReactions(model,'t_0082');

%% Introduce reactions:

t_0886 = 'ATP[c] + acetate[c] <=> ADP[c] + acetyl-phosphate[c]';% Acetate kinase
% ACKr (13383) from Rt_IFO0880 and iRhtoC

rxnsToAdd.equations = {t_0886}; 

% Define reaction names
t_0886 = 'Acetate kinase';
rxnsToAdd.rxnNames = {t_0886};
t_0886 = 't_0886';
rxnsToAdd.rxns = {t_0886};

% Define objective and bounds
rxnsToAdd.c  = [0];
rxnsToAdd.lb = [-1000];
rxnsToAdd.ub = [1000];

% Define EC numbers
t_0886 = '2.7.2.1';
rxnsToAdd.eccodes = {t_0886};

% Add genes
genesToAdd.genes          = {'RHTO_04461'};
genesToAdd.geneShortNames = {'RHTO_04461'};
rxnsToAdd.grRules         = {'RHTO_04461'};

%% Introduce changes to the model:

model = addGenesRaven(model,genesToAdd);
model = addRxns(model,rxnsToAdd,3);% add reactions in last position to avoid allowNewMets error

%% enhanceGEM pipeline:

%Clone the necessary repos:
%delete GECKO in case that a previous copy exists here
if isfolder('GECKO') 
    rmdir ('GECKO','s')
end
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git('fetch')
%
git('switch notReviewedDevel2')
cd ..

% Replace custom GECKO scripts
fileNames = struct2cell(dir('customGECKO'));
fileNames = fileNames(1,:);
fileNames(startsWith(fileNames,'.')) = [];
for i = 1:length(fileNames)
    GECKO_path = dir(['GECKO/**/' fileNames{i}]);
    copyfile(['customGECKO' filesep fileNames{i}],GECKO_path.folder)
    disp(['Replaced ' fileNames{i} ' at ' GECKO_path.folder '\'])
end
movefile('relative_proteomics.txt','GECKO/databases','f');
movefile('abs_proteomics.txt','GECKO/databases','f');
movefile('fermentationData.txt','GECKO/databases','f');
clear fileNames GECKO_path i

% Start GECKO/enhanceGEM pipeline
cd GECKO
delete databases/prot_abundance.txt;  % if not deleted,f factor will not be calculated;
GECKOver = git('describe --tags');
cd geckomat/get_enzyme_data
updateDatabases;

% Change lipid and protein conteint to Xexp condition
cd([root '/code'])
% First increase FA fraction by 25%, to prevent later problems in GexpUrea.
% scaleLipidProtein function compensates for the FA increase.
rxnIdx = find(contains(model.rxnNames,'lipid backbone pseudoreaction'));
metIdx = find(contains(model.metNames,'fatty acid backbone'));
model.S(metIdx,rxnIdx)=model.S(metIdx,rxnIdx)*1.25;
lipidData = loadLipidChainData(model,1);
model = scaleLipidProtein(model,lipidData,0.4385);
cd([root '/GECKO/geckomat/'])
[ecModel,ecModel_batch] = enhanceGEM(model,'RAVEN','ecRhtoGEM',modelVer);

% Ignore the sigma fitting, manually set sigma to 0.35, as specified in getModelParameters;
params = getModelParameters();
cd limit_proteins
f = measureAbundance(ecModel_batch.enzymes); % calculates f from average abundances of all conditions
ecModel_batch.ub(getIndexes(ecModel_batch,'prot_pool_exchange','rxns')) = params.Ptot*f*params.sigma;

% Overwrite the files exported by enhanceGEM, now with the new pool UB
cd ../../models
ecModel_batch = saveECmodel(ecModel_batch,'RAVEN','ecRhtoGEM_batch',modelVer);% RAVEN v2.4.1+
cd ecRhtoGEM
movefile('*','../../../models/')
cd ../../geckomat

%% generate_protModels pipeline:

cd([root '/GECKO/geckomat/utilities/integrate_proteomics'])
grouping = [2,2,2,2,2,2];
[~,~,fermParams] = load_Prot_Ferm_Data(grouping);

close all % Close all figures, to facilitate saving of new figures
% Add 
rxnsToAdd.rxns={'r_4062','r_4064'};
rxnsToAdd.equations={'s_3746 =>','s_3747 =>'};
rxnsToAdd.lb=[0,0]; rxnsToAdd.ub=[0,0]; %Blocked for now, should only be opened when scaling lipids
ecModel=addRxns(ecModel,rxnsToAdd,1);
ecModel_batch=addRxns(ecModel_batch,rxnsToAdd,1);
clear rxnsToAdd

generate_protModels(ecModel,grouping,'ecRhtoGEM',ecModel_batch);

% Save figure windows to files
figHandles = get(groot, 'Children'); 
[~,B] = sort([figHandles.Number]);
figHandles = figHandles(B);
for f = 1:numel(figHandles)
    if rem(f,2)
        figCond = fermParams.conds{round(f/2)};
        saveas(figHandles(f),[root, '/results/generate_protModels_pipeline/usage_', figCond, '.jpg']);
    else
        saveas(figHandles(f),[root, '/results/generate_protModels_pipeline/abundance_', figCond, '.jpg']);
    end
end

% Move the files outside of GECKO folder
cd([root '/GECKO/models/prot_constrained/ecRhtoGEM/'])
delete('ecRhtoGEM*.txt') % yml file more complete
delete('dependencies.txt')
movefile('*.txt',[root, '/results/generate_protModels_pipeline'])
movefile('*',[root '/models/'])

clear fileNames figHandles figCond B i

% Note: Models are saved with prot_pool_exchange as an objective function,
% constrained with GR (lb), carbon uptake (b, flexibilized), byproducts
% (ub, flexibilized). If 'chemostat constraints too stringent' models are saved with
% carbon uptake as objective

cd ([root]);

%% analyze_ecRhtoGEM_acetateKinase:
% prepare for random sampling:
%
%
% set measured byproduct constraints: (can be found at fermParams.byP_flux)
ecModelP_Xexp = setParam(ecModelP_Xexp,'eq','r_2104',0.223);    %xylitol exchange
ecModelP_Xexp = setParam(ecModelP_Xexp,'eq','r_4340',0.367);    %D-arabinitol

ecModelP_XNlim = setParam(ecModelP_XNlim,'eq','r_2104',0.004);    %xylitol exchange
ecModelP_XNlim = setParam(ecModelP_XNlim,'eq','r_4340',0.077);    %D-arabinitol

ecModelP_Aexp = setParam(ecModelP_Aexp,'eq','r_1687',0.122);    %citrate

%ANlim: citrate uptake 0.043 is constrained by UB and predicted by FBA;

ecModelP_GexpUrea = setParam(ecModelP_GexpUrea,'eq','r_1808',0.049);%glycerol

ecModelP_GNlimUrea = setParam(ecModelP_GNlimUrea,'eq','r_1808_REV',0.007);%glycerol uptake

% run FBA: with byProducts, optimize for prot_exchange
solveLP(ecModelP_GNlimUrea,1);%same for all models

% create temporary model and constrain it for random sampling:
ecModelP_XexpTmp = setParam(ecModelP_Xexp,'ub','r_2111',0.05399); %UB=flux+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'ub','r_4046',3.6865); %UB=min+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1672',2.957,10);   %CO2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1992_REV',2.288,10); %O2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','prot_pool_exchange',0.7413,10); % FBA flux with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_2104',0.223,10);  % measured xylitol
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_4340',0.367,10);  % measured D-arabinitol

ecModelP_XNlimTmp2 = setParam(ecModelP_XNlim,'lb','r_2111',0.01473); % different from 0.01782 w/o byproducts 
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_2111',0.01488); %UB=0.01473+1%
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_4046',2.8785); %UB=min+1%;
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1672',1.097,10);   %CO2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1992_REV',0.9532,10); %O2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','prot_pool_exchange',0.4792,10);
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_2104',0.004,10);  % measured xylitol
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_4340',0.077,10);  % measured D-arabinitol

ecModelP_AexpTmp = setParam(ecModelP_Aexp,'ub','r_2111',0.07299);%UB=flux+1%
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'ub','r_4046',3.636);%UB=min+1%
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1672',6.735,10);   %CO2
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1992_REV',6.667,10); %O2
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','prot_pool_exchange',0.9534,10);%flux
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1687',0.122,10);% measured citrate

ecModelP_ANlimTmp = setParam(ecModelP_ANlim,'ub','r_2111',0.01199);%UB=flux+1%
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'ub','r_4046',3.484);%UB=min+1%
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1672',1.39,10);   %CO2
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1992_REV',2.012,10); %O2
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','prot_pool_exchange',0.3002,10);%flux
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1687_REV',0.043,20);% citrate uptake (as in GECKO pipeline)

ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUrea,'ub','r_2111',0.17998);%0.1782+1%
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'ub','r_4046',1.814);% max NGAM GexpUrea
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1672',6.94,10)   %CO2
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1992_REV',5.812,10) %O2
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','prot_pool_exchange',1.41,10)
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1808',0.049,10)%glycerol

ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUrea,'ub','r_2111',0.02099);%0.02079+1%
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'ub','r_4046',2.7775);%UB=min+1%
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1672',1.444,10)   %CO2
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1992_REV',1.16,10) %O2
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','prot_pool_exchange',1.179,10)
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1808_REV',0.007,10)%glycerol uptake

% run random sampling:
%
% Note, goodRxns don't make sense because of different C sources.
%
goodRxns = [];
rsmatrix=rs(ecModelP_XexpTmp,2000,true,true,true,goodRxns,true);%MinFlux=True
rsmatrix=rs(ecModelP_XNlimTmp2,2000,true,true,true,goodRxns,true);
rsmatrix=rs(ecModelP_AexpTmp,2000,true,true,true,goodRxns,true);
rsmatrix=rs(ecModelP_ANlimTmp,2000,true,true,true,goodRxns,true);
rsmatrix=rs(ecModelP_GexpUreaTmp,2000,true,true,true,goodRxns,true);
rsmatrix=rs(ecModelP_GNlimUreaTmp,2000,true,true,true,goodRxns,true);

% output: ec
clear out
sol.median=full(median(rsmatrix,2));  %2 means for each row
sol.mean=full(mean(rsmatrix,2));
sol.std=full(std(rsmatrix,0,2));
clear out

%% convert fluxes to original, non-ecModel version:

% rs matrix: non-ec
rsXmapped = mapRxnsToOriginal(ecModelP,model,rsmatrix)
% (repeat for each condition)

% output: non-ec
solXmapped.median=full(median(rsXmapped,2));  %2 means for each row
solXmapped.mean=full(mean(rsXmapped,2));
solXmapped.std=full(std(rsXmapped,0,2));
% (repeat for each condition)

%% calculate ATP, NADPH, and NADH balances from non-ec fluxes:

% adjust for each compound and repeat for each condition
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,solXmapped.median,true);
    clear out
    out.rxns    = model.rxns(rxnIdx);
    out.rxnNames= model.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model,rxnIdx);
    out.fluxes  = num2cell(fluxes);
    out = [out.rxns out.rxnNames out.rxnEqns out.fluxes];
    end

%% I.b Phosphoketolase (XPK) - S6 Dataset:
% Gexp:
% new model - PK knockout:
ecModelP_Gexp = setParam(ecModelP_Gexp,'eq','t_0081No1',0)    % phosphoketolase

% GNlimUrea:
% new model - PK knockout:
ecModelP_GNlim = setParam(ecModelP_GNlim,'eq','t_0081No1',0)    % phosphoketolase

% Xexp:
% new model - PK knockout:
ecModelP_Xexp = setParam(ecModelP_Xexp,'eq','t_0081No1',0)    % phosphoketolase

% XNlim:
% new model - PK knockout:
ecModelP_XNlim = setParam(ecModelP_XNlim,'eq','t_0081No1',0)    % phosphoketolase

%% analyze_ecRhtoGEM_phosphoketolase:

% prepare for random sampling:
%
%

% run FBA:
solveLP(ecModelP,1);%repeat for all models

% create temporary model and constrain it for random sampling:
ecModelP_XexpTmp = setParam(ecModelP_Xexp,'ub','r_2111',0.05399); %UB=flux+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'ub','r_4046',3.6865); %UB=min+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1672',2.932,10);   %CO2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1992_REV',2.275,10); %O2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','prot_pool_exchange',1.692,10); % FBA flux with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_2104',0.223,10);  % measured xylitol
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_4340',0.367,10);  % measured D-arabinitol

ecModelP_XNlimTmp = setParam(ecModelP_XNlim,'lb','r_2111',0.01473); % different from 0.01782 w/o byproducts 
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp,'ub','r_2111',0.01488); %UB=0.01473+1%
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_4046',2.929); %UB=min+1%;
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1672',1.097,10);   %CO2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1992_REV',0.9543,10); %O2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','prot_pool_exchange',0.5338,10);
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_2104',0.004,10);  % measured xylitol
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_4340',0.077,10);  % measured D-arabinitol

ecModelP_AexpTmp = setParam(ecModelP_Aexp,'ub','r_2111',0.07299);%UB=flux+1%
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'ub','r_4046',3.686);%UB=min+1%
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1672',6.726,10);   %CO2
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1992_REV',6.668,10); %O2
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','prot_pool_exchange',3.209,10);%flux
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1687',0.122,10);% measured citrate

ecModelP_ANlimTmp = setParam(ecModelP_ANlim,'ub','r_2111',0.01199);%UB=flux+1%
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'ub','r_4046',3.484);%UB=min+1%
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1672',1.402,10);   %CO2
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1992_REV',2.011,10); %O2
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','prot_pool_exchange',0.5217,10);%flux
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1687_REV',0.043,20);% citrate uptake (as in GECKO pipeline)

ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUrea,'ub','r_2111',0.19099);%0.1891+1%
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'ub','r_4046',1.814);% max NGAM GexpUrea
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1672',6.943,10)   %CO2
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1992_REV',5.82,10) %O2
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','prot_pool_exchange',3.535,10)
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1808',0.049,10)%glycerol

ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUrea,'ub','r_2111',0.02099);%0.02079+1%
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'ub','r_4046',2.929);%UB=min+1%
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1672',1.445,10)   %CO2
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1992_REV',1.164,10) %O2
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','prot_pool_exchange',2.205,10)
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1808_REV',0.007,10)%glycerol uptake

ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'eq','t_0081No1',0) 
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'eq','t_0081No1',0)
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'eq','t_0081No1',0) 
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'eq','t_0081No1',0)
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'eq','t_0081No1',0) 
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'eq','t_0081No1',0)

% run random sampling:
%
% Note, goodRxns don't make sense because of different C sources.
%
goodRxns = [];
rsmatrix_Xexp=rs(ecModelP_XexpTmp,2000,true,true,true,goodRxns,true);%MinFlux=True
rsmatrix_XNlim=rs(ecModelP_XNlimTmp2,2000,true,true,true,goodRxns,true);
rsmatrix_Aexp=rs(ecModelP_AexpTmp,2000,true,true,true,goodRxns,true);
rsmatrix_ANlim=rs(ecModelP_ANlimTmp,2000,true,true,true,goodRxns,true);
rsmatrix_Gexp=rs(ecModelP_GexpUreaTmp,2000,true,true,true,goodRxns,true);
rsmatrix_GNlim=rs(ecModelP_GNlimUreaTmp,2000,true,true,true,goodRxns,true);

% output: ec
clear sol
sol.median=full(median(rsmatrix,2));  %2 means for each row
sol.mean=full(mean(rsmatrix,2));
sol.std=full(std(rsmatrix,0,2));
clear sol

%% convert fluxes to original, non-ecModel version:

% rs matrix: non-ec
rsXmapped = mapRxnsToOriginal(ecModelP,model,rsmatrix)
% (repeat for each condition)

% output: non-ec
solXmapped.median=full(median(rsXmapped,2));  %2 means for each row
solXmapped.mean=full(mean(rsXmapped,2));
solXmapped.std=full(std(rsXmapped,0,2));
% (repeat for each condition)

%% calculate ATP, NADPH, and NADH balances from non-ec fluxes:

% adjust for each compound and repeat for each condition
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,solXmapped.median,true);
    clear out
    out.rxns    = model.rxns(rxnIdx);
    out.rxnNames= model.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model,rxnIdx);
    out.fluxes  = num2cell(fluxes);
    out = [out.rxns out.rxnNames out.rxnEqns out.fluxes];
    end

%% II.a DAD-2 with NAD/NADH dependency test - S7 Dataset:

% Comment:
% the authors mentioned: “In support of this mechanism, the metabolic 
% model could only predict a feasible solution of LXR/DAD-2 only using 
% NADP(+) as the cofactor.” Was this simulation performed with enzyme 
% constraint (with kcat or kapp < kcat being active)? The authors need to 
% also perform simulation with FBA (i.e., without fluxes constrained by 
% enzyme levels) to assess if infeasibility is due to metabolic demand or 
% due to enzyme level constraining flux.

% Response:
% "The inability to grow on NADH was a problem with our analysis pipeline, 
% which we have now resolved (revision.m). Now, flux distributions with 
% NADH or NADPH are highly similar, except for NADPH regeneration."

% switch to master branch (without acetate kinase):
% load models:
load([root '/models/ecRhtoGEM_Xexp.mat'])
ecModelP_Xexp_masterBranch = model;
load([root '/models/ecRhtoGEM_XNlim.mat'])
ecModelP_XNlim_masterBranch = model;

model    = importModel(fullfile(root,'models','rhto.xml'));
model = changeRxnBounds(model, {'r_1714'},0,'l');      %D-glucose exchange
model = changeRxnBounds(model, {'r_1718'},-1.86,'l');  %D-xylose exchange

%check kcats:
[kcat,rxnName,rxn,MW] = getKcat(ecModelP_Xexp,'M7X791') %109 s-1 both fwd and rev;

% Remove reactions:
ecModelP_Xexp=removeReactions(ecModelP_Xexp_masterBranch,'t_0884_REVNo1');
ecModelP_Xexp=removeReactions(ecModelP_Xexp,'t_0884No1');

ecModelP_XNlim=removeReactions(ecModelP_XNlim_masterBranch,'t_0884_REVNo1');
ecModelP_XNlim=removeReactions(ecModelP_XNlim,'t_0884No1');

%% Introduce reactions:

t_0884No1 = 'H+[c] + NADH[c] + D-ribulose[c] + 0.0025484 prot_M7X791[c] => NAD[c] + D-arabinitol[c]';
t_0884_REVNo1 = 'NAD[c] + D-arabinitol[c] + 0.0025484 prot_M7X791[c] => H+[c] + NADH[c] + D-ribulose[c]';
rxnsToAdd.equations = {t_0884No1;t_0884_REVNo1}; 

t_0884No1 = 'D-arabinitol 2-dehydrogenase/D-ribulose reductase (No1)';
t_0884_REVNo1 = 'D-arabinitol 2-dehydrogenase/D-ribulose reductase (reversible) (No1)';
rxnsToAdd.rxnNames = {t_0884No1;t_0884_REVNo1};

t_0884No1 = 't_0884No1';
t_0884_REVNo1 = 't_0884_REVNo1';
rxnsToAdd.rxns = {t_0884No1;t_0884_REVNo1};

rxnsToAdd.c  = [0 0];
rxnsToAdd.lb = [0 0];
rxnsToAdd.ub = [1000 1000];

t_0884No1 = '1.1.1.10;1.1.1.138';
t_0884_REVNo1 = '1.1.1.10;1.1.1.138';

rxnsToAdd.eccodes = {t_0884No1;t_0884_REVNo1};

rxnsToAdd.grRules = {'RHTO_00373' 'RHTO_00373'};

% Add reactions:
ecModelP_Xexp = addRxns(ecModelP_Xexp,rxnsToAdd,3);
ecModelP_XNlim = addRxns(ecModelP_XNlim,rxnsToAdd,3);

% Check reactions:
printModel(ecModelP_Xexp,'t_0884No1');
printModel(ecModelP_Xexp,'t_0884_REVNo1');
printModel(ecModelP_XNlim,'t_0884No1');
printModel(ecModelP_XNlim,'t_0884_REVNo1');

clear kcat MW rxn rxnName
clear rxnsToAdd t_0884_REVNo1 t_0884No1

%% analyze_ecRhtoGEM:
% prepare for random sampling:
%
% set measured byproduct constraints: (can be found at fermParams.byP_flux)
ecModelP_Xexp = setParam(ecModelP_Xexp,'eq','r_2104',0.223);    %xylitol exchange
ecModelP_Xexp = setParam(ecModelP_Xexp,'eq','r_4340',0.367);    %D-arabinitol

ecModelP_XNlim = setParam(ecModelP_XNlim,'eq','r_2104',0.004);    %xylitol exchange
ecModelP_XNlim = setParam(ecModelP_XNlim,'eq','r_4340',0.077);    %D-arabinitol

% run FBA to get bounds: with byProducts, optimize for prot_exchange
checkObjective(ecModelP_XNlim);
%ecModelP_XexpTmp=changeObjective(ecModelP_XexpTmp,'r_2111');
solveLP(ecModelP_XNlimTmp2,1); %same for all models
printFluxVector(ecModelP_XNlimTmp2, ans.x, 'true', 'true');

% create temporary model and constrain it for random sampling:
ecModelP_XexpTmp = setParam(ecModelP_Xexp,'ub','r_2111',0.05399); %UB=flux+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'ub','r_4046',3.6865); %UB=min+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1672',2.961,10);   %CO2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1992_REV',2.296,10); %O2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','prot_pool_exchange',1.692,10); % FBA flux with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_2104',0.223,10);  % measured xylitol
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_4340',0.367,10);  % measured D-arabinitol

ecModelP_XNlimTmp2 = setParam(ecModelP_XNlim,'lb','r_2111',0.01473); % different from 0.01782 w/o byproducts 
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_2111',0.01488); %UB=0.01473+1%
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_4046',2.9); %UB=min+1%;
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1672',1.104,10);   %CO2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1992_REV',0.9615,10); %O2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','prot_pool_exchange',0.5338,10);
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_2104',0.004,10);  % measured xylitol
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_4340',0.077,10);  % measured D-arabinitol

% run random sampling:
%
% Note, goodRxns don't make sense because of different C sources.
%
goodRxns = [];
rsmatrix=rs(ecModelP_XexpTmp,2000,true,true,true,goodRxns,true);%MinFlux=True
rsmatrix=rs(ecModelP_XNlimTmp2,2000,true,true,true,goodRxns,true);

% output: ec
clear sol
sol.median=full(median(rsmatrix,2));  %2 means for each row
sol.mean=full(mean(rsmatrix,2));
sol.std=full(std(rsmatrix,0,2));
clear sol

%% convert fluxes to original, non-ecModel version:

% rs matrix: non-ec
rsXmapped = mapRxnsToOriginal(ecModelP,model,rsmatrix)
% (repeat for each condition)

% output: non-ec
solXmapped.median=full(median(rsXmapped,2));  %2 means for each row
solXmapped.mean=full(mean(rsXmapped,2));
solXmapped.std=full(std(rsXmapped,0,2));
% (repeat for each condition)

%% calculate ATP, NADPH, and NADH balances from non-ec fluxes:

% adjust for each compound and repeat for each condition
for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,solXmapped.median,true);
    clear out
    out.rxns    = model.rxns(rxnIdx);
    out.rxnNames= model.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model,rxnIdx);
    out.fluxes  = num2cell(fluxes);
    out = [out.rxns out.rxnNames out.rxnEqns out.fluxes];
    end

