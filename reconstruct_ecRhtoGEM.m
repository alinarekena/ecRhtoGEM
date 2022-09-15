%
% reconstruct_ecRhtoGEM
%
%   Reconstruction of enzyme-constrained genome-scale metabolic model for
%   R.toruloides (ecRhtoGEM). Firstly, creates ec-model with protein pool using
%   optimzed parameters (most importantly, manually curated kcat values).
%   Secondly, generates ec-model with individually constrained proteomics
%   measurements for the modelled enzymes. Thirdly, adds ribosomal subunits
%   to the ec-models by adding a translation pseudoreaction.
%
%   Last modified: 2021-05-03
%

% Prepare COBRA and set repo root path
initCobraToolbox
root = pwd; % Get the root directory of the folder

%% Load model
model    = importModel(fullfile(root,'models','rhto_edit.xml'));
modelVer = model.description(strfind(model.description,'_v')+1:end);

%% Define the model conditions and their "parameters"
% Expand on this conditions structure with more model specific information
% that is required to run the script. Define all of these parameters here
% in the beginning.
conditions.abbrev = {'Xexp','XNlim','Aexp','ANlim','GexpUrea','GNlimUrea'};
conditions.exch.rxns = {{'r_1718','r_2104','r_4340'},...    % D-xylose uptake, xylitol production, D-arabinitol production
    {'r_1718','r_2104','r_4340'},...                        % D-xylose uptake, xylitol production, D-arabinitol production
    {'r_1718','r_1634','r_1687'},...                        % block D-xylose uptake, allow acetate uptake, citrate(3-) production
    {'r_1718','r_1634','r_1687'},...                        % block D-xylose uptake, allow acetate uptake, citrate(3-) uptake
    {'r_1718','r_1654','r_1714','r_1808','r_2091'},...               % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake
    {'r_1718','r_1654','r_1714','r_1808','r_2091'}};                 % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake
conditions.exch.value = {[-1.86,0.223,0.367],... % Xexp
    [-0.4345,0.004,0.077],...                   % XNlim
    [0,-6.63,0.122],...                        % Aexp
    [0,-1.9706,-0.043],...                      % ANlim
    [0,0,-2.49,0.049,-1000],...                          % GexpUrea
    [0,0,-0.415,-0.007,-1000]};                        %GNlimUrea
conditions.exch.lbub = {{'lb','ub','ub'},...    % Xexp
    {'lb','ub','ub'},...                        % XNlim
    {'lb','lb','ub'},...                        % Aexp
    {'lb','lb','lb'},...                        % ANlim
    {'lb','lb','lb','ub','lb'},...                   % GexpUrea
    {'lb','lb','lb','lb','lb'}};                     % GNlimUrea
    
for i = 1:numel(conditions.abbrev) % Loop through the conditions
    modelTmp = model; % Work on a temporary model structure, leaving the original untouched for the next condition
    modelTmp = setParam(modelTmp, conditions.exch.lbub{i}, ...
        conditions.exch.rxns{i}, conditions.exch.value{i});
    fprintf(['\n=========== Results from model: ' conditions.abbrev{i}, ' ===========\n\n'])    
    printConstraints(modelTmp,-1000,1000);
    sol = solveLP(modelTmp,1);
    printFluxVector(modelTmp, sol.x, 'true', 'true');
    conditions.model{i} = modelTmp;
end

%% I.enhanceGEM pipeline

%Clone the necessary repos:
%delete GECKO in case that a previous copy exists here
if isfolder('GECKO') 
    rmdir ('GECKO','s')
end
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git('fetch')
% Switch GECKO to the 'notReviewedDevel' branch where relative changes to the
% GECKO code are tracked. This contains various necessary changes in
% e.g. measureAbundance, constrainEnzymes and generate_protModels.
    % files saveECmodel and sumBioMass were from older GECKO versions to avoid some errors.
    % TODO: updateDatabases, addProtein, getEnzymeCodes, convertToEnzymeModel
    % are in GECKO PR #122 and custom functions can be removed once PR is
    % merged.
git('switch notReviewedDevel2')
cd ..

% From here define a new loop that generates the condition-specific batch
% models.
% The folder 'customGECKO' contains identical files for all conditions.
    % files manualModifications, relative_proteomics, uniprot.tab,
    % abs_proteomics, fermentationData, ProtDatabase contain
    % R.toruloides-specific information.

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

% ecModel contains manually curated Kcat values, previously tested for each condition
% Check the solution of each condition by setting experimental data

conditions.abbrev = {'Xexp','XNlim','Aexp','ANlim','GexpUrea','GNlimUrea'};
conditions.exch.rxns = {{'r_1718_REV','r_2104','r_4340'},...    % D-xylose uptake, xylitol production, D-arabinitol production
    {'r_1718_REV','r_2104','r_4340'},...                        % D-xylose uptake, xylitol production, D-arabinitol production
    {'r_1718_REV','r_1634_REV','r_1687'},...                    % block D-xylose uptake, allow acetate uptake, citrate(3-) production
    {'r_1718_REV','r_1634_REV','r_1687_REV'},...                % block D-xylose uptake, allow acetate uptake, citrate(3-) uptake
    {'r_1718_REV','r_1654_REV','r_1714_REV','r_2091_REV','r_1808'},...  % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake, glycerol production
    {'r_1718_REV','r_1654_REV','r_1714_REV','r_2091_REV','r_1808_REV'}};             % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake, glycerol uptake
conditions.exch.value = {[1.86,0.223,0.367],...                 % Xexp
    [0.4345,0.004,0.077],...                                    % XNlim
    [0,6.63,0.122],...                                         % Aexp
    [0,1.9706,0.043],...                                        % ANlim
    [0,0,2.49,1000,0.049],...                                      % GexpUrea
    [0,0,0.415,1000,0.007]};                                          %GNlimUrea
conditions.exch.lbub = {{'ub','ub','ub'},...                    % Xexp
    {'ub','ub','ub'},...                                        % XNlim
    {'ub','ub','ub'},...                                        % Aexp
    {'ub','ub','ub'},...                                        % ANlim
    {'ub','ub','ub','ub','ub'},...                              % GexpUrea
    {'ub','ub','ub','ub','ub'}};                                     % GNlimUrea

cd utilities/integrate_proteomics
[pIDs,protData,fermParams,byProds] = load_Prot_Ferm_Data([2,2,2,2,2,2]);
cd ../../limit_proteins
% Quickly average the replicates, no filtering, just to be able to
% calculate conditon specific f-values.
for i=1:numel(fermParams.conds) 
    tmp = [protData{2*i-1} protData{2*i}];
    avgProtData(:,i) = mean(tmp,2);
    f(i) = measureAbundance(ecModel_batch.enzymes,pIDs,avgProtData(:,i));
end

for i = 1:numel(conditions.abbrev) % Loop through the conditions
    modelTmp_batch = ecModel_batch; % Work on a temporary model structure, leaving the original untouched for the next condition
    modelTmp_batch.ub(getIndexes(modelTmp_batch,'prot_pool_exchange','rxns')) = fermParams.Ptot(i)*f(i)*params.sigma;
    modelTmp_batch = setParam(modelTmp_batch, conditions.exch.lbub{i}, ...
        conditions.exch.rxns{i}, conditions.exch.value{i});
    fprintf(['\n=========== Results from model: ' conditions.abbrev{i}, ' ===========\n\n'])    
    printConstraints(modelTmp_batch,-1000,1000);
    sol = solveLP(modelTmp_batch,1);
    printFluxVector(modelTmp_batch, sol.x, 'true', 'true');
    topUsedEnzymes(sol.x,modelTmp_batch,{''},{''},false);
    conditions.ecModel_batch{i} = modelTmp_batch;
end

% If the solution shows that the ecRhtoGEM model can reach the experimental
% parameters and the majority of topUsedEnzymes are less than 10% of
% protein pool, proceed with the generate_protModels pipeline for the
% proteomics integration.

% sigma factor for each condition has been fitted and saved 'results/enhanceGEM_pipeline'

clear avgProtData tmp

%% II.generate_protModels pipeline
% Uncomment and run if the above code was not run in the same session
root = pwd; % Get the root directory of the folder
load([root '/models/ecRhtoGEM_batch.mat'])
ecModel_batch = model;
load([root '/models/ecRhtoGEM.mat'])
ecModel = model;
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
% TODO: save auto-flexibilized enzymes, now only shown in command window

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

% Models are stored in /models, not in the GECKO folder (as the GECKO
% folder would be deleted if you rerun this script)

clear fileNames figHandles figCond B i

% Note: Models are saved with prot_pool_exchange as an objective function,
% constrained with GR (lb), carbon uptake (b, flexibilized), byproducts
% (ub, flexibilized). If 'chemostat constraints too stringent' models are saved with
% carbon uptake as objective

cd ([root]);

%Remove the cloned repos:
rmdir('GECKO', 's')
