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
%   Last modified: 2021-02-20
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
    {'r_1718','r_1654','r_1714','r_2091'},...               % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake
    {'r_1718','r_1654','r_1714','r_2091'}};                 % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake
conditions.exch.value = {[-1.86,0.223,0.367],... % Xexp
    [-0.4345,0.004,0.077],...                   % XNlim
    [0,-6.941,0.127],...                        % Aexp
    [0,-1.9706,-0.033],...                      % ANlim
    [0,0,-2.45,-1000],...                          % GexpUrea
    [0,0,-0.415,-1000]};                        %GNlimUrea
conditions.exch.lbub = {{'lb','ub','ub'},...    % Xexp
    {'lb','ub','ub'},...                        % XNlim
    {'lb','lb','ub'},...                        % Aexp
    {'lb','lb','lb'},...                        % ANlim
    {'lb','lb','lb','lb'},...                   % GexpUrea
    {'lb','lb','lb','lb'}};                     % GNlimUrea
    
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
cd ..
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
    {'r_1718_REV','r_1654_REV','r_1714_REV','r_2091_REV'}};             % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake
conditions.exch.value = {[1.86,0.223,0.367],...                 % Xexp
    [0.4345,0.004,0.077],...                                    % XNlim
    [0,6.941,0.127],...                                         % Aexp
    [0,1.9706,0.033],...                                        % ANlim
    [0,0,2.45,1000,0.06],...                                      % GexpUrea
    [0,0,0.415,1000]};                                          %GNlimUrea
conditions.exch.lbub = {{'ub','ub','ub'},...                    % Xexp
    {'ub','ub','ub'},...                                        % XNlim
    {'ub','ub','ub'},...                                        % Aexp
    {'ub','ub','ub'},...                                        % ANlim
    {'ub','ub','ub','ub','ub'},...                              % GexpUrea
    {'ub','ub','ub','ub'}};                                     % GNlimUrea

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
% root = pwd; % Get the root directory of the folder
% load([root '/models/ecRhtoGEM_batch.mat'])
% ecModel_batch = model;
% load([root '/models/ecRhtoGEM.mat'])
% ecModel = model;
cd([root '/GECKO/geckomat/utilities/integrate_proteomics'])
grouping = [2,2,2,2,2,2];
[~,~,fermParams] = load_Prot_Ferm_Data(grouping);

close all % Close all figures, to facilitate saving of new figures
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

%% III.Add ribosome subunits
% Uncomment and run if the above code (step I and II) was not run in the
% same MATLAB session
%root = pwd; % Get the root directory of the repository
%grouping = [2,2,2,2,2,2];
cd([root, '/GECKO/geckomat/utilities/integrate_proteomics']);
[pIDs,protData,fermParams,byProds] = load_Prot_Ferm_Data(grouping);
protData = cell2mat(protData);

load([root '/models/ecRhtoGEM_Xexp.mat']); ecModels{1} = model;
load([root '/models/ecRhtoGEM_XNlim.mat']); ecModels{2} = model;
load([root '/models/ecRhtoGEM_Aexp.mat']); ecModels{3} = model;
load([root '/models/ecRhtoGEM_ANlim.mat']); ecModels{4} = model;
load([root '/models/ecRhtoGEM_GexpUrea.mat']); ecModels{5} = model;
load([root '/models/ecRhtoGEM_GNlimUrea.mat']); ecModels{6} = model;
% Order in ecModels matches protData and fermParams

%Get indexes for each replicate
clear repl
for i=1:length(grouping)
    try
        repl.first(i)=repl.last(end)+1;
        repl.last(i)=repl.last(end)+grouping(i);
    catch
        repl.first=1;
        repl.last=grouping(1);
    end
end

% Load ribosome information
fid=fopen([root '/data/ribosome.txt']);
data=textscan(fid,'%q %q %q %q %q %q','HeaderLines',1,'Delimiter','\t');
fclose(fid);
ribo=data{1};

%Keep only the subunits that are also in the proteomics data
[~,repRibo.dataIdx]=ismember(ribo,pIDs);
repRibo.protein=ribo(repRibo.dataIdx>0);
repRibo.dataIdx(repRibo.dataIdx==0)=[];
repRibo.avgLevel=mean(protData(repRibo.dataIdx,:),2);

%Density plot to see distribution of subunit abundances
[f,xi]=ksdensity(log10(repRibo.avgLevel),'Bandwidth',0.1);
plot(xi,f);
xlabel('Subunit abundance (log10(mmol/gDCW))');
ylabel('Density');
title('Distribution of average ribosomal subunit abundances');
saveas(gca,[root, '/results/ribosome_integration/average_riboSubunit_abundance.jpg']);

%Only include ribosomal subunite with abundance over 1e-5 mmol/gDCW, as
%these are likely essential subunits, while lower abundances are
%alternative subunits
rmRibo=repRibo.avgLevel<1e-5 | isnan(repRibo.avgLevel);
ribo=repRibo.protein(~rmRibo);

%Prepare subunit information to be added to each model
[~,idx]=ismember(ribo,data{1});
enzNames=data{2}(idx);
enzGenes=data{3}(idx);
MWs=str2double(data{5}(idx))/1000;
sequences=data{6}(idx);
pathways={'Ribosome'};

% Modify reactions
adjusted=cell.empty();
for j=1:numel(ecModels);
    disp(['Add ribosome subunits to condition: ' fermParams.conds{j}])
    model_P=ecModels{j};
    %Add enzymes and genes
    model_P.enzymes=[model_P.enzymes;ribo];
    model_P.enzNames=[model_P.enzNames;enzNames];
    model_P.enzGenes=[model_P.enzGenes;enzGenes];
    genesToAdd.genes=enzGenes;
    genesToAdd.geneShortNames=enzNames;
    model_P=addGenesRaven(model_P,genesToAdd);
    model_P.MWs=[model_P.MWs;MWs];
    model_P.sequences=[model_P.sequences;sequences];
    model_P.pathways=[model_P.pathways;repmat(pathways,numel(ribo),1)];

    %Include new amino acid pseudometabolite
    metsToAdd.mets={'aminoAcids'};
    metsToAdd.metNames={'amino acids for protein'};
    metsToAdd.compartments={'c'};
    model_P=addMets(model_P,metsToAdd);
    
    %Modify protein pseudoreaction to produce amino acid pseudometabolite
    protMetIdx=getIndexes(model_P,'protein','metnames');
    aaMetIdx=getIndexes(model_P,'aminoAcids','mets');
    protRxnIdx=getIndexes(model_P,'r_4047','rxns');
    model_P.S([protMetIdx,aaMetIdx],protRxnIdx)=[0,1];
    
    %Include ribosome subunits as pseudometabolites
    riboToAdd.mets=strcat('prot_',sort(ribo));
    riboToAdd.metNames=riboToAdd.mets;
    riboToAdd.compartments='c';
    model_P=addMets(model_P,riboToAdd);
    
    %Determine stoichiometric coefficient of ribosomes
    mmolAA=full(model_P.S(:,protRxnIdx));
    mmolAA=-sum(mmolAA(mmolAA<0)); %mmol amino acids in biomass
    riboKcat=10.5*3600; %10.5 aa/sec (doi:10.1042/bj1680409) -> aa/hour
    riboKcat=mmolAA/riboKcat; %compensate for the amount of amino acids elongated
    riboKcat=riboKcat*1000; % 1000-fold increase to prevent very low fluxes
    % Include new reaction representing ribosomes (=translation)
    rxnsToAdd.rxns={'translation'};
    rxnsToAdd.mets=[model_P.mets(protMetIdx),model_P.mets(aaMetIdx),riboToAdd.mets'];
    rxnsToAdd.stoichCoeffs=[1,-1,repmat(-riboKcat,1,numel(riboToAdd.mets))];
    rxnsToAdd.subSystem={'Ribosome'};
    rxnsToAdd.grRules={strjoin(enzGenes,' and ')};
    model_P=addRxns(model_P,rxnsToAdd);
    
    %Add exchange reactions or usage reactions
    riboExchId=numel(model_P.rxns)+1;
    model_P=addExchangeRxns(model_P,'in',riboToAdd.mets);
    riboExchId(2)=numel(model_P.rxns);
    %Add UB for enzyme exchange reactions based on measurements.
    %If UB is too low, then adjust to value predicted by model.
    %cd geckomat/utilities/integrate_proteomics
    abundances   = protData(:,repl.first(j):repl.last(j))*1000; % multiplied by 1000 to prevent very low fluxes
    [pIDsCond, filtAbundances] = filter_ProtData(pIDs,abundances,1.96,true);
    %cd ../../../..
    sol=solveLP(model_P);
    for i=riboExchId(1):riboExchId(2)
        protId=regexprep(model_P.rxnNames{i},'prot_(......).*','$1');
        k=find(strcmp(protId,pIDsCond));
        if ~isempty(k)
            model_P.rxns{i}=['prot_' protId '_exchange'];
            model_P.rxnNames{i}=['prot_' protId '_exchange'];
            model_P.concs(strcmp(protId,model_P.enzymes))=filtAbundances(k);
            if sol.x(i)<filtAbundances(k)
                model_P.ub(i)=filtAbundances(k);
            else
                model_P.ub(i)=sol.x(i)*1.01;
                adjusted{end+1,1}=protId;
                adjusted{end,2}=filtAbundances(k);
                adjusted{end,3}=sol.x(i)*1.01;
                adjusted{end,4}=fermParams.conds{j};
                fprintf('%s abundance adjusted. Measured: %e / Adjusted: %e\n',adjusted{end,1:3})
            end
        else
            model_P.ub(i)=Inf;
            model_P.rxns{i}=['draw_prot_' protId];
            model_P.rxnNames{i}=['draw_prot_' protId];
        end
    end
    model_P = rmfield(model_P,'subSystems');
    ecModels{j}=model_P;
    eval(['ecModelP_' fermParams.conds{j} ' = model_P;']);
    exportForGit(model_P,['ecRhtoGEM_' fermParams.conds{j}],[root, '/models'],{'xml','yml','mat'},false,false);
end
adjusted=adjusted';
fid=fopen([root, '/results/ribosome_integration/modifiedRibosomeSubunits.txt'],'w');
fprintf(fid,'protein_IDs previous_values modified_values condition\n');
fprintf(fid,'%s %f %f %s\n',adjusted{:});
fclose(fid);

clear abundances aaMetIdx adjusted ans enzGenes enzNames fid filtAbundances
clear genesToAdd i j k metsToAdd mmolAA model_P MWs pathways pIDsCond protId
clear protMetIdx protRxnIdx repl ribo riboExchId riboKcat riboToAdd 
clear rxnsToAdd sequences sol
cd ([root]);

%% Get fluxes and enzyme usages to each reaction

ecModelsRibo{1} = ecModelP_Xexp;
ecModelsRibo{2} = ecModelP_XNlim;
ecModelsRibo{3} = ecModelP_Aexp;
ecModelsRibo{4} = ecModelP_ANlim;
ecModelsRibo{5} = ecModelP_GexpUrea;
ecModelsRibo{6} = ecModelP_GNlimUrea;
% Order in ecModelsRibo matches fermParams.conds

for j = 1:numel(ecModelsRibo)
    disp(['======= Results from model: ' fermParams.conds{j}, ' ======='])
    modelTmp_proteome = ecModelsRibo{j};    
    printConstraints(modelTmp_proteome,-1000,1000);
    disp(['======= fluxes: ======='])
    sol = solveLP(modelTmp_proteome,1);
    printFluxVector(modelTmp_proteome, sol.x, 'true', 'true');
    disp(['======= enzyme usage: ======='])
    [absUsage,capUsage] = enzymeUsage(modelTmp_proteome,sol.x);
    ecModelsRibo{j} = modelTmp_proteome;
end