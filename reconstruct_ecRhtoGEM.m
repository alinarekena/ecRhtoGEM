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
%   Last modified: 2021-02-18
%

% Prepare COBRA and set repo root path
initCobraToolbox
root = pwd; % Get the root directory of the folder

%% Load model
model    = importModel(fullfile(root,'models','model_edit.xml'));
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
conditions.exch.value = {[-1.74,0.228,0.372],... % Xexp
    [-0.4345,0.004,0.077],...                   % XNlim
    [0,-6.941,0.127],...                        % Aexp
    [0,-1.9706,-0.033],...                      % ANlim
    [0,0,-3,-1000],...                          % GexpUrea
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
git('switch notReviewedDevel')
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
movefile('relative_proteomics.txt','GECKO/Databases','f');

% Start GECKO/enhanceGEM pipeline
cd GECKO
delete databases/prot_abundance.txt;  % if not deleted,f factor will not be calculated;
GECKOver = git('describe --tags');
cd geckomat/get_enzyme_data
updateDatabases;
cd ..
[ecModel,ecModel_batch] = enhanceGEM(model,'RAVEN','ecRhtoGEM',modelVer);%error fitGAM line 54 column 10

% Ignore the sigma fitting, manually set sigma to 1;
params = getModelParameters();
cd limit_proteins
f = measureAbundance(ecModel_batch.enzymes); % calculates f from average abundances of all conditions
ecModel_batch = updateProtPool(ecModel_batch,params.Ptot,f*params.sigma); % f should be multiplied with params.sigma

% Overwrite the files exported by enhanceGEM, now with the new pool UB
cd ../../models
ecModel_batch = saveECmodel(ecModel_batch,'RAVEN','ecRhtoGEM_batch',modelVer);%error using exportForGit: too many input arguments, saveECmodel line60
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
conditions.exch.value = {[1.74,0.228,0.372],...                 % Xexp
    [0.4345,0.004,0.077],...                                    % XNlim
    [0,6.941,0.127],...                                         % Aexp
    [0,1.9706,0.033],...                                        % ANlim
    [0,0,3,1000,0.071],...                                      % GexpUrea
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
    modelTmp_batch = updateProtPool(modelTmp_batch,fermParams.Ptot(i),f(i)*0.5); % Sigma = 0.5
    modelTmp_batch = setParam(modelTmp_batch, conditions.exch.lbub{i}, ...
        conditions.exch.rxns{i}, conditions.exch.value{i});
    fprintf(['\n=========== Results from model: ' conditions.abbrev{i}, ' ===========\n\n'])    
    printConstraints(modelTmp_batch,-1000,1000);
    sol = solveLP(modelTmp_batch,1);
    printFluxVector(modelTmp_batch, sol.x, 'true', 'true');
    topUsedEnzymes(sol.x,modelTmp_batch,{''},{''},false); % currently shows in 'ans' variable, can it be written to command window?
    conditions.ecModel_batch{i} = modelTmp_batch;
end

% If the solution shows that the ecRhtoGEM model can reach the experimental
% parameters and the majority of topUsedEnzymes are less than 10% of
% protein pool, proceed with the generate_protModels pipeline for the
% proteomics integration.

% sigma factor for each condition has been fitted and saved 'results/enhanceGEM_pipeline'

clear avgProtData tmp

cd ../..

%% II.generate_protModels pipeline

grouping   = [2 2 2 2 2 2];   %Number represents replicates per condition, count of numbers - how many conditions
flexFactor = 1.05;  %Allowable flexibilization factor for fixing carbon uptake rate

cd([root '/GECKO/geckomat/utilities/integrate_proteomics'])
generate_protModels(ecModel,grouping,'ecRhtoGEM',ecModel_batch);


% The idea would be at this point to load models that are saved on GECKO
% folder, because they are not automatically loaded in comparison with
% 'enhanceGEM' pipeline. Can they maybe be loaded automatically?
load('GECKO/models/prot_constrained/ecYeastGEM/ecYeastGEM_Xexp.mat');
% ..and so on for all others, as they would be required for the next step


clear fileNames flexFactor i

% Some questions in regard to 'generate_protModels' pipeline:
% should graphical output of abundances and usage be saved?
     %saveas(gca,fullfile('../../../..','results','generate_protModels_pipeline','abundances.jpg'));
     %saveas(gca,fullfile('../../../..','results','generate_protModels_pipeline','usages.jpg'));
% how to: save auto-flexibilized enzymes from 1st round? currently appears only in command window, "Limiting abundance for..."
% how to: save modifiedEnzymes_.txt (2nd round of auto flexibilization) at
%    '../../../../results/generate_protModels_pipeline'? currently it is
%    saved at same place where ecModelP


% Note: Models are saved with prot_pool_exchange as an objective function,
% constrained with GR (lb), carbon uptake (b, flexibilized), byproducts
% (ub). If 'chemostat constraints too stringent' models are saved with
% carbon uptake as objective

% At this point don't test for growth optimization yet?

%% III.Add ribosome subunits

%Load proteomics data
fID       = fopen('GECKO/Databases/abs_proteomics.txt');
prot.cond = textscan(fID,['%s' repmat(' %s',1,14)],1); %n-1
prot.data = textscan(fID,['%s %s' repmat(' %f',1,12)],'TreatAsEmpty',{'NA','na','NaN'}); %n-2
prot.cond = [prot.cond{3:end}];
prot.IDs  = prot.data{1};
prot.data = cell2mat(prot.data(3:end));
fclose(fID);

%Set some additional parameters
oxPhos = ecModel.rxns(startsWith(ecModel.rxns,params.oxPhos));
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

%Get indexes for carbon source uptake and biomass pseudoreactions
positionsEC(1) = find(strcmpi(ecModel.rxnNames,params.c_source));
positionsEC(2) = find(strcmpi(ecModel.rxns,params.bioRxn));

% Load ribosome information
fid=fopen('data/ribosome.txt');
data=textscan(fid,'%q %q %q %q %q %q','HeaderLines',1,'Delimiter','\t');
fclose(fid);
ribo=data{1};

%Keep only the subunits that are also in the proteomics data
[~,repRibo.dataIdx]=ismember(ribo,prot.IDs);
repRibo.protein=ribo(repRibo.dataIdx>0);
repRibo.dataIdx(repRibo.dataIdx==0)=[];
repRibo.avgLevel=mean(prot.data(repRibo.dataIdx,:),2);

%Density plot to see distribution of subunit abundances
[f,xi]=ksdensity(log10(repRibo.avgLevel),'Bandwidth',0.1);
plot(xi,f);
xlabel('Subunit abundance (log10(mmol/gDCW))');
ylabel('Density');
title('Distribution of average ribosomal subunit abundances');
saveas(gca,fullfile('results','ribosome_integration','average_riboSubunit_abundance.jpg'));

%Only include ribosomal subunite with abundance over 1e-5 mmol/gDCW 
rmRibo=repRibo.avgLevel<1e-5 | isnan(repRibo.avgLevel);
ribo=repRibo.protein(~rmRibo);

%Prepare subunit information to be added to each model
[~,idx]=ismember(ribo,data{1});
enzNames=data{2}(idx);
enzGenes=data{3}(idx);
MWs=str2double(data{5}(idx))/1000;
sequences=data{6}(idx);
pathways={'Ribosome'};

clear repRibo rmRibo data idx fid ans f xi

ecModels{1}=ecModelP_Xexp;
ecModels{2}=ecModelP_XNlim;
ecModels{3}=ecModelP_Aexp;
ecModels{4}=ecModelP_ANlim;
ecModels{5}=ecModelP_GexpUrea;
ecModels{6}=ecModelP_GNlimUrea;

cd([root '/GECKO/geckomat/utilities/integrate_proteomics'])

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
    abundances   = prot.data(:,repl.first(j):repl.last(j));
    [pIDs, filtAbundances] = filter_ProtData(prot.IDs,abundances,1.96,true);
    %cd ../../../..
    sol=solveLP(model_P);
    for i=riboExchId(1):riboExchId(2)
        protId=regexprep(model_P.rxnNames{i},'prot_(......).*','$1');
        k=find(strcmp(protId,pIDs));
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
    ecModels{j}=model_P;
    eval(['ecModelP_' fermParams.conds{j} ' = model_P;']);
    save(fullfile('../../../../models',['ecRhtoGEM_P_' fermParams.conds{j}]), ['ecModelP_' fermParams.conds{j}]);
end
adjusted=adjusted';
fid=fopen(fullfile('results','ribosome_integration','modifiedRibosomeSubunits.txt'),'w');
fprintf(fid,'protein_IDs previous_values modified_values condition\n');
fprintf(fid,'%s %f %f %s\n',adjusted{:});
fclose(fid);

clear abundances aaMetIdx filtAbundances i j genesToAdd metsToAdd mmolAA
clear protId pathways MWs enzGenes enzNames enzAdjust protMetIdx
clear protRxnIdx sequences riboExchId riboKcat riboToAdd rxnsToAdd ribo
clear ecModels model_P
clear pIDs protData fermParams byProds
clear sol fID adjusted GECKO_path k kegg
clear params grouping flexFactor prot oxPhos repl positionsEC

cd root

%% Get fluxes and enzyme usages to each reaction

% optimize for growth rate
conditions.abbrev = {'Xexp','XNlim','Aexp','ANlim','GexpUrea','GNlimUrea'};
conditions.exch.rxns = {{'r_1718_REV','r_2104','r_4340','r_1718_REV','r_2111'},...    % D-xylose uptake, xylitol production, D-arabinitol production
    {'r_1718_REV','r_2104','r_4340','r_1718_REV','r_2111'},...                        % D-xylose uptake, xylitol production, D-arabinitol production
    {'r_1718_REV','r_1634_REV','r_1687','r_2111'},...                    % block D-xylose uptake, allow acetate uptake, citrate(3-) production
    {'r_1718_REV','r_1634_REV','r_1687_REV','r_2111'},...                % block D-xylose uptake, allow acetate uptake, citrate(3-) uptake
    {'r_1718_REV','r_1654_REV','r_1714_REV','r_2091_REV','r_1808','r_1714_REV','r_2111'},...  % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake, glycerol production
    {'r_1718_REV','r_1654_REV','r_1714_REV','r_2091_REV','r_1714_REV','r_2111'}};             % block D-xylose uptake, ammonium uptake, allow D-glucose uptake, urea uptake
conditions.exch.value = {[1.74,0.228,0.372,0,0],...                 % Xexp
    [0.4345,0.004,0.077,0,0],...                                    % XNlim
    [0,6.941,0.127,0],...                                         % Aexp
    [0,1.9706,0.033,0],...                                        % ANlim
    [0,0,3,1000,0.071,0,0],...                                      % GexpUrea
    [0,0,0.415,1000,0,0]};                                          %GNlimUrea
conditions.exch.lbub = {{'ub','ub','ub','lb','lb'},...                    % Xexp
    {'ub','ub','ub','lb','lb'},...                                        % XNlim
    {'ub','ub','ub','lb'},...                                        % Aexp
    {'ub','ub','ub','lb'},...                                        % ANlim
    {'ub','ub','ub','ub','ub','lb','lb'},...                              % GexpUrea
    {'ub','ub','ub','ub','lb','lb'}};                                     % GNlimUrea


for i = 1:numel(conditions.abbrev) % Loop through the conditions
    ecModelP_Tmp = ecModelP; % each ecModelP % Work on a temporary model structure, leaving the original untouched for the next condition
    ecModelP_Tmp = setParam(ecModelP_Tmp, conditions.exch.lbub{i}, ...
        conditions.exch.rxns{i}, conditions.exch.value{i});
    ecModelP_Tmp=changeObjective(ecModelP_Tmp,'r_2111');
    fprintf(['\n=========== Results from model: ' conditions.abbrev{i}, ' ===========\n\n'])    
    printConstraints(modelTmp_batch,-1000,1000);
    sol = solveLP(ecModelP_Tmp,1);
    printFluxVector(ecModelP_Tmp, sol.x, 'true', 'true');
    [absUsage,capUsage] = enzymeUsage(ecModelP_Tmp,sol.x)
    conditions.ecModelP{i} = ecModelP_Tmp;
end
