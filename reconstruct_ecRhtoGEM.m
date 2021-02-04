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
%   Last modified: 2021-02-03
%

% Prepare COBRA and set repo root path
initCobraToolbox
root = pwd; % Get the root directory of the folder

%% Load model
model    = importModel(fullfile(root,'models','model_edit.xml'));
modelVer = model.description(strfind(model.description,'_v')+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% Switch GECKO to the last commit (4480b07) that is part of GECKO pull
% request #130: https://github.com/SysBioChalmers/GECKO/pull/130
% This contains the necessary changes in measureAbundance, constrainEnzymes
% and generate_protModels.
git('switch -d 4408b07be48be5acaabe2acbf8c5724ede4013f2')
cd ..

% From here define a new loop that generates the condition-specific batch
% models.
% The folder 'customGECKO' contains identical files for all conditions,
% folders 'customGECKO_' contain condition-specific files

% files saveECmodel and sumBioMass are from older GECKO versions to avoid some errors.

% TODO: updateDatabases, addProtein, getEnzymeCodes, convertToEnzymeModel
% are in GECKO PR #122 and custom functions can be removed once PR is merged

% files manualModifications, relative_proteomics, uniprot.tab contain
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

delete relative_proteomics.txt;
copyfile('customGECKO/relative_proteomics.txt','GECKO/Databases','f');


% Start GECKO/enhanceGEM pipeline
cd GECKO
delete databases/prot_abundance.txt;  % if not deleted,f factor will not be calculated;
GECKOver = git('describe --tags');
cd geckomat/get_enzyme_data
updateDatabases;
cd ..
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA','ecRhtoGEM',modelVer);

% Ignore the sigma fitting, manually set sigma to 1
params = getModelParameters();
cd limit_proteins
f = measureAbundance(ecModel_batch.enzymes);
ecModel_batch = updateProtPool(ecModel_batch,params.Ptot,f); % Should f not be multiplied with params.sigma?

% Overwrite the files exported by enhanceGEM, now with the new pool UB
cd ../../models
ecModel_batch = saveECmodel(ecModel_batch,'COBRA','ecRhtoGEM_batch',modelVer);
cd ../geckomat


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
    conditions.ecModel_batch{i} = modelTmp_batch;
end



%if tempModel_batch results are successful, constrain ecModels accordingly:
ecModel_batch = setParam(ecModel_batch,'eq','r_2104',0.004);       %xylitol production
ecModel_batch = setParam(ecModel_batch,'eq','r_4340',0.077);       %arabitol production
ecModel = setParam(ecModel,'ub','r_1718_REV',0.4345);              %xylose uptake
ecModel = setParam(ecModel,'eq','r_2104',0.004);                   %xylitol uptake
ecModel = setParam(ecModel,'eq','r_4340',0.077);                   %arabitol uptake

%if tempModel_batch results are successful, constrain ecModels accordingly:
%     remove byproduct constraints:
%     (with citrate uptake 'r_1687_REV'0.033 modelP can't be generated
%     probably because in that case it would be the 2nd uptake constraint)
ecModel_batch = setParam(ecModel_batch,'lb','r_1687',0);
ecModel_batch = setParam(ecModel_batch,'ub','r_1687',Inf);
%     constrain ecModel:
ecModel = setParam(ecModel,'ub','r_1634_REV',1.9706);                   % acetate uptake
ecModel = setParam(ecModel,'lb','r_1687',0);
ecModel = setParam(ecModel,'ub','r_1687',Inf);


%%% Gexp_in_urea:
printConstraints(ecModel_batch,-1000,1000);
tempModel_batch = setParam(ecModel_batch,'ub','r_1714_REV',3);          % glucose uptake
tempModel_batch = setParam(tempModel_batch,'eq','r_1808',0.071);        % glycerol production
tempModel_batch = setParam(tempModel_batch,'ub','r_1654_REV',0);        % block ammonium uptake
tempModel_batch = setParam(tempModel_batch,'ub','r_2091_REV',Inf);      % allow urea uptake
tempModel_batch = setParam(tempModel_batch,'ub','prot_pool_exchange',0.2);

ecModel_batch = setParam(ecModel_batch,'ub','r_1654_REV',0);
ecModel_batch = setParam(ecModel_batch,'ub','r_2091_REV',Inf);
ecModel_batch = setParam(ecModel_batch,'ub','r_1808',0,.071);   % ?

ecModel = setParam(ecModel,'ub','r_1714_REV',3.0);

ecModel = setParam(ecModel,'ub','r_1654_REV',0);
ecModel = setParam(ecModel,'ub','r_2091_REV',Inf);

%%% GNim_in_urea:
tempModel_batch = setParam(tempModel_batch,'ub','r_1714_REV',0.415);

ecModel = setParam(ecModel,'ub','r_1714_REV',0.415);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printConstraints(tempModel_batch,-1000,1000);
solveLP(tempModel_batch,1);
printFluxVector(tempModel_batch, ans.x, 'true', 'true');
topUsedEnzymes(ans.x,tempModel_batch,{''},{''},false);    %how to save and write a table?


% If the solution shows that the ecRhtoGEM model can reach the experimental
% parameters and the majority of topUsedEnzymes are less than 10% of
% protein pool, proceed with the generate_protModels pipeline for the
% proteomics integration.

% sigma factor for each condition has been fitted and saved 'results/enhanceGEM_pipeline'


cd ../..

%% II.generate_protModels pipeline

% Replace custom GECKO scripts
fileNames = struct2cell(dir('customGECKO_proteomics'));
fileNames = fileNames(1,:);
fileNames(startsWith(fileNames,'.')) = [];
for i = 1:length(fileNames)
    GECKO_path = dir(['GECKO/**/' fileNames{i}]);
    copyfile(['customGECKO_proteomics' filesep fileNames{i}],GECKO_path.folder)
    disp(['Replaced ' fileNames{i} ' at ' GECKO_path.folder '\'])
end

% Incorporate proteomics

grouping   = [2];   %Number represents replicates per condition, count of numbers - how many conditions
flexFactor = 1.0;  %Allowable flexibilization factor for fixing carbon uptake rate

cd GECKO/geckomat/utilities/integrate_proteomics
generate_protModels(ecModel,grouping,'ecYeastGEM',ecModel_batch);


% how to: save at '../../../../results/generate_protModels_pipeline' graphical output?
     %saveas(gca,fullfile('../../../..','results','generate_protModels_pipeline','abundances.jpg'));
     %saveas(gca,fullfile('../../../..','results','generate_protModels_pipeline','usages.jpg'));
% how to: save auto-flexibilized enzymes from 1st round? currently appears only in command window, "Limiting abundance for..."
% how to: save modifiedEnzymes_.txt (2nd round of auto flexibilization) at
%    '../../../../results/generate_protModels_pipeline'? currently it is
%    saved at same place where ecModelP

% is it possible to have ecModelP on workspace automatically? otherwise:
load('../../../models/prot_constrained/ecYeastGEM/ecYeastGEM_Xexp.mat');
load('../../../models/prot_constrained/ecYeastGEM/ecYeastGEM_XNlim.mat');
load('../../../models/prot_constrained/ecYeastGEM/ecYeastGEM_Aexp.mat');
load('../../../models/prot_constrained/ecYeastGEM/ecYeastGEM_ANlim.mat');
load('../../../models/prot_constrained/ecYeastGEM/ecYeastGEM_Gexp.mat');
load('../../../models/prot_constrained/ecYeastGEM/ecYeastGEM_GNlim.mat');



%%% Xexp:
printConstraints(ecModelP,-1000,1000);
tempModelP = setParam(ecModelP,'ub','r_1718_REV',1.74)   %xylose uptake
tempModelP = setParam(tempModelP,'lb','r_1718_REV',0)
tempModelP = setParam(tempModelP,'lb','r_2111',0)        %growth
tempModelP = setParam(tempModelP,'lb','r_4041',0)        %biomass pseudoreaction
tempModelP = setParam(tempModelP,'eq','r_2104',0.228)    %xylitol exchange
tempModelP = setParam(tempModelP,'eq','r_4340',0.372)    %D-arabinitol

checkObjective(tempModelP);
tempModelP=changeObjective(tempModelP,'r_2111');


%%% XNlim:
printConstraints(ecModelP,-1000,1000);
tempModelP = setParam(ecModelP,'ub','r_1718_REV',0.4345)%xylose uptake
tempModelP = setParam(tempModelP,'lb','r_1718_REV',0)
tempModelP = setParam(tempModelP,'lb','r_2111',0)        %growth
tempModelP = setParam(tempModelP,'lb','r_4041',0)        %biomass pseudoreaction
tempModelP = setParam(tempModelP,'eq','r_2104',0.004)    %xylitol exchange
tempModelP = setParam(tempModelP,'eq','r_4340',0.077)    %D-arabinitol

checkObjective(tempModelP);
tempModelP=changeObjective(tempModelP,'r_2111');


%%% Aexp:
printConstraints(ecModelP,-1000,1000);
tempModelP = setParam(ecModelP,'ub','r_1634_REV',6.941);  %acetate uptake
tempModelP = setParam(tempModelP,'lb','r_1634_REV',6.941);%if C balance < 1
tempModelP = setParam(tempModelP,'lb','r_2111',0);        %growth
tempModelP = setParam(tempModelP,'lb','r_4041',0);        %biomass pseudoreaction
tempModelP = setParam(tempModelP,'eq','r_1687',0.127);    %citrate(3-) secretion

checkObjective(tempModelP);
tempModelP=changeObjective(tempModelP,'r_2111');


%%% ANlim:
printConstraints(ecModelP,-1000,1000);
tempModelP = setParam(ecModelP,'ub','r_1634_REV',1.9706) %carbon uptake
tempModelP = setParam(tempModelP,'lb','r_1634_REV',1.9706)
tempModelP = setParam(tempModelP,'lb','r_2111',0)        %growth
tempModelP = setParam(tempModelP,'lb','r_4041',0)        %biomass pseudoreaction 
tempModelP = setParam(tempModelP,'eq','r_1687_REV',0.033)%citrate(3-) uptake 

checkObjective(tempModelP);
tempModelP=changeObjective(tempModelP,'r_2111');


%%% Gexp_in_urea:
printConstraints(ecModelP,-1000,1000);
tempModelP = setParam(ecModelP,'ub','r_1714_REV',3.0) %carbon uptake
tempModelP = setParam(tempModelP,'lb','r_1714_REV',0)
tempModelP = setParam(tempModelP,'lb','r_2111',0)        %growth
tempModelP = setParam(tempModelP,'lb','r_4041',0)        %biomass pseudoreaction
tempModelP = setParam(tempModelP,'eq','r_1808',0.071)    %glycerol production

checkObjective(tempModelP);
tempModelP=changeObjective(tempModelP,'r_2111');


%%% GNlim_in_urea:
printConstraints(ecModelP,-1000,1000);
tempModelP = setParam(ecModelP,'ub','r_1714_REV',0.415) %carbon uptake
tempModelP = setParam(tempModelP,'lb','r_1714_REV',0)
tempModelP = setParam(tempModelP,'lb','r_2111',0)        %growth
tempModelP = setParam(tempModelP,'lb','r_4041',0)        %biomass pseudoreaction

checkObjective(tempModelP);
tempModelP=changeObjective(tempModelP,'r_2111');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printConstraints(tempModelP,-1000,1000);
solveLP(tempModelP,1);
printFluxVector(tempModelP, ans.x, 'true', 'true');
%topUsedEnzymes(ans.x,tempModelP,{''},{''},false)
%[absUsage,capUsage] = enzymeUsage(tempModelP,ans.x)
%exportToExcelFormat(tempModelP,'ecModelP_Gexp_v2.xlsx');

clear fileNames flexFactor grouping i


%% III.Add ribosome subunits

%Load proteomics data
fID       = fopen('../../../Databases/abs_proteomics.txt');
prot.cond = textscan(fID,['%s' repmat(' %s',1,3)],1); %n-1
prot.data = textscan(fID,['%s %s' repmat(' %f',1,2)],'TreatAsEmpty',{'NA','na','NaN'}); %n-2
prot.cond = [prot.cond{3:end}];
prot.IDs  = prot.data{1};
prot.data = cell2mat(prot.data(3:end));
fclose(fID);

%Load total protein content and fermentation data
fID       = fopen('../../../Databases/fermentationData.txt');
% xylose (2 byproducts)
byProds   = textscan(fID,['%s' repmat(' %s',1,7)],1,'Delimiter','\t'); %n-1
data      = textscan(fID,['%s' repmat(' %f',1,7)],'TreatAsEmpty',{'NA','na','NaN'}); %n-1
% acetate & glucose (1 byproduct)
byProds   = textscan(fID,['%s' repmat(' %s',1,6)],1,'Delimiter','\t'); %n-1
data      = textscan(fID,['%s' repmat(' %f',1,6)],'TreatAsEmpty',{'NA','na','NaN'}); %n-1
fclose(fID);

flux.conds     = data{1};
flux.Ptot      = data{2};
flux.Drate     = data{3};
flux.GUR       = data{4};
flux.CO2prod   = data{5};
flux.OxyUptake = data{6};
flux.byP_flux  = [data{7:end}];
flux.byP_flux(isnan(flux.byP_flux))=0;
flux.byProds   = [byProds{7:end}];

% Load model and additional parameters
cd ../..
params = getModelParameters;
cd utilities/integrate_proteomics

%Set some additional parameters
oxPhos = ecModelP.rxns(startsWith(ecModelP.rxns,params.oxPhos));
grouping=[2];
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
positionsEC(1) = find(strcmpi(ecModelP.rxnNames,params.c_source));
positionsEC(2) = find(strcmpi(ecModelP.rxns,params.bioRxn));
clear ans fID data byProds fileNames GECKO_path i fID

ecModels{1}=tempModelP;
%ecModels{2}=tempModelP;

% Load ribosome information
fid=fopen('../../../Databases/ribosome.txt');
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
saveas(gca,fullfile('../../../..','results','ribosome_integration','average_riboSubunit_abundance.jpg'));

%Only include ribosomal subunite with abundance over 1e-5 mmol/gDCW 
rmRibo=repRibo.avgLevel<1e-5 | isnan(repRibo.avgLevel);
ribo=repRibo.protein(~rmRibo);

%Prepare subunit information to be added to each model
[~,idx]=ismember(ribo,data{1});
enzNames=data{2}(idx);
enzGenes=data{3}(idx);
MWs=str2double(data{5}(idx))/1000;
sequences=data{6}(idx);
pathways={'sce03010  Ribosome'};

clear repRibo rmRibo idx data fid ans f xi

% Modify reactions
adjusted=cell.empty();
for j=1:numel(ecModels);
    disp(['Add ribosome subunits to condition: ' flux.conds{j}])
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
    rxnsToAdd.subSystem={'sce03010  Ribosome'};
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
                adjusted{end,4}=flux.conds{j};
                fprintf('%s abundance adjusted. Measured: %e / Adjusted: %e\n',adjusted{end,1:3})
            end
        else
            model_P.ub(i)=Inf;
            model_P.rxns{i}=['draw_prot_' protId];
            model_P.rxnNames{i}=['draw_prot_' protId];
        end
    end
    ecModels{j}=model_P;
    eval(['ecModelP_' flux.conds{j} ' = model_P;']);
    save(fullfile('../../../..','models',['ecModel_P_' flux.conds{j}]), ['ecModelP_' flux.conds{j}]);
end
adjusted=adjusted';
fid=fopen(fullfile('../../../..','results','ribosome_integration','modifiedRibosomeSubunits.txt'),'w');
fprintf(fid,'protein_IDs previous_values modified_values condition\n');
fprintf(fid,'%s %f %f %s\n',adjusted{:});
fclose(fid);

clear cond abundances aaMetIdx filtAbundances i j genesToAdd metsToAdd mmolAA
clear protId prot pIDs pathways MWs enzGenes enzNames enzAdjust protMetIdx
clear protRxnIdx repl sequences riboExchId riboKcat riboToAdd rxnsToAdd
clear sequences sol fid adjusted
clear oxPhos params positionsEC grouping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get fluxes and enzyme usages to each reaction

cd ../../../..


%%% Xexp:
printConstraints(ecModelP_Xexp,-1000,1000);
checkObjective(ecModelP_Xexp);

sol = solveLP(ecModelP_Xexp,1);
printFluxVector(ecModelP_Xexp, sol.x, 'true', 'true');
[absUsage,capUsage] = enzymeUsage(ecModelP_Xexp,sol.x)
% for now, manually copied ans.x, model.rxns, model.subSystems,
% model.enzymes, absUsage and capUsage in excel and merged with names
% geckopy_data_prepare.ipynb

exportToExcelFormat(ecModelP_Xexp,'models/ecModel_P_Xexp.xlsx');


%%% XNlim:
printConstraints(ecModelP_XNlim,-1000,1000);
checkObjective(ecModelP_XNlim);

sol = solveLP(ecModelP_XNlim,1);
printFluxVector(ecModelP_XNlim, sol.x, 'true', 'true');
[absUsage,capUsage] = enzymeUsage(ecModelP_XNlim,sol.x)

exportToExcelFormat(ecModelP_XNlim,'models/ecModel_P_XNlim.xlsx');


%%% Aexp:
printConstraints(ecModelP_Aexp,-1000,1000);
checkObjective(ecModelP_Aexp);

sol = solveLP(ecModelP_Aexp,1);
printFluxVector(ecModelP_Aexp, sol.x, 'true', 'true');
[absUsage,capUsage] = enzymeUsage(ecModelP_Aexp,sol.x)

exportToExcelFormat(ecModelP_Aexp,'models/ecModel_P_Aexp.xlsx');


%%% ANlim:
printConstraints(ecModelP_ANlim,-1000,1000);
checkObjective(ecModelP_ANlim);

sol = solveLP(ecModelP_ANlim,1);
printFluxVector(ecModelP_ANlim, sol.x, 'true', 'true');
[absUsage,capUsage] = enzymeUsage(ecModelP_ANlim,sol.x)

exportToExcelFormat(ecModelP_ANlim,'models/ecModel_P_ANlim.xlsx');


%%% Gexp:
printConstraints(ecModelP_Gexp,-1000,1000);
checkObjective(ecModelP_Gexp);

sol = solveLP(ecModelP_Gexp,1);
printFluxVector(ecModelP_Gexp, sol.x, 'true', 'true');
[absUsage,capUsage] = enzymeUsage(ecModelP_Gexp,sol.x)

exportToExcelFormat(ecModelP_Gexp,'models/ecModel_P_Gexp.xlsx');


%%% GNlim:
printConstraints(ecModelP_GNlim,-1000,1000);
checkObjective(ecModelP_GNlim);

sol = solveLP(ecModelP_GNlim,1);
printFluxVector(ecModelP_GNlim, sol.x, 'true', 'true');
[absUsage,capUsage] = enzymeUsage(ecModelP_GNlim,sol.x)

exportToExcelFormat(ecModelP_GNlim,'models/ecModel_P_GNlim.xlsx');

