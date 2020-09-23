%
% reconstruct_ecRhtoGEM
%
%   Reconstruction of enzyme-constrained R.toruloides genome-scale metabolic
%   model (ecRhtoGEM). Firstly, creates ec-model with protein pool using
%   optimzed parameters (most importantly, manually curated kcat values).
%   Secondly, generates ec-model with individually constrained proteomics
%   measurements for the modelled enzymes. Thirdly, adds ribosomal subunits
%   to the ec-models by adding a translation pseudoreaction.
%
%   Last modified: 2020-09-22
%

%% Prepare software
code = pwd();

%Clone the necessary repos:
cd ..
%delete GECKO in case that a previous copy exists here
if isfolder('GECKO') 
    rmdir ('GECKO','s')
end
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git('pull')
cd ..

% Replace custom GECKO scripts
fileNames = struct2cell(dir('customGECKO_XP1_v2_C'));
fileNames = fileNames(1,:);
fileNames(startsWith(fileNames,'.')) = [];
for i = 1:length(fileNames)
    GECKO_path = dir(['GECKO/**/' fileNames{i}]);
    copyfile(['customGECKO_XP1_v2_C' filesep fileNames{i}],GECKO_path.folder)
    disp(['Replaced ' fileNames{i} ' at ' GECKO_path.folder '\'])
end

delete relative_proteomics.txt
copyfile('customGECKO_XP1_v2_C/relative_proteomics.txt','GECKO/Databases','f')
delete GECKO/Databases/prot_abundance.txt

%% Load model
model    = load('models/rhto_edit_v1_XP1.mat');
model    = model.model;
modelVer = model.description(strfind(model.description,'_v')+1:end);

%% I.enhanceGEM pipeline: checklist for GECKO
%
%   /databases: uniprot.tab, relative_proteomics.txt, chemostatData.tsv;
%   /geckomat/getModelParameters.m: sigma, Ptot, gR_exp, c_source,
%   exch_name{2};
%   /geckomat/enhanceGEM.m: disp(['Sigma factor (not fitted)...
%   /geckomat/limit_proteins/getConstrainedModel.m: OptSigma=sigma;
%   ecModel_batch=...currentEnzymeUB;
%   /geckomat/kcat_sensitivity_analysis/changeMedia_batch.m: %block xylose
%   and oxygen production;
%   /geckomat/change_model/manualModifications.m: %growth limiting Kcats
%   section;
%
%% Run GECKO/enhanceGEM pipeline
cd GECKO
GECKOver = git('describe --tags');
cd geckomat/get_enzyme_data
updateDatabases;
cd ..
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA','ecYeastGEM',modelVer);
save('../../models/ecModel_XP1_v2_C.mat','ecModel')
save('../../models/ecModel_batch_XP1_v2_C.mat','ecModel_batch')

% Check model solutions: set experimental conditions
printConstraints(ecModel_batch,-1000,1000);
tempModel = setParam(ecModel_batch,'ub','r_1718_REV',1.74)  %xylose uptake
tempModel = setParam(tempModel,'eq','r_2104',0.228)         %xylitol exchange
tempModel = setParam(tempModel,'eq','r_4340',0.372)         %D-arabinitol
printConstraints(tempModel,-1000,1000);

% Run flux balance analysis and save results
solveLP(tempModel,1)
% 1. fluxes
%printFluxVector(tempModel, ans.x, 'true', 'true');
printFluxes(tempModel,ans.x,false,0,fullfile('../..','results','enhanceGEM_pipeline',['allFluxes_XP1_v2_C.csv']),'%rxnID\t%rxnName\t%eqn\t%flux\n');
% 2. top used enzymes (how to save?)
%topUsedEnzymes(ans.x,tempModel,{''},{''},false)
topUsedEnzymes(ans.x,tempModel,{''},{''})
% 3. sigma fitted from previous optimization (already saved in
% results/enhanceGEM_pipeline)

%    If the solution shows that the ecRhtoGEM model can reach the experimental
% parameters and the majority of topUsedEnzymes are less than 10% of
% protein pool, proceed with the generate_protModels pipeline for the
% proteomics integration.

save('../../models/ecModel_batch_XP1_v2_C_tempModel.mat','tempModel')
exportToExcelFormat(tempModel,'../../models/ecModel_batch_XP1_v2_C_tempModel.xlsx')
cd ../..

%% II.generate_protModels pipeline: checklist for GECKO
%
%   /databases: abs_proteomics.txt, fermentationData.txt,
%   chemostatData.tsv and relative_proteomics.txt (for the protein pool);
%   /geckomat/getModelParameters.m: sigma=1, Ptot, gR_exp, c_source,
%   exch_names{2};
%   /geckomat/kcat_sensitivity_analysis/changeMedia_batch.m: %block glucose
%   and oxygen production;
%   
% Replace custom GECKO scripts
fileNames = struct2cell(dir('customGECKO_XP1_v2_C_proteomics'));
fileNames = fileNames(1,:);
fileNames(startsWith(fileNames,'.')) = [];
for i = 1:length(fileNames)
    GECKO_path = dir(['GECKO/**/' fileNames{i}]);
    copyfile(['customGECKO_XP1_v2_C_proteomics' filesep fileNames{i}],GECKO_path.folder)
    disp(['Replaced ' fileNames{i} ' at ' GECKO_path.folder '\'])
end

%% generate_protModels pipeline: incorporate proteomics

grouping   = [3];   %Our dataset contains three replicates per condition
flexFactor = 1.05;  %Allowable flexibilization factor for fixing carbon uptake rate

cd GECKO/geckomat/utilities/integrate_proteomics
generate_protModels(ecModel,grouping,'ecYeastGEM',ecModel_batch);

load('../../../models/prot_constrained/ecYeastGEM/ecYeastGEM_XP1.mat');
save('../../../../models/ecModelP_XP1_v2_C.mat','ecModelP')

% Run flux balance analysis and save results
% 1. Model information, modified enzymes 2nd round of auto flexibilization,
% etc.
movefile ../../../models/prot_constrained/ecYeastGEM/ ../../../../results/generate_protModels_pipeline
delete ../../../../results/generate_protModels_pipeline/ecYeastGEM_XP1.mat
% 2. Auto-flexibilized enzymes from 1st round (how to save? currently
% appears only in command window, "Limiting abundance for...")

checkObjective(ecModelP);
ecModelP=changeObjective(ecModelP,'r_2111');
printConstraints(ecModelP,-1000,1000);
% set constraints according to experimental results (change from flexibilized
% rates when model was generated)
tempModel = setParam(ecModelP,'ub','r_1718_REV',1.74)  %xylose uptake
tempModel = setParam(tempModel,'lb','r_1718_REV',0)
tempModel = setParam(tempModel,'lb','r_2111',0)        %growth
tempModel = setParam(tempModel,'lb','r_4041',0)        %biomass pseudoreaction
tempModel = setParam(tempModel,'eq','r_2104',0.228)    %xylitol exchange
tempModel = setParam(tempModel,'eq','r_4340',0.372)    %D-arabinitol

printConstraints(tempModel,-1000,1000);

solveLP(tempModel,1)
% 3. fluxes
%printFluxVector(tempModel, ans.x, 'true', 'true');
printFluxes(tempModel,ans.x,false,0,fullfile('../../../..','results','generate_protModels_pipeline',['allFluxes_XP1_v2_C.csv']),'%rxnID\t%rxnName\t%eqn\t%flux\n');

save('../../../../models/ecModelP_XP1_v2_C_tempModel.mat','tempModel')
exportToExcelFormat(tempModel,'../../../../models/ecModelP_XP1_v2_C_tempModel.xlsx')
cd
clear fileNames flexFactor grouping i

%% III.Add ribosome subunits

%Load proteomics data
fID       = fopen('../../../databases/abs_proteomics.txt');
prot.cond = textscan(fID,['%s' repmat(' %s',1,4)],1);
prot.data = textscan(fID,['%s %s' repmat(' %f',1,3)],'TreatAsEmpty',{'NA','na','NaN'});
prot.cond = [prot.cond{3:end}];
prot.IDs  = prot.data{1};
prot.data = cell2mat(prot.data(3:end));
fclose(fID);

%Load total protein content and fermentation data
fID       = fopen('../../../databases/fermentationData.txt');
byProds   = textscan(fID,['%s' repmat(' %s',1,7)],1,'Delimiter','\t');
data      = textscan(fID,['%s' repmat(' %f',1,7)],'TreatAsEmpty',{'NA','na','NaN'});
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
grouping=[3];
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

ecModels{1}=tempModel;

% Load ribosome information
fid=fopen('../../../../data/ribosome.txt');
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
% why f is NaN?
plot(xi,f);
xlabel('Subunit abundance (log10(mmol/gDCW))');
ylabel('Density');
title('Distribution of average ribosomal subunit abundances');
saveas(gca,fullfile('..','results','modelGeneration','average_riboSubunit_abundance.pdf'));

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
    model=ecModels{j};
    %Add enzymes and genes
    model.enzymes=[model.enzymes;ribo];
    model.enzNames=[model.enzNames;enzNames];
    model.enzGenes=[model.enzGenes;enzGenes];
    genesToAdd.genes=enzGenes;
    genesToAdd.geneShortNames=enzNames;
    model=addGenesRaven(model,genesToAdd);
    model.MWs=[model.MWs;MWs];
    model.sequences=[model.sequences;sequences];
    model.pathways=[model.pathways;repmat(pathways,numel(ribo),1)];

    %Include new amino acid pseudometabolite
    metsToAdd.mets={'aminoAcids'};
    metsToAdd.metNames={'amino acids for protein'};
    metsToAdd.compartments={'c'};
    model=addMets(model,metsToAdd);
    
    %Modify protein pseudoreaction to produce amino acid pseudometabolite
    protMetIdx=getIndexes(model,'protein','metnames');
    aaMetIdx=getIndexes(model,'aminoAcids','mets');
    protRxnIdx=getIndexes(model,'r_4047','rxns');
    model.S([protMetIdx,aaMetIdx],protRxnIdx)=[0,1];
    
    %Include ribosome subunits as pseudometabolites
    riboToAdd.mets=strcat('prot_',sort(ribo));
    riboToAdd.metNames=riboToAdd.mets;
    riboToAdd.compartments='c';
    model=addMets(model,riboToAdd);
    
    %Determine stoichiometric coefficient of ribosomes
    mmolAA=full(model.S(:,protRxnIdx));
    mmolAA=-sum(mmolAA(mmolAA<0)); %mmol amino acids in biomass
    riboKcat=10.5*3600; %10.5 aa/sec (doi:10.1042/bj1680409) -> aa/hour
    riboKcat=mmolAA/riboKcat; %compensate for the amount of amino acids elongated
    % Include new reaction representing ribosomes (=translation)
    rxnsToAdd.rxns={'translation'};
    rxnsToAdd.mets=[model.mets(protMetIdx),model.mets(aaMetIdx),riboToAdd.mets'];
    rxnsToAdd.stoichCoeffs=[1,-1,repmat(-riboKcat,1,numel(riboToAdd.mets))];
    rxnsToAdd.subSystem={'sce03010  Ribosome'};
    rxnsToAdd.grRules={strjoin(enzGenes,' and ')};
    model=addRxns(model,rxnsToAdd);
    
    %Add exchange reactions or usage reactions
    riboExchId=numel(model.rxns)+1;
    model=addExchangeRxns(model,'in',riboToAdd.mets);
    riboExchId(2)=numel(model.rxns);
    %Add UB for enzyme exchange reactions based on measurements.
    %If UB is too low, then adjust to value predicted by model.
    %cd geckomat/utilities/integrate_proteomics
    abundances   = prot.data(:,repl.first(j):repl.last(j));
    [pIDs, filtAbundances] = filter_ProtData(prot.IDs,abundances,1.96,true);
    %cd ../../../..
    sol=solveLP(model);
    for i=riboExchId(1):riboExchId(2)
        protId=regexprep(model.rxnNames{i},'prot_(......).*','$1');
        k=find(strcmp(protId,pIDs));
        if ~isempty(k)
            model.rxns{i}=['prot_' protId '_exchange'];
            model.rxnNames{i}=['prot_' protId '_exchange'];
            model.concs(strcmp(protId,model.enzymes))=filtAbundances(k);
            if sol.x(i)<filtAbundances(k)
                model.ub(i)=filtAbundances(k);
            else
                model.ub(i)=sol.x(i)*1.01;
                adjusted{end+1,1}=protId;
                adjusted{end,2}=filtAbundances(k);
                adjusted{end,3}=sol.x(i)*1.01;
                adjusted{end,4}=flux.conds{j};
                fprintf('%s abundance adjusted. Measured: %e / Adjusted: %e\n',adjusted{end,1:3})
            end
        else
            model.ub(i)=Inf;
            model.rxns{i}=['draw_prot_' protId];
            model.rxnNames{i}=['draw_prot_' protId];
        end
    end
    ecModels{j}=model;
    eval(['ecModelP_' flux.conds{j} ' = model;']);
    save(fullfile('../../../..','models',['ecModel_P_' flux.conds{j}]), ['ecModelP_' flux.conds{j}]);
end
adjusted=adjusted';
fid=fopen(fullfile('../../../..','results','ribosome_integration','modifiedRibosomeSubunits.txt'),'w');
fprintf(fid,'protein_IDs previous_values modified_values condition\n');
fprintf(fid,'%s %f %f %s\n',adjusted{:});
fclose(fid);

exportToExcelFormat(ecModelP_XP1,'../../../../models/ecModel_P_XP1.xlsx')

clear cond abundances aaMetIdx filtAbundances i j genesToAdd metsToAdd mmolAA
clear protId prot pIDs pathways MWs model enzGenes enzNames enzAdjust protMetIdx
clear protRxnIdx repl sequences riboExchId riboKcat riboToAdd rxnsToAdd
clear sequences sol fid adjusted

cd ../../../..

