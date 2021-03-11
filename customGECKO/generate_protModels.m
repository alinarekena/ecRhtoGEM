function generate_protModels(ecModel,grouping,name,ecModel_batch)
% generate_protModels
%
% Function that takes an ecModel and constraints it with absolute proteomics 
% datasets for the experimental conditions provided by the user.
%
%   ecModel       (structure) ecModel MATLAB structure (without the total protein pool constraint)
%   grouping      (vector) Number of biological replicates for each
%                 experimental condition in the dataset.
%   name          (string) name string for the ecModel + proteomics and its
%                 container folder (GECKO/models/prot_constrained/name/name...)
%   ecModel_batch (structure) ecModel MATLAB structure with total protein pool
%                 constraint, if provided the process of fitting GAM for
%                 each new protein content is speed-up.
%
% Usage:  generate_protModels(ecModel,grouping,name,ecModel_batch)
%
% Last modified.  Ivan Domenzain 2020-04-06
close all
current = pwd;

if nargin<4
    ecModel_batch = [];
end
%Flexibilization factor for carbon source uptake rate (needed for
%flexibilizeProteins step in constrainEnzymes).
flexFactor = 1.05;
%get model parameters
cd ../..
parameters = getModelParameters;
Ptot_model = parameters.Ptot;
growthRxn  = parameters.exch_names{1};
NGAM       = parameters.NGAM;
if isfield(parameters,'GAM')
    GAM = parameters.GAM;
else
    cd limit_proteins
    GAM = splitGAEC(ecModel_batch);
    cd ..
end
%Get oxPhos related rxn IDs
oxPhos = getOxPhosRxnIDs(ecModel,parameters);
%create subfolder for ecModelProts output files
mkdir(['../models/prot_constrained/' name])
%Get indexes for carbon source uptake and growth reaction
positionsEC(2) = find(strcmpi(ecModel.rxnNames,growthRxn));

%Calculate f-factor from paxDB file, to be used for all conditions
if isfile('../databases/prot_abundance.txt')
    cd limit_proteins
    f = measureAbundance(ecModel.enzymes);
    cd ..
else
	warning('prot_abundance.txt file not found in Databases folder, f factor will be calculated from each condition-specific abundance dataset instead')
    f = [];
end

%Load absolute proteomics dataset [mmol/gDw]
%and fermentation data (GUR, OUR, CO2 production, byProducts, Ptot, Drate)
cd utilities/integrate_proteomics/
[initialProts,absValues,fermData,byProducts] = load_Prot_Ferm_Data(grouping);
conditions = fermData.conds;
Ptot       = fermData.Ptot;
Drate      = fermData.Drate;
GUR        = fermData.GUR;
CO2prod    = fermData.CO2prod;
OxyUptake  = fermData.OxyUptake;
byP_flux   = fermData.byP_flux;
c_source   = fermData.c_source;

%Increase enzyme usage fluxes 1000-fold, to prevent very low fluxes.
models={ecModel,ecModel_batch};
for m=1:2
    protIdx = find(startsWith(models{m}.mets,'prot_'));
    pExcIdx = find(contains(models{m}.rxns,'prot_'));
    for n=1:length(protIdx)
        enzRxnIdx = find(models{m}.S(protIdx(n),:)); % All reactions involving protein
        enzRxnIdx = setdiff(enzRxnIdx,pExcIdx); % Leave out protein exchange
        coeffs = models{m}.S(protIdx(n),enzRxnIdx)*1000;
        models{m}.S(protIdx(n),enzRxnIdx) = coeffs;
    end
end
ecModel=models{1};
ecModel_batch=models{2};
ecModel_batch.ub(end)=ecModel_batch.ub(end)*1000;

%For each condition create a protein constrained model
for i=1:length(conditions)
    cd (current)
    disp(conditions{i})
    c_source_exch = c_source{i};
    positionsEC(1) = find(strcmpi(ecModel.rxnNames,c_source_exch));
    %Extract data for the i-th condition
    abundances   = cell2mat(absValues(1:grouping(i)))*1000;
    absValues    = absValues(grouping(i)+1:end);
    %Filter proteomics data, to only keep high quality measurements
    [pIDs, abundances] = filter_ProtData(initialProts,abundances,1.96,true);
    cd ..
    %correct oxPhos complexes abundances
    if ~isempty(oxPhos)
        for j=1:length(oxPhos)
            [abundances,pIDs] = fixComplex(oxPhos{j},ecModel,abundances,pIDs);
        end
    end
    filteredProts = pIDs;
    %Set minimal medium
    cd ../kcat_sensitivity_analysis
    ecModelP  = changeMedia_batch(ecModel,c_source_exch);
    tempModel = changeMedia_batch(ecModel_batch,c_source_exch);
    cd ../limit_proteins
    if isempty(f) | f==0
        f = measureAbundance(ecModel.enzymes,pIDs,abundances);
    end
    %If the ecModel's protein content is not the same as the Ptot for i-th
    %condition then biomass should be rescaled and GAM refitted to this condition.
    %For fitting GAM a functional model is needed therefore an ecModel with
    if Ptot_model ~= Ptot(i)
        ecModelP = scaleBioMass(ecModelP,Ptot(i),GAM,true);
        disp(' ')
    end 
    %Block production of non-observed metabolites before data incorporation
    %and flexibilization
    expData  = [GUR(i),CO2prod(i),OxyUptake(i)];
    flexGUR  = flexFactor*GUR(i);
    ecModelP = setParam(ecModelP,'ub',positionsEC(1),flexGUR);
    ecModelP = DataConstrains(ecModelP,byProducts,byP_flux(i,:),1.1);
    %Get a temporary model structure with the same constraints to be used
    %for minimal enzyme requirements analysis. For all measured enzymes
    %(present in the dataset) a minimal usage is obtained from a FBA
    %simulation with the ecModel_batch, then 
    tempModel       = DataConstrains(tempModel,byProducts,byP_flux(i,:),1.1);
    tempModel       = setParam(tempModel,'ub',positionsEC(1),flexGUR);
    [matchedEnz,iA] = intersect(pIDs,tempModel.enzymes);
    enzModel        = setParam(tempModel,'lb',positionsEC(2),Drate(i));
    for j=1:length(matchedEnz)
        rxnIndex  = find(contains(tempModel.rxnNames,matchedEnz{j}));
        tempModel = setParam(enzModel,'obj',rxnIndex,-1);
        tempSol   = solveLP(tempModel);
        %Compare enzyme minimum usage with abundance value
        if (tempSol.x(rxnIndex)-abundances(iA(j)))>0
            %Flexibilize limiting values
            disp(['Limiting abundance found for: ' matchedEnz{j} '/Previous value: ' num2str(abundances(iA(j))) ' /New value: ' num2str(tempSol.x(rxnIndex))])
            abundances(iA(j)) = 1.01*tempSol.x(rxnIndex);
        end
        enzIndex = find(contains(tempModel.enzymes,matchedEnz{j}));
    end
    %Get model with proteomics
    disp(['Incorporation of proteomics constraints for ' conditions{i} ' condition'])
    %Get sum of measured protein after filter, adding flexFactor and setting minimum value. 
    %If this is higher than the sum of raw measured protein (sumP), then increase the total 
    %protein content by the same ratio, so that the protein pool is receiving the similar 
    %flexibilization as applied to the measured proteins.
    [ecModelP,usagesT,modificationsT,~,coverage] = constrainEnzymes(ecModelP,f,GAM,Ptot(i)*1000,pIDs,abundances,Drate(i));
    matchedProteins = usagesT.prot_IDs;
    prot_input = {initialProts filteredProts matchedProteins ecModel.enzymes coverage};
    writeProtCounts(conditions{i},prot_input,name); 
    cd (current)
    %NGAM interval for fitting
    interval = [0 5];
    ecModelP = setStressConditions(ecModelP,Drate(i),positionsEC,expData,NGAM,interval);
    %Fix experimental Glucose uptake rate and save models
    cd ..
    ecModelP = setChemostatConstraints(ecModelP,positionsEC,Drate(i),true,0.01,GUR(i));
    %Get optimal flux distribution and display exchange fluxes
    solution = solveLP(ecModelP,1);
    if ~isempty(solution.f)
        fileFluxes = ['../../models/prot_constrained/' name '/fluxes_Exch_' conditions{i} '.txt'];
        printFluxes(ecModelP,solution.x,true,1E-4,fileFluxes)
    end
    cd ../change_model
    [~,~,version] = preprocessModel(ecModelP,'','');
    cd ../../models/prot_constrained
    addpath('..')
    saveECmodel(ecModelP,'RAVEN',name,version);
    rmpath('..')
    %Rename model file names to include conditions{i}
    cd(name)
    fileNames = struct2cell(dir('.'));
    fileNames = fileNames(1,:);
    fileNames(cellfun(@isempty, regexp(fileNames, [name, '(\.\w{3})']))) = [];
    for j = 1:length(fileNames)
        newFileName = regexprep(fileNames{j},[name, '(\.\w{3})'], [name, '_', conditions{i}, '$1']);
        movefile(fileNames{j}, newFileName);
    end
    cd ../../../geckomat/limit_proteins
    %save .txt file
    writetable(usagesT,['../../models/prot_constrained/' name '/enzymeUsages_' conditions{i} '.txt'],'Delimiter','\t')
    writetable(modificationsT,['../../models/prot_constrained/' name '/modifiedEnzymes_' conditions{i} '.txt'],'Delimiter','\t')
    f=[];
end
cd (current)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function condModel = setStressConditions(model,Drate,pos,expData,NGAM,interval)
minProt = true;
ecFlag  = true;
%Rescale biomass set fitted GAM for std conditions and fit non-growth
%associated maintenance for stress conditions
cd ..
condModel = setChemostatConstraints(model,pos,Drate,minProt,0.01);
cd integrate_proteomics
disp(' ')
condModel = fitNGAM(condModel,NGAM,expData,interval,ecFlag);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = DataConstrains(model,compounds,bounds,flexBounds)
if ~isempty(compounds)
    disp('Constraining byproducts exchange fluxes with fermentation data')
    for i=1:length(compounds)
        %Get exchange rxn index
        if le(0,bounds(i))
            inFlux = 0; % Excretion
        else
            inFlux = 1; % Uptake
        end
        rxnNames = strcat(compounds{i}, {' exchange',' exchange (reversible)'});
        [~, BPindex] = ismember(lower(rxnNames),lower(model.rxnNames));
        if ~isempty(BPindex(inFlux+1))
            disp([compounds{i} ' exchange has been constrained to: ' num2str(bounds(i)) ' [mmol/gDw h]'])
            %Allow some flexibility 
            model = setParam(model,'ub',BPindex(inFlux+1),flexBounds(1)*abs(bounds(i)));
            if numel(flexBounds)>1
                model = setParam(model,'lb',BPindex(inFlux+1),flexBounds(2)*abs(bounds(i)));
            end
            model = setParam(model,'eq',BPindex(~inFlux+1),0); % Block exchange reaction in other direction
        else
            disp(['No exchange rxn for ' compounds{i} ' was found in ecModel'])
        end
    end
end
disp(' ')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function oxPhos = getOxPhosRxnIDs(model,parameters)
if isfield(parameters,'oxPhos')
    oxPhos = [];
    for i=1:length(parameters.oxPhos)
        ID = parameters.oxPhos{i};
        isoEnzymes = model.rxns(find(contains(model.rxns,ID)));
        isoEnzymes = isoEnzymes(~contains(isoEnzymes,'arm_'));
        oxPhos = [oxPhos; isoEnzymes];
    end
else
    oxPhos = [];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeProtCounts(condition,prot_parameters,name)
initial_prots  = length(prot_parameters{1});
filtered_prots = length(prot_parameters{2});
matched_prots  = length(prot_parameters{3});
model_prots    = length(prot_parameters{4});
mass_coverage  = prot_parameters{5};
disp(['The total number of proteins in the dataset is:                ' num2str(initial_prots)])
disp(['The total number of proteins in the filtered dataset is:       ' num2str(filtered_prots)])
disp(['The total number of filtered proteins present in the model is: ' num2str(matched_prots)])
disp(['The mass ratio between measured and unmeasured protein is:     ' num2str(mass_coverage)])

T = table(initial_prots,filtered_prots,matched_prots,model_prots,mass_coverage);
writetable(T,['../../models/prot_constrained/' name '/prot_counts_' condition '.txt'],'Delimiter','\t')
end
