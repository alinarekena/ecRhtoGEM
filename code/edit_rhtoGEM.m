%
% edit_rhtoGEM
%
%   This script prepares the original rhtoGEM (from Tiukova et al. 2019)
%   for the integration of enzymatic constraints. Sets the carbon source uptake,
%   adds D-arabinitol pathway and makes some other minor changes.
%
%   Last modified: 2020-09-28

code = pwd();

%% Load the model MATLAB structure
load('../models/rhto.mat')
disp(model)
printConstraints(model,-1000,1000);

%% Introduce a set of reactions for D-arabinitol
% Define reactions equations
t_0883 = 'D-xylulose[c] + NADH[c] + H+[c] <=> D-arabinitol[c] + NAD[c]';
r_4339 = 'D-arabinitol[c] <=> D-arabinitol[e]';
r_4340 = 'D-arabinitol[e] => ';
rxnsToAdd.equations = {t_0883; r_4339; r_4340}; 
% Define reaction names
t_0883 = 'D-arabinitol dehydrogenase';
r_4339 = 'D-arabinitol transport';
r_4340 = 'D-arabinitol exchange';
rxnsToAdd.rxnNames = {t_0883; r_4339; r_4340};
t_0883 = 't_0883';
r_4339 = 'r_4339';
r_4340 = 'r_4340';
rxnsToAdd.rxns = {t_0883; r_4339; r_4340};
% Define objective and bounds
rxnsToAdd.c  = [0 0 0];
rxnsToAdd.lb = [-1000 -1000 0];
rxnsToAdd.ub = [1000 1000 1000];

% Metabolites to Add
metsToAdd.mets          = {'s_D-arabinitol_c' 's_D-arabinitol_e'};
metsToAdd.metNames      = {'D-arabinitol' 'D-arabinitol'};
metsToAdd.compartments  = {'c' 'e'};

%genes to add
genesToAdd.genes          = {'RHTO_07844'};
genesToAdd.geneShortNames = {'RHTO_07844'};
rxnsToAdd.grRules         = {'RHTO_07844' '' ''};
% Introduce changes to the model
model = addGenesRaven(model,genesToAdd);
model = addMets(model,metsToAdd);
model = addRxns(model,rxnsToAdd,3);
%Standardize gene related fields
[grRules, rxnGeneMat] = standardizeGrRules(model,true);
model.grRules     = grRules;
model.rxnGeneMat  = rxnGeneMat;

%% Change gene associations for the following reactions
%String or cell array of reaction IDs
rxnID = 'r_0697';
%String of additional or replacement gene association
geneAssoc = 'RHTO_02522 or RHTO_06289';
% Introduce changes to the model
model = changeGeneAssoc(model,rxnID,geneAssoc);

rxnID = {'t_0240' 't_0241'};
geneAssoc = 'RHTO_03646';
model = changeGeneAssoc(model,rxnID,geneAssoc);

%% Save models

exportToExcelFormat(model,'../models/model_edit.xlsx')
exportModel(model,'../models/model_edit.xml')
save('../models/model_edit.mat','model')

clear geneAssoc genesToAdd grRules metsToAdd
clear r_4339 r_4340 rxnGeneMat rxnID rxnsToAdd t_0883

