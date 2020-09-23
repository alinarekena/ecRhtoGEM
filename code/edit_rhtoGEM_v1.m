%
% edit_rhtoGEM
%
%   This script prepares the original rhtoGEM (from Tiukova et al. 2019)
%   for the integration of enzymatic constraints. Sets the carbon source,
%   adds D-arabinitol pathway and makes some other minor changes.
%
%   Last modified: 2020-09-22

code = pwd();

%% Load the model MATLAB structure
load('../models/rhto.mat')
disp(model)

%% Set glucose bound to zero and xylose as a carbon source
printConstraints(model,-1000,1000);
model = changeRxnBounds(model, {'r_1714'},0,'l');      %D-glucose exchange

model = changeRxnBounds(model, {'r_1718'},-1.74,'l');  %D-xylose exchange

%% If separate base model for xylose P3 is to be created, then change bounds
%model = changeRxnBounds(model, {'r_2104'},-0.039,'l'); %xylitol exchange
%model = changeRxnBounds(model, {'r_4340'},-0.142,'l');%D-arabinitol
%exchange
%printConstraints(model,-1000,1000);

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


%model = changeRxnBounds(model, {'r_4340'},-0.142,'l');%D-arabinitol
printConstraints(model,-1000,1000);

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

%% Save new models
exportToExcelFormat(model,'../models/rhto_edit_v1_XP1.xlsx')
exportModel(model,'../models/rhto_edit_v1_XP1.xml')
save('../models/rhto_edit_v1_XP1.mat','model')

clear geneAssoc genesToAdd grRules metsToAdd
clear r_4339 r_4340 rxnGeneMat rxnID rxnsToAdd t_0883

