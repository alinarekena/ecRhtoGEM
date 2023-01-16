%
% edit_rhtoGEM
%
%   Preparing original rhtoGEM (https://github.com/SysBioChalmers/rhto-GEM)
%   to integrate enzymatic constraints. This script includes curation of
%   xylose assimilation pathway and setting carbon uptake for xylose.
%
%   Last modified: 2021-02-11
%

% Prepare COBRA and set repo root path
initCobraToolbox
root = pwd; % Get the root directory of the folder

%% Load model
model    = importModel(fullfile(root,'models','rhto.xml'));
model = changeRxnBounds(model, {'r_1714'},0,'l');      %D-glucose exchange
model = changeRxnBounds(model, {'r_1718'},-1.86,'l');  %D-xylose exchange

%% Remove reactions
model=removeReactions(model,'t_0082');

%% Introduce reactions

% Define reaction equations
t_0883 = 'D-xylulose[c] + NADH[c] + H+[c] <=> D-arabinitol[c] + NAD[c]';% D-arabinitol 4-dehydrogenase (DAD-4)
r_4339 = 'D-arabinitol[c] <=> D-arabinitol[e]';
r_4340 = 'D-arabinitol[e] => ';

t_0884 = 'D-ribulose[c] + NADPH[c] + H+[c] <=> D-arabinitol[c] + NADP(+)[c]';% D-arabinitol 2-dehydrogenase (DAD-2)/D-ribulose reductase (RiR)
% RiR as analogue to L-xylulose reductase (LXR):L-xylulose[c] + H+[c] + NADPH[c] <=> NADP(+)[c] + xylitol[c]
% https://doi.org/10.1074/jbc.M312533200
t_0885 = 'ATP[c] + D-ribulose[c] => ADP[c] + D-ribulose 5-phosphate[c] + H+[c]';% D-ribulokinase (RiK)
% analogue to ATP:D-xylulose 5-phosphotransferase from Kluyveromyces
% marxianus GEM (https://github.com/SysBioChalmers/Kluyveromyces_marxianus-GEM) and Pichia stipitis GEM
% (http://biomet-toolbox.chalmers.se/index.php?page=models-Stipitis)
t_0886 = 'ATP[c] + acetate[c] <=> ADP[c] + acetyl-phosphate[c]';% Acetate kinase
% ACKr (13383) from Rt_IFO0880 and iRhtoC

rxnsToAdd.equations = {t_0883; r_4339; r_4340; t_0884; t_0885; t_0886}; 

% Define reaction names
t_0883 = 'D-arabinitol 4-dehydrogenase';
r_4339 = 'D-arabinitol transport';
r_4340 = 'D-arabinitol exchange';
t_0884 = 'D-arabinitol 2-dehydrogenase/D-ribulose reductase';
t_0885 = 'D-ribulokinase';
t_0886 = 'Acetate kinase';
rxnsToAdd.rxnNames = {t_0883; r_4339; r_4340; t_0884; t_0885; t_0886};
t_0883 = 't_0883';
r_4339 = 'r_4339';
r_4340 = 'r_4340';
t_0884 = 't_0884';
t_0885 = 't_0885';
t_0886 = 't_0886';
rxnsToAdd.rxns = {t_0883; r_4339; r_4340; t_0884; t_0885; t_0886};

% Define objective and bounds
rxnsToAdd.c  = [0 0 0 0 0 0];
rxnsToAdd.lb = [-1000 -1000 0 -1000 0 -1000];
rxnsToAdd.ub = [1000 1000 1000 1000 1000 1000];

% Define EC numbers
t_0883 = '1.1.1.11;1.1.1.138';
r_4339 = '';
r_4340 = '';
t_0884 = '1.1.1.10;1.1.1.138';
t_0885 = '2.7.1.16;2.7.1.47';
t_0886 = '2.7.2.1';
rxnsToAdd.eccodes = {t_0883; r_4339; r_4340; t_0884; t_0885; t_0886};

% Add metabolites
metsToAdd.mets          = {'s_D-arabinitol_c' 's_D-arabinitol_e' 's_D-ribulose'};
metsToAdd.metNames      = {'D-arabinitol' 'D-arabinitol' 'D-ribulose'};
metsToAdd.compartments  = {'c' 'e' 'c'};

% Add genes
genesToAdd.genes          = {'RHTO_07844' 'RHTO_00950' 'RHTO_04461'};
genesToAdd.geneShortNames = {'RHTO_07844' 'RHTO_00950' 'RHTO_04461'};
rxnsToAdd.grRules         = {'RHTO_07844' '' '' 'RHTO_00373' 'RHTO_00950' 'RHTO_04461'};

%% Introduce changes to the model

model = addMets(model,metsToAdd);
model = addGenesRaven(model,genesToAdd);
model = addRxns(model,rxnsToAdd,3);% add reactions in last position to avoid allowNewMets error

%Standardize gene related fields
[grRules, rxnGeneMat] = standardizeGrRules(model,true);
model.grRules     = grRules;
model.rxnGeneMat  = rxnGeneMat;

%% Change gene associations for the following reactions

rxnID = 'r_0697'; % lactoylglutathione lyase (EC 4.4.1.5)
geneAssoc = 'RHTO_02522 or RHTO_06289';
model = changeGeneAssoc(model,rxnID,geneAssoc);

rxnID = {'t_0240' 't_0241'}; % glycerol-3-phosphate acyltransferases (16:0 and 18:0)
geneAssoc = 'RHTO_03646';
model = changeGeneAssoc(model,rxnID,geneAssoc);

%% Save models

exportToExcelFormat(model,'models/model_edit.xlsx')
exportModel(model,'models/rhto_edit.xml')
save('models/rhto_edit.mat','model')

clear geneAssoc genesToAdd grRules metsToAdd
clear r_4339 r_4340 rxnGeneMat rxnID rxnsToAdd t_0883 t_0884 t_0885 t_0886

