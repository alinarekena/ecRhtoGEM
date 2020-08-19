%
% ecRhtoGEM_generation
%
%   Generation of enzyme-constrained R.toruloides genome-scale metabolic
%   model (ecRhtoGEM). Firstly, enzyme-constrained model for growth on
%   xylose at 0.064 [h-1] with total protein pool is generated. Then
%   secondly, absolute proteomics data are incorporated, allowing the
%   modelled but not individually measured enzymes be drawn from the
%   enzymatic protein pool. Thus, generation of ecRhtoGEM employs two GECKO
%   pipelines: enhanceGEM and generate_protModels
%
%   Alina Rekena.   Last modified: 2020-08-19
%

%% Load model
model    = load('rhto_edit_v1.mat');
model    = model.model;
modelVer = model.description(strfind(model.description,'_v')+1:end);

%% Checklist for GECKO scripts
