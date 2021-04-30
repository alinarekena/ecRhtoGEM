function [absUsage,capUsage,UB,protId] = enzymeUsage(ecModel,fluxes,zero,contrastRevRxns)
% enzymeUsage
%
%   Gives enzyme usages based on a provided flux distribution. It can give:
%   1)  absolute usage: the specific enzyme usage in mmol/gDCW/h, which can
%       be given for enzymes both with- and without abundance information;
%   2)  capacity usage: the ratio of available enzyme that is used, calcuted
%       by (absUsage/UB) (note that capacity usage is 0 if an enzyme
%       abundance was not constrained in the model);
%   3)  UB: the upper bound of each enzyme exchange reaction;
%   4)  protId: the protein identifiers for each enzyme (if the model has an
%       enzymes field than this order is used, otherwise it is given
%       alphabetically.
%
%   Input:
%   ecModel         enzyme-constrained model
%   fluxes          vector of fluxes, for instance sol.x
%   zero            logical whether also zero absolute enzyme usages should be
%                   included (opt, default true)
%   contrastRevRxns logical whether usage of enzymes involved in reversible
%                   reactions should be subtracted from either direction,
%                   to yield "net" usage. If set to false, usage is defined
%                   from the enzyme exchange reactions instead (opt,
%                   default false)
%
%   Output:
%   capUsage        vector of enzyme capacity usages
%   absUsage        vector of absolute enzyme usages
%   UB              vector of enzyme exchange reaction upper bounds
%   protId          string array of protein IDs matching the other output
%
% Usage: [absUsage,capUsage,UB,protId] = enzymeUsage(ecModel,fluxes,zero,contrastRevRxns)

if nargin<3
    zero=true;
end
if nargin<4
    contrastRevRxns=false;
end

%Get enzyme exchange reactions (exclude protein pool)
expression  = '^(draw_prot_.*)|(prot_[^(pool)].*_exchange$)';
enzExch     = regexp(ecModel.rxnNames,expression);
enzExch     = find(~cellfun(@isempty,enzExch));

%Get list of Uniprot protein IDs
if ~isfield(ecModel,'enzymes')
    expression  = '(^draw_prot_)|(^prot_)|(_exchange$)';
    protId      = regexprep(ecModel.rxnNames(enzExch),expression,'');
else
    protId      = ecModel.enzymes;
end

%Map reversible reactions to their original flux value
if contrastRevRxns == true
    revRxns     = find(contains(ecModel.rxns,'_REV'));
    revRxnsOri  = regexprep(ecModel.rxns(revRxns),'_REV','');
    [revMatch,oriRxns] = ismember(revRxnsOri,ecModel.rxns);
    oriRxns     = oriRxns(revMatch);
    revRxns     = revRxns(revMatch);
    fluxes(oriRxns) = fluxes(oriRxns) - fluxes(revRxns);
else
    revRxns = false(numel(ecModel.rxns),1);
end

%Index of enzymes in S-matrix
enzMets     = strcat('prot_',protId);
[~,enzMets] = ismember(enzMets,ecModel.mets);
%Index of reactions that use enzymes
[~,enzRxns,~]   = find(ecModel.S(enzMets,:));
enzRxns   = unique(enzRxns);
%Leave out exchange and reversed reactions
enzRxns   = setdiff(enzRxns,enzExch,'sorted');
enzRxns   = setdiff(enzRxns,revRxns,'sorted');

%Calculate absolute usage
enzFlux = -full(ecModel.S(enzMets,enzRxns)) * fluxes(enzRxns);
%Take absolute value, as net negative flux through reversible reaction
%should still yield a positive enzyme usage
absUsage = abs(enzFlux);

%For capacity usage
UB = ecModel.ub(enzExch);
capUsage = absUsage ./ UB;

if true(~zero)
    nonzero     = absUsage>0;
    absUsage    = absUsage(nonzero);
    capUsage    = capUsage(nonzero);
    UB          = UB(nonzero);
    protId      = protId(nonzero);
end
end