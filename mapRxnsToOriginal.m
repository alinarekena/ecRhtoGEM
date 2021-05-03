function solXmapped = mapRxnsToOriginal(ecModel,model,solX)
% mapRxnsToOriginal
%
%   Maps fluxes as obtained from an ecModel to the reactions as they were
%   defined in the non-ecModel that forms the basis for the ecModel. If
%   used for further summarizing of the data, use model to provide the
%   model structure.
%
%   Input:
%   ecModel     model structure of enzyme-constrained model that was
%               simulated
%   model       model structure of the classical, non-ecModel that was used
%               to construct the ecModel
%   solX        vector with flux solution from simulation with ecModel
%
%   Output:
%   solXmapped  vector with flux solution after mapping to the original
%               model
%
%   Example:
%   sol = solveLP(ecModel);
%   solX = sol.x;
%   solXmapped = mapsRxnsToOriginal(ecModel,model.solX);
%   printFluxes(model,solXmapped);
%
%   Usage: solXmapped = mapRxnsToOriginal(ecModel,model,solX)

rxns = ecModel.rxns;
%Unchanged reactions
oriRxns = cellfun(@isempty,regexp(rxns,'(^arm_)|(^(draw_)?prot)|(No\d{1,2}$)'));

%Remove protein exchange reactions
protEx = ~cellfun(@isempty,regexp(rxns,'^(draw_)?prot_'));
rxns(protEx) = [];
solX(protEx) = [];

%Keep arm reactions, remove arm_ prefix
armRxns         = find(~cellfun(@isempty,regexp(rxns,'^arm_.+')));
armRxnIds       = regexprep(rxns(armRxns),'^arm_','');
rxns(armRxns)   = armRxnIds;

%Remove arm-extension reactions
rmArmRxns = contains(rxns,armRxnIds);
rxns(rmArmRxns) = [];
solX(rmArmRxns) = [];

%Merge rev reactions
revRxns     = find(contains(rxns,'_REV'));
revRxnsOri  = regexprep(rxns(revRxns),'_REV','');
[revMatch,fwdRxns] = ismember(revRxnsOri,rxns);
fwdRxns     = fwdRxns(revMatch);
revRxns     = revRxns(revMatch);
solX(fwdRxns) = solX(fwdRxns) - solX(revRxns);
solX(revRxns) = [];
rxns(revRxns) = [];

%Remove No1 suffix from single-enzyme reactions
rxns = regexprep(rxns,'No1$','');

%Match to original model
solXmapped = zeros(numel(model.rxns),1);
[a,b] = ismember(rxns,model.rxns);
solXmapped(b(a)) = solX(a);
