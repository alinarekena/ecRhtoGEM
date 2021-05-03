function [fluxes, rxnIdx] = getMetProduction(model,metName,solVector,nonZero)
% Which fluxes produce/consume a metabolite (excluding transport
% reactions).
% nonZero   only include non zero fluxes, default false

metLoc      = find(ismember(model.metNames,metName));
subS        = full(model.S(metLoc,:));
[~,rxnIdx]  = find(subS);
coeff       = sum(subS(:,rxnIdx),1);

fluxes = solVector(rxnIdx,:) .* transpose(coeff);

if nonZero
    nonZero     = sum(fluxes,2)~=0;
    fluxes      = fluxes(nonZero,:);
    rxnIdx      = rxnIdx(nonZero);
end
end
