function model = changeMedia_batch(model,c_source)
%changeMedia_batch
%
% Function that modifies the ecModel and makes it suitable for batch growth
% simulations on the carbon source of choice. You can add more changes in
% this function if you need to do so for your model.
%
% INPUT:
%   model       An enzyme constrained model.
%	c_source	The name of the exchange reaction that supplies the model
%				with carbon.
%
% OUTPUT:
%   model       The enzyme constrained model with modified boundaries.
%
% Usage: model = changeMedia_batch(model,c_source)
%
% Benjamin J. Sanchez	2018-12-11
% Ivan Domenzain        2020-02-07

%first block any uptake
[rxnIDs,exchange]  = getExchangeRxns(model);
exchange           = exchange(find(contains(rxnIDs,'_REV')));
model.ub(exchange) = 0;
%Allow main carbon source uptake
c_id  = model.rxns(strcmp(model.rxnNames,c_source));
model = setParam(model,'ub',c_id,Inf);
%block xylose and oxygen production
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-xylose exchange')) = 0;
%Allow uptake of essential components
model = setParam(model, 'ub', 'r_1654_REV', Inf); % 'ammonium exchange';
model = setParam(model, 'ub', 'r_2100_REV', Inf); % 'water exchange' ;
model = setParam(model, 'ub', 'r_1992_REV', Inf); % 'oxygen exchange';
model = setParam(model, 'ub', 'r_2049_REV', Inf); % sodium exchange
model = setParam(model, 'ub', 'r_1861_REV', Inf); % 'iron(2+) exchange';
model = setParam(model, 'ub', 'r_2005_REV', Inf); % 'phosphate exchange';
model = setParam(model, 'ub', 'r_2060_REV', Inf); % 'sulphate exchange';
model = setParam(model, 'ub', 'r_1832_REV', Inf); % 'H+ exchange' ;
model = setParam(model, 'ub', 'r_2020_REV', Inf); % potassium exchange
%Allow biomass production 
model = setParam(model, 'ub', 'r_2111', Inf); % growth
end