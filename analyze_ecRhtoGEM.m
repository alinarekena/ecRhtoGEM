%
% analyze_ecRhtoGEM
%
%   run random sampling, calculate flux mean value and standard deviation;
%   generate ATP, NADPH, and NADH production balances;
%   convert fluxes to original, non-ecModel version.
%

% set 'gurobi' as RAVEN solver. With 'cobra' it gives warnings.
getpref('RAVEN')
setRavenSolver('gurobi')

% prepare for random sampling:
%
%
% set measured byproduct constraints: (can be found at fermParams.byP_flux)
ecModelP_Xexp = setParam(ecModelP_Xexp,'eq','r_2104',0.223);    %xylitol exchange
ecModelP_Xexp = setParam(ecModelP_Xexp,'eq','r_4340',0.367);    %D-arabinitol

ecModelP_XNlim = setParam(ecModelP_XNlim,'eq','r_2104',0.004);    %xylitol exchange
ecModelP_XNlim = setParam(ecModelP_XNlim,'eq','r_4340',0.077);    %D-arabinitol

ecModelP_Aexp = setParam(ecModelP_Aexp,'eq','r_1687',0.122);    %citrate

%ANlim: citrate uptake 0.043 is constrained by UB and predicted by FBA;

ecModelP_GexpUrea = setParam(ecModelP_GexpUrea,'eq','r_1808',0.049);%glycerol

ecModelP_GNlimUrea = setParam(ecModelP_GNlimUrea,'eq','r_1808_REV',0.007);%glycerol uptake

% run FBA: with byProducts, optimize for prot_exchange
solveLP(ecModelP_GNlimUrea,1);%same for all models

% create temporary model and constrain it for random sampling:
ecModelP_XexpTmp = setParam(ecModelP_Xexp,'ub','r_2111',0.05399); %UB=flux+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'ub','r_4046',3.6865); %UB=min+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1672',2.957,10);   %CO2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1992_REV',2.288,10); %O2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','prot_pool_exchange',0.7413,10); % FBA flux with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_2104',0.223,10);  % measured xylitol
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_4340',0.367,10);  % measured D-arabinitol

ecModelP_XNlimTmp2 = setParam(ecModelP_XNlim,'lb','r_2111',0.01473); % different from 0.01782 w/o byproducts 
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_2111',0.01488); %UB=0.01473+1%
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_4046',2.8785); %UB=min+1%;
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1672',1.097,10);   %CO2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1992_REV',0.9532,10); %O2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','prot_pool_exchange',0.4792,10);
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_2104',0.004,10);  % measured xylitol
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_4340',0.077,10);  % measured D-arabinitol

ecModelP_AexpTmp = setParam(ecModelP_Aexp,'ub','r_2111',0.07299);%UB=flux+1%
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'ub','r_4046',3.636);%UB=min+1%
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1672',6.735,10);   %CO2
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1992_REV',6.667,10); %O2
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','prot_pool_exchange',0.9534,10);%flux
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1687',0.122,10);% measured citrate

ecModelP_ANlimTmp = setParam(ecModelP_ANlim,'ub','r_2111',0.01199);%UB=flux+1%
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'ub','r_4046',3.484);%UB=min+1%
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1672',1.39,10);   %CO2
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1992_REV',2.012,10); %O2
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','prot_pool_exchange',0.3002,10);%flux
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1687_REV',0.043,20);% citrate uptake (as in GECKO pipeline)

ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUrea,'ub','r_2111',0.17998);%0.1782+1%
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'ub','r_4046',1.814);% max NGAM GexpUrea
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1672',6.94,10)   %CO2
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1992_REV',5.812,10) %O2
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','prot_pool_exchange',1.41,10)
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1808',0.049,10)%glycerol

ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUrea,'ub','r_2111',0.02099);%0.02079+1%
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'ub','r_4046',2.7775);%UB=min+1%
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1672',1.444,10)   %CO2
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1992_REV',1.16,10) %O2
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','prot_pool_exchange',1.179,10)
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1808_REV',0.007,10)%glycerol uptake

% run random sampling:
%
% Note, goodRxns don't make sense because of different C sources.
%
goodRxns = [];
rsmatrix=rs(ecModelP_XexpTmp,2000,true,true,true,goodRxns,true);%MinFlux=True
rsmatrix=rs(ecModelP_XNlimTmp2,2000,true,true,true,goodRxns,true);
rsmatrix=rs(ecModelP_AexpTmp,2000,true,true,true,goodRxns,true);
rsmatrix=rs(ecModelP_ANlimTmp,2000,true,true,true,goodRxns,true);
rsmatrix=rs(ecModelP_GexpUreaTmp,2000,true,true,true,goodRxns,true);
rsmatrix=rs(ecModelP_GNlimUreaTmp,2000,true,true,true,goodRxns,true);

% output: ec
clear out
sol.median=full(median(rsmatrix,2));  %2 means for each row
sol.mean=full(mean(rsmatrix,2));
sol.std=full(std(rsmatrix,0,2));
clear out

%% convert fluxes to original, non-ecModel version:

% rs matrix: non-ec
rsXmapped = mapRxnsToOriginal(ecModelP,model,rsmatrix)
% (repeat for each condition)

% output: non-ec
solXmapped.median=full(median(rsXmapped,2));  %2 means for each row
solXmapped.mean=full(mean(rsXmapped,2));
solXmapped.std=full(std(rsXmapped,0,2));
% (repeat for each condition)

%% calculate ATP, NADPH, and NADH balances from non-ec fluxes:

% adjust for each compound and repeat for each condition
for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model,i,solXmapped.median,true);
    clear out
    out.rxns    = model.rxns(rxnIdx);
    out.rxnNames= model.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model,rxnIdx);
    out.fluxes  = num2cell(fluxes);
    out = [out.rxns out.rxnNames out.rxnEqns out.fluxes];
    end

