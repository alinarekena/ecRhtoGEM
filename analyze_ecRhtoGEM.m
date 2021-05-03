%
% analyze_ecRhtoGEM
%
%   run random sampling, calculate flux mean value and standard deviation;
%   calculate enzyme usage; generate ATP, NADPH, and NADH production
%   balances; calculate energy turnover per compartment; convert fluxes to
%   original, non-ecModel version.
%
%

% better set 'gurobi' as RAVEN solver. With 'cobra' it gives warnings.
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

% run FBA:
solveLP(ecModelP_Xexp,1);%same for all models

% create temporary model and constrain it for random sampling:
ecModelP_XexpTmp = setParam(ecModelP_Xexp,'ub','r_2111',0.05399); %UB=flux+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'ub','r_4046',3.6865); %UB=min+1%
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1672',2.932,10);   %CO2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_1992_REV',2.275,10); %O2 flux from FBA with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','prot_pool_exchange',1.97,10); % FBA flux with byProducts
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_2104',0.223,10);  % measured xylitol
ecModelP_XexpTmp = setParam(ecModelP_XexpTmp,'var','r_4340',0.367,10);  % measured D-arabinitol

ecModelP_XNlimTmp2 = setParam(ecModelP_XNlim,'lb','r_2111',0.01473); % different from 0.01782 w/o byproducts 
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_2111',0.01488); %UB=0.01473+1%
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'ub','r_4046',3.131); %UB=min+1%;
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1672',1.124,10);   %CO2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_1992_REV',0.9867,10); %O2
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','prot_pool_exchange',0.5185,10);
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_2104',0.004,10);  % measured xylitol
ecModelP_XNlimTmp2 = setParam(ecModelP_XNlimTmp2,'var','r_4340',0.077,10);  % measured D-arabinitol

ecModelP_AexpTmp = setParam(ecModelP_Aexp,'ub','r_2111',0.07299);%UB=flux+1%
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'ub','r_4046',3.686);%UB=min+1%
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1672',7.334,10);   %CO2
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1992_REV',7.27,10); %O2
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','prot_pool_exchange',3.249,10);%flux
ecModelP_AexpTmp = setParam(ecModelP_AexpTmp,'var','r_1687',0.122,10);% measured citrate

ecModelP_ANlimTmp = setParam(ecModelP_ANlim,'ub','r_2111',0.01199);%UB=flux+1%
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'ub','r_4046',3.484);%UB=min+1%
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1672',1.604,10);   %CO2
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1992_REV',2.096,10); %O2
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','prot_pool_exchange',0.64,10);%flux
ecModelP_ANlimTmp = setParam(ecModelP_ANlimTmp,'var','r_1687_REV',0.043,20);% citrate uptake (as in GECKO pipeline)

ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUrea,'ub','r_2111',0.19099);%0.1891+1%
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'ub','r_4046',1.814);% max NGAM GexpUrea
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1672',7.387,10)   %CO2
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1992_REV',6.101,10) %O2
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','prot_pool_exchange',4.104,10)
ecModelP_GexpUreaTmp = setParam(ecModelP_GexpUreaTmp,'var','r_1808',0.049,10)%glycerol

ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUrea,'ub','r_2111',0.018998);%0.01881+1%
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'ub','r_4046',3.333);%UB=min+1%
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1672',1.421,10)   %CO2
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1992_REV',1.179,10) %O2
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','prot_pool_exchange',2.028,10)
ecModelP_GNlimUreaTmp = setParam(ecModelP_GNlimUreaTmp,'var','r_1808_REV',0.007,10)%glycerol uptake

% run random sampling:
%
% Note, goodRxns don't make sense because of different C sources.
%
goodRxns = [];
sol=rs(ecModelP_XexpTmp,2000,true,true,true,goodRxns,true);%MinFlux=True
sol=rs(ecModelP_XNlimTmp2,2000,true,true,true,goodRxns,true);
sol=rs(ecModelP_AexpTmp,2000,true,true,true,goodRxns,true);
sol=rs(ecModelP_ANlimTmp,2000,true,true,true,goodRxns,true);
sol=rs(ecModelP_GexpUreaTmp,2000,true,true,true,goodRxns,true);
sol=rs(ecModelP_GNlimUreaTmp,2000,true,true,true,goodRxns,true);

% output:
clear out
out.mean=full(mean(sol,2)); %2 means for each row
out.std=full(std(sol,0,2));
clear out

%% enzyme usage:

[absUsage,capUsage,UB,protID] = enzymeUsage(ecModelP_XexpTmp,out_Xexp.mean,true,true);
[absUsage,capUsage,UB,protID] = enzymeUsage(ecModelP_XNlimTmp2,out_XNlim.mean,true,true);
[absUsage,capUsage,UB,protID] = enzymeUsage(ecModelP_AexpTmp,out_Aexp.mean,true,true);
[absUsage,capUsage,UB,protID] = enzymeUsage(ecModelP_ANlimTmp,out_ANlim.mean,true,true);
[absUsage,capUsage,UB,protID] = enzymeUsage(ecModelP_GexpUreaTmp,out_GexpUrea.mean,true,true);
[absUsage,capUsage,UB,protID] = enzymeUsage(ecModelP_GNlimUreaTmp,out_GNlimUrea.mean,true,true);

% prepare output:
% repeat for each condition
clear out
out(:,1)=ecModelP_Xexp.enzymes;
out(:,2)=ecModelP_Xexp.enzGenes;
out(:,3)=ecModelP_Xexp.enzNames;
for i=1:1
    out(:,3+i)=strtrim(cellstr(num2str(capUsage,3)));
end
for i=1:1
    out(:,4+i)=strtrim(cellstr(num2str(capUsage,3)));
end
for i=1:1
    out(:,5+i)=strtrim(cellstr(num2str(absUsage,3)));
end
for i=1:1
    out(:,6+i)=strtrim(cellstr(num2str(UB,3)));
end

head={'protID','geneID','protName','capUse_Xexp','capUse_Xexp','absUse_Xexp','UB_Xexp'};
out=cell2table(out,'VariableNames',head);
writetable(out,fullfile('results','model_simulation','enzymeUsages_Xexp.txt'),'Delimiter','\t')

%% calculate ATP, NADPH, and NADH balances:

for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(ecModelP_GNlimUrea,i,out_GNlimUrea.mean,true);
    clear out
    out.rxns    = ecModelP_GNlimUrea.rxns(rxnIdx);
    out.rxnNames= ecModelP_GNlimUrea.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(ecModelP_GNlimUrea,rxnIdx);
    out.fluxes  = num2cell(fluxes);
    out = [out.rxns out.rxnNames out.rxnEqns out.fluxes];
    %fid = fopen([root '/results/model_simulation/rs_' i{1} '_productionFluxes.tsv'],'w');
    %fprintf(fid,['%s' repmat('\t%s',1,8) '\n'],["rxnID" "rxnName" ...
    %    "rxnEqn" string(conditions.abbrev')]);
    %for j=1:length(out)
    %    fprintf(fid,['%s\t%s\t%s' repmat('\t%f',1,6) '\n'],out{j,:});
    end
    fclose(fid);

%% calculate energy turnover per compartment:

% Define flux indexes:
fluxRxns={'r_0226','r_1022No1','r_0892No1','r_0962No1','r_0886No1',...
    'arm_r_0534','r_4046','r_4041','r_1718_REV','r_0961No1','arm_r_0658',...
    'r_0714No1','r_0714_REVNo1','r_0713No1','r_0713_REVNo1','arm_r_0770'};
idxs=getIndexes(ecModelP_Xexp,fluxRxns,'rxns');
rXexp=out_Xexp.mean(idxs(9),:);
fluxes(1,:)=rXexp;

% calculate flux rates related to NAD turnover, per compartment:
% repeat for each condition
NADmets=find(strcmp(ecModelP_GNlimUrea.metNames,'NAD'));
[~,NADrxns]=find(ecModelP_GNlimUrea.S(NADmets,:)); clear NAD
for i=1:numel(ecModelP_GNlimUrea)
    NADtmp=full(ecModelP_GNlimUrea.S(NADmets,NADrxns).*out_GNlimUrea.mean(NADrxns,i)');
    NADtmp(NADtmp<0)=0;
    NAD(:,i)=sum(NADtmp,2);
end
zeroFlux=sum(NAD,2)==0;
NAD=[NAD(~zeroFlux,:);sum(NAD,1)];
NADmets=[strcat('rNAD[',ecModelP_GNlimUrea.comps(ecModelP_GNlimUrea.metComps(NADmets(~zeroFlux))),']'); ...
    'rNAD[tot]']';
fluxes=[fluxes;NAD];

% NADPH:
NADPmets=find(strcmp(ecModelP_GNlimUrea.metNames,'NADP(+)'));
[~,NADPrxns]=find(ecModelP_GNlimUrea.S(NADPmets,:)); clear NADP
for i=1:numel(ecModelP_GNlimUrea)
    NADPtmp=full(ecModelP_GNlimUrea.S(NADPmets,NADPrxns).*out_GNlimUrea.mean(NADPrxns,i)');
    NADPtmp(NADPtmp<0)=0;
    NADP(:,i)=sum(NADPtmp,2);
end
zeroFlux=sum(NADP,2)==0;
NADP=[NADP(~zeroFlux,:);sum(NADP,1)];
NADPmets=[strcat('rNADP[',ecModelP_GNlimUrea.comps(ecModelP_GNlimUrea.metComps(NADPmets(~zeroFlux))),']'); ...
    'rNADP[tot]']';
fluxes=[fluxes;NADP];

% Define row and variable names:
rowNames=[NADmets, NADPmets];
%varNames={'Xexp','XNlim','Aexp','ANlim','GexpUrea','GNlimUrea'};

%% convert fluxes to original, non-ecModel version:

% repeat for each condition
solXmapped = mapRxnsToOriginal(ecModelP_Xexp,model,out_Xexp.mean)
