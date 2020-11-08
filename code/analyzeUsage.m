%
% analyzeUsage
% 
%   Analyze enzyme usage assuming that models have been reconstructed by
%   generate_protModels and ribosomal subunits have been added.
%
%   Last modified: 2020-09-23
%

%% Load models

load('../models/ecModel_P_Xexp.mat')
load('../models/ecModel_P_XNlim.mat')

printConstraints(ecModelP_Xexp,-1000,1000); % check that model is constrained by experimental data

ecModels{1}=ecModelP_Xexp;
ecModels{2}=ecModelP_XNlim;

%% Get enzymes usages to each reaction
for i=1:2
    disp(['Now testing: ' flux.conds{i}])
    sol{i} = solveLP(ecModels{i});
    [absUsage{i}, capUsage{i}, UB{i}, protName{i}]=enzymeUsage(ecModels{i},sol{i}.x,true);
    printFluxes(ecModels{i},sol{i}.x,false,0,fullfile('..','results','modelSimulation',['allFluxes_' flux.conds{i},'.txt']),'%rxnID\t%rxnName\t%eqn\t%flux\n');
end

%% Prepare output
clear out
out(:,1)=ecModels{1}.enzymes;
out(:,2)=ecModels{1}.enzGenes;
out(:,3)=ecModels{1}.enzNames;
for i=1:2
    out(:,3+i)=strtrim(cellstr(num2str(capUsage{i},3)));
end
for i=1:2
    out(:,5+i)=strtrim(cellstr(num2str(absUsage{i},3)));
end
for i=1:2
    out(:,7+i)=strtrim(cellstr(num2str(UB{i},3)));
end

%% All usage per subSystem
head={'protID','geneID','protName','capUse_Xexp','capUse_XNlim','absUse_Xexp','absUse_XNlim','UB_Xexp','UB_XNlim'};
out=cell2table(out,'VariableNames',head);
writetable(out,fullfile('..','results','modelSimulation','enzymeUsages.txt'),'Delimiter','\t')

%Remove the cloned repos:
%rmdir('GECKO', 's')
