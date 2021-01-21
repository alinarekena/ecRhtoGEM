%
% analyzeUsage
% 
%   Analyze enzyme usage assuming that models have been reconstructed by
%   generate_protModels and ribosomal subunits have been added.
%
%   Last modified: 2021-01-21
%

%% Load models

%load('models/ecModel_P_Xexp.mat')
%load('models/ecModel_P_XNlim.mat')
%load('models/ecModel_P_Aexp.mat')
%load('models/ecModel_P_ANlim.mat')
%load('models/ecModel_P_Gexp.mat')
%load('models/ecModel_P_GNlim.mat')

getpref('RAVEN')
setRavenSolver('cobra')

ecModels{1}=ecModelP_XNlim;
disp(['Now testing: ' flux.conds{1}])


%% Get enzymes usages to each reaction
sol = solveLP(ecModelP_XNlim,1);
[absUsage, capUsage, UB, protName]=enzymeUsage(ecModelP_XNlim,sol.x,true);


%% Prepare output
clear out
out(:,1)=ecModels{1}.enzymes;
out(:,2)=ecModels{1}.enzGenes;
out(:,3)=ecModels{1}.enzNames;
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

%% All usage per subSystem
head={'protID','geneID','protName','capUse_XNlim','capUse_XNim2','absUse_XNlim','UB_XNlim'};
out=cell2table(out,'VariableNames',head);
writetable(out,fullfile('results','model_simulation','enzymeUsages_XNlim.txt'),'Delimiter','\t')

%Remove the cloned repos:
rmdir('GECKO', 's')
