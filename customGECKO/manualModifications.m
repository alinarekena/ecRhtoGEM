function [model,modifications] = manualModifications(model)
%
% Ivan Domenzain.      Last edited: 2020-02-07

modifications{1} = cell(0,1);
modifications{2} = cell(0,1);

for i = 1:length(model.rxns)
    reaction = model.rxnNames{i};
    %Find set of proteins present in rxn:
    S        = full(model.S);
    subs_pos = find(S(:,i) < 0);
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos);
    prot_set = cell(size(int_pos));
    MW_set   = 0;
    for j = 1:length(int_pos)
        met_name    = model.mets{int_pos(j)};
        prot_set{j} = met_name(6:end);
        MW_set      = MW_set + model.MWs(strcmp(model.enzymes,prot_set{j}));
    end
    %Update int_pos:
    S        = full(model.S);
    subs_pos = find(S(:,i) < 0);
    %Get the proteins that are part of the i-th rxn
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%  Individual Changes:  %%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(int_pos)
        enzName = model.mets(int_pos(j));
        %%%%%%%%%%%%%%%%%% MANUAL CURATION FOR TOP GROWTH LIMITING ENZYMES:
        [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications);
        if ~isempty(newValue)
            model.S(int_pos(j),i) = newValue;
        else
            %%%%%%%%%%%%%%% MANUAL CURATION FOR TOP USED ENZYMES:
            [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW_set,modifications);
            if ~isempty(newValue)
                model.S(int_pos(j),i) = newValue;
            end
        end
    end
    if rem(i,100) == 0 || i == length(model.rxns)
        disp(['Improving model with curated data: Ready with rxn ' num2str(i)])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%% Other manual changes: %%%%%%%%%%%%%%%%%%%%%%%%%%%
%model = otherChanges(model);
% Remove repeated reactions (2017-01-16):
rem_rxn = false(size(model.rxns));
for i = 1:length(model.rxns)-1
    for j = i+1:length(model.rxns)
        if isequal(model.S(:,i),model.S(:,j)) && model.lb(i) == model.lb(j) && ...
                model.ub(i) == model.ub(j)
            rem_rxn(j) = true;
            disp(['Removing repeated rxn: ' model.rxns{i} ' & ' model.rxns{j}])
        end
    end
end
model = removeRxns(model,model.rxns(rem_rxn));
% Merge arm reactions to reactions with only one isozyme (2017-01-17):
arm_pos = zeros(size(model.rxns));
p       = 0;
for i = 1:length(model.rxns)
    rxn_id = model.rxns{i};
    if contains(rxn_id,'arm_')
        rxn_code  = rxn_id(5:end);
        k         = 0;
        for j = 1:length(model.rxns)
            if contains(model.rxns{j},[rxn_code 'No'])
                k   = k + 1;
                pos = j;
            end
        end
        if k == 1
            %Condense both reactions in one:
            new_id     = model.rxns{pos};
            new_name   = model.rxnNames{pos};
            stoich     = model.S(:,i) + model.S(:,pos);
            model      = addReaction(model,{new_id,new_name},model.mets,stoich,true,0,1000);
            p          = p + 1;
            arm_pos(p) = i;
            disp(['Merging reactions: ' model.rxns{i} ' & ' model.rxns{pos}])
        end
    end
end
% Remove saved arm reactions:
model = removeRxns(model,model.rxns(arm_pos(1:p)));
%Change gene rules:
if isfield(model,'rules')
    model.rules = cell(length(model.grRules),1);
end
% Remove unused enzymes after manual curation (2017-01-16):
rem_enz = false(size(model.enzymes));
for i = 1:length(model.enzymes)
    pos_met = strcmp(model.mets,['prot_' model.enzymes{i}]);
    if sum(model.S(pos_met,:)~=0) == 1
        rem_enz(i) = true;
    end
end
rem_enz = model.enzymes(rem_enz);
for i = 1:length(rem_enz)
    model = deleteProtein(model,rem_enz{i});
    disp(['Removing unused protein: ' rem_enz{i}])
end
% Block O2 and glucose production (for avoiding multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-xylose exchange')) = 0;
% Map the index of the modified Kcat values to the new model (after rxns
% removals).
modifications = mapModifiedRxns(modifications,model);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the top growth limiting enzymes that were detected by the
% modifyKcats.m script in a preliminary run.
function [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications)
        newValue = [];
        reaction = string(reaction);
        % It was glucan synthase [M7WFB4, EC2.4.1.34], which is curated in
        % topUsedEnz function below.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify those kcats involved in extreme misspredictions for growth on
% several carbon sources. This values were obtained by specific searches on
% the involved pathways for the identification of the ec numbers and then
% its associated Kcat values were gotten from BRENDA.
function [newValue,modifications] = curation_carbonSources(reaction,enzName,MW_set,modifications)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After the growth limiting Kcats analysis and the curation for several
% carbon sources, a simulation for the model growing on minimal xylose
% media yielded a list of the top used enzymes (mass-wise), those that were
% taking more than 10% of the total proteome are chosen for manual curation
function [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW_set,modifications)
        newValue = [];
        reaction = string(reaction);
        % 1. [M7WT33//EC2.4.1.109] dolichyl-phosphate-mannose--protein
        % mannosyltransferase (No2):
        % prev_Kcat:0.0053666 new_Kcat:0.0053667 CC:6.6589 Err:-81.3496%
        % No Kcats available in BRENDA; S.A. values for EC2.4.1.109 were
        % very low;
        % New Kcat chosen according to a wild card:
        % EC2.4.1.-/any/'substrate';
        % Required for Xexp condition;
          if contains(reaction,'dolichyl-phosphate-mannose--protein mannosyltransferase')
              if strcmpi('prot_M7WT33',enzName)
                 % S.cerevisiae EC2.4.1.83 Kcat [2020-04-28];
                 % https://www.brenda-enzymes.org/literature.php?e=2.4.1.83&r=659501;
                 newValue     = -(70.9*3600)^-1;
                 modifications{1} = [modifications{1}; 'M7WT33'];
                 modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7WZY6',enzName)
                  newValue      = -(70.9*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7WZY6'];
                  modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7XBY2',enzName)
                  newValue      = -(70.9*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7XBY2'];
                  modifications{2} = [modifications{2}; reaction];
              end
          end
        % 2. [M7XI04/EC1.1.1.34] HMG-CoA reductase/ 3-hydroxy-3-methylglutaryl
        % coenzyme A reductase:
        % 1.7122e-01 of total protein pool at round A of ecRhtoGEM generation;
        % Higher than prev_Kcat not available in BRENDA (Rattus
        % Norvegicus).
        % New value of S.A. chosen according to: EC number/any organism;
        % Required for Xexp condition;
          if (strcmpi('prot_M7XI04',enzName) && (contains(reaction,'hydroxymethylglutaryl CoA reductase')))
             %Homarus americanus S.A. [2020-05-14];
             %https://www.brenda-enzymes.org/literature.php?e=1.1.1.34&r=654546;
             %newValue         = -(60*140*60)^-1;            %140 [1/s]
             % Streptomyces sp. S.A. [2020-04-28];
             % https://www.brenda-enzymes.org/literature.php?e=1.1.1.34&r=286581;
             newValue         = -(1.6*100*60)^-1;          %2.67
             modifications{1} = [modifications{1}; 'M7XI04'];
             modifications{2} = [modifications{2}; reaction];
          end
         % 3. [M7XKV7//EC#:6.3.5.3] 5'-phosphoribosylformyl glycinamidine
         % synthetase:
         % 8.4551e-02 of total protein pool at round A of ecRhtoGEM generation;
         % Higher than prev_Kcat not available in BRENDA (E.coli).
         % New value of S.A. chosen according to: EC number/any organism;
         % Required for Xexp condition;
          if (strcmpi('prot_M7XKV7',enzName)) && (contains(reaction,'phosphoribosylformyl'))
             % E.coli S.A. [2020-04-28];
             % https://www.brenda-enzymes.org/literature.php?e=6.3.5.3&r=1670;
             newValue         = -(2.15*141.418*60)^-1;      %5.0675
             modifications{1} = [modifications{1}; 'M7XKV7'];
             modifications{2} = [modifications{2}; reaction];
          end
        % 4. [M7WFB4//EC#:2.4.1.34] 1,3-beta-glucan synthase (No1):
        % 7.0144e-02 of total protein pool at round A of ecRhtoGEM generation;
        % No higher than prev_Kcat available in BRENDA (S.aureus).
        % New value of S.A. chosen according to: EC number/any organism;
        % Required for Xexp condition;
          if (strcmpi('prot_M7WFB4',enzName)) && (contains(reaction,'1,3-beta-glucan synthase'))
             % S.cerevisiae S.A. [2020-04-28];
             % https://www.brenda-enzymes.org/literature.php?e=2.4.1.34&r=489076;
             newValue         = -(4*200*60)^-1;             %13.3:
             %isn't enough for X20 model;
             % 2.0556e-02 of total protein pool; prev_Kcat 13.3;
             % S.A. of 27.9 reported for P.sativum, so Kcat increased 7 times:
             % https://www.brenda-enzymes.org/literature.php?e=2.4.1.34&r=489091;
             %newValue         = -(7*4*200*60)^-1;             %93.1
             modifications{1} = [modifications{1}; 'M7WFB4'];
             modifications{2} = [modifications{2}; reaction];
          end
          % [M7X0R7/EC2.2.1.2] transaldolase (No1):
          % 2.4464e-02 of total protein pool (ecModel version C);
          % One higher than prev_Kcat value available in BRENDA (Thermotoga
          % maritima) with the value of 22.3 [1/s];
          % But new value of S.A. chosen according to: EC number/any organism;
          % and the highest S.A. of the list was chosen;
          %if strcmpi('prot_M7X0R7',enzName) && contains(reaction,'transaldolase')
              % E.coli S.A. [2020-04-28];
              % https://www.brenda-enzymes.org/literature.php?e=2.2.1.2&r=486041;
              %newValue      = -(60*75*60)^-1;               %75
              %modifications{1} = [modifications{1}; 'M7X0R7'];
              %modifications{2} = [modifications{2}; reaction];
          %end
        % 5. [M7XI95//EC#:1.14.19.1] Delta-9 fatty acid desaturase ER membrane
        % (No1); stearoyl- or oleoyl- or palmitoyl-CoA desaturase
        % (n-C18:0CoA -> n-C18:1CoA),(n-C18:1CoA -n-C18:2CoA),(n-C16:0CoA
        % -> n-C16:1CoA);
        % 2.6441e-02  of total protein pool (version B of ecModel);
        % Kcat values available in BRENDA were for Rattus norvegicus and
        % spinach. Further artificially increased. 10 times.
        % Required for Xexp condition;
          if (strcmpi('prot_M7XI95',enzName)) && (contains(reaction,'-CoA desaturase'))
             % Spinacia oleracea [2020-05-01]
             % https://www.brenda-enzymes.org/literature.php?e=1.14.19.1&r=437670;
             %newValue         = -(0.5*3600)^-1;            %prev_Kcat
             newValue         = -(10*0.5*3600)^-1;
             modifications{1} = [modifications{1}; 'M7XI95'];
             modifications{2} = [modifications{2}; reaction];
          end
          % 6. [M7WPW0/EC2.1.1.41] sterol 24-C-methyltransferase:
          % 2.3058e-02 of total protein pool (version B of ecModel);
          % prev_Kcat:0.013;
          % Higher than prev_Kcat not available in BRENDA (S.cerevisiae);
          % New value of S.A. chosen according to: EC number/any organism;
          % Required for Xexp condition;
          if strcmpi('prot_M7WPW0',enzName) && contains(reaction,'24-sterol-c-methyltransferase')
              % S.cerevisiae S.A. [2020-05-01];
              % https://www.brenda-enzymes.org/literature.php?e=2.1.1.41&r=485321;
              newValue      = -(0.53*172*60)^-1;             %1.52
              modifications{1} = [modifications{1}; 'M7WPW0'];
              modifications{2} = [modifications{2}; reaction];
          end
          % [M7XNL9/EC2.2.1.1] transketolase (No1):
          % 5.3160e-02 of total protein pool; prev_Kcat:121.8003;
          % The highest Kcats available in BRENDA for the correct
          % substrates were for S.cerevisiae up to 69 [1/s].
          % New value of S.A. chosen according to: EC number/any organism;
          % the highest was chosen; Further artificially increased 2 times;
          %if strcmpi('prot_M7XNL9',enzName) && contains(reaction,'transketolase')
              % E.coli S.A. [2020-05-01];
              % https://www.brenda-enzymes.org/literature.php?e=2.2.1.1&r=486023;
              %newValue      = -(2*50.4*145*60)^-1;
              %newValue      = -(50.4*145*60)^-1;             %121.8 prev_Kcat
              %newValue      = -(69*3600)^-1;
              %modifications{1} = [modifications{1}; 'M7XNL9'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7X5F4/EC4.1.2.13] fructose-bisphosphate aldolase:
          % 4.4339e-01 of total protein pool (C);
          % prev_Kcat: 64.5005 (forward) and 0.002 (reverse);
          % Many Kcats available in BRENDA according to different
          % substrates; by curating, some Kcats calculated from S.A., some
          % increased for the respective substrate (if found):
          %if strcmpi('prot_M7X5F4',enzName) && contains(reaction,'fructose-bisphosphate aldolase (No1)')
              % T. aquaticus S.A. [2020-05-24];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.13&r=653734;
              %newValue      = -(46*165*60)^-1;             % 126.5 1/s
              %modifications{1} = [modifications{1}; 'M7X5F4'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          %if strcmpi('prot_M7X5F4',enzName) && contains(reaction,'fructose-bisphosphate aldolase (reversible) (No1)')
              % B. methanolicus [2020-05-24];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.13&r=728260;
              %newValue      = -(12.6*3600)^-1;             % D-glyceraldehyde 3-phosphate
              %modifications{1} = [modifications{1}; 'M7X5F4'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % prev_Kcat: 6.08 (forward) and 0.002 (reverse);
          % D-glyceraldehyde and dihydroxyacetone phosphate (reversible
          % reaction) no higher Kcat values available in BRENDA;
          %if strcmpi('prot_M7X5F4',enzName) && contains(reaction,'D-fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase')
              % B. subtilis S.A. [2020-05-24];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.13&r=4923;
              %newValue      = -(5.3*150*60)^-1;             % 13.25 [1/s]
              %modifications{1} = [modifications{1}; 'M7X5F4'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % Sedoheptulose 1,7-diphosphate highest Kcat in BRENDA 25.2 [1/s]
          % for D.carota, but D-erythrose 4-phosphate and dihydroxyacetone
          % phosphate have no Kcats. Therefore, S.A. of B.subtilis used for
          % both reactions as a middle way for both FWD and REV rxns:
          %if strcmpi('prot_M7X5F4',enzName) && contains(reaction,'sedoheptulose 1,7-bisphosphate D-glyceraldehyde-3-phosphate-lyase')
              % B. subtilis S.A. [2020-05-24];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.13&r=4923;
              %newValue      = -(5.3*150*60)^-1;             % 13.25 [1/s]
              %modifications{1} = [modifications{1}; 'M7X5F4'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % 7. [M7WFJ0/EC2.3.1.35] Glutamate N-acetyltransferase/ ornithine
          % transacetylase (No1):
          % 3.1109e-02 of total protein pool (version B of ecModel); prev_Kcat:0.22;
          % No Kcats available in BRENDA.
          % New value of S.A. chosen according to: EC number/organism;
          % Required for Xexp condition;
          if strcmpi('prot_M7WFJ0',enzName) && contains(reaction,'ornithine transacetylase (No1)')
             % S.cerevisiae S.A. [2020-05-01];
             % https://www.brenda-enzymes.org/literature.php?e=2.3.1.35&r=486801;
             newValue      = -(22*57*60)^-1;                %20.9
             modifications{1} = [modifications{1}; 'M7WFJ0'];
             modifications{2} = [modifications{2}; reaction];
          end
        % [M7X749/EC4.2.1.11] Enolase/ phosphopyruvate hydratase:
        % prev_Kcat:230.0056;
        % Kcats found in BRENDA: 71.4 (1/s) is reported for
        % 2-phospho-D-glycerate; 230 (1/s) is reported for
        % 2-phosphoglycerate; both measurements are for native enzymes.
        % New value of S.A. chosen according to: EC number/organism;
          %%if strcmpi('prot_M7X749',enzName)  && contains(reaction,'enolase')
               % S.cerevisiae S.A. [2020-05-01];
               % https://www.brenda-enzymes.org/literature.php?e=4.2.1.11&r=653043;
               %%newValue      = -1/(442*96*60);               %707.2
               %%modifications{1} = [modifications{1}; 'M7X749'];
               %%modifications{2} = [modifications{2}; reaction];
          %%end
          % 8. [M7X6X3/EC4.2.1.3] aconitase/ aconitate hydratase:
          % 1.1823e-01 of total protein pool (version B of ecModel); prev_Kcat:5.2;
          % Higher than prev_Kcat not available in BRENDA (Salmonella
          % enterica subsp. enterica serovar Typhimurium;
          % New S.A. chosen;
          % Required for Xexp condition;
          if strcmpi('prot_M7X6X3',enzName)
              if contains(reaction,'citrate to cis-aconitate(3-)')
                 % M.tuberculosis S.A. [2020-05-13];
                 % https://www.brenda-enzymes.org/literature.php?e=4.2.1.3&r=680491;
                 newValue     = -(118*102*60)^-1;               %200 [1/s]
                 modifications{1} = [modifications{1}; 'M7X6X3'];
                 modifications{2} = [modifications{2}; reaction];
              elseif contains(reaction,'cis-aconitate(3-) to isocitrate')
                  newValue     = -(118*102*60)^-1;
                  modifications{1} = [modifications{1}; 'M7X6X3'];
                  modifications{2} = [modifications{2}; reaction];
              end
          end
          if strcmpi('prot_M7WQ73',enzName)
              if contains(reaction,'citrate to cis-aconitate(3-)')
                 newValue     = -(118*102*60)^-1;             %200 [1/s]
                 modifications{1} = [modifications{1}; 'M7WQ73'];
                 modifications{2} = [modifications{2}; reaction];
              elseif contains(reaction,'cis-aconitate(3-) to isocitrate')
                  newValue     = -(118*102*60)^-1;
                  modifications{1} = [modifications{1}; 'M7WQ73'];
                  modifications{2} = [modifications{2}; reaction];
              end
          end
          % [M7X6S3/EC2.7.1.11] phosphofructokinase:
          % prev_Kcat:356.9995;
          % Higher than prev_Kcat not available in BRENDA;
          % New value of S.A. chosen according to: EC number/any organism;
          %%if strcmpi('prot_M7X6S3',enzName) && contains(reaction,'phosphofructokinase')
              % S.cerevisiae S.A. [2020-05-02];
              % https://www.brenda-enzymes.org/literature.php?e=2.7.1.11&r=640462;
              %%newValue      = -(60*835*60)^-1;             %835 [1/s]
              %%modifications{1} = [modifications{1}; 'M7X6S3'];
              %%modifications{2} = [modifications{2}; reaction];
          %%end
          % [M7XKF0/EC2.7.6.1] ribose-phosphate diphosphokinase/
          % phosphoribosylpyrophosphate synthetase;
          % 2.1884e-02 of total protein pool; prev_Kcat:0.65999;
          % Highest Kcats available in BRENDA were for Mycobacterium
          % tuberculosis and Pyrobaculum calidifontis;
          % New value of Kcat chosen according to: EC number/any/D-ribose
          % 5-phosphate;
          %if contains(reaction,'phosphoribosylpyrophosphate synthetase')
              %if strcmpi('prot_M7XKF0',enzName)
                 % Mycobacterium tuberculosis [2020-05-02];
                 % https://www.brenda-enzymes.org/literature.php?e=2.7.6.1&r=723571;
                 %newValue     = -(60.68*3600)^-1;
                 %modifications{1} = [modifications{1}; 'M7XKF0'];
                 %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7XHB5',enzName)
                  %newValue      = -(60.68*3600)^-1;
                  %modifications{1} = [modifications{1}; 'M7XHB5'];
                  %modifications{2} = [modifications{2}; reaction];
              %end
          %end
        % [M7XVW7/EC1.1.1.86] Ketol-acid Reductoisomerase/ acetohydroxy
        % acid isomeroreductase (No1):
        % 2.0775e-02  of total protein pool (version C); prev_Kcat:0.97999; (new_Kcat:78.3);
        % New highest Kcat chosen according to: EC number/any/(any)substrate;
        % Substrate: 2-acetyllactic acid - as in the model;
          %if strcmpi('prot_M7XVW7',enzName) 
              %if contains(reaction,'acetohydroxy acid isomeroreductase (')
                 % Salmonella enterica subsp. enterica serovar Typhimurium
                 % Kcat [2020-05-14];
                 % https://www.brenda-enzymes.org/literature.php?e=1.1.1.86&r=639174
                 %newValue         = -(78.3*3600)^-1;    %BRENDA: 2-acetolactate                 
                 % Salmonella enterica subsp. enterica serovar Typhimurium
                 % Kcat [2015-08-28]
                 %newValue         = -(18.3*3600)^-1;    %BRENDA: 2-acetolactate
                 %modifications{1} = [modifications{1}; 'M7XVW7'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7X5E3/EC2.6.1.52] Phosphoserine aminotransferase/transaminase:
          % 1.9110e-02 of total protein pool; prev_Kcat:2.4899;
          % Higher than prev_Kcat not available in BRENDA;
          % New value of S.A. chosen according to: EC number/any organism;
          %if strcmpi('prot_M7X5E3',enzName) && contains(reaction,'phosphoserine transaminase')
              % A.thaliana S.A. [2020-05-13];
              % https://www.brenda-enzymes.org/literature.php?e=2.6.1.52&r=675866;
              %newValue      = -(21.22*41*60)^-1;             %14.5 [1/s]
              %modifications{1} = [modifications{1}; 'M7X5E3'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7X3Z4/EC1.1.1.44] 6-phosphogluconate dehydrogenase
          % (NADP+-dependent, decarboxylating):
          % prev_Kcat:72.9995 new_Kcat:198.445 CC:0.11226 Err:-7.212%;
          % Highest available Kcat in BRENDA for T.maritima 325 [1/s];
          % Difficult to calculate from S.A. anything higher, because of no
          % MWs;
          % New value chosen according to: EC number/any organism;
          %if strcmpi('prot_M7X3Z4',enzName) && contains(reaction,'phosphogluconate dehydrogenase')
              % T.maritima Kcat [2020-05-13];
              % https://www.brenda-enzymes.org/literature.php?e=1.1.1.44&r=748598;
              %newValue      = -(325*3600)^-1;
              %modifications{1} = [modifications{1}; 'M7X3Z4'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7WW40/EC6.3.4.13] phosphoribosylamine-glycine ligase:
          % 1.8040e-02 of total protein pool; prev_Kcat:0.999;
          % Kcats in BRENDA reported only for G.gallus, max 7 [1/s];
          % New value chosen according to: EC number/any organism;
          %if strcmpi('prot_M7WW40',enzName)
              %if contains(reaction,'phopshoribosylaminoimidazole synthetase')
                 % E.coli S.A. [2020-05-13];
                 % https://www.brenda-enzymes.org/literature.php?e=6.3.4.13&r=649840;
                 %newValue     = -(19*49*60)^-1;               %15.5 [1/s] 
                 %modifications{1} = [modifications{1}; 'M7WW40'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          %if strcmpi('prot_M7WW40',enzName)
              %if contains(reaction,'phosphoribosylglycinamidine synthetase')
                 %newValue     = -(19*49*60)^-1;               %15.5 [1/s] 
                 %modifications{1} = [modifications{1}; 'M7WW40'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7WXM9/EC2.1.1.14]
          % 5-methyltetrahydropteroyltriglutamate-homocysteine
          % S-methyltransferase:
          % 4.7892e-02 of total protein pool; prev_Kcat:0.45;
          % Highest reported Kcat in BRENDA for S.cerevisiae 0.6 [1/s];
          % New value chosen according to: EC number/any organism;
          %if strcmpi('prot_M7WXM9',enzName)
              %if contains(reaction,'5-methyltetrahydropteroyltriglutamate-homocysteine S-methyltransferase')
                 % E.coli S.A. [2020-05-13];
                 % https://www.brenda-enzymes.org/literature.php?e=2.1.1.14&r=441323;
                 %newValue     = -(2.5*84*60)^-1;               %3.5 [1/s] 
                 %modifications{1} = [modifications{1}; 'M7WXM9'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          %if strcmpi('prot_M7WXM9',enzName)
              %if contains(reaction,'methionine synthase')
                 %newValue     = -(2.5*84*60)^-1;               %3.5 [1/s] 
                 %modifications{1} = [modifications{1}; 'M7WXM9'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7WL03/EC1.3.1.70] DELTA14-sterol reductase:
          % 1.9823e-02 of total protein pool (topUsedEnzymes);
          % prev_Kcat:0.064;
          % No Kcats available in BRENDA; the highest S.A. reported for
          % R.norvegicus already used as prev_Kcat;
          % New value artificially increased 10 times;
          %if strcmpi('prot_M7WL03',enzName) && contains(reaction,'C-14 sterol reductase')
              %newValue      = -(10*0.064*3600)^-1;
              %modifications{1} = [modifications{1}; 'M7WL03'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7XBZ9/EC1.14.14.17]squalene epoxidase:
          % 1.9289e-02 of total protein pool; prev_Kcat:0.076;
          % Highest reported Kcat in BRENDA for R.norvegicus 0.076 [1/s];
          % New value chosen according to: EC number/any organism;
          %if strcmpi('prot_M7XBZ9',enzName)
              %if contains(reaction,'squalene epoxidase (NAD)')
                 % R.norvegicus S.A. [2020-05-13];
                 % https://www.brenda-enzymes.org/literature.php?e=1.14.14.17&r=438206;
                 %newValue     = -(0.17*100*60)^-1;               %0.283 [1/s] 
                 %modifications{1} = [modifications{1}; 'M7XBZ9'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          %if strcmpi('prot_M7XBZ9',enzName)
              %if contains(reaction,'squalene epoxidase (NADP)')
                 %newValue     = -(0.17*100*60)^-1;               %0.283 [1/s]
                 %modifications{1} = [modifications{1}; 'M7XBZ9'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7WKX0/EC4.2.3.5] chorismate synthase:
          % 3.5680e-02 of total protein pool; prev_Kcat:0.87;
          % Highest Kcat in BRENDA for N.crassa 0.87 [1/s];
          % New value artificially increased 10 times;
          %if strcmpi('prot_M7WKX0',enzName) && contains(reaction,'chorismate synthase')
              % N.crassa S.A. [2020-05-13];
              % https://www.brenda-enzymes.org/literature.php?e=4.2.3.5&r=137645;
              %newValue      = -(32.1*198*60)^-1;            %105.93 [1/s]
              %modifications{1} = [modifications{1}; 'M7WKX0'];
              %modifications{2} = [modifications{2}; reaction];
          %end
        % [M7WMT2/EC4.2.1.36] Homoaconitase, mitochondrial:
        % 2.7797e-02 of total protein pool; prev_Kcat 2.5 [1/s];
        % For (Z)-but-1-ene-1,2,4-tricarboxylate highest Kcat available in
        % BRENDA 2.5 [1/s]; no S.A. available in BRENDA;
        % New Kcat chosen according to: EC number/any/any substrate;
        % In the 4th round, increased artificially 2 times;
          %if contains(reaction,'homoacontinate hydratase')
              %if strcmpi('prot_M7X6X3',enzName)
                 % M. jannaschii Kcat [2020-05-14];
                 % https://www.brenda-enzymes.org/literature.php?e=4.2.1.36&r=702368;
                 %newValue     = -(2*6.6*3600)^-1;         %(Z)-pent-1-ene-1,2,5-tricarboxylate;
                 %modifications{1} = [modifications{1}; 'M7X6X3'];
                 %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7WMT2',enzName)
                  %newValue      = -(2*6.6*3600)^-1;
                  %modifications{1} = [modifications{1}; 'M7WMT2'];
                  %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7WQ73',enzName)
                  %newValue      = -(2*6.6*3600)^-1;
                  %modifications{1} = [modifications{1}; 'M7WQ73'];
                  %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          %if strcmpi('prot_M7WMT2',enzName) && contains(reaction,'2-methylcitrate dehydratase')
              %newValue      = -(2*6.6*3600)^-1;
              %modifications{1} = [modifications{1}; 'M7WMT2'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7XAJ6/EC2.6.1.13] Ornithine aminotransferase :
          % 3.2100e-01 of total protein pool (version C); prev_Kcat: 0.033 and 4.3;
          % Highest reported Kcat in BRENDA for P.sativum 4.3 [1/s] and
          % M.sexta 0.033 [1/s];
          % New value chosen according to: EC number/any organism;
          %if strcmpi('prot_M7XAJ6',enzName)
              %if contains(reaction,'adenosylmethionine-8-amino-7-oxononanoate transaminase')
                 % R.norvegicus S.A. [2020-05-13];
                 %https://www.brenda-enzymes.org/literature.php?e=2.6.1.13&r=637122;
                 %newValue     = -(16*161*60)^-1;               %42.9 [1/s] 
                 %modifications{1} = [modifications{1}; 'M7XAJ6'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          %if strcmpi('prot_M7XAJ6',enzName)
              %if contains(reaction,'ornithine transaminase')
                 %newValue     = -(16*161*60)^-1;               %42.9 [1/s]
                 %modifications{1} = [modifications{1}; 'M7XAJ6'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7WHZ8/EC2.3.1.15] Glycerol-3-phosphate-acyltransferase/
          % glycerol-3-phosphate 1-O-acyltransferase:
          % 7.0512e-02 of total protein pool (version C); prev_Kcat: 0.0918;
          % Highest reported Kcat in BRENDA for H.annuus 0.0918 [1/s];
          % New value chosen according to: EC number/any organism;
          %if strcmpi('prot_M7WHZ8',enzName)
              %if contains(reaction,'glycerol-3-phosphate acyltransferase (16:0)')
                 % E.coli S.A. [2020-05-13];
                 %https://www.brenda-enzymes.org/literature.php?e=2.3.1.15&r=486357;
                 %newValue     = -(6.8*88*60)^-1;               %9.97 [1/s] 
                 %modifications{1} = [modifications{1}; 'M7WHZ8'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          %if strcmpi('prot_M7WHZ8',enzName)
              %if contains(reaction,'glycerol-3-phosphate acyltransferase (18:0)')
                 %newValue     = -(6.8*88*60)^-1;               %9.97 [1/s]
                 %modifications{1} = [modifications{1}; 'M7WHZ8'];
                 %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7XIT2/EC6.1.1.17] Glutamyl-tRNA synthetase:
          % 2.7554e-02 of total protein pool; prev_Kcat:2.2;
          % Highest Kcat available in BRENDA for Plasmodium falciparum 3.6
          % [1/s]; No higher S.A. available in BRENDA;
          % In 4th round, artificially increased 10 times;
          %if strcmpi('prot_M7XIT2',enzName) && contains(reaction,'glutamyl-tRNA synthetase')
              % P.falciparum Kcat [2020-05-14];
              % https://www.brenda-enzymes.org/literature.php?e=6.1.1.17&r=727963;
              %newValue      = -(10*3.6*3600)^-1;
              %newValue      = -(3.6*3600)^-1;
              %modifications{1} = [modifications{1}; 'M7XIT2'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7XJ70/EC2.6.1.11] Acetylornithine aminotransferase:
          % 3.7934e-02 of total protein pool; prev_Kcat:1.55;
          % Highest Kcat available in BRENDA for PSalmonella enterica
          % subsp. enterica serovar Typhimurium 1.55[1/s];
          %if strcmpi('prot_M7XJ70',enzName) && contains(reaction,'acteylornithine transaminase')
              % Klebsiella aerogenes S.A. [2020-05-14];
              %https://www.brenda-enzymes.org/literature.php?e=2.6.1.11&r=644949;
              %newValue      = -(34.6*56*60)^-1;
              %modifications{1} = [modifications{1}; 'M7XJ70'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7WNB5/EC1.14.19.41] C-22 sterol desaturase:
          % 3.7515e-02 of total protein pool; prev_Kcat:0.17;
          % Highest Kcat available in BRENDA for A.thaliana 0.17[1/s]; No
          % S.A. data available on BRENDA; Artificially increased 10 times;
          %if strcmpi('prot_M7WNB5',enzName) && contains(reaction,'C-22 sterol desaturase (NADP)')
              %newValue      = -(10*0.17*3600)^-1;
              %modifications{1} = [modifications{1}; 'M7WNB5'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % 9. [M7WG08/EC1.1.1.302]
          % 2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone 5'-phosphate
          % reductase (NADPH):
          % 2.0551e-02 of total protein pool; prev_Kcat:0.0006;
          % No Kcats reported on BRENDA; highest S.A. available in BRENDA
          % for Methanocaldococcus jannaschii 0.0008 (prev_Kcat calculated
          % from that);
          % Artificially increased 10 times;
          % Required for Xexp condition;
          if strcmpi('prot_M7WG08',enzName) && contains(reaction,'2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone')
              newValue      = -(10*0.0008*50*60)^-1;
              modifications{1} = [modifications{1}; 'M7WG08'];
              modifications{2} = [modifications{2}; reaction];
          end
          % [M7WGA7/EC4.1.2.9] Phosphoketolase:
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat:39.2;
          % Previous Kcat is the highest S.A. for Bifidobacterium animalis
          % subsp. lactis;
          % Higher Kcats reported in BRENDA for different substrates or
          % organisms;
          %if strcmpi('prot_M7WGA7',enzName) && contains(reaction,'phosphoketolase')
              % Lactococcus lactis subsp. lactis [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.9&r=746881;
              %newValue      = -(117.8*3600)^-1;             %D-fructose 6-phosphate
              %modifications{1} = [modifications{1}; 'M7WGA7'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7WHA9/EC1.3.3.6] Acyl-coenzyme A oxidase:
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat:from 1.53 (decanoyl-CoA) to 77.99 (palmitoyl-CoA);
          % Highest Kcats reported on BRENDA used as prev_Kcat;
          % New Kcat chosen from S.A.;
          %if strcmpi('prot_M7WHA9',enzName) && contains(reaction,'acyl-CoA oxidase')
              % Paenarthrobacter ureafaciens S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=1.3.3.6&r=672429;
              %newValue      = -(60.9*190*60)^-1;             %192.85
              %modifications{1} = [modifications{1}; 'M7WHA9'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          %if strcmpi('prot_M7WHA9',enzName) && contains(reaction,'fatty acid oxidation')
              %newValue      = -(60.9*190*60)^-1;             %192.85
              %modifications{1} = [modifications{1}; 'M7WHA9'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7WHZ8/EC2.3.1.15] Glycerol-3-phosphate-acyltransferase:
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat:0.0918;
          % No higher Kcats reported on BRENDA;
          %if strcmpi('prot_M7WHZ8',enzName) && contains(reaction,'glycerol-3-phosphate acyltransferase')
              % E.coli S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=2.3.1.15&r=486357;
              %newValue      = -(6.8*88*60)^-1;             %9.97 [1/s]
              %modifications{1} = [modifications{1}; 'M7WHZ8'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % 10. [M7X6B9/EC6.2.1.3] fatty-acid--CoA ligase:
          % 6.5097e-02 of total protein pool (version B of ecModel); prev_Kcat: 0.029;
          % and [M7WNY8/EC6.2.1.3]; prev_Kcat:1.2;
          % No higher Kcats reported on BRENDA;
          % Required for Xexp condition;
          if strcmpi('prot_M7WNY8',enzName) && contains(reaction,'fatty-acid--CoA ligase')
              % Y.lipolytica S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=6.2.1.3&r=759;
              newValue      = -(10.8*84*60)^-1;             %15.12 [1/s]
              modifications{1} = [modifications{1}; 'M7WNY8'];
              modifications{2} = [modifications{2}; reaction];
          end
          if strcmpi('prot_M7X6B9',enzName) && contains(reaction,'fatty-acid--CoA ligase')
              % Y.lipolytica S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=6.2.1.3&r=759;
              newValue      = -(10.8*84*60)^-1;             %15.12 [1/s]
              modifications{1} = [modifications{1}; 'M7X6B9'];
              modifications{2} = [modifications{2}; reaction];
          end
          % [M7WS17/EC6.4.1.1]: Pyruvate carboxylase;
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat:60.0 (highest reported Kcat in BRENDA);
          %if contains(reaction,'pyruvate carboxylase')
              %if strcmpi('prot_M7XNE0',enzName)
                 % R.norvegicus S.A. [2020-06-09];
                 % https://www.brenda-enzymes.org/literature.php?e=6.4.1.1&r=1758;
                 %newValue     = -(34*501*60)^-1;            % 283.9 [1/s]
                 %modifications{1} = [modifications{1}; 'M7XNE0'];
                 %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7WMD7',enzName)
                  %newValue      = -(34*501*60)^-1;
                  %modifications{1} = [modifications{1}; 'M7WMD7'];
                  %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7WS17',enzName)
                  %newValue      = -(34*501*60)^-1;
                  %modifications{1} = [modifications{1}; 'M7WS17'];
                  %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7WUV4/EC4.4.1.1] Cystathionine gamma-lyase:
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat: from 0.67 to 7.5 (L-cystathionine);
          % No higher Kcats reported on BRENDA;
          %if contains(reaction,'cystathionine g-lyase')
              %if strcmpi('prot_M7WUV4',enzName)
                 % R.norvegicus S.A. [2020-06-09];
                 % https://www.brenda-enzymes.org/literature.php?e=4.4.1.1&r=721432;
                 %newValue      = -(24*160*60)^-1;             %64 [1/s]
                 %modifications{1} = [modifications{1}; 'M7WUV4'];
                 %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7WT61',enzName)
                  %newValue      = -(24*160*60)^-1;             %64 [1/s]
                  %modifications{1} = [modifications{1}; 'M7WT61'];
                  %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % Same enzyme M7WUV4 in different reaction, where isoenzymes had
          % higher Kcat, so for all Kcat increased to:
          %if contains(reaction,'O-succinylhomoserine lyase')
              %if strcmpi('prot_M7WWQ2',enzName)
                 % Aspergillus carneus S.A. [2020-06-09];
                 % https://www.brenda-enzymes.org/literature.php?e=4.4.1.1&r=748430;
                 %newValue     = -(48.71*204*60)^-1;            % 165.61 [1/s]
                 %modifications{1} = [modifications{1}; 'M7WWQ2'];
                 %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7WUV4',enzName)
                  %newValue     = -(48.71*204*60)^-1;            % 165.61 [1/s]
                  %modifications{1} = [modifications{1}; 'M7WUV4'];
                  %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7WT61',enzName)
                  %newValue     = -(48.71*204*60)^-1;            % 165.61 [1/s]
                  %modifications{1} = [modifications{1}; 'M7WT61'];
                  %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7WWQ2/EC2.5.1.48] cystathionine gamma-synthase:
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat:121;
          % No higher Kcats reported on BRENDA; only 130 [1/s] for
          % O-succinyl-L-homoserine as substrate;
          %if strcmpi('prot_M7WWQ2',enzName) && contains(reaction,'cystathionine gamma-synthase')
              % Lysinibacillus sphaericus S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=2.5.1.48&r=637439;
              %newValue      = -(118*165*60)^-1;             %324.5 [1/s]
              %modifications{1} = [modifications{1}; 'M7WWQ2'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          %if strcmpi('prot_M7WWQ2',enzName) && contains(reaction,'O-Succinyl-L-homoserine succinate-lyase')
              % Lysinibacillus sphaericus S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=2.5.1.48&r=637439;
              %newValue      = -(118*165*60)^-1;             %324.5 [1/s]
              %modifications{1} = [modifications{1}; 'M7WWQ2'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7WY92/EC1.4.1.14] GOGAT, glutamate synthase:
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat:70; highest in BRENDA (Lupinus angustifolius);
          % No higher Kcats reported on BRENDA;
          %if strcmpi('prot_M7WY92',enzName) && contains(reaction,'glutamate synthase (NADH2)')
              % S.cerevisiae S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=1.4.1.14&r=391460;
              %newValue      = -(41.7*265*60)^-1;             %184.18[1/s]
              %modifications{1} = [modifications{1}; 'M7WY92'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7X5E3/EC2.6.1.52] Phosphoserine aminotransferase:
          % 3.2274e-02 of total protein pool;
          % prev_Kcat:2.49; highest in BRENDA (Synechocystis sp.);
          % No higher Kcats reported on BRENDA;
          %if strcmpi('prot_M7X5E3',enzName) && contains(reaction,'phosphoserine transaminase')
              % Bos taurus S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=2.6.1.52&r=640160;
              %newValue      = -(15*90.7*60)^-1;             %22.68[1/s]
              %modifications{1} = [modifications{1}; 'M7X5E3'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7XHV2/EC1.1.1.38] Malic enzyme:
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat:134; highest in BRENDA (Escherichia coli K-12);
          % No higher Kcats reported on BRENDA;
          %if strcmpi('prot_M7XHV2',enzName) && contains(reaction,'malic enzyme')
              % E.coli W S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=1.1.1.38&r=286697;
              %newValue      = -(177*220*60)^-1;             %649[1/s]
              %modifications{1} = [modifications{1}; 'M7XHV2'];
              %modifications{2} = [modifications{2}; reaction];
          %end
          % [M7XHV2/EC5.4.2.11] phosphoglycerate mutase:
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat:320; incorrect name in Uniprot, correct EC number;
          %if contains(reaction,'phosphoglycerate mutase')
              %if strcmpi('prot_M7XL58',enzName)
                 % S.cerevisiae [2020-06-09];
                 % https://www.brenda-enzymes.org/literature.php?e=5.4.2.11&r=652799;
                 %newValue     = -(530*3600)^-1;
                 %modifications{1} = [modifications{1}; 'M7XL58'];
                 %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7XSI3',enzName)
                  %newValue     = -(530*3600)^-1;
                  %modifications{1} = [modifications{1}; 'M7XSI3'];
                  %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7XRQ3',enzName)
                  %newValue     = -(530*3600)^-1;
                  %modifications{1} = [modifications{1}; 'M7XRQ3'];
                  %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7XV41/EC2.4.1.15] alpha,alpha-trehalose-phosphate synthase
          % (UDP-forming):
          % 100% enzyme usage for growth on xylose;
          % prev_Kcat: 5.1;
          %if contains(reaction,'alpha,alpha-trehalose-phosphate synthase (UDP-forming)')
              %if strcmpi('prot_M7XCG8',enzName)
                 % S.cerevisiae S.A. [2020-06-09];
                 % https://www.brenda-enzymes.org/literature.php?e=2.4.1.15&r=134590;
                 %newValue      = -(15*300*60)^-1;             %75 [1/s]
                 %modifications{1} = [modifications{1}; 'M7XCG8'];
                 %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7XV41',enzName)
                  %newValue      = -(15*300*60)^-1;             %75 [1/s]
                  %modifications{1} = [modifications{1}; 'M7XV41'];
                  %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % [M7X383/EC6.2.1.3] Long-chain acyl-CoA synthetase:
          % 2.8999e-02 of total protein pool; prev_Kcat 0.029 (on ATP as
          % substrate);
          %if contains(reaction,'fatty-acid--CoA ligase')
              %if strcmpi('prot_M7X383',enzName)
                 % Trypanosoma brucei [2020-07-12];
                 % https://www.brenda-enzymes.org/literature.php?e=6.2.1.3&r=649717;
                 %newValue     = -(1.2*3600)^-1;          % tetradecanoate
                 %modifications{1} = [modifications{1}; 'M7X383'];
                 %modifications{2} = [modifications{2}; reaction];
              %elseif strcmpi('prot_M7WW26',enzName)
                  %newValue     = -(1.2*3600)^-1;
                  %modifications{1} = [modifications{1}; 'M7WW26'];
                  %modifications{2} = [modifications{2}; reaction];
              %end
          %end
          % 11. [M7WT53/EC2.3.1.158] Phospholipid:diacylglycerol acyltransferase:
          % 5.1163e-01 of total protein pool (B version of ecModel);
          % 3.0163e-01 of total protein pool (C version of ecModel);
          % prev_Kcat:0.00346; it's a specific activity of S.cerevisiae as
          % reported on BRENDA; No Kcats available on BRENDA;
          % Therefore, Kcat value artificially increased 10 and 100 times;
          % Required for Xexp condition;
          if strcmpi('prot_M7WT53',enzName) && contains(reaction,'diacylglycerol')
              % S.cerevisiae S.A. [2020-07-12];
              % https://www.brenda-enzymes.org/literature.php?e=2.3.1.158&r=685341;
              newValue      = -(10*0.002597*80*60)^-1;             % 10*0.00346[1/s]
              %newValue      = -(100*0.002597*80*60)^-1;             % 100*0.00346[1/s]
              modifications{1} = [modifications{1}; 'M7WT53'];
              modifications{2} = [modifications{2}; reaction];
          end
           % 12. EC1.9.3.1 - cytochrome-c oxidase;
           % 0.034 to 0.086 % of total protein pool at version A;
           % prev_Kcat for M7WUI0 0.15634;
           % The highest available S.A. in BRENDA gives Kcat of 712.5 [1/s]
           % for S.cerevisiae;
           % Required for Xexp condition;
          if contains(reaction,'ferrocytochrome-c:oxygen oxidoreductase')
              if strcmpi('prot_M7WUI0',enzName)
                 % S.cerevisiae Kcat [2020-05-08];
                 % https://www.nature.com/articles/srep22264;
                 newValue     = -(693.5*3600)^-1;
                 modifications{1} = [modifications{1}; 'M7WUI0'];
                 modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7XMJ4',enzName)
                  newValue      = -(693.5*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7XMJ4'];
                  modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7XKJ8',enzName)
                  newValue      = -(693.5*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7XKJ8'];
                  modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7WI20',enzName)
                  newValue      = -(693.5*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7WI20'];
                  modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7WMJ7',enzName)
                  newValue      = -(693.5*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7WMJ7'];
                  modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7WD99',enzName)
                  newValue      = -(693.5*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7WD99'];
                  modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7X3R7',enzName)
                  newValue      = -(693.5*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7X3R7'];
                  modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7WWJ6',enzName)
                  newValue      = -(693.5*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7WWJ6'];
                  modifications{2} = [modifications{2}; reaction];
              elseif strcmpi('prot_M7XKA8',enzName)
                  newValue      = -(693.5*3600)^-1;
                  modifications{1} = [modifications{1}; 'M7XKA8'];
                  modifications{2} = [modifications{2}; reaction];
              end
          end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = otherChanges(model)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modified = mapModifiedRxns(modifications,model)
modified = [];
for i=1:length(modifications{1})
    rxnIndex = find(strcmp(model.rxnNames,modifications{2}(i)),1);
    str      = {horzcat(modifications{1}{i},'_',num2str(rxnIndex))};
    modified = [modified; str];
end
end