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
        % xylose-->prev_Kcat:0.0053666 new_Kcat:0.0053667 CC:6.6589 Err:-81.3496%
        % acetate-->prev_Kcat:0.0053666 new_Kcat:0.0053667 CC:3.4618 Err:-94.1527%
        % No Kcats available in BRENDA; S.A. values for EC2.4.1.109 were
        % very low;
        % New Kcat chosen according to a wild card:
        % EC2.4.1.-/any/'substrate';
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
          % [M7X5F4/EC4.1.2.13] fructose-bisphosphate aldolase:
          % acetate-->prev_Kcat:0.002 new_Kcat:1562 CC:138.5014 Err:344.9769%
          % Many Kcats available in BRENDA according to different
          % substrates; by curating, some Kcats calculated from S.A., some
          % increased for the respective substrate (if found):
          if strcmpi('prot_M7X5F4',enzName) && contains(reaction,'fructose-bisphosphate aldolase (No1)')
              % T. aquaticus S.A. [2020-05-24];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.13&r=653734;
              newValue      = -(46*165*60)^-1;             % 126.5 1/s
              modifications{1} = [modifications{1}; 'M7X5F4'];
              modifications{2} = [modifications{2}; reaction];
          end
          if strcmpi('prot_M7X5F4',enzName) && contains(reaction,'fructose-bisphosphate aldolase (reversible) (No1)')
              % B. methanolicus [2020-05-24];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.13&r=728260;
              newValue      = -(12.6*3600)^-1;             % D-glyceraldehyde 3-phosphate
              modifications{1} = [modifications{1}; 'M7X5F4'];
              modifications{2} = [modifications{2}; reaction];
          end
          % prev_Kcat: 6.08 (forward) and 0.002 (reverse);
          % D-glyceraldehyde and dihydroxyacetone phosphate (reversible
          % reaction) no higher Kcat values available in BRENDA;
          if strcmpi('prot_M7X5F4',enzName) && contains(reaction,'D-fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase')
              % B. subtilis S.A. [2020-05-24];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.13&r=4923;
              newValue      = -(5.3*150*60)^-1;             % 13.25 [1/s]
              modifications{1} = [modifications{1}; 'M7X5F4'];
              modifications{2} = [modifications{2}; reaction];
          end
          % Sedoheptulose 1,7-diphosphate highest Kcat in BRENDA 25.2 [1/s]
          % for D.carota, but D-erythrose 4-phosphate and dihydroxyacetone
          % phosphate have no Kcats. Therefore, S.A. of B.subtilis used for
          % both reactions as a middle way for both FWD and REV rxns:
          if strcmpi('prot_M7X5F4',enzName) && contains(reaction,'sedoheptulose 1,7-bisphosphate D-glyceraldehyde-3-phosphate-lyase')
              % B. subtilis S.A. [2020-05-24];
              % https://www.brenda-enzymes.org/literature.php?e=4.1.2.13&r=4923;
              newValue      = -(5.3*150*60)^-1;             % 13.25 [1/s]
              modifications{1} = [modifications{1}; 'M7X5F4'];
              modifications{2} = [modifications{2}; reaction];
          end
        % 2. [M7XI04/EC1.1.1.34] HMG-CoA reductase/ 3-hydroxy-3-methylglutaryl
        % coenzyme A reductase:
        % xylose-->1.7122e-01 of total protein pool at round A of ecRhtoGEM generation;
        % acetate-->3.6845e-01 of total pool at round A;
        % Higher than prev_Kcat not available in BRENDA (Rattus
        % Norvegicus).
        % New value of S.A. chosen according to: EC number/any organism;
          if (strcmpi('prot_M7XI04',enzName) && (contains(reaction,'hydroxymethylglutaryl CoA reductase')))
             %Homarus americanus S.A. [2020-05-14];
             %https://www.brenda-enzymes.org/literature.php?e=1.1.1.34&r=654546;
             %newValue         = -(60*140*60)^-1;            %140 [1/s]
             % Streptomyces sp. S.A. [2020-04-28];
             % https://www.brenda-enzymes.org/literature.php?e=1.1.1.34&r=286581;
             newValue         = -(1.6*100*60)^-1;          %2.67 (it works)
             modifications{1} = [modifications{1}; 'M7XI04'];
             modifications{2} = [modifications{2}; reaction];
          end
         % 3. [M7XKV7//EC#:6.3.5.3] 5'-phosphoribosylformyl glycinamidine
         % synthetase:
         % xylose-->8.4551e-02 of total protein pool at round A of ecRhtoGEM generation;
         % acetate-->1.3833e-01 of total pool at round A
         % Higher than prev_Kcat not available in BRENDA (E.coli).
         % New value of S.A. chosen according to: EC number/any organism;
          if (strcmpi('prot_M7XKV7',enzName)) && (contains(reaction,'phosphoribosylformyl'))
             % E.coli S.A. [2020-04-28];
             % https://www.brenda-enzymes.org/literature.php?e=6.3.5.3&r=1670;
             newValue         = -(2.15*141.418*60)^-1;      %5.0675
             modifications{1} = [modifications{1}; 'M7XKV7'];
             modifications{2} = [modifications{2}; reaction];
          end
        % 4. [M7WFB4//EC#:2.4.1.34] 1,3-beta-glucan synthase (No1):
        % xylose-->7.0144e-02 of total protein pool at round A of ecRhtoGEM generation;
        % acetate-->1.5094e-01 of total pool at round A
        % No higher than prev_Kcat available in BRENDA (S.aureus).
        % New value of S.A. chosen according to: EC number/any organism;
          if (strcmpi('prot_M7WFB4',enzName)) && (contains(reaction,'1,3-beta-glucan synthase'))
             % S.cerevisiae S.A. [2020-04-28];
             % https://www.brenda-enzymes.org/literature.php?e=2.4.1.34&r=489076;
             newValue         = -(4*200*60)^-1;             %13.3 (works!)
             % 2.0556e-02 of total protein pool for X20 model; prev_Kcat 13.3;
             % S.A. of 27.9 reported for P.sativum, so Kcat increased 7 times:
             % https://www.brenda-enzymes.org/literature.php?e=2.4.1.34&r=489091;
             %newValue         = -(7*4*200*60)^-1;             %93.1
             modifications{1} = [modifications{1}; 'M7WFB4'];
             modifications{2} = [modifications{2}; reaction];
          end
        % 5. [M7XI95//EC#:1.14.19.1] Delta-9 fatty acid desaturase ER membrane
        % (No1); stearoyl- or oleoyl- or palmitoyl-CoA desaturase
        % (n-C18:0CoA -> n-C18:1CoA),(n-C18:1CoA -n-C18:2CoA),(n-C16:0CoA
        % -> n-C16:1CoA);
        % xylose-->2.6441e-02  of total protein pool at round B of ecModel generation;
        % acetate-->3.5704e-02 of total pool at round A
        % Kcat values available in BRENDA were for Rattus norvegicus and
        % spinach. Further artificially increased.
          if (strcmpi('prot_M7XI95',enzName)) && (contains(reaction,'-CoA desaturase'))
             % Spinacia oleracea [2020-05-01]
             % https://www.brenda-enzymes.org/literature.php?e=1.14.19.1&r=437670;
             %newValue         = -(0.5*3600)^-1;            %prev_Kcat
             newValue         = -(10*0.5*3600)^-1;
             modifications{1} = [modifications{1}; 'M7XI95'];
             modifications{2} = [modifications{2}; reaction];
          end
          % 6. [M7WPW0/EC2.1.1.41] sterol 24-C-methyltransferase:
          % xylose-->2.3058e-02 of total protein pool (version B of ecModel);
          % acetate-->3.1136e-02 of total pool at round A
          % prev_Kcat:0.013;
          % Higher than prev_Kcat not available in BRENDA (S.cerevisiae);
          % New value of S.A. chosen according to: EC number/any organism;
          if strcmpi('prot_M7WPW0',enzName) && contains(reaction,'24-sterol-c-methyltransferase')
              % S.cerevisiae S.A. [2020-05-01];
              % https://www.brenda-enzymes.org/literature.php?e=2.1.1.41&r=485321;
              newValue      = -(0.53*172*60)^-1;             %1.52
              modifications{1} = [modifications{1}; 'M7WPW0'];
              modifications{2} = [modifications{2}; reaction];
          end
          % 7. [M7WFJ0/EC2.3.1.35] Glutamate N-acetyltransferase/ ornithine
          % transacetylase (No1):
          % xylose-->3.1109e-02 of total protein pool (version B of ecModel);
          % acetate-->3.3393e-02 of total pool at round B
          % No Kcats available in BRENDA; prev_Kcat:0.22;
          % New value of S.A. chosen according to: EC number/organism;
          if strcmpi('prot_M7WFJ0',enzName) && contains(reaction,'ornithine transacetylase (No1)')
             % S.cerevisiae S.A. [2020-05-01];
             % https://www.brenda-enzymes.org/literature.php?e=2.3.1.35&r=486801;
             newValue      = -(22*57*60)^-1;                %20.9
             modifications{1} = [modifications{1}; 'M7WFJ0'];
             modifications{2} = [modifications{2}; reaction];
          end
          % 8. [M7X6X3/EC4.2.1.3] aconitase/ aconitate hydratase:
          % xylose-->2.8584e-01 of total protein pool )topUsedEnzymes); prev_Kcat:5.2;
          % acetate-->prev_Kcat:5.2 new_Kcat:5.2 CC:0.097124 Err:-94.1527%
          % Higher than prev_Kcat not available in BRENDA (Salmonella
          % enterica subsp. enterica serovar Typhimurium;
          % New S.A. chosen;
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
          % 9. [M7WG08/EC1.1.1.302]
          % 2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone 5'-phosphate
          % reductase (NADPH):
          % xylose-->2.0551e-02 of total protein pool; prev_Kcat:0.0006;
          % acetate-->2.1099e-02 of total pool at round A
          % No Kcats reported on BRENDA; highest S.A. available in BRENDA
          % for Methanocaldococcus jannaschii 0.0008 (prev_Kcat calculated
          % from that);
          % Artificially increased 10 times;
          if strcmpi('prot_M7WG08',enzName) && contains(reaction,'2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone')
              newValue      = -(10*0.0008*50*60)^-1;
              modifications{1} = [modifications{1}; 'M7WG08'];
              modifications{2} = [modifications{2}; reaction];
          end
          % 10. [M7X6B9/M7WNY8/M7X383/M7WW26/EC6.2.1.3] fatty-acid--CoA
          % ligase/long-chain acyl-coa synthetase: prev_Kcat: 0.029;
          % xylose-->6.5097e-02 of total protein pool (version B of ecModel);
          % acetate-->1.0845e-01 of total pool at round B
          % and [M7WNY8/EC6.2.1.3]; prev_Kcat:1.2;
          % No higher Kcats reported on BRENDA;
          if strcmpi('prot_M7WNY8',enzName) && contains(reaction,'fatty-acid--CoA ligase')
              % Y.lipolytica S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=6.2.1.3&r=759;
              %newValue      = -(10.8*84*60)^-1;             %15.12 [1/s]
              newValue      = -(10*10.8*84*60)^-1;
              modifications{1} = [modifications{1}; 'M7WNY8'];
              modifications{2} = [modifications{2}; reaction];
          end
          if strcmpi('prot_M7X6B9',enzName) && contains(reaction,'fatty-acid--CoA ligase')
              % Y.lipolytica S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=6.2.1.3&r=759;
              %newValue      = -(10.8*84*60)^-1;             %15.12 [1/s]
              newValue      = -(10*10.8*84*60)^-1;
              modifications{1} = [modifications{1}; 'M7X6B9'];
              modifications{2} = [modifications{2}; reaction];
          end
          if strcmpi('prot_M7X383',enzName) && contains(reaction,'fatty-acid--CoA ligase')
              % Y.lipolytica S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=6.2.1.3&r=759;
              %newValue      = -(10.8*84*60)^-1;             %15.12 [1/s]
              newValue      = -(10*10.8*84*60)^-1;
              modifications{1} = [modifications{1}; 'M7X383'];
              modifications{2} = [modifications{2}; reaction];
          end
          if strcmpi('prot_M7WW26',enzName) && contains(reaction,'fatty-acid--CoA ligase')
              % Y.lipolytica S.A. [2020-06-09];
              % https://www.brenda-enzymes.org/literature.php?e=6.2.1.3&r=759;
              %newValue      = -(10.8*84*60)^-1;             %15.12 [1/s]
              newValue      = -(10*10.8*84*60)^-1;
              modifications{1} = [modifications{1}; 'M7WW26'];
              modifications{2} = [modifications{2}; reaction];
          end
          % 11. [M7WT53/EC2.3.1.158] Phospholipid:diacylglycerol acyltransferase:
          % xylose-->5.1163e-01 of total protein pool (B version of ecModel);
          % xylose-->3.0163e-01 of total protein pool (C version of ecModel);
          % acetate-->4.3341e-01 of total pool at round B
          % prev_Kcat:0.00346; it's a specific activity of S.cerevisiae as
          % reported on BRENDA; No Kcats available on BRENDA;
          % Therefore, Kcat value artificially increased;
          if strcmpi('prot_M7WT53',enzName) && contains(reaction,'diacylglycerol')
              % S.cerevisiae S.A. [2020-07-12];
              % https://www.brenda-enzymes.org/literature.php?e=2.3.1.158&r=685341;
              %newValue      = -(10*0.002597*80*60)^-1;             % 10*0.00346[1/s]
              newValue      = -(100*0.002597*80*60)^-1;             % 100*0.00346[1/s]
              %newValue      = -(1000*0.002597*80*60)^-1;
              modifications{1} = [modifications{1}; 'M7WT53'];
              modifications{2} = [modifications{2}; reaction];
          end
          % [M7WXM9/EC2.1.1.14]
          % 5-methyltetrahydropteroyltriglutamate-homocysteine
          % S-methyltransferase:
          % acetate-->2.4333e-02 of total pool at round A; prev_Kcat:0.45;
          % Highest reported Kcat in BRENDA for S.cerevisiae 0.6 [1/s];
          % New value chosen according to: EC number/any organism;
          if strcmpi('prot_M7WXM9',enzName)
              if contains(reaction,'5-methyltetrahydropteroyltriglutamate-homocysteine S-methyltransferase')
                 % E.coli S.A. [2020-05-13];
                 % https://www.brenda-enzymes.org/literature.php?e=2.1.1.14&r=441323;
                 newValue     = -(2.5*84*60)^-1;               %3.5 [1/s] 
                 modifications{1} = [modifications{1}; 'M7WXM9'];
                 modifications{2} = [modifications{2}; reaction];
              end
          end
          if strcmpi('prot_M7WXM9',enzName)
              if contains(reaction,'methionine synthase')
                 newValue     = -(2.5*84*60)^-1;               %3.5 [1/s] 
                 modifications{1} = [modifications{1}; 'M7WXM9'];
                 modifications{2} = [modifications{2}; reaction];
              end
          end
          % [M7WSW5/M7XM89/EC2.3.1.86] fatty acid synthase, subunits beta and alpha:
          % acetate-->limiting enzymes at Aexp condition
          % prev_Kcat:120 is specific activity of S.cerevisiae,as
          % reported on BRENDA; No Kcats available on BRENDA;
          % Therefore, Kcat value artificially increased;
          if strcmpi('prot_M7WSW5',enzName) && contains(reaction,'fatty-acyl-CoA synthase')
              % S.cerevisiae S.A. [2020-10-29];
              % https://www.brenda-enzymes.org/literature.php?e=2.3.1.86&r=487603;
              %newValue      = -(3*2400*60)^-1;             % 120[1/s]
              newValue      = -(10*3*2400*60)^-1;
              modifications{1} = [modifications{1}; 'M7WSW5'];
              modifications{2} = [modifications{2}; reaction];
          end
          if strcmpi('prot_M7XM89',enzName) && contains(reaction,'fatty-acyl-CoA synthase')
              %newValue      = -(3*2400*60)^-1;             % 120[1/s]
              newValue      = -(10*3*2400*60)^-1;
              modifications{1} = [modifications{1}; 'M7XM89'];
              modifications{2} = [modifications{2}; reaction];
          end
          % [M7WLQ0/EC2.3.1.7] Carnitine O-acetyltransferase:
          % acetate-->limiting enzymes at Aexp condition
          % prev_Kcat:97.5994 Mus Musculus (acetyl-Coa) as reported on
          % BRENDA;
          if strcmpi('prot_M7WLQ0',enzName) && contains(reaction,'carnitine O-acetyltransferase')
              % S.cerevisiae S.A. [2021-01-07];
              % https://www.brenda-enzymes.org/literature.php?e=2.3.1.7&r=487418;
              newValue      = -(200*65*60)^-1;             % 217[1/s]
              modifications{1} = [modifications{1}; 'M7WLQ0'];
              modifications{2} = [modifications{2}; reaction];
          end
          if strcmpi('prot_M7WLQ0',enzName) && contains(reaction,'cytoplasmatic carnitine acyltransferase')
              newValue      = -(200*65*60)^-1;             % 217[1/s]
              modifications{1} = [modifications{1}; 'M7WLQ0'];
              modifications{2} = [modifications{2}; reaction];
          end
           % 13. EC1.9.3.1 - cytochrome-c oxidase;
           % xylose-->0.043 to 0.137 of total protein pool at round A;
           % acetate-->prev_Kcat:12.0001 new_Kcat:2000 CC:0.085716 Err:-93.654%
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