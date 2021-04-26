function model = scaleLipidProtein(model,lipidData,protein,GAMnonPol)
% scaleLipidProtein
%
% Function that takes an ecModel and constraints it with absolute proteomics 
% datasets for the experimental conditions provided by the user.
%
%   model       model structure
%   lipid       lipid content in g/gDW
%   protein     protein content in g/gDW
%   GAMnonPol   non-polymerization part of GAM (e.g. as defined in
%               parameters.GAM). (optional, is otherwise determined from
%               model
%
% Usage:  scaleLipidProtein(model,lipid,protein,GAMnonPol)

%Components of biomass:
%        id             MW [g/mol]  class     name
comps = {'s_0404'	89.09       'P'     % A     Alanine         ala
         's_0542'    121.16      'P'     % C     Cysteine        cys
         's_0432'    133.11      'P'     % D     Aspartic acid   asp
         's_0748'    147.13      'P'     % E     Glutamic acid   glu
         's_1314'    165.19      'P'     % F     Phenylalanine   phe
         's_0757'    75.07       'P'     % G     Glycine         gly
         's_0832'    155.15      'P'     % H     Histidine       his
         's_0847'    131.17      'P'     % I     Isoleucine      ile
         's_1099'    146.19      'P'     % K     Lysine          lys
         's_1077'    131.17      'P'     % L     Leucine         leu
         's_1148'    149.21      'P'     % M     Methionine      met
         's_0430'    132.12      'P'     % N     Asparagine      asn
         's_1379'    115.13      'P'     % P     Proline         pro
         's_0747'    146.14      'P'     % Q     Glutamine       gln
         's_0428'    174.2       'P'     % R     Arginine        arg
         's_1428'    105.09      'P'     % S     Serine          ser
         's_1491'    119.12      'P'     % T     Threonine       thr
         's_1561'    117.15      'P'     % V     Valine          val
         's_1527'    204.23      'P'     % W     Tryptophan      trp
         's_1533'    181.19      'P'     % Y     Tyrosine        tyr
         's_0001'	 180.16      'C'     % (1->3)-beta-D-glucan
         's_0004'	 180.16      'C'     % (1->6)-beta-D-glucan
         's_0773'    180.16      'C'     % glycogen
         's_1107'    180.16      'C'     % mannan
         's_1520'    342.296 	 'C'     % trehalose
         's_0423'    347.22      'R'     % AMP
         's_0526'    323.2       'R'     % CMP
         's_0782'    363.22      'R'     % GMP
         's_1545'    324.18      'R'     % UMP
         's_0584'    331.22      'D'     % dAMP
         's_0589'    307.2       'D'     % dCMP
         's_0615'    345.21      'D'     % dGMP
         's_0649'    322.21      'D'     % dTMP
         's_3714'    852.83      'N'     % heme a
         's_1405'    376.36      'N'     % riboflavin
         's_1467'    96.06       'N'};   % sulphate

% Browse to location of sumBioMass
cd ../GECKO/geckomat/limit_proteins/

% Get current non-polymerization part of GAM
if nargin<4
    GAMnonPol = splitGAEC(model);
end

% Scale protein
[~,P,~,~,~,~]   = sumBioMass(model);
protPos         = strcmp(model.rxnNames,'protein pseudoreaction');
fP              = protein/P;
isAA            = contains(model.metNames,'tRNA');	%protein components
model.S(isAA,protPos) = model.S(isAA,protPos)*fP;

% Scale lipids
lipidPos        = find(contains(model.rxnNames,{'lipid backbone pseudoreaction','lipid chain pseudoreaction'}));
[scaleLp,~]     = find(model.S(:,lipidPos)<0);
fL = 1.01;
while abs(fL-1) > 1e-12
    fL0 = fL;
    model.S(scaleLp,lipidPos) = model.S(scaleLp,lipidPos)*fL;
    [~,~,~,~,~,L] = sumBioMass(model);
    fL = lipidData.lipidContent/L;
end

model = adjustLipidChain(model,lipidData);
cd ../../../code
model=scaleLipidsRhto(model,'tails');
cd ../GECKO/geckomat/limit_proteins/

% Scale carbohydrate to correct for protein and lipid changes
[X,~,C,~,~,~] = sumBioMass(model);
delta         = X - 1;                      %difference to balance

%Balance out mass with all sugars:
mets     = comps(strcmp(comps(:,3),'C'),1);
massPre  = C;
massPost = massPre - delta;
if massPost<0
    error('Adjusted biomass sums to > 1 g/gDW')
end
fC       = massPost/massPre;
carbPos  = find(strcmp(model.rxnNames,'carbohydrate pseudoreaction'));
for i = 1:length(mets)
    modelPos = strcmp(model.mets,mets{i});
    model.S(modelPos,carbPos) = model.S(modelPos,carbPos)*fC;
end

%Estimate maintenance belonging to polymerization:
[~,P,C,R,D,~] = sumBioMass(model);
GAMpol        = P*37.7 + C*12.8 + R*26.0 + D*26.0;    %Forster 2003 (sup table 8)
GAM           = GAMpol + GAMnonPol;
model         = setGAM(model,GAM);
cd ../../../code


end

function model = setGAM(model,GAM)
xr_pos = getIndexes(model','r_4041','rxns');
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,xr_pos) = sign(S_ix) * GAM;
    end
end
end

function model = adjustLipidChain(model,lipidData)
% Change lipid chain length distribution
rxnIdx                  =   getIndexes(model,'r_4065','rxns');%lipid chain pseudoreaction
metIdx                  =   cell2mat(getIndexes(model,strcat('C',lipidData.chainData.chain,' chain'),'metnames'));
bbIdx                   =   getIndexes(model,'s_3747','mets');%lipid chain
model.S(:,rxnIdx)       =   0;
% Normalize lipid chain data to lipid backbone level, to get reasonable
% estimate before final scaling
scaling = lipidData.lipidContent / sum([lipidData.chainData.abundance;lipidData.chainData.std],'omitnan');
scaling = (lipidData.chainData.abundance + lipidData.chainData.std) * scaling;
model.S(metIdx,rxnIdx)  =   -scaling;
model.S(bbIdx,rxnIdx)   =   1;

chainExIdx  = getIndexes(model,'r_4064','rxns');
backbExIdx  = getIndexes(model,'r_4062','rxns');
model = setParam(model,'ub',[chainExIdx,backbExIdx],[1000,1000]);
end
