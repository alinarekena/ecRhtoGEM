
%Load chain data:

i=1;
fid = fopen([root '/data/lipidChain.csv']);
chainData = textscan(fid,['%q %q %q %f32 %q' repmat('%f32 ',[1,12])],'Delimiter',',','HeaderLines',1);
data.chainData.chain     = chainData{1};
data.chainData.FA        = chainData{2};
data.chainData.formulas  = chainData{3};
data.chainData.MW        = chainData{4};
data.chainData.comp      = chainData{5};
data.chainData.abundance = chainData{4+2*i};%first sample is Xexp
data.chainData.std       = chainData{5+2*i};
for j=1:length(data.chainData.chain)
    metIds = ecModelP_Xexp.mets(getIndexes(ecModelP_Xexp,['C' data.chainData.chain{j} ...
        ' chain'],'metnames'));
    data.chainData.metIds(j,1) = metIds;
end
fclose(fid);


% Change lipid chain length distribution
rxnIdx                  =   getIndexes(ecModelP_Xexp,'r_4065','rxns');%lipid chain pseudoreaction
metIdx                  =   cell2mat(getIndexes(ecModelP_Xexp,strcat('C',data.chainData.chain,' chain'),'metnames'));
bbIdx                   =   getIndexes(model,'s_3747','mets');%lipid chain
ecModelP_Xexptest.S(:,rxnIdx)       =   0;
% Normalize lipid chain data to lipid backbone level, to get reasonable
% estimate before final scaling
scaling = sum([data.lipidData.abundance;data.lipidData.std]);
scaling = scaling / sum([data.chainData.abundance;data.chainData.std],'omitnan');
scaling = (data.chainData.abundance + data.chainData.std) * scaling;
ecModelP_Xexp.S(metIdx,rxnIdx)  =   -scaling;
ecModelP_Xexp.S(bbIdx,rxnIdx)   =   1;

chainExIdx  = getIndexes(ecModelP_Xexp,'r_4064','rxns');
backbExIdx  = getIndexes(ecModelP_Xexp,'r_4062','rxns');
ecModelP_Xexp = setParam(ecModelP_Xexp,'ub',[chainExIdx,backbExIdx],1000);

[newModel,k] = scaleLipidsRhto(ecModelP_Xexp,'tails');% ??
