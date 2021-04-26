function lipidData = loadLipidChainData(model,i)
root = regexprep(pwd,[filesep 'code'],'');

%Load chain data:
fid = fopen([root '/data/lipidChain.csv']);
chainData = textscan(fid,['%q %q %q %f32 %q' repmat('%f32 ',[1,12])],'Delimiter',',','HeaderLines',1);
lipidData.chainData.chain     = chainData{1};
lipidData.chainData.FA        = chainData{2};
lipidData.chainData.formulas  = chainData{3};
lipidData.chainData.MW        = chainData{4};
lipidData.chainData.comp      = chainData{5};
lipidData.chainData.abundance = chainData{4+2*i};%first sample is Xexp
lipidData.chainData.std       = chainData{5+2*i};
for j=1:length(lipidData.chainData.chain)
    metIds = model.mets(getIndexes(model,['C' lipidData.chainData.chain{j} ...
        ' chain'],'metnames'));
    lipidData.chainData.metIds(j,1) = metIds;
end
fclose(fid);

%Load lipid data
fID      = fopen([root '/data/lipidContent.txt']);
lipidData.lipidContent = textscan(fID,'%s %f','HeaderLines',1);
lipidData.lipidContent = lipidData.lipidContent{2}(i); fclose(fID);
end
