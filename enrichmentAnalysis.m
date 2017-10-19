%%% Enrichment analysis

clear all;

modFolder = 'modules';
modFiles = dir(modFolder);

for i = 3 : length(modFiles)
    
    modFile = modFiles(i);
    temp = textread([modFolder '\' modFile.name], '%s', 'delimiter', '\n');
    
    modNames{i-2} = strrep(modFile.name, 'GeneNames-', '');
    modNames{i-2} = strrep(modNames{i-2}, '.txt', '');
    
    for j = 1 : length(temp)
        
        temp(j) = strrep(temp(j), '"', '');
        
    end
    
    modD{i-2} = temp;
    
end

dName = 'ND';
[num txt] = xlsread(['Data\Donors\' dName '.xls']);
clear num;
dgNames = txt(:,1);
clear txt;
for i = 1 : length(modD)
    
    dGenesInMod{i} = intersect(dgNames, modD{i});
    
end

%%% calculate RR scores
[num2 txt2] = xlsread(['Data\Donors\All.xls']);
clear num2;
gNames = txt2(2:end,1);
clear txt2;
N = 13563;
for i = 1: length(modD)
    
    nModDG = length(dGenesInMod{i});
    nMod = length(modD{i});
    nDG = length(dgNames);
    RR(i) = (nModDG/nMod) / ((nDG - nModDG)/(N - nMod));
    
%     %%% permutation testing
%     for k = 1 : 10000
%         rs = randsample(N, nMod);
%         randModG = gNames(rs);
%         R_dGenesInMod = intersect(dgNames, randModG);
%         
%         R_nModDG = length(R_dGenesInMod);
%         R_nMod = nMod;
%         R_nDG = nDG;
%         PT(k) = (R_nModDG/R_nMod) / ((R_nDG - R_nModDG)/(N - R_nMod));
%     end
%     
%     [h, p] = ttest(PT, RR(i));
end

figure, 
bar(RR, 'r', 'EdgeColor', 'black');
grid on
title([dName ' Genes Enrichment in Modules'], 'fontweight', 'bold')
xlabel('Modules', 'fontweight', 'bold')
ylabel('Enrichment Score', 'fontweight', 'bold')
set(gca, 'XTick', [1:length(modD)])






