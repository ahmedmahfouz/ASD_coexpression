function enrichAnalyze(clusts, dList, outF)

load('files/genesStatus_5RPKM.mat');
load('files/gNames.mat');
load('files/enIDs.mat');
gIND = find(genesStatus_5RPKM == 1);
gNames_5RPKM = gNames(gIND);
enIDs_5RPKM = enIDs(gIND);
N = length(gIND);
clear gNames; clear enIDs;
clear genesStatus_5RPKM;

%%% extract each cluster genes
for j = 1 : max(clusts)
    clust = find(clusts == j);
    tempList = gNames_5RPKM(clust);
    tempList2 = enIDs_5RPKM(clust);
    clustNames{j} = tempList;
    clustIDs{j} = tempList2;
    clear clust; clear tempList; clear tempList2;
end

%%% Disease Genes Enrichment
for dL = 1 : length(dList)
    clear dName; clear dgNames;
    dName = dList{dL};
    [num txt] = xlsread(['Data\Donors\' dName 'List.xls']);
    dgNames = txt(2:end,1);
    dgIDs = txt(2:end,2);
    clear num; clear txt;
    %%% calculate RR score for each cluster    
    nDG = length(dgNames);
    for k = 1 : length(clustNames)
        clear dGenesInMod; clear nModDG; clear nMod;
        % find the disease genes in each cluster
%         dGenesInMod = intersect(dgNames, clustNames{k});
        dGenesInMod = intersect(dgIDs, clustIDs{k});
        nModDG = length(dGenesInMod);
        nMod = length(clustNames{k});
        RR(dL, k) = hygepdf(nModDG, N, nDG, nMod);
        clustSize(k) = length(clustNames{k});
        interG(k) = nModDG;
    end
%     tempP = RR(dL,:); 
%     pVals = mafdr(tempP);
% %     pVals = mafdr(tempP, 'Method', 'polynomial');
%     RR(dL,:) = pVals;
%     clear tempP; clear pVals;
end

% save([outF '_enrichmentScores.mat'], 'RR');
% csvwrite([outF '_enrichmentScores.csv'], RR);

xlswrite([outF dName '_enrichmentScores.xls'], RR);
xlswrite([outF dName '_clustSize.xls'], clustSize);
xlswrite([outF dName '_asdGenesInModules.xls'], interG);

f = figure, bar([interG; clustSize]'); grid on
legend({'ASD genes in Module', 'All genes in Modules'});
saveas(f, [outF dName '_ModulesSize.fig']);
saveas(f, [outF dName '_ModulesSize.jpg']);
clear f;

f = figure;
bar((-1*log(RR')), 'EdgeColor', 'black');
colormap(autumn)
grid on
% line([0 length(RR)+1], [(-1*log(0.0001)) (-1*log(0.0001))])
% line([0 length(RR)+1], [(-1*log(0.05/(length(clustNames)))) (-1*log(0.05/(length(clustNames))))])
title('Clusters enrichment of disease-related genes', 'fontweight', 'bold')
xlabel('Clusters', 'fontweight', 'bold')
ylabel('P-values', 'fontweight', 'bold')
set(gca, 'XTick', [1:length(clustNames)])
axis([0 size(RR,2)+1 0 20])
legend(dList)
% legend({'ASDs Only', 'SZ Only', 'Common'})
saveas(f, [outF dName '_p-values.fig']);
saveas(f, [outF dName '_p-values.jpg']);
% close(f);
    
% f2 = figure;
% bar(L', 'EdgeColor', 'black');
% colormap(autumn)
% grid on
% % line([0 length(RR)+1], [(-1*log(0.05)) (-1*log(0.05))])
% title('Distribution of disease-related genes over clusters', 'fontweight', 'bold')
% xlabel('Clusters', 'fontweight', 'bold')
% ylabel('Number of Genes', 'fontweight', 'bold')
% set(gca, 'XTick', [1:length(clustNames)])
% % axis([0 size(RR,2)+1 0 max(max(RR))+0.001])
% % legend(dList)
% legend({'ASDs Only', 'SZ Only', 'Common'})











