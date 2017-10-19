function enrichAnalyze_cellType(clustIDs, ctList, D, outF)

load('files/genesStatus_5RPKM.mat');
gIND = find(genesStatus_5RPKM == 1);
N = length(gIND);
clear genesStatus_5RPKM;

[num txt] = xlsread('files\DataFiles\cellTypeEnrichedGenes.xlsx');
neo = num(:,1); n1 = find(isnan(neo) == 1); neo(n1) = [];
oligo = num(:,2); n2 = find(isnan(oligo) == 1); oligo(n2) = [];
astro = num(:,3); n3 = find(isnan(astro) == 1); astro(n3) = [];
ctEIDs{1} = neo;
ctEIDs{2} = oligo;
ctEIDs{3} = astro;
clear num; clear txt;

%%% Cell Type Enrichment
for dL = 1 : length(ctList)
    clear ctName; clear current_ctEIDs; 
    ctName = ctList{dL};
    current_ctEIDs = ctEIDs{dL};
    %%% calculate RR score for each cluster    
    for k = 1 : length(clustIDs)
        % find the cell-type enriched genes in each cluster
        ctGenesInMod = intersect(current_ctEIDs, clustIDs{k});
        nModDG = length(ctGenesInMod);
        nMod = length(clustIDs{k});
        nDG = length(current_ctEIDs);
        RR(dL, k) = hygepdf(nModDG, N, nDG, nMod);
        clustSize(k) = length(clustIDs{k});
    end
    tempP = RR(dL,:); 
    pVals = tempP;
%     pVals = mafdr(tempP, 'Method', 'polynomial');
    RR(dL,:) = pVals;
    clear tempP; clear pVals;
end

csvwrite([outF 'cellType_enrichmentScores.csv'], RR);

f = figure;
bar((-1*log(RR')), 'EdgeColor', 'black');
colormap(autumn)
grid on
% line([0 length(RR)+1], [(-1*log(0.0001)) (-1*log(0.0001))])
% line([0 length(RR)+1], [(-1*log(0.05/(length(clustNames)))) (-1*log(0.05/(length(clustNames))))])
title([D ' Clusters enrichment of cell-type enriched genes'], 'fontweight', 'bold')
xlabel('Clusters', 'fontweight', 'bold')
ylabel('P-values', 'fontweight', 'bold')
set(gca, 'XTick', [1:length(clustIDs)])
axis([0 size(RR,2)+1 0 40])
legend(ctList)
legend({'Neurons', 'Oligodendrytes', 'Astrocytes'})
saveas(f, ['results\corrChange\'  D '_cellType_p-values.fig']);
saveas(f, ['results\corrChange\'  D '_cellType_p-values.jpg']);
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











