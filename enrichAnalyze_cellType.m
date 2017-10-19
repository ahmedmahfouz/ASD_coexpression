function enrichAnalyze_cellType(clusts, ctList, outF)

load('results\hClustResults_NEW\exp2\allSamples_corrDist.mat'); 
gCorr = squareform(corrDist);

load('files/gNames.mat');
load('files/genesStatus_5RPKM.mat');
load('files/entrezIDs.mat');
load('files/enIDs.mat');
gIND = find(genesStatus_5RPKM == 1);
gNames_5RPKM = gNames(gIND);
entrezIDs_5RPKM = entrezIDs(gIND);
enIDs_5RPKM = enIDs(gIND);

N = length(gIND);
clear gNames;
clear genesStatus_5RPKM;

%%% extract each cluster genes' entrezIDs
for j = 1 : max(clusts)
    clust = find(clusts == j);
    clustIDs{j} = entrezIDs_5RPKM(clust);
    clear clust;
end

[num txt] = xlsread('files\DataFiles\cellTypeEnrichedGenes.xlsx');
neo = num(:,1); n1 = find(isnan(neo) == 1); neo(n1) = [];
oligo = num(:,2); n2 = find(isnan(oligo) == 1); oligo(n2) = [];
astro = num(:,3); n3 = find(isnan(astro) == 1); astro(n3) = [];
[neo_2 ia ib] = intersect(neo, entrezIDs_5RPKM);
[oligo_2 ja jb] = intersect(oligo, entrezIDs_5RPKM);
[astro_2 ka kb] = intersect(astro, entrezIDs_5RPKM);
ctEIDs{1} = neo_2;
ctEIDs{2} = oligo_2;
ctEIDs{3} = astro_2;
clear num; clear txt;

Nnames = gNames_5RPKM(ib);
Nen = enIDs_5RPKM(ib);
Onames = gNames_5RPKM(jb);
Oen = enIDs_5RPKM(jb);
Anames = gNames_5RPKM(kb);
Aen = enIDs_5RPKM(kb);
fname = 'files\DataFiles\cellTypeGenes.xls';
xlswrite(fname, Nnames, 1, 'A2');
xlswrite(fname, Nen, 1, 'B2');
xlswrite(fname, Onames, 1, 'C2');
xlswrite(fname, Oen, 1, 'D2');
xlswrite(fname, Anames, 1, 'E2');
xlswrite(fname, Aen, 1, 'F2');

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
%     pVals = tempP;
    pVals = mafdr(tempP);%, 'Method', 'polynomial');
    RR(dL,:) = pVals;
    clear tempP; clear pVals;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [num txt] = xlsread(['Data\Donors\ASDs.xls']);
% asd = txt(:,1);
% clear num; clear txt;
% [num txt] = xlsread(['Data\Donors\SZ.xls']);
% sz = txt(:,1);
% clear num; clear txt;
% com = intersect(asd, sz);
% asdO = setdiff(asd, sz);
% szO = setdiff(sz, asd);
% dgNames{1} = asdO;
% dgNames{2} = szO;
% dgNames{3} = com;
% 
% for dL = 1 : 3
%     
%     %%% calculate RR score for each cluster    
%     for k = 1 : length(clustNames)
%         clear dGenesInMod; clear nModDG; clear nMod; clear nDG; 
%         % find the disease genes in each cluster
%         dGenesInMod = intersect(dgNames{dL}, clustNames{k});
%         nModDG = length(dGenesInMod);
%         nMod = length(clustNames{k});
%         nDG = length(dgNames{dL});
%         RR(dL, k) = hygepdf(nModDG, N, nDG, nMod);
%         clustSize(k) = length(clustNames{k});
%         L(dL,k) = length(dGenesInMod);
%     end
% %     tempP = RR(dL,:); 
% % %     pVals = mafdr(tempP);
% %     pVals = mafdr(tempP, 'Method', 'polynomial');
% %     RR(dL,:) = pVals;
% %     clear tempP; clear pVals;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save([outF '_enrichmentScores.mat'], 'RR');
csvwrite([outF 'cellType_enrichmentScores.csv'], RR);

f = figure;
bar((-1*log(RR')), 'EdgeColor', 'black');
colormap(autumn)
grid on
% line([0 length(RR)+1], [(-1*log(0.0001)) (-1*log(0.0001))])
% line([0 length(RR)+1], [(-1*log(0.05/(length(clustNames)))) (-1*log(0.05/(length(clustNames))))])
title('Clusters enrichment of cell-type enriched genes', 'fontweight', 'bold')
xlabel('Clusters', 'fontweight', 'bold')
ylabel('P-values', 'fontweight', 'bold')
set(gca, 'XTick', [1:length(clustIDs)])
axis([0 size(RR,2)+1 0 220])
legend(ctList)
legend({'Neurons', 'Oligodendrytes', 'Astrocytes'})
saveas(f, 'results\hClustResults_NEW\exp2\clustersEnrich\cellType_p-values.fig');
saveas(f, 'results\hClustResults_NEW\exp2\clustersEnrich\cellType_p-values.jpg');
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











